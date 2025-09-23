
#include <poll.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <arpa/inet.h>
#include <errno.h>
#include <netdb.h>
#include <netinet/in.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <x86intrin.h>

#include "client_loop.h"
#include "inner_loop.h"
#include "lcg.h"
#include "proxy_loop.h"
#include "setup.h"
#include "tlv.h"

static char server[128];
static unsigned short server_port;

#define STATE_UNUSED 0
#define STATE_PENDING 1
#define STATE_RUNNING 2
#define STATE_DONE 3
#define STATE_DEAD 4

#define MAX_BLOCK (MAX_CID * 32)          // power of 2
#define BLOCK_TIME (120000000000ull)      // approx 10 seconds
#define RATE (10000ull)                   // approx 1 check per 10000 ticks
#define BLOCK_TIMEOUT (2ull * BLOCK_TIME) // approx 20 seconds
#define INIT_SEED (1)

typedef struct connection_s
{
    uint64_t rate;
    uint64_t t_connect;
    uint64_t progress;
    int s;
    unsigned state;
    unsigned cid;
} connection_t;

typedef struct progress_s
{
    uint128_t seed;
    uint64_t count;
    uint64_t t_start;
    uint64_t expected_t_end;
    unsigned cid;
    unsigned state;
} progress_t;

static connection_t connections[MAX_CID];
static progress_t progress[MAX_BLOCK];
static uint64_t head = MAX_CID;
static uint64_t tail = MAX_CID;
static uint128_t done_count = 0;

static Lcg lll;

#if 0
#if __clang__
#if __clang_major__ < 4
static uint64_t __rdtsc(void)
{
	unsigned a, d;
	__asm volatile ("rdtsc":"=a" (a), "=d"(d));
	return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}
#endif
#endif
#endif

static void report(const char *s, uint128_t v)
{
    // display
    print128(s, v);

    // append to log file
    FILE *l = fopen("lnrc.log", "at");
    if (l)
    {
        print128(s, v, l);
        fclose(l);
    }
}

static void connections_init(void)
{
    int i;
    for (i = 0; i < MAX_CID; i++)
    {
        connections[i].state = STATE_UNUSED;
        connections[i].progress = 0;
        connections[i].cid = i;
        connections[i].s = -1;
        connections[i].rate = 20 * RATE;
    }
    head = MAX_CID;
    tail = MAX_CID;
}

// check if a socket is unexpectedly closed, This can happen from network disruption, tcp timeout
// Enter reality : the impact of a broken socket can include SIG_PIPE, read return -1, or return 0 and errno is set, or
// .....
static bool is_closed_socket(int fd)
{
    struct pollfd pfd;
    int ret;

    pfd.fd = fd;
    pfd.events = POLLRDHUP;
    pfd.revents = 0;
    ret = poll(&pfd, 1, 0);
    if (ret > 0 && (pfd.revents & (POLLHUP | POLLERR | POLLRDHUP)))
    {
        return true;
    }
    else
    {
        return false;
    }
}

// close all connections related to socket s
static void set_broken_socket(int s)
{
    int i;
    for (i = 0; i < MAX_CID; i++)
    {
        if (connections[i].s == s && connections[i].state != STATE_UNUSED)
        {
            if (connections[i].state == STATE_RUNNING)
            {
                progress_t *pg = &progress[connections[i].progress % MAX_BLOCK];
                if (pg->state == STATE_RUNNING)
                {
                    pg->state = STATE_DEAD;
                }
            }
            connections[i].state = STATE_UNUSED;
            printf("End cid %d\n", i);
        }
    }
}

static void set_timeout(void)
{
    uint64_t current_t = __rdtsc();

    while (tail < head)
    {
        uint64_t i = tail % MAX_BLOCK;
        if (progress[i].state != STATE_DONE)
            break;
        tail++;
    }

    // printf("%lu %lu\n", tail, head);
    if (tail < head)
    {
        uint64_t i = tail % MAX_BLOCK;
        if (progress[i].state == STATE_RUNNING)
        {
            if (current_t > progress[i].expected_t_end + BLOCK_TIMEOUT)
            {
                unsigned cid = progress[i].cid;
                progress[i].state = STATE_DEAD;
                connections[cid].state = STATE_DEAD;
                printf("Server Timeout Cid %d\n", cid);
                fflush(stdout);
            }
        }
    }
}

static uint64_t get_count_from_rate(uint64_t rate)
{
    uint64_t ncount = (BLOCK_TIME + rate - 1) / rate;
    ncount |= 1;
    ncount += 3;
    return ncount;
}

static void get_next(uint128_t *seed, uint64_t *count, uint64_t rate)
{
    uint64_t ncount;
    int j, i;
    for (j = tail; j < head; j++)
    {
        i = j % MAX_BLOCK;
        if (progress[i].state == STATE_DEAD && progress[i].count > 0)
        {
            *seed = progress[i].seed;
            ncount = get_count_from_rate(rate);
            if (ncount < progress[i].count)
            {
                Lcg u;
                u.set_seed(progress[i].seed);
                progress[i].seed = u.get_seed(ncount);
                progress[i].count -= ncount;
            }
            else
            {
                ncount = progress[i].count;
                progress[i].count = 0;
                progress[i].state = STATE_DONE;
            }
            *count = ncount;
            return;
        }
    }

    *seed = lll.get_seed(0);
    ncount = get_count_from_rate(rate);
    lll.get_seed(ncount);
    *count = ncount;
    return;
}

static void display_progress()
{
    if (lll.sequential())
    {
        bool update = false;
        int j, i;
        uint128_t min_seed = lll.get_seed(0);
        for (j = tail; j < head; j++)
        {
            i = j % MAX_BLOCK;
            if (progress[i].seed < min_seed)
            {
                min_seed = progress[i].seed;
                update = true;
            }
        }

        time_t now = time(NULL);
        static time_t last_display = 0;
        static uint128_t last_seed = 0;
        if (update && now != last_display && last_seed != min_seed)
        {
            last_display = now;
            last_seed = min_seed;
            char buff[60];
            time_t rawtime = time(NULL);
            ctime_r(&rawtime, buff);
            int len = strlen(buff);
            while (len > 0 && isspace(buff[len - 1]))
            {
                len--;
            }
            buff[len++] = 0;
            min_seed = convert_seed_to_number(min_seed);
            strcat(buff, " Completed");
            print128(buff, min_seed);
            fflush(stdout);
        }
    }
    else
    {
        time_t now = time(NULL);
        static time_t last_display = 0;
        if (now - last_display > 5)
        {
            last_display = now;
            char buff[60];
            time_t rawtime = time(NULL);
            ctime_r(&rawtime, buff);
            int len = strlen(buff);
            while (len > 0 && isspace(buff[len - 1]))
            {
                len--;
            }
            buff[len++] = 0;
            strcat(buff, " Completed count");
            print128(buff, done_count);
            fflush(stdout);
        }
    }
}

void *server_thread(void *arg)
{
    int max_sd = 0, listen_sd = 0, conn_sd = 0;
    int rc, i, j;
    struct sockaddr_in serv_addr;
    fd_set read_sds, master_sds, except_sds;
    struct timeval timeout;
    uint64_t current_t;
    uint64_t current_p;

    connections_init();

    listen_sd = socket(AF_INET, SOCK_STREAM, 0);
    memset(&serv_addr, '0', sizeof(serv_addr));

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
    serv_addr.sin_port = htons(server_port);

    rc = bind(listen_sd, (struct sockaddr *)&serv_addr, sizeof(serv_addr));
    if (rc < 0)
    {
        perror("bind");
    }
    rc = listen(listen_sd, 10);
    if (rc < 0)
    {
        perror("listen");
    }

    FD_ZERO(&master_sds);
    FD_SET(listen_sd, &master_sds);
    max_sd = listen_sd;

    while (1)
    {
        set_timeout();
        memcpy(&read_sds, &master_sds, sizeof(read_sds));
        memcpy(&except_sds, &master_sds, sizeof(except_sds));
        timeout.tv_sec = 3 * 60;
        timeout.tv_usec = 0;
        rc = select(max_sd + 1, &read_sds, NULL, &except_sds, &timeout);
        if (rc < 0)
        {
            perror("select");
            exit(1);
        }
        if (rc > 0)
        {
            for (i = 0; i <= max_sd; i++)
            {
                if (FD_ISSET(i, &except_sds))
                {
                    if (i == listen_sd)
                    {
                        printf("Exception listen socket %d\n", i);
                        fflush(stdout);
                        exit(1);
                    }
                    printf("Exception connection %d\n", i);
                    fflush(stdout);
                    set_broken_socket(i);
                    FD_CLR(i, &master_sds);
                    FD_CLR(i, &read_sds);
                    close(i);
                }
            }
            for (i = 0; i <= max_sd; i++)
            {
                if (FD_ISSET(i, &read_sds))
                {
                    if (i == listen_sd)
                    {
                        struct sockaddr_in peer_addr;
                        socklen_t peer_addr_size;
                        char peer_name[INET_ADDRSTRLEN];

                        peer_addr_size = sizeof(struct sockaddr_in);
                        memset(&peer_addr, 0, peer_addr_size);
                        conn_sd = accept(listen_sd, (struct sockaddr *)&peer_addr, &peer_addr_size);
                        if (conn_sd < 0)
                        {
                            if (errno != EWOULDBLOCK)
                            {
                                goto err;
                            }
                            continue;
                        }
                        FD_SET(conn_sd, &master_sds);
                        if (conn_sd > max_sd)
                        {
                            max_sd = conn_sd;
                        }
                        inet_ntop(AF_INET, &peer_addr.sin_addr, peer_name, INET_ADDRSTRLEN);
                        printf("New connection %d from %s\n", conn_sd, peer_name);
                        fflush(stdout);
                    }
                    else
                    {
                        uint8_t t;
                        uint128_t v, v_seed;
                        uint64_t v_count;
                        uint16_t cid;
                        int n = 0;
                        rc = ioctl(i, FIONREAD, &n);
                        if (rc < 0)
                        {
                            printf("Broken connection %d\n", i);
                            fflush(stdout);
                        }
                        else if (n == 0)
                        {
                            if (is_closed_socket(i))
                            {
                                printf("Closed connection %d\n", i);
                                fflush(stdout);
                                rc = -1;
                            }
                            else
                            {
                                continue;
                            }
                        }
                        else
                        {
                            rc = tlv_read(i, &cid, &t, &v);
                            if (rc < 0)
                            {
                                printf("Invalid read on connection %d\n", i);
                                fflush(stdout);
                            }
                        }
                        if (rc < 0)
                        {
                            set_broken_socket(i);
                            FD_CLR(i, &master_sds);
                            close(i);
                            continue;
                        }
                        switch (t)
                        {
                        case TLV_NEW:
                            for (j = 0; j < MAX_CID; j++)
                            {
                                if (connections[j].state == STATE_UNUSED)
                                {
                                    current_t = __rdtsc();
                                    connections[j].state = STATE_PENDING;
                                    connections[j].s = i;
                                    connections[j].rate = RATE * 20;
                                    connections[j].progress = 0;
                                    connections[j].t_connect = current_t;
                                    rc = tlv_write(i, j, TLV_NEW, 0);
                                    printf("New cid %d\n", connections[j].cid);
                                    if (rc < 0)
                                    {
                                        set_broken_socket(i);
                                        FD_CLR(i, &master_sds);
                                        close(i);
                                    }
                                    break;
                                }
                            }
                            break;
                        case TLV_STOP:
                            if (connections[cid].state == STATE_RUNNING)
                            {
                                progress[connections[cid].progress % MAX_BLOCK].state = STATE_DEAD;
                            }
                            connections[cid].state = STATE_UNUSED;
                            connections[cid].progress = 0;
                            break;
                        case TLV_READY:
                            if (connections[cid].progress && connections[cid].state == STATE_RUNNING)
                            {
                                progress_t *pg = &progress[connections[cid].progress % MAX_BLOCK];
                                current_t = __rdtsc();
                                connections[cid].rate = (current_t - pg->t_start) / pg->count;
                                pg->state = STATE_DONE;
                                done_count += pg->count;
                            }
                            else
                            {
                                connections[cid].rate = RATE * 20;
                            }
                            get_next(&v_seed, &v_count, connections[cid].rate);
                            current_t = __rdtsc();
                            current_p = head % MAX_BLOCK;
                            progress[current_p].seed = v_seed;
                            progress[current_p].count = v_count;
                            progress[current_p].t_start = current_t;
                            progress[current_p].expected_t_end = v_count * connections[cid].rate + current_t;
                            progress[current_p].cid = cid;
                            progress[current_p].state = STATE_RUNNING;
                            connections[cid].state = STATE_RUNNING;
                            connections[cid].progress = head;
                            rc = tlv_write(i, cid, TLV_SEED, progress[current_p].seed);
                            rc = rc < 0 ? rc : tlv_write(i, cid, TLV_COUNT, progress[current_p].count);
                            rc = rc < 0 ? rc : tlv_write(i, cid, TLV_GO, 0);
                            if (rc < 0)
                            {
                                set_broken_socket(i);
                                FD_CLR(i, &master_sds);
                                close(i);
                                continue;
                            }
                            head++;
                            display_progress();
                            break;
                        case TLV_PSEUDOPRIME:
                            report("Pseudoprime", v);
                            break;
                        case TLV_PSEUDOCOMPOSITE:
                            report("Pseudocomposite", v);
                            break;
                        case TLV_B1:
                            report("B == 1", v);
                            break;
                        case TLV_SEED:
                        case TLV_COUNT:
                        default:
                            break;
                        }
                    }
                }
            }
        }
    }

err:
    for (i = 0; i < max_sd; i++)
    {
        if (FD_ISSET(i, &master_sds))
        {
            if (i != listen_sd)
                tlv_write(i, 0, TLV_STOP, 0);
            FD_CLR(i, &master_sds);
            close(i);
        }
    }
    max_sd = 0;
    return 0;
}

int main(int argc, char **argv)
{
    pthread_t tid;
    long i;
    long t = 1;
    bool start_server = false;
    bool start_proxy = false;
    bool start_client = false;
    unsigned short port = SERVER_PORT;

    strcpy(server, "127.0.0.1");

    for (i = 1; i < (long)argc; i++)
    {
        if (!strcmp(argv[i], "-t"))
        {
            t = atol(argv[++i]);
            start_client = true;
            continue;
        }
        if (!strcmp(argv[i], "-p"))
        {
            port = atol(argv[++i]);
            continue;
        }
        if (!strcmp(argv[i], "-s"))
        {
            strcpy(server, argv[++i]);
            continue;
        }
        if (!strcmp(argv[i], "-e"))
        {
            scan128(argv[++i], &done_count);
            lll.set_seed(INIT_SEED);
            lll.get_seed(convert_number_to_seed(done_count));
            continue;
        }
        if (!strcmp(argv[i], "-server"))
        {
            start_server = true;
            continue;
        }
        if (!strcmp(argv[i], "-proxy"))
        {
            start_proxy = true;
            continue;
        }
        if (!strcmp(argv[i], "-st"))
        {
            printf("self-test %s\n", (inner_self_test() == 0) ? "passed" : "failed");
            exit(1);
        }
        if (!strcmp(argv[i], "-h"))
        {
            printf("Usage %s -s remote_server -p remote_port [-t client_thread_count] [-e next] [-server] [-proxy]\n",
                   argv[0]);
            exit(1);
        }
    }

    if (start_server)
    {
        printf("Starting server ....\n");
        server_port = port;
        pthread_create(&tid, 0, server_thread, 0);
        sleep(1);
    }
    if (start_proxy)
    {

        proxy_setup("127.0.0.1", port + 1, server, port);

        printf("Starting proxy ....\n");
        pthread_create(&tid, 0, proxy_thread, 0);
        sleep(1);
    }

    if (start_client)
    {
        if (start_proxy)
        {
            client_setup("127.0.0.1", port + 1);
        }
        else
        {
            client_setup(server, port);
        }

        for (i = 0; i < t - 1; i++)
        {
            pthread_create(&tid, 0, client_thread, 0);
            if (i < t - 4)
            {
                i += 3;
                pthread_create(&tid, 0, client_thread, 0);
                pthread_create(&tid, 0, client_thread, 0);
                pthread_create(&tid, 0, client_thread, 0);
            }
            sleep(1);
        }
        client_thread(0);
    }

    while (1)
        sleep(1);
}
