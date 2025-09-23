
#include <poll.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <arpa/inet.h>
#include <errno.h>
#include <netdb.h>
#include <netinet/in.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "proxy_loop.h"
#include "setup.h"
#include "tlv.h"

static char server_addr[128];
static unsigned short server_port = SERVER_PORT;

static char proxy_addr[128] = "127.0.0.1";
static unsigned short proxy_port = PROXY_PORT;

void proxy_setup(const char *local_addr, unsigned short local_port, const char *remote_addr, unsigned short remote_port)
{
    strcpy(proxy_addr, local_addr);
    proxy_port = local_port;
    strcpy(server_addr, remote_addr);
    server_port = remote_port;
}

static bool closed_socket(int fd)
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

static int proxy_connect(void)
{
    int server_sd;
    struct sockaddr_in serv_addr;

    /* connect to server */
    if ((server_sd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
        printf("Error : Could not create socket \n");
        goto done;
    }

    memset(&serv_addr, '0', sizeof(serv_addr));

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(server_port);

    if (inet_pton(AF_INET, server_addr, &serv_addr.sin_addr) <= 0)
    {
        printf("Error : inet_pton error occured\n");
        goto err;
    }

    if (connect(server_sd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
    {
        printf("Error : Connect Failed \n");
        goto err;
    }
    printf("Proxy : New connection with server %d\n", server_sd);
    fflush(stdout);
    return server_sd;

err:
    close(server_sd);
done:
    return -1;
}

typedef struct proxy_s
{
    int s;
    bool in_use;
} proxy_t;

static proxy_t proxy[MAX_CID];

typedef struct new_conn_s
{
    int s;
    bool in_use;
} new_conn_t;

#define MAX_NEW_CONN 300
static new_conn_t new_conn[MAX_NEW_CONN];

void *proxy_thread(void *arg)
{
    int max_sd = 0, listen_sd = 0, conn_sd = 0, server_sd = -1;
    int rc, i, j;
    struct sockaddr_in proxy_addr;
    fd_set read_sds, master_sds, except_sds;
    struct timeval timeout;
    uint8_t t;
    uint128_t v;
    uint16_t cid;

    while (1)
    {
        for (j = 0; j < MAX_NEW_CONN; j++)
        {
            new_conn[j].in_use = false;
        }

        /* connect to server */
        server_sd = proxy_connect();
        if (server_sd < 0)
        {
            printf("Proxy : No connection with server\n");
            fflush(stdout);
            goto done;
        }

        /* setup as a server */
        listen_sd = socket(AF_INET, SOCK_STREAM, 0);
        if (listen_sd < 0)
        {
            goto done;
        }
        memset(&proxy_addr, '0', sizeof(proxy_addr));

        proxy_addr.sin_family = AF_INET;
        proxy_addr.sin_addr.s_addr = htonl(INADDR_ANY);
        proxy_addr.sin_port = htons(proxy_port);

        rc = bind(listen_sd, (struct sockaddr *)&proxy_addr, sizeof(proxy_addr));
        if (rc < 0)
        {
            perror("bind");
            goto done;
        }
        rc = listen(listen_sd, 10);
        if (rc < 0)
        {
            perror("listen");
            goto done;
        }

        FD_ZERO(&master_sds);
        FD_SET(listen_sd, &master_sds);
        FD_SET(server_sd, &master_sds);
        max_sd = listen_sd;
        if (server_sd > max_sd)
        {
            max_sd = server_sd;
        }

        while (1)
        {
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
                        if (i == server_sd)
                        {
                            printf("Proxy: Exception server socket %d\n", i);
                            fflush(stdout);
                            goto err;
                        }
                        if (i == listen_sd)
                        {
                            printf("Proxy: Exception listen socket %d\n", i);
                            fflush(stdout);
                            goto err;
                        }
                        printf("Proxy: Exception connection %d\n", i);
                        fflush(stdout);
                        for (j = 0; j < MAX_CID; j++)
                        {
                            if (proxy[j].in_use == true && proxy[j].s == i)
                            {
                                rc = tlv_write(server_sd, j, TLV_STOP, 0);
                                proxy[j].in_use = false;
                            }
                        }
                        close(i);
                        FD_CLR(i, &master_sds);
                        FD_CLR(i, &read_sds);
                    }
                }

                for (i = 0; i <= max_sd; i++)
                {
                    if (FD_ISSET(i, &read_sds))
                    {
                        if (i == listen_sd)
                        {
                            conn_sd = accept(listen_sd, (struct sockaddr *)NULL, NULL);
                            if (conn_sd < 0)
                            {
                                if (errno != EWOULDBLOCK)
                                {
                                    goto err;
                                }
                                continue;
                            }
                            for (j = 0; j < MAX_NEW_CONN; j++)
                            {
                                if (new_conn[j].in_use == false)
                                {
                                    new_conn[j].in_use = true;
                                    new_conn[j].s = conn_sd;
                                    break;
                                }
                            }
                            if (j == MAX_NEW_CONN)
                            {
                                close(conn_sd);
                                continue;
                            }

                            FD_SET(conn_sd, &master_sds);
                            if (conn_sd > max_sd)
                            {
                                max_sd = conn_sd;
                            }
                            printf("Proxy : New connection %d\n", conn_sd);
                            fflush(stdout);
                        }
                        else if (i == server_sd)
                        {
                            int n = 0;
                            rc = ioctl(i, FIONREAD, &n);
                            if (rc < 0)
                            {
                                rc = -1;
                                printf("Proxy : End connection %d\n", i);
                                fflush(stdout);
                            }
                            else if (n == 0)
                            {
                                continue;
                            }
                            else
                            {
                                rc = tlv_read(i, &cid, &t, &v);
                            }
                            if (rc < 0)
                            {
                                goto err;
                            }
                            if (proxy[cid].in_use == false)
                            {
                                for (j = 0; j < MAX_NEW_CONN; j++)
                                {
                                    if (new_conn[j].in_use == true)
                                    {
                                        proxy[cid].s = new_conn[j].s;
                                        proxy[cid].in_use = true;
                                        new_conn[j].in_use = false;
                                        break;
                                    }
                                }
                                if (j == MAX_NEW_CONN)
                                {
                                    rc = tlv_write(i, cid, TLV_STOP, 0);
                                    goto err;
                                }
                            }
                            rc = tlv_write(proxy[cid].s, cid, t, v);
                            if (rc < 0)
                            {
                                printf("Proxy : write to client error on connection %d\n", proxy[cid].s);
                                fflush(stdout);
                                for (j = 0; j < MAX_CID; j++)
                                {
                                    if (proxy[j].in_use == true && proxy[j].s == proxy[cid].s)
                                    {
                                        rc = tlv_write(i, j, TLV_STOP, 0);
                                        proxy[j].in_use = false;
                                    }
                                }
                                close(proxy[cid].s);
                                FD_CLR(proxy[cid].s, &master_sds);
                            }
                        }
                        else
                        {
                            int n = 0;
                            rc = ioctl(i, FIONREAD, &n);
                            if (rc < 0)
                            {
                                rc = -1;
                                printf("Proxy : broken connection %d\n", i);
                                fflush(stdout);
                            }
                            else if (n == 0)
                            {
                                if (closed_socket(i))
                                {
                                    rc = -1;
                                    printf("Proxy : closed connection %d\n", i);
                                    fflush(stdout);
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
                                    printf("Proxy : read error on connection %d\n", i);
                                    fflush(stdout);
                                }
                            }
                            if (rc < 0)
                            {
                                for (j = 0; j < MAX_CID; j++)
                                {
                                    if (proxy[j].in_use == true && proxy[j].s == i)
                                    {
                                        rc = tlv_write(server_sd, j, TLV_STOP, 0);
                                        proxy[j].in_use = false;
                                    }
                                }
                                close(i);
                                FD_CLR(i, &master_sds);
                                continue;
                            }
                            rc = tlv_write(server_sd, cid, t, v);
                            if (rc < 0)
                            {
                                printf("Proxy : write to server error on connection %d\n", server_sd);
                                fflush(stdout);
                                goto err;
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
                FD_CLR(i, &master_sds);
                if (i != listen_sd && i != server_sd)
                {
                    tlv_write(i, 0, TLV_STOP, 0);
                }
                close(i);
            }
        }
    done:
        max_sd = 0;
        if (listen_sd > 0)
            close(listen_sd);
        listen_sd = -1;
        if (server_sd > 0)
            close(server_sd);
        server_sd = -1;
        conn_sd = 0;

        sleep(1);
    }

    return 0;
}
