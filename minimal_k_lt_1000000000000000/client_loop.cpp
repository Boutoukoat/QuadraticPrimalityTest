
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "client_loop.h"
#include "inner_loop.h"
#include "setup.h"
#include "tlv.h"

static int client_outer_loop(int s)
{
    int r;
    uint8_t t;
    uint16_t cid = 0;
    uint64_t count;
    uint128_t v, seed;

    r = tlv_write(s, cid, TLV_NEW, 0);
    if (r < 0)
        return -1;
    r = tlv_read(s, &cid, &t, &v);
    if (r < 0)
        return -1;
    if (t == TLV_STOP)
        return 0;
    printf("Client new cid %d\n", (int)cid);
    fflush(stdout);
    r = tlv_write(s, cid, TLV_READY, 0);
    if (r < 0)
        return -1;
    while (1)
    {
        r = tlv_read(s, &cid, &t, &v);
        if (r < 0)
            return -1;
        switch (t)
        {
        case TLV_STOP:
            return 0;
        case TLV_SEED:
            seed = v;
            break;
        case TLV_COUNT:
            count = v;
            break;
        case TLV_GO:
            inner_loop(s, cid, seed, count);
            r = tlv_write(s, cid, TLV_READY, 0);
            if (r < 0)
                return -1;
            break;
        default:
            printf("Worker Error : Unknown tlv\n");
            return -1;
        }
    }
}

static char partner_addr[128] = "127.0.0.1";
static unsigned short partner_port = PROXY_PORT;

void client_setup(const char *remote_addr, unsigned short remote_port)
{
    strcpy(partner_addr, remote_addr);
    partner_port = remote_port;
}

void *client_thread(void *arg)
{
    int conn_sd = -1;
    int ret_code = 0;
    struct sockaddr_in serv_addr;
    int back_off = 0;

    while (1)
    {

        ret_code = 0;

        /* a socket is created through call to socket() function */
        if ((conn_sd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
        {
            printf("Worker Error : Could not create socket \n");
            goto done;
        }

        memset(&serv_addr, '0', sizeof(serv_addr));

        serv_addr.sin_family = AF_INET;
        serv_addr.sin_port = htons(partner_port);

        if (inet_pton(AF_INET, partner_addr, &serv_addr.sin_addr) <= 0)
        {
            printf("Worker Error : inet_pton error occured\n");
            goto err;
        }

        if (connect(conn_sd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
        {
            printf("Worker Error : Connect Failed \n");
            goto err;
        }

        ret_code = client_outer_loop(conn_sd);

    err:
        if (conn_sd > 0)
            close(conn_sd);
        conn_sd = -1;
    done:
        if (ret_code != 0)
        {
            // a protocol/disconnection error occured, must retry
            sleep(back_off % 60);
            back_off = (back_off * 3) / 2 + 1;
        }
        else
        {
            // end of thread
            return 0;
        }
    }
}
