
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>

#include "tlv.h"

int block_read(int s, unsigned len, uint8_t *p)
{
    unsigned l = 0;
    while (l < len)
    {
        int r = read(s, &p[l], len - l);
        if (r < 0)
        {
            perror("read");
            return -1;
        }
        l += r;
    }
    return len;
}

int block_write(int s, unsigned len, uint8_t *p)
{
    unsigned l = 0;
    while (l < len)
    {
        int r = write(s, &p[l], len - l);
        if (r < 0)
        {
            perror("read");
            return -1;
        }
        l += r;
    }
    return len;
}

int tlv_read(int s, uint16_t *cid, uint8_t *type, uint128_t *value)
{
    uint8_t *p;
    uint128_t v;
    uint16_t l;
    p = (uint8_t *)alloca(5);
    int len = block_read(s, 5, p);
    if (len < 5)
        return -1;
    *type = p[0];
    *cid = p[2] * 256 + p[1];
    l = p[4] * 256 + p[3];
    if (l < 1)
        return -1;
    p = (uint8_t *)alloca(l);
    len = block_read(s, l, p);
    if (len < (int)l)
        return -1;
    v = 0;
    while (l--)
    {
        v = v * 256 + p[l];
    }
    *value = v;
    return 0;
}

int tlv_write(int s, uint16_t cid, uint8_t type, uint128_t value)
{
    uint8_t *p;
    uint128_t cmp, v = value;
    uint16_t i, l;
    l = 1;
    cmp = 0xff;
    while (l < 16 && v > cmp)
    {
        cmp = (cmp << 8) + 0xff;
        l++;
    }
    p = (uint8_t *)alloca(l + 4);
    p[0] = type;
    p[1] = cid % 256;
    p[2] = cid / 256;
    p[3] = l % 256;
    p[4] = l / 256;
    int len = block_write(s, 5, p);
    if (len < 5)
        return -1;
    for (i = 0; i < l; i++)
    {
        p[i] = v & 0xff;
        v >>= 8;
    }
    len = block_write(s, l, p);
    if (len < (int)l)
        return -1;
    return 0;
}
