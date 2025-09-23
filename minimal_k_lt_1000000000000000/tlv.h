
#ifndef TLV_H_H
#define TLV_H_H

#include <ctype.h>
#include <stdio.h>

#define TLV_SEED 1
#define TLV_COUNT 2
#define TLV_STOP 3
#define TLV_GO 4
#define TLV_PSEUDOCOMPOSITE 10
#define TLV_PSEUDOPRIME 11
#define TLV_B1 20
#define TLV_READY 12
#define TLV_NEW 13

typedef unsigned __int128 uint128_t;

static inline void print128(const char *s, const uint128_t v, FILE *fd = stdout)
{
    uint64_t lo = (uint64_t)v, hi = (uint64_t)(v >> 64);
    fprintf(fd, "%s 0x%16.16llx%16.16llx\n", s, (unsigned long long)hi, (unsigned long long)lo);
}

static inline void scan128(const char *s, uint128_t *value)
{
    uint128_t v = 0;
    char c;
    const char *t = s;
    unsigned basis = 10;
    if (*t == '0' && *(t + 1) == 'x')
    {
        t += 2;
        basis = 16;
    }
    while (*t)
    {
        v *= basis;
        c = tolower(*t++);
        if (c >= 0 && c <= '9')
        {
            v += c - '0';
        }
        else if (c >= 'a' and c < 'a' + basis - 10)
        {
            v += c - 'a' + 10;
        }
    }
    *value = v;
}

int block_read(int s, unsigned len, uint8_t *p);
int block_write(int s, unsigned len, uint8_t *p);
int tlv_read(int s, uint16_t *cid, uint8_t *type, uint128_t *value);
int tlv_write(int s, uint16_t cid, uint8_t type, uint128_t value);

#endif
