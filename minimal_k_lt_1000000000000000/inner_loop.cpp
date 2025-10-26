
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <x86intrin.h>

#include "inner_loop.h"
#include "lcg.h"
#include "tlv.h"

template <class T, class TT> static T square_mod(const T &u, const T &q)
{
    TT r = u;
    r *= u;
    r %= q;
    return (T)r;
}

template <class T, class TT> static T mul_mod(const T &u, const T &v, const T &q)
{
    TT r = u;
    r *= v;
    r %= q;
    return (T)r;
}

template <class T, class TT> static T long_mod(const T &u, const T &v, const T &q)
{
    TT r = u;
    r <<= sizeof(T) * 8;
    r += v;
    r %= q;
    return (T)r;
}

template <class T, class TT> static T shift_mod(const T &u, const T &v, const T &q)
{
    TT r = u;
    r <<= v;
    r %= q;
    return (T)r;
}

template <class T, class TT> static T longlong_mod(const TT &u, const T &q)
{
    TT r = u;
    r %= q;
    return (T)r;
}

// (u << 64 + v) mod n
template <> uint64_t long_mod<uint64_t, uint128_t>(const uint64_t &u, const uint64_t &v, const uint64_t &n)
{
#ifdef __x86_64__
    uint64_t r, a;
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(u), "1"(v), "r"(n));
    return r;
#else
    uint128_t t = ((uint128_t)u << 64) + v;
    return t % n;
#endif
}

// (u) mod n
template <> uint64_t longlong_mod(const uint128_t &u, const uint64_t &n)
{
#ifdef __x86_64__
    uint64_t r = (uint64_t)(u >> 64), a = (uint64_t)u;
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#else
    return u % n;
#endif
}

// (u * v) mod n
template <> uint64_t mul_mod<uint64_t, uint128_t>(const uint64_t &u, const uint64_t &v, const uint64_t &n)
{
#ifdef __x86_64__
    uint64_t r, a;
    asm("mulq %3" : "=d"(r), "=a"(a) : "1"(u), "r"(v));
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#else
    uint128_t t = (uint128_t)u * v;
    return t % n;
#endif
}

// (u << v) mod n
template <> uint64_t shift_mod<uint64_t, uint128_t>(const uint64_t &u, const uint64_t &v, const uint64_t &n)
{
#ifdef __x86_64__
    uint64_t r = 0, a = u;
    asm("shldq %b3, %2, %0" : "=d"(r) : "0"(r), "a"(a), "c"(v));
    asm("shlxq %2, %1, %0" : "=a"(a) : "0"(a), "c"(v));
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#else
    uint128_t t = (uint128_t)u << v;
    return t % n;
#endif
}

// (u * u) mod n
template <> uint64_t square_mod<uint64_t, uint128_t>(const uint64_t &u, const uint64_t &n)
{
#ifdef __x86_64__
    uint64_t r, a;
    asm("mulq %3" : "=d"(r), "=a"(a) : "1"(u), "r"(u));
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#else
    uint128_t t = (uint128_t)u * u;
    return t % n;
#endif
}

// count trailing zeroed bits
static inline uint64_t tzcnt(uint64_t a)
{
#ifdef __x86_64__
    uint64_t r;
    asm("tzcntq %1,%0" : "=r"(r) : "r"(a));
    return r;
#else
    return __builtin_ctzll(a);
#endif
}

// count leading zeroed bits
static inline uint64_t lzcnt(uint64_t a)
{
#ifdef __x86_64__
    uint64_t r;
    asm("lzcntq %1,%0" : "=r"(r) : "r"(a));
    return r;
#else
    return __builtin_clzll(a);
#endif
}

template <class T> unsigned log_2(const T &x)
{
    if (x)
    {
        return 63 - lzcnt(x);
    }
    else
    {
        return 0;
    }
}

template <> unsigned log_2<uint128_t>(const uint128_t &x)
{
    uint64_t lo = (uint64_t)x;
    uint64_t hi = (uint64_t)(x >> 64);
    if (hi)
    {
        return 127 - lzcnt(hi);
    }
    return 63 - lzcnt(lo);
}

template <> unsigned log_2<uint64_t>(const uint64_t &x)
{
    return 63 - lzcnt(x);
}

template <class T> static T sub_mod(const T &u, const T &v, const T &q)
{
    T s = u - v;
    if (s <= u)
    {
        s %= q;
    }
    else
    {
        s = -s;
        s %= q;
        s = q - s;
    }
    return s;
}

template <class T, class TT> static T add_mod(const T &u, const T &v, const T &q)
{
    TT r = u;
    r += v;
    r %= q;
    return (T)r;
}

template <class T, class TT> static T div2_mod(const T &u, const T &q)
{
    if (u & 1)
    {
        TT tu = u;
        tu += q;
        tu >>= 1;
        return (T)tu;
    }
    else
    {
        return u >> 1;
    }
}

template <class T, class TT> static T mod(const TT &u, const T &q)
{
    TT r = u;
    r %= q;
    return (T)r;
}

template <class T> static T mod(const T &u, const T &q)
{
    return u % q;
}

template <class T, class TT> static TT square(const T &u)
{
    TT r = u;
    r *= u;
    return (T)r;
}

template <class T, class TT> static TT mul(const T &u, const T &v)
{
    TT r = u;
    r *= v;
    return (T)r;
}

template <class T, class TT> static TT mul2(const T &u, const T &v)
{
    TT r = u;
    r *= v;
    r <<= 1;
    return (T)r;
}

// gcd (x, y) based on Stein's algorithm
template <class T> static T gcd(const T &x, const T &y)
{
    if (x == 0)
        return y;
    if (y == 0)
        return x;
    unsigned tu = tzcnt(x);
    unsigned tv = tzcnt(y);
    unsigned h = tu > tv ? tv : tu;
    T u = x >> tu;
    T v = y >> tv;
    while (1)
    {
        if (u > v)
        {
            T t = u;
            u = v;
            v = t;
        }
        v -= u;
        if (v == 0)
        {
            return u << h;
        }
        v >>= tzcnt(v);
    }
}

// Jacobi symbol (x, y) based on Stein's algorithm
template <class T> static int jacobi(const T &x, const T &y)
{
    // assert((x & 1) == 0);
    // assert((y & 1) == 1);
    if (__builtin_constant_p(x) && x == 2)
    {
        return ((y + 2) & 4) ? -1 : 1;
    }
    if (y == 1 || x == 1)
    {
        return 1;
    }

    if (x == 2)
    {
        // char j[4] = { -1,-1,1,1};
        // return j[((y - 3) >> 1) % 4];
        return ((y + 2) & 4) ? -1 : 1;
    }
    if (x == 3)
    {
        char j[6] = {0, -1, -1, 0, 1, 1};
        return j[((y - 3) >> 1) % 6];
    }
    if (x == 5)
    {
        char j[5] = {-1, 0, -1, 1, 1};
        return j[((y - 3) >> 1) % 5];
    }
    if (x == 7)
    {
        char j[14] = {1, -1, 0, 1, -1, -1, -1, -1, 1, 0, -1, 1, 1, 1};
        return j[((y - 3) >> 1) % 14];
    }
    if (x == 11)
    {
        char j[22] = {-1, 1, 1, 1, 0, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 0, 1, 1, 1, -1, 1, 1};
        return j[((y - 3) >> 1) % 22];
    }
    if (x == 13)
    {
        char j[13] = {1, -1, -1, 1, -1, 0, -1, 1, -1, -1, 1, 1, 1};
        return j[((y - 3) >> 1) % 13];
    }
    if (x == 17)
    {
        char j[17] = {-1, -1, -1, 1, -1, 1, 1, 0, 1, 1, -1, 1, -1, -1, -1, 1, 1};
        return j[((y - 3) >> 1) % 17];
    }
    if (x == 19)
    {
        char j[19] = {1, 1, -1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, -1, 1, -1, -1, -1, -1};
        unsigned t = ((y - 3) >> 1) % 38;
        return t >= 19 ? -j[t - 19] : j[t];
    }
    if (x == 23)
    {
        char j[23] = {-1, -1, 1, 1, 1, 1, 1, -1, 1, -1, 0, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1};
        unsigned t = ((y - 3) >> 1) % 46;
        return t >= 23 ? -j[t - 23] : j[t];
    }
    if (x == 29)
    {
        char j[29] = {-1, 1, 1,  1,  -1, 1,  -1, -1, -1, -1, 1, 1,  -1, 0, -1,
                      1,  1, -1, -1, -1, -1, 1,  -1, 1,  1,  1, -1, 1,  1};
        return j[((y - 3) >> 1) % 29];
    }
    if (x == 31)
    {
        char j[31] = {1,  1,  -1, 1, 1, -1, 1,  -1, -1, -1, 1, 1,  1,  -1, 0, 1,
                      -1, -1, -1, 1, 1, 1,  -1, 1,  -1, -1, 1, -1, -1, -1, -1};
        unsigned t = ((y - 3) >> 1) % 62;
        return t >= 31 ? -j[t - 31] : j[t];
    }

    int t = 1;
    T a = x;
    T n = y;
    unsigned v = n & 7;
    unsigned c = (v == 3) || (v == 5);
    while (a)
    {
        v = tzcnt(a);
        a >>= v;
        t = (c & (v & 1)) ? -t : t;

        if (a < n)
        {
            T tmp = a;
            a = n;
            n = tmp;
            t = ((a & n & 3) == 3) ? -t : t;
            v = n & 7;
            c = (v == 3) || (v == 5);
        }

        a -= n;
    }

    return (n == 1) ? t : 0;
}

// Kronecker symbol (x, y) based on Stein's algorithm
template <class T> static int kronecker(const T &x, const T &y)
{
    unsigned x1 = x & 1;
    unsigned y1 = y & 1;
    if (y1 == 1)
    {
        // bit wizardry for K(+/- 2,y), K(+/- 3,y) when y odd

        if (x == 0)
        {
            return (y == 1) ? 1 : 0;
        }
        if (x == 1)
        {
            return 1;
        }
        if (x == -1)
        {
            return (y & 2) ? -1 : 1;
        }
        if (x == 2)
        {
            return ((y + 2) & 4) ? -1 : 1;
        }
        if (x == -2)
        {
            return (y & 4) ? -1 : 1;
        }
        if (x1 == 1)
        {
            if (x >= 0)
            {
                return jacobi<T>(x, y);
            }
            else
            {
                int j = jacobi<T>(-x, y);
                return (y & 2) ? -j : j;
            }
        }
    }
    else
    {
        // bit wizardry when y even
        if ((x1 == 0 && y1 == 0) || (y == 0 && !(x == 1 || x == -1)))
        {
            return 0;
        }
        if (y == 0)
        {
            return 1;
        }
    }

    // generic case
    T a = x;
    T n = y;
    T a0 = a;
    unsigned v, cn, ca;
    int t = 1;
    if (a < 0)
    {
        if (n < 0)
        {
            n = -n;
            a = -a;
            t = -t;
        }
        else
        {
            a = -a;
        }
    }
    else if (n < 0)
    {
        n = -n;
    }
    // gulp trailing zeroes from n
    v = a & 7;
    ca = (v == 3) || (v == 5);
    v = tzcnt(n);
    n >>= v;
    t = (ca & (v & 1)) ? -t : t;

    if (a0 < 0 && (n & 3) == 3)
    {
        t = -t;
    }
    v = n & 7;
    cn = (v == 3) || (v == 5);

    // gulp trailing zeroes from a within Stein's algorithm.
    while (a > 0)
    {
        v = tzcnt(a);
        a >>= v;
        t = (cn & (v & 1)) ? -t : t;
        if (a == 1)
        {
            n = 1;
            break;
        }
        if (a < n)
        {
            T r = a;
            a = n;
            n = r;
            t = ((a & n & 3) == 3) ? -t : t;
            v = n & 7;
            cn = (v == 3) || (v == 5);
        }
        a = a - n;
    }
    return (n == 1) ? t : 0;
}

template <class T, class TT> static T mod_inv(const T &x, const T &m)
{
    if (m < 3)
        return 0;
    if (x < 2)
        return x;
    T a = x, b = m, u = 1, v = 0;
    while (a != 0)
    {
        unsigned ta = tzcnt(a);
        a >>= ta;
        while (ta--)
        {
            u = div2_mod<T, TT>(u, m);
        }
        if (a < b)
        {
            T t = a;
            T s = u;
            a = b;
            u = v;
            b = t;
            v = s;
        }
        a -= b;
        u = sub_mod<T>(u, v, m);
    }
    return b == 1 ? v : 0;
}

template <class T, class TT> static bool is_perfect_square(const T &a)
{
    if (0xffedfdfefdecull & (1ull << (a % 48)))
        return false;
    if (0xfdfdfdedfdfcfdecull & (1ull << (a % 64)))
        return false;
    if (0x7bfdb7cfedbafd6cull & (1ull << (a % 63)))
        return false;
    if (0x7dcfeb79ee35ccull & (1ull << (a % 55)))
        return false;
    if (0x8ec196bf5a60dc4ull & (1ull << (a % 61)))
        return false;
    if (0x5d49de7c1846d44ull & (1ull << (a % 59)))
        return false;
    if (0xd228fccfc512cull & (1ull << (a % 53)))
        return false;
    if (0x7bcae4d8ac20ull & (1ull << (a % 47)))
        return false;
    if (0x4a77c5c11acull & (1ull << (a % 43)))
        return false;
    if (0x4c7d4af8c8ull & (1ull << (a % 41)))
        return false;
    if (0x9a1dee164ull & (1ull << (a % 37)))
        return false;
    if (0x6de2b848ull & (1ull << (a % 31)))
        return false;
    if (0xc2edd0cull & (1ull << (a % 29)))
        return false;
    if (0x7acca0ull & (1ull << (a % 23)))
        return false;
    if (0x4f50cull & (1ull << (a % 19)))
        return false;
    if (0x5ce8ull & (1ull << (a % 17)))
        return false;
    if (0x9e4ull & (1ull << (a % 13)))
        return false;

    // approximation of square root with floating point accuracy
    double d = (double)a;
    d = exp(log(d) / 2.0); // square root
    double dl = d * 0.999999;
    double dh = d * 1.000001;
    T c, m;
    // binary search (1 more bit of square root per iteration)
    T r = (T)d;
    T l = (T)dl;
    T h = (T)dh;
    while (l <= h)
    {
        m = (l + h) >> 1;
        c = m * m;
        if (c == a)
        {
            return true; // perfect square
        }
        if (c < a)
        {
            l = m + 1;
            r = m;
        }
        else
        {
            h = m - 1;
        }
    }
    c = r * r; // check perfect square
    return (c == a);
}

template <class T, class TT> static bool is_perfect_cube(const T &a)
{
    if (0x3f7fffe7e7fffefcull & (1ull << (a % 63)))
        return false;
    if (0x1fafd7e3f5fafcull & (1ull << (a % 54)))
        return false;
    if (0xbcbfd99e66ff4f4ull & (1ull << (a % 61)))
        return false;
    if (0x176f79ef6e8ull & (1ull << (a % 43)))
        return false;
    if (0xf537fb2bcull & (1ull << (a % 37)))
        return false;
    if (0x177e7ee8 & (1 << (a % 31)))
        return false;
    if (0x3e67c & (1 << (a % 19)))
        return false;
    if (0xedc & (1 << (a % 13)))
        return false;

    // approximation of cubic root with floating point accuracy
    double d = (double)a;
    d = exp(log(d) / 3.0); // cubic root
    double dl = d * 0.999999;
    double dh = d * 1.000001;
    T c, m;
    // binary search (1 more bit of cube root per iteration)
    T r = (T)d;
    T l = (T)dl;
    T h = (T)dh;
    while (l <= h)
    {
        m = (l + h) >> 1;
        c = m * m * m;
        if (c == a)
        {
            return true; // perfect cube
        }
        if (c < a)
        {
            l = m + 1;
            r = m;
        }
        else
        {
            h = m - 1;
        }
    }
    c = r * r * r; // check perfect cube
    return (c == a);
}

template <class T, class TT> static bool is_perfect_sursolid(const T &a)
{
    if (0x1f7fef8fbff7cull & (1ull << (a % 50)))
        return false;
    if (0x7fcff5fe7fcull & (1ull << (a % 44)))
        return false;
    if (0xffa7efedfdf97fcull & (1ull << (a % 61)))
        return false;
    if (0xbef7ffbdf4ull & (1ull << (a % 41)))
        return false;
    if (0x39ffff9cull & (1ull << (a % 31)))
        return false;
    if (0x1249248ull & (1ull << (a % 27)))
        return false;
    if (0x40810204080ull & (1ull << (a % 49)))
        return false;

    // approximation of fifth root with floating point accuracy
    double d = (double)a;
    d = exp(log(d) / 5.0); // fifth root
    double dl = d * 0.999999;
    double dh = d * 1.000001;
    T c, m;
    // binary search (1 more bit of fifth root per iteration)
    T r = (T)d;
    T l = (T)dl;
    T h = (T)dh;
    while (l <= h)
    {
        m = (l + h) >> 1;
        c = m * m * m * m * m;
        if (c == a)
        {
            return true; // perfect sursolid
        }
        if (c < a)
        {
            l = m + 1;
            r = m;
        }
        else
        {
            h = m - 1;
        }
    }
    c = r * r * r * r * r; // check perfect sursolid
    return (c == a);
}

// 2^e mod m
template <class T, class TT> static T pow2_mod(const T &d, const T &m)
{
    // 2^e mod m : hardcode 6 first iterations
    T n = log_2(d);
    n = (n > 5) ? n - 5 : 0;
    T e = (d >> n) & 0x3f;
    T result = shift_mod<T, TT>(1ull, e, m);

    while (n >= 6)
    {
        n -= 6;
        result = square_mod<T, TT>(result, m);
        result = square_mod<T, TT>(result, m);
        result = square_mod<T, TT>(result, m);
        result = square_mod<T, TT>(result, m);
        result = square_mod<T, TT>(result, m);
        result = square_mod<T, TT>(result, m);
        e = (d >> n) & 0x3f;
        if (e)
        {
            result = shift_mod<T, TT>(result, e, m);
        }
    }
    while (n--)
    {
        result = square_mod<T, TT>(result, m);
        e = (d >> n) & 1;
        if (e)
        {
            result <<= 1;
            result -= result >= m ? m : 0;
        }
    }
    return result;
}

// a^e mod m
template <class T, class TT> static T pow_mod(const T &a, const T &d, const T &m)
{
    if (a == 2)
    {
        return pow2_mod<T, TT>(d, m);
    }

    T n = d;
    T s = a;
    T result = 1;
    while (n)
    {
        if (n & 1)
            result = mul_mod<T, TT>(result, s, m);
        s = square_mod<T, TT>(s, m);
        n >>= 1;
    }
    return result;
}

// MR strong test
template <class T, class TT> static bool witness(const T &n, int s, const T &d, const T &a)
{
    T x, y;
    if (n == a)
        return true;
    x = pow_mod<T, TT>(a, d, n);
    while (s)
    {
        y = square_mod<T, TT>(x, n);
        if (y == 1 && x != 1 && x != n - 1)
            return false;
        x = y;
        --s;
    }
    if (y != 1)
        return false;
    return true;
}

// sieve small factors <= 151
template <class T> static bool sieve(const T &n)
{
    if (n <= 152)
    {
        // return false for small numbers which are composite for sure , without checking further.
        bool stooopid_prime_table[] = {
            true,  true,  true,  true,  false, true,  false, true,  false, false, false, true,  false, true,  false,
            false, false, true,  false, true,  false, false, false, true,  false, false, false, false, false, true,
            false, true,  false, false, false, false, false, true,  false, false, false, true,  false, true,  false,
            false, false, true,  false, false, false, false, false, true,  false, false, false, false, false, true,
            false, true,  false, false, false, false, false, true,  false, false, false, true,  false, true,  false,
            false, false, false, false, true,  false, false, false, true,  false, false, false, false, false, true,
            false, false, false, false, false, false, false, true,  false, false, false, true,  false, true,  false,
            false, false, true,  false, true,  false, false, false, true,  false, false, false, false, false, false,
            false, false, false, false, false, false, false, true,  false, false, false, true,  false, false, false,
            false, false, true,  false, true,  false, false, false, false, false, false, false, false, false, true,
            false, true,  false, false, false, false, false, true,  false, false, false, false, false, true,  false,
            false, false, true,  false, false, false, false, false, true,  false, false, false, false, false, true,
            false, true,  false, false, false, false, false, false, false, false, false, true,  false, true,  false,
            false, false, true,  false, true,  false};
        return stooopid_prime_table[n];
    }
    if ((uint64_t)(n * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
        return false; // divisible by 3
    if ((uint64_t)(n * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
        return false; // divisible by 5
    if ((uint64_t)(n * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
        return false; // divisible by 7
    if ((uint64_t)(n * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
        return false; // divisible by 11
    if ((uint64_t)(n * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
        return false; // divisible by 13
    if ((uint64_t)(n * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
        return false; // divisible by 17
    if ((uint64_t)(n * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
        return false; // divisible by 19
    if ((uint64_t)(n * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
        return false; // divisible by 23
    if ((uint64_t)(n * 0x34f72c234f72c235ull) <= 0x08d3dcb08d3dcb08ull)
        return false; // divisible by 29
    if ((uint64_t)(n * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
        return false; // divisible by 31
    if (n < 37 * 37)
        return true; // prime
    if ((uint64_t)(n * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
        return false; // divisible by 37
    if ((uint64_t)(n * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
        return false; // divisible by 41
    if ((uint64_t)(n * 0x82fa0be82fa0be83ull) <= 0x05f417d05f417d05ull)
        return false; // divisible by 43
    if ((uint64_t)(n * 0x51b3bea3677d46cfull) <= 0x0572620ae4c415c9ull)
        return false; // divisible by 47
    if ((uint64_t)(n * 0x21cfb2b78c13521dull) <= 0x04d4873ecade304dull)
        return false; // divisible by 53
    if ((uint64_t)(n * 0xcbeea4e1a08ad8f3ull) <= 0x0456c797dd49c341ull)
        return false; // divisible by 59
    if ((uint64_t)(n * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
        return false; // divisible by 61
    if ((uint64_t)(n * 0xf0b7672a07a44c6bull) <= 0x03d226357e16ece5ull)
        return false; // divisible by 67
    if ((uint64_t)(n * 0x193d4bb7e327a977ull) <= 0x039b0ad12073615aull)
        return false; // divisible by 71
    if ((uint64_t)(n * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
        return false; // divisible by 73
    if ((uint64_t)(n * 0x9b8b577e613716afull) <= 0x033d91d2a2067b23ull)
        return false; // divisible by 79
    if ((uint64_t)(n * 0xa3784a062b2e43dbull) <= 0x03159721ed7e7534ull)
        return false; // divisible by 83
    if ((uint64_t)(n * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
        return false; // divisible by 89
    if ((uint64_t)(n * 0xa3a0fd5c5f02a3a1ull) <= 0x02a3a0fd5c5f02a3ull)
        return false; // divisible by 97
    if (n < 101 * 101)
        return true; // prime
    if ((uint64_t)(n * 0x3a4c0a237c32b16dull) <= 0x0288df0cac5b3f5dull)
        return false; // divisible by 101
    if ((uint64_t)(n * 0xdab7ec1dd3431b57ull) <= 0x027c45979c95204full)
        return false; // divisible by 103
    if ((uint64_t)(n * 0x77a04c8f8d28ac43ull) <= 0x02647c69456217ecull)
        return false; // divisible by 107
    if ((uint64_t)(n * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
        return false; // divisible by 109
    if ((uint64_t)(n * 0x90fdbc090fdbc091ull) <= 0x0243f6f0243f6f02ull)
        return false; // divisible by 113
    if ((uint64_t)(n * 0x7efdfbf7efdfbf7full) <= 0x0204081020408102ull)
        return false; // divisible by 127
    if ((uint64_t)(n * 0x03e88cb3c9484e2bull) <= 0x01f44659e4a42715ull)
        return false; // divisible by 131
    if ((uint64_t)(n * 0xe21a291c077975b9ull) <= 0x01de5d6e3f8868a4ull)
        return false; // divisible by 137
    if ((uint64_t)(n * 0x3aef6ca970586723ull) <= 0x01d77b654b82c339ull)
        return false; // divisible by 139
    if ((uint64_t)(n * 0xdf5b0f768ce2cabdull) <= 0x01b7d6c3dda338b2ull)
        return false; // divisible by 149
    if ((uint64_t)(n * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
        return false; // divisible by 151
    return true;      // might be prime
}

template <class T, class TT> static bool isprime(const T &n)
{
    if (!sieve<T>(n))
        return false; // composite
    if (n < 157 * 157)
        return true; // prime
    T d = n / 2;
    int s = 1;
    while (!(d & 1))
    {
        d /= 2;
        ++s;
    }

    if (n < 1373653)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3);
    if (n < 9080191)
        return witness<T, TT>(n, s, d, 31) && witness<T, TT>(n, s, d, 73);
    if (n < 4759123141)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 61);
    if (n < 1122004669633)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 13) && witness<T, TT>(n, s, d, 23) &&
               witness<T, TT>(n, s, d, 1662803);
    if (n < 2152302898747)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3) && witness<T, TT>(n, s, d, 5) &&
               witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 11);
    if (n < 3474749660383)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3) && witness<T, TT>(n, s, d, 5) &&
               witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 11) && witness<T, TT>(n, s, d, 13);
    if (n < 341550071728321)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3) && witness<T, TT>(n, s, d, 5) &&
               witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 11) && witness<T, TT>(n, s, d, 13) &&
               witness<T, TT>(n, s, d, 17);
    if (n < 3825123056546413051)
        return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3) && witness<T, TT>(n, s, d, 5) &&
               witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 11) && witness<T, TT>(n, s, d, 13) &&
               witness<T, TT>(n, s, d, 17) && witness<T, TT>(n, s, d, 19) && witness<T, TT>(n, s, d, 23);
    // n < 318665857834031151167461
    return witness<T, TT>(n, s, d, 2) && witness<T, TT>(n, s, d, 3) && witness<T, TT>(n, s, d, 5) &&
           witness<T, TT>(n, s, d, 7) && witness<T, TT>(n, s, d, 11) && witness<T, TT>(n, s, d, 13) &&
           witness<T, TT>(n, s, d, 17) && witness<T, TT>(n, s, d, 19) && witness<T, TT>(n, s, d, 23) &&
           witness<T, TT>(n, s, d, 29) && witness<T, TT>(n, s, d, 31) && witness<T, TT>(n, s, d, 37);
}

//  Mod(Mod(x+t,n),x^2-(sgn*a))^e
//
//  if T is uint64_t, assume t,a,n,e are 61 bit numbers   (require 3 guard bits)
template <class T, class TT> inline void exponentiate2(T &s, T &t, const T e, const T n, const int sgn, const T a)
{
    T t0 = t;
    T t2, s2, ss, tt;
    TT tmp, ss2, tt2;
    unsigned bit = log_2<T>(e);
    while (bit--)
    {
        if (__builtin_constant_p(sgn) && sgn < 0 && __builtin_constant_p(a) && a == 1)
        {
            tt = mul_mod<T, TT>(t + s, t + n - s, n); // f bits
            ss = mul_mod<T, TT>(s, t, n);             // f bits
            ss += ss;                                 // f+1 bits
        }
        else
        {

            t2 = square_mod<T, TT>(t, n); // f bits
            s2 = square_mod<T, TT>(s, n); // f bits
            ss = mul_mod<T, TT>(s, t, n); // f bits
            ss += ss;                     // f+1 bits
            if (__builtin_constant_p(sgn) && sgn == 1)
            {
                if (__builtin_constant_p(a) && a == 1)
                {
                    tt = s2 + t2; // f+1 bits
                }
                else if (__builtin_constant_p(a) && a == 2)
                {
                    tt = s2 + s2 + t2; // f+2 bits
                }
                else
                {
                    tt = mul_mod<T, TT>(s2, a, n); // f bits
                    tt += t2;                      // f+1 bits
                }
            }
            else if (__builtin_constant_p(sgn) && sgn == -1)
            {
                if (__builtin_constant_p(a) && a == 1)
                {
                    tt = n - s2 + t2; // f+1 bits
                }
                else if (__builtin_constant_p(a) && a == 2)
                {
                    tt = n + n - (s2 + s2) + t2; // f+2 bits
                }
                else
                {
                    tt = mul_mod<T, TT>(s2, n - a, n); // f bits
                    tt += t2;                          // f+1 bits
                }
            }
            else if (sgn == 1)
            {
                tt = mul_mod<T, TT>(s2, a, n); // f bits
                tt += t2;                      // f+1 bits
            }
            else if (sgn == -1)
            {
                tt = mul_mod<T, TT>(s2, n - a, n); // f bits
                tt += t2;                          // f+1 bits
            }
            else
            {
                assert(0);
            }
        }

        if (e & ((T)1 << bit))
        {
            if (__builtin_constant_p(sgn) && sgn == 1)
            {
                if (__builtin_constant_p(a) && a == 1)
                {
                    tmp = ss; // 2f + 1 bits
                }
                else if (__builtin_constant_p(a) && a == 2)
                {
                    tmp = ss + ss; // 2f + 1 bits
                }
                else
                {
                    tmp = (TT)ss * a; // 2f + 1 bits
                }
            }
            else if (__builtin_constant_p(sgn) && sgn == -1)
            {
                if (__builtin_constant_p(a) && a == 1)
                {
                    tmp = n + n - ss; // 2f + 1 bits
                }
                else if (__builtin_constant_p(a) && a == 2)
                {
                    tmp = n + n - ss; // 2f + 1 bits
                    tmp <<= 1;
                }
                else
                {
                    tmp = (TT)ss * (n - a); // 2f + 1 bits
                }
            }
            else if (sgn == 1)
            {
                tmp = (TT)ss * a; // 2f + 1 bits
            }
            else if (sgn == -1)
            {
                tmp = (TT)ss * (n - a); // 2f + 1 bits
            }
            else
            {
                assert(0);
            }

            if (__builtin_constant_p(t0) && t0 == 0)
            {
                ss2 = tt;
                tt2 = tmp;
            }
            else if (__builtin_constant_p(t0) && t0 == 1)
            {
                ss2 = ss + tt;
                tt2 = tt + tmp;
            }
            else if (__builtin_constant_p(t0) && t0 == 2)
            {
                ss2 = ss + ss + tt;
                tt2 = tt + tt + tmp;
            }
            else
            {
                ss2 = (TT)ss * t0 + tt;  // 3f + 2 bits
                tt2 = (TT)tt * t0 + tmp; // 4f + 2 bits
            }
            s = mod<T, TT>(ss2, n);
            t = mod<T, TT>(tt2, n);
        }
        else
        {
            s = mod<T>(ss, n);
            t = mod<T>(tt, n);
        }
    }
}

//  Mod(Mod(x,n),x^2-2*x-1)^e
//
//  if T is uint64_t, assume t,a,n,e are 63 bit numbers   (require 1 guard bit)
template <class T, class TT> void exponentiate_half(T &s, T &t, const T e, const T n)
{
    unsigned bit = log_2<T>(e);

    while (bit--)
    {
        // (s*x+t)^2 = 2*(s+t)*s*x + t^2+s^2
        T s2 = square_mod<T, TT>(s, n);
        T t2 = square_mod<T, TT>(t, n);
        s = mul_mod<T, TT>(s + t, 2 * s, n);
        t = add_mod<T, TT>(s2, t2, n);

        if (e & ((T)1 << bit))
        {
            // (s*x+t)*x = (2*s+t)*x + s
            T tmp = s;
            s = add_mod<T, TT>(s + s, t, n);
            t = tmp;
        }
    }
}

/*
If n==3 mod 4 test Mod(Mod(x+2,n),x^2+1)^(n+1)==5.
If n==5 mod 8 test Mod(Mod(x+2,n),x^2+2)^(n+1)==6.
If n==1 mod 8 test Mod(Mod(x+2,n),x^2-a)^(n+1)==4-a and Mod(Mod(x+2,n),x^2+a)^(n+1)==4+a for kronecker(a,n)==-1
*/

template <class T, class TT> static bool islnrc2prime(const T &n, int s = 0, int cid = 0)
{
    if (n < 23)
        return (n == 1 || n == 2 || n == 3 || n == 5 || n == 7 || n == 11 || n == 13 || n == 17 ||
                n == 19); // prime for sure

    T mod8 = n & 7;
    if (mod8 == 3 || mod8 == 7)
    {
        // (x+2)^(n+1) mod (n, x^2+1) == 5
        T bs = 1;
        T bt = 2;
        exponentiate2<T, TT>(bs, bt, n + 1, n, -1, 1);
        // printf("3 mod 4 : %lx %lx %lx\n", bs, bt, n);
        return (bs == 0 && bt == 5); // ?? n prime ? n composite for sure ?
    }
    if (mod8 == 5)
    {
        // (x+2)^(n+1) mod (n, x^2+2) == 6
        T bs = 1;
        T bt = 2;
        exponentiate2<T, TT>(bs, bt, n + 1, n, -1, 2);
        // printf("5 mod 8 : %lx %lx %lx\n", bs, bt, n);
        return (bs == 0 && bt == 6); // ?? n prime ? n composite for sure ?
    }

    if (is_perfect_square<T, TT>(n))
        return false; // n composite perfect square, for any x, kronecker(x, n)==1 always

    // search minimal a where Kronecker(a, n) == -1
    T a = 3;
    T da = 2;
    int j = jacobi<T>(a, n);
    if (j == 0)
        return false; // composite for sure
    if (j == 1)
    {
        T da = 2;
        for (a = 5;; a += da, da = 6 - da)
        {
            if (!isprime<T, TT>(a))
                continue;

            j = jacobi<T>(a, n);
            if (j == 0)
                return false; // composite for sure
            if (j == -1)
                break;
        }
    }
    // (x+2)^(n+1) mod (n, x^2+a) == 4+a
    T bs = 1;
    T bt = 2;
    exponentiate2<T, TT>(bs, bt, n + 1, n, -1, a);
    if (!(bs == 0 && bt == add_mod<T, TT>(4, a, n)))
        return false; // composite for sure

    // (x+2)^(n+1) mod (n, x^2-a) == 4-a
    bs = 1;
    bt = 2;
    exponentiate2<T, TT>(bs, bt, n + 1, n, 1, a);
    if (!(bs == 0 && bt == add_mod<T, TT>(4, n - a, n)))
        return false; // composite for sure
    return true;      // ?? n prime ?
}

int inner_loop(int s, uint16_t cid, uint128_t seed, uint64_t count)
{
    bool r, rl;
    uint128_t v = 0;
    char buff[60];
    fflush(stdout);
    Lcg u;
    u.set_seed(seed);
    v = convert_seed_to_number(seed);
    sprintf(buff, "Client %4u Started ....", (unsigned)cid);
    print128(buff, v);
    fflush(stdout);
    while (count--)
    {
        v = u.get_seed(1);
        v = convert_seed_to_number(v);
        if (v >> 61)
        {
            assert(0);
        }
        else
        {
            r = isprime<uint64_t, uint128_t>(v);
            rl = islnrc2prime<uint64_t, uint128_t>(v, s, cid);
        }
        if (r)
        {
            if (!rl)
            {
                // printf("pseudocomposite %lx\n", (uint64_t)v);
                if (s)
                {
                    tlv_write(s, cid, TLV_PSEUDOCOMPOSITE, v);
                }
            }
            else
            {
                // printf("prime %lx\n", (uint64_t)v);
            }
        }
        else
        {
            if (rl)
            {
                // printf("pseudoprime %lx\n", (uint64_t)v);
                if (s)
                {
                    tlv_write(s, cid, TLV_PSEUDOPRIME, v);
                }
            }
            else
            {
                // printf("composite %lx\n", (uint64_t)v);
            }
        }
    }

    sprintf(buff, "Client %4u Completed ..", (unsigned)cid);
    print128(buff, v);
    fflush(stdout);
    return 0;
}

static int inner_self_test_64(void)
{
    uint64_t s, r, t;
    bool b;
    int j;

    printf("Modular operations ...\n");
    s = 10103;
    t = 10101;
    r = square_mod<uint64_t, uint128_t>(s, t);
    assert(r == 4);
    s = 10103;
    t = 10101;
    r = mul_mod<uint64_t, uint128_t>(s, s, t);
    assert(r == 4);

    s = 1;
    t = 65535;
    r = shift_mod<uint64_t, uint128_t>(s, 16, t);
    assert(r == 1);
    r = shift_mod<uint64_t, uint128_t>(s, 32, t);
    assert(r == 1);
    r = shift_mod<uint64_t, uint128_t>(s, 48, t);
    assert(r == 1);

    s = 3ull << 47;
    t = (1ull << 60) - 1;
    r = shift_mod<uint64_t, uint128_t>(s, 13, t);
    assert(r == 3);
    r = shift_mod<uint64_t, uint128_t>(s, 23, t);
    assert(r == 3072);
    r = shift_mod<uint64_t, uint128_t>(s, 33, t);
    assert(r == 3145728);

    printf("Perfect square ...\n");
    b = is_perfect_square<uint64_t, uint128_t>(6);
    if (b)
        return -1;
    b = is_perfect_square<uint64_t, uint128_t>(64);
    if (!b)
        return -1;
    b = is_perfect_square<uint64_t, uint128_t>(27);
    if (b)
        return -1;
    b = is_perfect_square<uint64_t, uint128_t>(0x1002001);
    if (!b)
        return -1;
    b = is_perfect_square<uint64_t, uint128_t>(0x1002000);
    if (b)
        return -1;
    b = is_perfect_square<uint64_t, uint128_t>(0x1002002);
    if (b)
        return -1;

    printf("Perfect cube ...\n");
    b = is_perfect_cube<uint64_t, uint128_t>(6);
    if (b)
        return -1;
    b = is_perfect_cube<uint64_t, uint128_t>(64);
    if (!b)
        return -1;
    b = is_perfect_cube<uint64_t, uint128_t>(81);
    if (b)
        return -1;
    b = is_perfect_cube<uint64_t, uint128_t>(0x1003003001);
    if (!b)
        return -1;
    b = is_perfect_cube<uint64_t, uint128_t>(0x1003003000);
    if (b)
        return -1;
    b = is_perfect_cube<uint64_t, uint128_t>(0x1003003002);
    if (b)
        return -1;

    printf("Perfect sursolid ...\n");
    b = is_perfect_sursolid<uint64_t, uint128_t>(6);
    if (b)
        return -1;
    b = is_perfect_sursolid<uint64_t, uint128_t>(64 * 16);
    if (!b)
        return -1;
    b = is_perfect_sursolid<uint64_t, uint128_t>(81);
    if (b)
        return -1;
    b = is_perfect_sursolid<uint64_t, uint128_t>(0x100500A00A005001ull);
    if (!b)
        return -1;
    b = is_perfect_sursolid<uint64_t, uint128_t>(0x100500A00A005000ull);
    if (b)
        return -1;
    b = is_perfect_sursolid<uint64_t, uint128_t>(0x100500A00A005002ull);
    if (b)
        return -1;

    printf("Gcd ...\n");
    s = 12;
    t = 15;
    r = gcd<uint64_t>(s, t);
    if (r != 3)
        return -1;
    s = 12;
    t = 30;
    r = gcd<uint64_t>(s, t);
    if (r != 6)
        return -1;

    printf("Jacobi ...\n");
    s = 33;
    t = 9999;
    j = jacobi<uint64_t>(s, t);
    if (j != 0)
        return -1;
    s = 34;
    t = 9999;
    j = jacobi<uint64_t>(s, t);
    if (j != -1)
        return -1;
    s = 35;
    t = 9999;
    j = jacobi<uint64_t>(s, t);
    if (j != 1)
        return -1;

    printf("Kronecker ...\n");
    s = 33;
    t = 9999;
    j = kronecker<int64_t>(s, t);
    if (j != 0)
        return -1;
    s = 34;
    t = 9999;
    j = kronecker<int64_t>(s, t);
    if (j != -1)
        return -1;
    s = 35;
    t = 9999;
    j = kronecker<int64_t>(s, t);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(11, 101);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(-11, 101);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(13, 101);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(-13, 101);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(-1, 101);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(0, 101);
    if (j != 0)
        return -1;
    j = kronecker<int64_t>(1, 101);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(1, 0);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(2, 0);
    if (j != 0)
        return -1;
    j = kronecker<int64_t>(13, -101);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(-13, -101);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(-2, -11);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(-2, -9);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(-2, -7);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(-2, -5);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(-2, -3);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(-2, -1);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(-2, 1);
    if (j != 1)
        return 1;
    j = kronecker<int64_t>(-2, 3);
    if (j != 1)
        return 1;
    j = kronecker<int64_t>(-2, 5);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(-2, 7);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(-2, 9);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(-2, 11);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(2, 9);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(2, -9);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(2, 11);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(2, -11);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(3, 11);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(-3, 11);
    if (j != -1)
        return -1;
    j = kronecker<int64_t>(3, 13);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(-3, 13);
    if (j != 1)
        return -1;
    j = kronecker<int64_t>(3, 15);
    if (j != 0)
        return -1;
    j = kronecker<int64_t>(-3, 15);
    if (j != 0)
        return -1;

    printf("Modinv ...\n");
    s = 11;
    t = 15;
    r = mod_inv<uint64_t, uint128_t>(s, t);
    if (r != 11)
        return -1;
    s = 12;
    t = 31;
    r = mod_inv<uint64_t, uint128_t>(s, t);
    if (r != 13)
        return -1;
    s = 1234567;
    t = 87654321;
    r = mod_inv<uint64_t, uint128_t>(s, t);
    if (r != 75327931)
        return -1;
    s = 11;
    t = 99;
    r = mod_inv<uint64_t, uint128_t>(s, t);
    if (r != 0)
        return -1;

    printf("Power ...\n");
    r = pow_mod<uint64_t, uint128_t>(2, 0xfedc, 197);
    s = pow2_mod<uint64_t, uint128_t>(0xfedc, 197);
    if (r != 182 || r != s)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x8765, 197);
    s = pow2_mod<uint64_t, uint128_t>(0x8765, 197);
    if (r != 103 || r != s)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x81, 197);
    if (r != 153)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x80, 197);
    s = pow2_mod<uint64_t, uint128_t>(0x80, 197);
    if (r != 175 || r != s)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x41, 197);
    if (r != 122)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x40, 197);
    s = pow2_mod<uint64_t, uint128_t>(0x40, 197);
    if (r != 61 || r != s)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x22, 197);
    if (r != 155)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x21, 197);
    if (r != 176)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x20, 197);
    s = pow2_mod<uint64_t, uint128_t>(0x20, 197);
    if (r != 88 || r != s)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x1f, 197);
    if (r != 44)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x1e, 197);
    if (r != 22)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x1d, 197);
    s = pow2_mod<uint64_t, uint128_t>(0x1d, 197);
    if (r != 11 || r != s)
        return -1;
    r = pow_mod<uint64_t, uint128_t>(2, 0x1c, 197);
    s = pow2_mod<uint64_t, uint128_t>(0x1c, 197);
    if (r != 104 || r != s)
        return -1;

    r = pow_mod<uint64_t, uint128_t>(3, 0xaa55, 197);
    if (r != 0xa7)
        return -1;

    printf("Known primes ...\n");
    b = isprime<uint64_t, uint128_t>(200003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = isprime<uint64_t, uint128_t>(2000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = isprime<uint64_t, uint128_t>(20000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = isprime<uint64_t, uint128_t>(2000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = isprime<uint64_t, uint128_t>(20000000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }
    b = isprime<uint64_t, uint128_t>(200000000000000003ull);
    if (!b)
    {
        printf("expected prime failed\n");
        return (-1);
    }

    printf("Linear recurrence second order ...\n");
    if (1)
    {
        uint64_t s, t;
        s = 1;
        t = 2;
        exponentiate2<uint64_t, uint128_t>(s, t, 2, 101, 1, 13);
        if (s != 4 || t != 17)
        {
            printf("linear recurrence 2 failed _2_\n");
            return (-1);
        }
        s = 1;
        t = 2;
        exponentiate2<uint64_t, uint128_t>(s, t, 4, 101, 1, 13);
        if (s != 35 || t != 93)
        {
            printf("linear recurrence 2 failed _4_\n");
            return (-1);
        }
        s = 1;
        t = 2;
        exponentiate2<uint64_t, uint128_t>(s, t, 3, 101, 1, 13);
        if (s != 25 || t != 86)
        {
            printf("linear recurrence 2 failed _3_\n");
            return (-1);
        }
        s = 1;
        t = 2;
        exponentiate2<uint64_t, uint128_t>(s, t, 5, 101, 1, 13);
        if (s != 62 || t != 35)
        {
            printf("linear recurrence 2 failed _5_\n");
            return (-1);
        }
        s = 1;
        t = 0;
        exponentiate2<uint64_t, uint128_t>(s, t, 12, 101, 1, 13);
        if (s != 0 || t != 19)
        {
            printf("linear recurrence 2 failed _12_\n");
            return (-1);
        }
    }

    printf("Isprime ...\n");
    t = 1;
    t <<= 3;
    t -= 1;
    b = isprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("isprime M(3) failed\n");
        return -1;
    }
    b = islnrc2prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc2 M(3) failed\n");
        return -1;
    }

    t = 101;
    b = isprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("isprime 101 failed\n");
        return -1;
    }
    b = islnrc2prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc2 101 failed\n");
        return -1;
    }

    t = 4493;
    b = isprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("isprime 4493 failed\n");
        return -1;
    }
    b = islnrc2prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc2 4493 failed\n");
        return -1;
    }

    t = 1;
    t <<= 31;
    t -= 1;
    b = isprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("isprime M(31) failed\n");
        return -1;
    }
    b = islnrc2prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc2 M(31) failed\n");
        return -1;
    }

    t = 1;
    t <<= 61;
    t -= 1;
    b = isprime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("isprime M(61) failed\n");
        return -1;
    }
    b = islnrc2prime<uint64_t, uint128_t>(t);
    if (!b)
    {
        printf("islnrc2 M(61) failed\n");
        return -1;
    }

#define COUNT 5000
    volatile uint64_t t0, t1;
    t0 = __rdtsc();
    inner_loop(0, 2222, 0x4000000000, COUNT);
    t1 = __rdtsc();
    double d = (double)(t1 - t0);
    d /= COUNT;
    printf("Average: %8.1f ticks/iteration\n", d);
    printf("Self-test completed\n");
    // pass
    return 0;
}

int inner_self_test(void)
{
    int rc = 0;
    rc = inner_self_test_64();
    if (rc)
    {
        printf("64/128 bit failed\n");
        return rc;
    }
    return 0;
}
