
class Lcg
{
  private:
    uint128_t seed;

    const uint64_t A;
    const uint64_t C;
    const uint64_t M;

    uint64_t modpow(uint64_t b, uint64_t e)
    {
        uint64_t result = e & 1 ? b : 1;
        uint128_t t;
        while (e)
        {
            e >>= 1;
            t = b;
            t *= b;
            t %= M;
            b = (uint64_t)t;
            if (e & 1)
            {
                t = b;
                t *= result;
                t %= M;
                result = (uint64_t)t;
            }
        }
        return result;
    }

    uint128_t next(void)
    {
        uint128_t t = seed;
        t *= A;
        t += C;
        seed = t % M;
        return seed;
    }

  public:
    Lcg(void)
        : A(1), C(1), M((1ul << 60)-1)
        // : A(137), C(13), M((1ul << 60) - 1)
    {
        seed = 0x8765432187654321ull;
        seed %= M;
    }

    void set_seed(uint128_t s)
    {
        seed = s;
    }

    uint128_t get_seed(uint64_t n = 1)
    {
        uint128_t xn;
        if (n == 0)
        {
            xn = seed;
        }
        else if (n == 1)
        {
            if (A == 1)
            {
                // sequential numbering
                uint128_t t = C;
                t += seed;
                seed = t % M;
                xn = seed;
            }
            else
            {
                // simple case
                xn = next();
            }
        }
        else
        {
            /// corner case where A-1 has no inverse
            if (A == 1)
            {
                uint128_t t = n;
                t *= C;
                t += seed;
                seed = t % M;
                xn = seed;
            }
            else
            {
                // general case when A-1 and M are coprime
                uint128_t x0 = seed;
                uint128_t x1 = next();
                uint128_t t = modpow(A, n) - 1;
                t *= modpow(A - 1, M - 2);
                t %= M;
                t *= (x1 + M) - x0;
                t += x0;
                seed = t % M;
                xn = seed;
            }
        }
        return xn;
    }

    bool sequential(void)
    {
        return (A == 1 && C == 1);
    }
};

static inline uint128_t convert_seed_to_number(uint128_t v)
{
    v <<= 1;
    v |= 1;
    return v;
}

static inline uint128_t convert_number_to_seed(uint128_t v)
{
    return v >> 1;
}
