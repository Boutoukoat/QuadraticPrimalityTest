// Version 0.3.0
// Compile: gcc -03 -o Lucas-and-Fermats Lucas-and-Fermats.c -lgmp

#include <gmp.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define FERMAT_THREAD_COUNT 3

FILE *input_file;
int jac;
uint64_t r, shft, sz1, sz2;
pthread_t f2, fa, fc, l1;
time_t next_time = 0;
mpz_t n, a, b, c, Q, g, np1, s, t, tmp, tmp2, tmp3, thou, barrett;
mpz_t fermat_args[FERMAT_THREAD_COUNT];
bool composite = false;

void barrett_precalc()
{
    sz1 = mpz_sizeinbase(n, 2);
    sz2 = sz1 >> 1;
    shft = sz1 + sz2;
    mpz_set_ui(thou, 1);
    mpz_mul_2exp(thou, thou, shft);
    mpz_div(barrett, thou, n);
    mpz_mod(thou, thou, n);
}

void barrett_mod(mpz_t x)
{
    // First step :
    // efficiently strip the most significant quarter of the number x
    // In this file, the input number can be immensely larger than twice the
    // modulus size x = x % (2^(1.5 log2(n))) + (x >> (2^(1.5 log2(n)))) *
    // ((2^(1.5 log2(n))) % n)
    while (mpz_sizeinbase(x, 2) > 2 * sz1 + 2)
    {
        // really large input numbers > 2 * log2(n) + guard bits
        mpz_mod_2exp(tmp, x, shft + sz2);
        mpz_div_2exp(tmp2, x, shft + sz2);
        mpz_mul(tmp2, tmp2, thou);
        mpz_mul_2exp(tmp2, tmp2, sz2);
        mpz_add(x, tmp, tmp2);
    }
    if (mpz_sizeinbase(x, 2) > shft + 2)
    {
        // large input numbers > 1.5 * log2(n) + guard bits
        mpz_div_2exp(tmp2, x, shft);
        mpz_mod_2exp(tmp, x, shft);
        mpz_mul(tmp2, tmp2, thou);
        mpz_add(x, tmp, tmp2);
    }

    // Second step :
    // strip the most significant part of the number x (regular Barrett
    // subtraction) x -= ((x >> log2(n)*(2^(1.5 log2(n))/n)) >> (0.5 log2(n)))*n
    mpz_div_2exp(tmp, x, sz1);
    mpz_mul(tmp, tmp, barrett);
    mpz_div_2exp(tmp, tmp, sz2);
    mpz_mul(tmp, tmp, n);
    mpz_sub(x, x, tmp);

    // Third step :
    // make sure the remainder mod n is < n  (at most 4 subtractions)
    // There is a trade-off between the number of subtractions and the guard bits
    // in the first step.
    // x -= n
    while (mpz_cmp(x, n) >= 0)
        mpz_sub(x, x, n);

    // this happens when the input number is negative
    // (the calculations in this file do not prevent negative numbers)
    // x += n
    while (mpz_cmp_ui(x, 0) < 0)
        mpz_add(x, x, n);

    // now 0 <= x < n
}

void *fermat(void *arg)
{
    mpz_srcptr fermat_arg = (mpz_srcptr)arg;
    mpz_t result;
    mpz_init(result);
    mpz_powm(result, fermat_arg, n, n);
    if (mpz_cmp(result, fermat_arg) != 0)
    {
        gmp_printf("Failed Fermat %Zd-PRP test.\n", fermat_arg);
        composite = true;
    }
    else
    {
        gmp_printf("Passed Fermat %Zd-PRP test.\n", fermat_arg);
    }
    mpz_clear(result);
    return 0;
}

void *lucas(void *)
{ 
    // Mod(Mod(x+b),x^2-a)^(n+1) == b^2-a where b=2^r
    mpz_set_ui(s, 1);
    mpz_set(t, b);
    for (int bit = mpz_sizeinbase(np1, 2) - 2; bit >= 0; bit--)
    {
        if (time(NULL) > next_time)
        {
            printf("%d iterations to go...\r", bit);
            fflush(stdout);
            next_time = time(NULL) + 5;
            if (composite)
            {
                return 0;
            }
        }
        // (s*x+t)^2 = 2*s*t*x + (t+s)*(t+a*s)-s*t*(a+1)
        // (s*x+t)^2 = ((s+t)^2-s^2-t^2)*x + t^2 + s^2 * a
        mpz_add(tmp, s, t);
        mpz_mul(tmp, tmp, tmp);
        mpz_mul(tmp2, s, s);
        mpz_mul(tmp3, t, t);
        mpz_add(s, tmp2, tmp3);
        mpz_sub(s, tmp, s);
        mpz_mul(t, tmp2, a);
        mpz_add(t, t, tmp3);
        if (mpz_tstbit(np1, bit))
        {
            // (s*x+t)*(x+b)=(t+b*s)*x + a*s+b*t
            mpz_mul(tmp, b, s);
            mpz_add(tmp, tmp, t);
            mpz_mul(tmp2, a, s);
            mpz_mul(tmp3, b, t);
            mpz_add(t, tmp2, tmp3);
            mpz_set(s, tmp);
        }
        barrett_mod(s);
        barrett_mod(t);
    }
    mpz_mul(Q, b, b);
    mpz_sub(Q, Q, a);
    mpz_mod(Q, Q, n);
    if (mpz_cmp_ui(s, 0) != 0 || mpz_cmp(t, Q) != 0)
    {
        printf("Failed Lucas PRP test.    \n");
        composite = true;
    }
    else
    {
        printf("Passed Lucas PRP test.    \n");
    }
    return 0;
}

int main(int argc, char **argv)
{
    //	mpz_inits(n, a, b, c, Q, g, np1, s, t, 0);
    //	mpz_inits(tmp, tmp2, tmp3, thou, barrett);

    mpz_init(n);
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
    mpz_init(Q);
    mpz_init(g);
    mpz_init(np1);
    mpz_init(s);
    mpz_init(t);
    mpz_init(tmp);
    mpz_init(tmp2);
    mpz_init(tmp3);
    mpz_init(thou);
    mpz_init(barrett);
    for (int i = 0; i++; i < FERMAT_THREAD_COUNT)
        mpz_init(fermat_args[i]);
    if (argc != 2)
    {
        printf("Usage: ./Lucas-and-Fermats <file-name>\n");
        exit(-1);
    }
    input_file = fopen(argv[1], "r");
    mpz_inp_str(n, input_file, 10);
    fclose(input_file);

    printf("This threaded program performs key Fermat PRP tests\n");
    printf("and Lucas (x+b)^(n+1)==b^2-a (mod n,x^2-a) where b=a-1.\n");

    gmp_printf("Input number has approximately %ld digits.\n", mpz_sizeinbase(n, 10));

    if (mpz_cmp_ui(n, 3) != 1)
    {
        printf("Too negative to be tested!\n");
        exit(1);
    }
    if (mpz_even_p(n))
    {
        printf("This program is for odd integers!\n");
        exit(1);
    }
    if (mpz_root(tmp, n, 2))
    {
        printf("Input is a square number.\n");
        exit(1);
    }

    mpz_set_ui(fermat_args[0], 2);
    mpz_add_ui(np1, n, 1);
    barrett_precalc();

    if (mpz_ui_kronecker(2, n) == -1)
    {
        printf("Parameter a=1/2.\n");
        mpz_gcd_ui(g, n, 3);
        if (mpz_cmp_ui(g, 1) != 0)
        {
            printf("Has factor!\n");
            exit(1);
        }
        mpz_set_ui(a, 2);
        mpz_set_ui(b, 1);
        pthread_create(&l1, NULL, lucas, NULL);
        pthread_create(&f2, NULL, fermat, &fermat_args[0]);
        pthread_join(f2, NULL);
        pthread_join(l1, NULL);
    }
    else
    {
        mpz_mod_ui(tmp, n, 4);
        if (mpz_cmp_ui(tmp, 3) == 0)
        {
            mpz_gcd_ui(g, n, 33);
            if (mpz_cmp_ui(g, 1) != 0)
            {
                printf("Has factor!\n");
                exit(1);
            }
            printf("Parameter a=-1.\n");
            mpz_set_si(a, -1);
            mpz_set_ui(b, 2);
            mpz_set_ui(fermat_args[1], 3);
            pthread_create(&l1, NULL, lucas, NULL);
            pthread_create(&f2, NULL, fermat, &fermat_args[0]);
            pthread_create(&fc, NULL, fermat, &fermat_args[1]);
            pthread_join(f2, NULL);
            pthread_join(fc, NULL);
            pthread_join(l1, NULL);
        }
        else
        {
            for (r = 1;; r++)
            {
                mpz_set_ui(a, 1);
                mpz_mul_2exp(b, a, r);
                mpz_add_ui(a, b, 1);
                mpz_mod(a, a, n);
                jac = mpz_jacobi(a, n);
                if (jac == 0 && mpz_cmp_ui(a, 0) != 0)
                {
                    printf("Has factor!\n");
                    exit(1);
                }
                if (jac == 1)
                    continue;
                printf("Parameter a=");
                mpz_out_str(NULL, 10, a);
                printf(".\n");
                // gcd(a^2-a+1,n)
                mpz_mul(tmp, b, b);
                mpz_add(tmp, tmp, a);
                mpz_gcd(g, tmp, n);
                if (mpz_cmp_ui(g, 1) != 0)
                {
                    printf("Has factor!\n");
                    exit(1);
                }
                // gcd(3*b^2+a,n)
                mpz_mul(g, b, b);
                mpz_mul_ui(g, g, 3);
                mpz_sub(g, g, a);
                mpz_gcd(g, g, n);
                if (mpz_cmp_ui(g, 1) != 0)
                {
                    printf("Has factor!\n");
                    exit(1);
                }
                mpz_set(fermat_args[1], a);
                mpz_mul_2exp(c, a, 1);
                mpz_sub_ui(fermat_args[2], c, 1);
                pthread_create(&l1, NULL, lucas, NULL);
                pthread_create(&f2, NULL, fermat, &fermat_args[0]);
                pthread_create(&fa, NULL, fermat, &fermat_args[1]);
                pthread_create(&fc, NULL, fermat, &fermat_args[2]);
                pthread_join(f2, NULL);
                pthread_join(fa, NULL);
                pthread_join(fc, NULL);
                pthread_join(l1, NULL);
                break;
            }
        }
    }
    if (composite)
    {
        printf("Composite!\n");
        exit(1);
    }
    else
    {
        printf("Likely Prime!\n");
        exit(0);
    }
}
