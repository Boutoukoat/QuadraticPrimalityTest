#pragma once

// -----------------------------------------------------------------------
// Cubic primality test
//
// mpz_quadratic_primality():
//    true: composite for sure
//    false: might be prime
//
// quadratic_primality_self_test()
//    simplified unit tests to detect a possible compiler/platform issue.
//    assert when fail (this should not happen).
// -----------------------------------------------------------------------

#include "gmp.h"
#include <stdbool.h>

bool mpz_quadratic_primality(mpz_t v, bool verbose = false);
void quadratic_primality_self_test(void);
