#pragma once

// -----------------------------------------------------------------------
// Cubic primality test
//
// interface to a multithreaded version for heavy computations
// -----------------------------------------------------------------------

#include <stdint.h>

void *quadratic_allocate_function(size_t alloc_size);
void *quadratic_reallocate_function(void *ptr, size_t old_size, size_t new_size);
void quadratic_free_function(void *ptr, size_t size);
