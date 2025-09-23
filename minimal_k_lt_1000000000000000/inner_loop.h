
#ifndef INNER_LOOP_H_H
#define INNER_LOOP_H_H

#include "tlv.h"

int inner_loop(int s, uint16_t cid, uint128_t seed, uint64_t count);
int inner_self_test(void);

#endif
