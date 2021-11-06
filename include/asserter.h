#ifndef ASSERTER_H
#define ASSERTER_H

#ifdef DEBUG
#include "assert.h"
inline __attribute__((always_inline)) void abort_here(int x) {throw x;}
#else
#define assert(...) { }
#define abort_here(...) { }
#endif

#endif
