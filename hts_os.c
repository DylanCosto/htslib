/// @file hts_os.c
/// Operating System specific tweaks, for compatibility with POSIX.
/*
   Copyright (C) 2017, 2019-2020 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>
#include "htslib/hts_defs.h"

#include <stddef.h>
#include <stdint.h>
#include "htslib/hts_os.h"

// Optional per-thread state for hts_drand48(), see hts_drand48_set_thread_state()
#if defined(_MSC_VER)
#define HTS_THREAD_LOCAL __declspec(thread)
#else
#define HTS_THREAD_LOCAL __thread
#endif
static HTS_THREAD_LOCAL unsigned short *hts_tls_xseed = NULL;
static HTS_THREAD_LOCAL uint64_t *hts_tls_ndraw = NULL;

HTSLIB_EXPORT
void hts_drand48_set_thread_state(unsigned short *xseed, uint64_t *ndraw)
{
    hts_tls_xseed = xseed;
    hts_tls_ndraw = xseed ? ndraw : NULL;
}

HTSLIB_EXPORT
void hts_rand48_skip(unsigned short xseed[3], uint64_t n)
{
    // x -> a*x + c (mod 2^48), applied n times by binary exponentiation
    const uint64_t mask = (((uint64_t)1)<<48) - 1;
    uint64_t a = 0x5DEECE66DULL, c = 0xB, A = 1, C = 0;
    while ( n )
    {
        if ( n & 1 ) { A = (A*a) & mask; C = (C*a + c) & mask; }
        c = (c*a + c) & mask;
        a = (a*a) & mask;
        n >>= 1;
    }
    uint64_t x = (uint64_t)xseed[0] | (uint64_t)xseed[1]<<16 | (uint64_t)xseed[2]<<32;
    x = (A*x + C) & mask;
    xseed[0] = x & 0xffff;
    xseed[1] = (x >> 16) & 0xffff;
    xseed[2] = (x >> 32) & 0xffff;
}

// Windows (maybe more) lack a drand48 implementation.
#ifndef HAVE_DRAND48
#include "os/rand.c"
#else
#include <stdlib.h>
#include <string.h>
HTSLIB_EXPORT
void hts_srand48(long seed)
{
#ifdef HAVE_SRAND48_DETERMINISTIC
    srand48_deterministic(seed);
#else
    srand48(seed);
#endif
}

HTSLIB_EXPORT
double hts_erand48(unsigned short xseed[3]) { return erand48(xseed); }

HTSLIB_EXPORT
double hts_drand48(void)
{
    if ( hts_tls_xseed )
    {
        if ( hts_tls_ndraw ) (*hts_tls_ndraw)++;
        return erand48(hts_tls_xseed);
    }
    return drand48();
}

HTSLIB_EXPORT
int hts_drand48_get_state(unsigned short xseed[3])
{
    // seed48() returns the previous state, put it back immediately
    unsigned short tmp[3] = {0,0,0}, cur[3];
    unsigned short *prev = seed48(tmp);
    memcpy(cur, prev, sizeof(cur));
    seed48(cur);
    memcpy(xseed, cur, sizeof(cur));
    return 0;
}

HTSLIB_EXPORT
long hts_lrand48(void) { return lrand48(); }
#endif

// // On Windows when using the MSYS or Cygwin terminals, isatty fails
// #ifdef _WIN32
// #define USE_FILEEXTD
// #include "os/iscygpty.c"
// #endif
