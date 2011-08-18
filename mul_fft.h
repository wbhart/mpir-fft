/* mul_fft -- split-radix fft routines for MPIR.

Copyright 2009, 2011 William Hart. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY William Hart ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL William Hart OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of William Hart.

*/

#include <stdio.h>
#include <stdlib.h>
#include "mpir.h"
#include "gmp-impl.h"

#ifndef MUL_FFT_H
#define MUL_FFT_H

/*
   Add the signed limb c to the value r which is an integer 
   modulo 2^GMP_LIMB_BITS*l + 1. We assume that the generic case
   is that c is very small and allow the compiler to inline that.
*/
static inline
void mpn_addmod_2expp1_1(mp_limb_t * r, mp_size_t l, mp_limb_signed_t c)
{
   mp_limb_t sum = r[0] + c;

   // check if adding c would cause a carry to propagate
   if ((mp_limb_signed_t)(sum ^ r[0]) >= 0)
      r[0] = sum;
   else
   {
      if (c >= 0) mpn_add_1(r, r, l + 1, c);
      else mpn_sub_1(r, r, l + 1, -c);
   }
}

void mpn_to_mpz(mpz_t m, mp_limb_t * i, mp_size_t limbs);

void set_p(mpz_t p, mp_size_t n, mp_bitcnt_t w);

void rand_n(mp_limb_t * n, gmp_randstate_t state, mp_size_t limbs);

void ref_mul_2expmod(mpz_t m, mpz_t i2, mpz_t p, mp_size_t n, mp_bitcnt_t w, mp_bitcnt_t d);

void ref_norm(mpz_t m, mpz_t p);

void ref_sumdiff_rshBmod(mpz_t t, mpz_t u, mpz_t i1,
                      mpz_t i2, mpz_t p, mp_size_t n, mp_bitcnt_t w, mp_bitcnt_t x, mp_bitcnt_t y);

void FFT_radix2_twiddle(mp_limb_t ** ii, mp_size_t is,
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
      mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs);

void IFFT_radix2_twiddle(mp_limb_t ** ii, mp_size_t is,
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
              mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs);

#endif

