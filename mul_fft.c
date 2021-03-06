/* mul_fft -- radix 2 fft routines for MPIR.

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

#define TEST 0
#define TIME 1

#include <stdio.h>
#include <stdlib.h>
#include "mpir.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mul_fft.h"

mp_limb_t new_mpn_mulmod_2expp1(mp_limb_t * r, mp_limb_t * i1, mp_limb_t * i2, 
                           mp_limb_t c, mp_limb_t bits, mp_limb_t * tt);

/*
   NOTES: throughout the following we use the following notation:
   
   * convolution length 2*n where n is a power of 2
   * p = 2^{wn} + 1 with wn divisible by GMP_LIMB_BITS
   * l = {wn}/GMP_LIMB_BITS (number of limbs)
*/

const mp_limb_t revtab0[1] = { 0 };
const mp_limb_t revtab1[2] = { 0, 1 };
const mp_limb_t revtab2[4] = { 0, 2, 1, 3 };
const mp_limb_t revtab3[8] = { 0, 4, 2, 6, 1, 5, 3, 7 };
const mp_limb_t revtab4[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

const mp_limb_t * revtab[5] = { revtab0, revtab1, revtab2, revtab3, revtab4 };

/*
   computes the reverse binary of a binary number of the given number of bits
 */
mp_limb_t mpir_revbin(mp_limb_t in, mp_bitcnt_t bits)
{
    mp_bitcnt_t i;
    mp_limb_t out = 0;
    
    if (bits <= 4)
        return revtab[bits][in];

    for (i = 0; i < bits; i++)
    {   
       out <<= 1;
       out += (in & 1);
       in >>= 1;
    }

    return out;
}

/*
   Splits an mpn into segments of length coeff_limbs and stores in 
   zero padded coefficients of length output_limbs, for use in FFT 
   convolution code. Assumes that the input is total_limbs in length. 
   The total number of coefficients written is returned.
*/
mp_size_t FFT_split(mp_limb_t ** poly, mp_limb_t * limbs, 
                mp_size_t total_limbs, mp_size_t coeff_limbs, mp_size_t output_limbs)
{
   mp_size_t length = (total_limbs - 1)/coeff_limbs + 1;
   mp_size_t i, j, skip;
   
   for (skip = 0, i = 0; skip + coeff_limbs <= total_limbs; skip += coeff_limbs, i++)
   {
      //if (i + 1 < length)
		   //for (j = 0; j + 8 < output_limbs; j += 8) PREFETCH(poly[i+1], j);
      
      MPN_ZERO(poly[i], output_limbs + 1);
      // convert a coefficient
      MPN_COPY(poly[i], limbs + skip, coeff_limbs);
   }
   if (i < length) MPN_ZERO(poly[i], output_limbs + 1);
   if (total_limbs > skip) MPN_COPY(poly[i], limbs + skip, total_limbs - skip);
   
   return length;
}

/*
   Splits an mpn into segments of length _bits_ and stores in
   zero padded coefficients of length output_limbs, for use in FFT 
   convolution code. Assumes that the input is total_limbs in length. 
   Returns the total number of coefficient written. 
*/

mp_size_t FFT_split_bits(mp_limb_t ** poly, mp_limb_t * limbs, 
               mp_size_t total_limbs, mp_size_t bits, mp_size_t output_limbs)
{
   mp_size_t length = (GMP_LIMB_BITS*total_limbs - 1)/bits + 1;
   mp_size_t i, j;
   
   mp_bitcnt_t top_bits = ((GMP_LIMB_BITS - 1) & bits);
   if (top_bits == 0)
   {
      return FFT_split(poly, limbs, total_limbs, bits/GMP_LIMB_BITS, output_limbs);
   }

   mp_size_t coeff_limbs = (bits/GMP_LIMB_BITS) + 1;
   mp_limb_t mask = (1L<<top_bits) - 1L;
   mp_bitcnt_t shift_bits = 0L;
   mp_limb_t * limb_ptr = limbs;                      
    
   for (i = 0; i < length - 1; i++)
   {
      //for (j = 0; j + 8 < output_limbs; j += 8) PREFETCH(poly->coeffs[i+1], j);
      
      MPN_ZERO(poly[i], output_limbs + 1);
      // convert a coefficient
      if (!shift_bits)
      {
         MPN_COPY(poly[i], limb_ptr, coeff_limbs);
         poly[i][coeff_limbs - 1] &= mask;
         limb_ptr += (coeff_limbs - 1);
         shift_bits += top_bits;
      } else
      {
         mpn_rshift(poly[i], limb_ptr, coeff_limbs, shift_bits);
         limb_ptr += (coeff_limbs - 1);
         shift_bits += top_bits;
         if (shift_bits >= GMP_LIMB_BITS)
         {
            limb_ptr++;
            poly[i][coeff_limbs - 1] += (limb_ptr[0] << (GMP_LIMB_BITS - (shift_bits - top_bits)));
            shift_bits -= GMP_LIMB_BITS; 
         }
         poly[i][coeff_limbs - 1] &= mask;
      }                      
   }
   
   MPN_ZERO(poly[i], output_limbs + 1);
   mp_size_t limbs_left = total_limbs - (limb_ptr - limbs);
   if (!shift_bits)
   {
      MPN_COPY(poly[i], limb_ptr, limbs_left);
   } else
   {
      mpn_rshift(poly[i], limb_ptr, limbs_left, shift_bits);
   }                      
      
   return length;
}

/*
   Recombines coefficients after doing a convolution. Assumes each of the 
   coefficients of the poly of the given length is output_limbs long, that each 
   of the coefficients is being shifted by a multiple of coeff_limbs and added
   to an mpn which is total_limbs long. It is assumed that the mpn has been 
   zeroed in advance.
*/

void FFT_combine(mp_limb_t * res, mp_limb_t ** poly, mp_size_t length, 
            mp_size_t coeff_limbs, mp_size_t output_limbs, mp_size_t total_limbs)
{
   mp_size_t skip, i, j;
   
   for (skip = 0, i = 0; (i < length) && (skip + output_limbs <= total_limbs); i++, skip += coeff_limbs)
   { 
      //for (j = 0; j < output_limbs; j += 8) PREFETCH(poly->coeffs[i+1], j);
      mpn_add(res + skip, res + skip, output_limbs + 1, poly[i], output_limbs);      
   } 

   while ((skip < total_limbs) && (i < length))
   {
      mpn_add(res + skip, res + skip, total_limbs - skip, poly[i], MIN(total_limbs - skip, output_limbs));
      i++;
      skip += coeff_limbs;
   }  
}

/*
   Recombines coefficients of a poly after doing a convolution. Assumes 
   each of the coefficients of the poly of the given length is output_limbs 
   long, that each is being shifted by a multiple of _bits_ and added
   to an mpn which is total_limbs long. It is assumed that the mpn has been 
   zeroed in advance.
*/

void FFT_combine_bits(mp_limb_t * res, mp_limb_t ** poly, mp_size_t length, 
                  mp_size_t bits, mp_size_t output_limbs, mp_size_t total_limbs)
{
   mp_bitcnt_t top_bits = ((GMP_LIMB_BITS - 1) & bits);
   if (top_bits == 0)
   {
      FFT_combine(res, poly, length, bits/GMP_LIMB_BITS, output_limbs, total_limbs);
      return;
   }
   TMP_DECL;
   TMP_MARK;

   mp_size_t coeff_limbs = (bits/GMP_LIMB_BITS) + 1;
   mp_size_t i, j;
   mp_limb_t * temp = (mp_limb_t *) TMP_BALLOC_LIMBS(output_limbs + 1);
   mp_bitcnt_t shift_bits = 0;
   mp_limb_t * limb_ptr = res;
   mp_limb_t * end = res + total_limbs;
   
   for (i = 0; (i < length) && (limb_ptr + output_limbs < end); i++)
   { 
      //for (j = 0; j < output_limbs; j += 8) PREFETCH(poly->coeffs[i+1], j);
      if (shift_bits)
      {
         mpn_lshift(temp, poly[i], output_limbs + 1, shift_bits);
         mpn_add_n(limb_ptr, limb_ptr, temp, output_limbs + 1);
      } else
      {
         mpn_add(limb_ptr, limb_ptr, output_limbs + 1, poly[i], output_limbs);
      }
      shift_bits += top_bits;
      limb_ptr += (coeff_limbs - 1);
      if (shift_bits >= GMP_LIMB_BITS)
      {
         limb_ptr++;
         shift_bits -= GMP_LIMB_BITS;
      }      
   } 

   while ((limb_ptr < end) && (i < length))
   {
      if (shift_bits)
      {
         mpn_lshift(temp, poly[i], output_limbs + 1, shift_bits);
         mpn_add_n(limb_ptr, limb_ptr, temp, end - limb_ptr);
      } else
      {
         mpn_add_n(limb_ptr, limb_ptr, poly[i], end - limb_ptr);
      }
      shift_bits += top_bits;
      limb_ptr += (coeff_limbs - 1);
      if (shift_bits >= GMP_LIMB_BITS)
      {
         limb_ptr++;
         shift_bits -= GMP_LIMB_BITS;
      }  
      i++;    
   }
   
   TMP_FREE;     
}

/*
   Normalise t to be in the range [0, 2^nw]
*/
void mpn_normmod_2expp1(mp_limb_t * t, mp_size_t l)
{
   mp_limb_signed_t hi = t[l];
   
   if (hi)
   {
      t[l] = CNST_LIMB(0);
      mpn_addmod_2expp1_1(t, l, -hi);
      hi = t[l];

      // hi will be in [-1,1]
      if (t[l])
      {
         t[l] = CNST_LIMB(0);
         mpn_addmod_2expp1_1(t, l, -hi);
         if (t[l] == ~CNST_LIMB(0)) // if we now have -1 (very unlikely)
         {
            t[l] = CNST_LIMB(0);
            mpn_addmod_2expp1_1(t, l, 1);
         }
      }
   }
}

/*
   We are given two integers modulo 2^wn+1, i1 and i2, which are 
   not necessarily normalised and are given n and w. We compute 
   t = (i1 + i2)*B^x, u = (i1 - i2)*B^y. Aliasing between inputs and 
   outputs is not permitted. We require x and y to be less than the 
   number of limbs of i1 and i2.
*/
void mpn_lshB_sumdiffmod_2expp1(mp_limb_t * t, mp_limb_t * u, mp_limb_t * i1, 
                      mp_limb_t * i2, mp_size_t limbs, mp_size_t x, mp_size_t y)
{
   mp_limb_t cy, cy1, cy2;

   if (x == 0)
   {
      if (y == 0)
      {
         mpn_sumdiff_n(t + x, u + y, i1, i2, limbs + 1);
      } else
      {
         cy = mpn_sumdiff_n(t, u + y, i1, i2, limbs - y);
         u[limbs] = -(cy&1);
         cy1 = cy>>1;
         cy = mpn_sumdiff_n(t + limbs - y, u, i2 + limbs - y, i1 + limbs - y, y);
         t[limbs] = cy>>1;
         mpn_add_1(t + limbs - y, t + limbs - y, y + 1, cy1);
         cy1 = -(cy&1) + (i2[limbs] - i1[limbs]);
         mpn_addmod_2expp1_1(u + y, limbs - y, cy1);
         cy1 = -(i1[limbs] + i2[limbs]);
         mpn_addmod_2expp1_1(t, limbs, cy1);
      }
   } else if (y == 0)
   {
      cy = mpn_sumdiff_n(t + x, u, i1, i2, limbs - x);
      t[limbs] = cy>>1;
      cy1 = cy&1;
      cy = mpn_sumdiff_n(t, u + limbs - x, i1 + limbs - x, i2 + limbs - x, x);
      cy2 = mpn_neg_n(t, t, x);
      u[limbs] = -(cy&1);
      mpn_sub_1(u + limbs - x, u + limbs - x, x + 1, cy1);
      cy1 = -(cy>>1) - cy2;
      cy1 -= (i1[limbs] + i2[limbs]);
      mpn_addmod_2expp1_1(t + x, limbs - x, cy1);
      cy1 = (i2[limbs] - i1[limbs]);
      mpn_addmod_2expp1_1(u, limbs, cy1);
   } else if (x > y)
   {
      cy = mpn_sumdiff_n(t + x, u + y, i1, i2, limbs - x);
      t[limbs] = cy>>1;
      cy1 = cy&1;
      cy = mpn_sumdiff_n(t, u + y + limbs - x, i1 + limbs - x, i2 + limbs - x, x - y);
      cy2 = mpn_neg_n(t, t, x - y);
      u[limbs] = -(cy&1);
      mpn_sub_1(u + y + limbs - x, u + y + limbs - x, x - y + 1, cy1);
      cy1 = (cy>>1) + cy2;
      cy = mpn_sumdiff_n(t + x - y, u, i2 + limbs - y, i1 + limbs - y, y);
      cy2 = mpn_neg_n(t + x - y, t + x - y, y);
      cy1 = -(cy>>1) - mpn_sub_1(t + x - y, t + x - y, y, cy1) - cy2;
      cy1 -= (i1[limbs] + i2[limbs]);
      mpn_addmod_2expp1_1(t + x, limbs - x, cy1);
      cy1 = -(cy&1) + (i2[limbs] - i1[limbs]);
      mpn_addmod_2expp1_1(u + y, limbs - y, cy1);
   } else if (x < y)
   {
      cy = mpn_sumdiff_n(t + x, u + y, i1, i2, limbs - y);
      u[limbs] = -(cy&1);
      cy1 = cy>>1;
      cy = mpn_sumdiff_n(t + x + limbs - y, u, i2 + limbs - y, i1 + limbs - y, y - x);
      t[limbs] = cy>>1;
      mpn_add_1(t + x + limbs - y, t + x + limbs - y, y - x + 1, cy1);
      cy1 = cy&1;
      cy = mpn_sumdiff_n(t, u + y - x, i2 + limbs - x, i1 + limbs - x, x);
      cy1 = -(cy&1) - mpn_sub_1(u + y - x, u + y - x, x, cy1);
      cy1 += (i2[limbs] - i1[limbs]);
      mpn_addmod_2expp1_1(u + y, limbs - y, cy1);
      cy2 = mpn_neg_n(t, t, x);
      cy1 = -(cy>>1) - (i1[limbs] + i2[limbs]) - cy2;
      mpn_addmod_2expp1_1(t + x, limbs - x, cy1);
   } else // x == y
   {
      cy = mpn_sumdiff_n(t + x, u + x, i1, i2, limbs - x);
      t[limbs] = cy>>1;
      u[limbs] = -(cy&1);
      cy = mpn_sumdiff_n(t, u, i2 + limbs - x, i1 + limbs - x, x);
      cy2 = mpn_neg_n(t, t, x);
      cy1 = -(cy>>1) - (i1[limbs] + i2[limbs]) - cy2;
      mpn_addmod_2expp1_1(t + x, limbs - x, cy1);
      cy1 = -(cy&1) + i2[limbs] - i1[limbs];
      mpn_addmod_2expp1_1(u + x, limbs - x, cy1);
  }
}

/*
   We are given two integers modulo 2^wn+1, i1 and i2, which are 
   not necessarily normalised and are given n and w. We compute 
   t = i1 + i2/B^x, u = i1 - i2/B^x. Aliasing between inputs and 
   outputs is not permitted. We require x be less than the 
   number of limbs of i1 and i2.
*/
void mpn_sumdiff_rshBmod_2expp1(mp_limb_t * t, mp_limb_t * u, mp_limb_t * i1, 
                      mp_limb_t * i2, mp_size_t limbs, mp_size_t x, mp_size_t y)
{
   mp_limb_t cy, cy1, cy2, cy3;

   if (x == 0)
   {
      if (y == 0)
      {
         mpn_sumdiff_n(t, u, i1, i2, limbs + 1);     
      } else // y != 0
      {
         cy = mpn_sumdiff_n(t, u, i1, i2 + y, limbs - y);
         cy1 = (cy>>1);
         cy2 = -(cy&1);
         cy = mpn_sumdiff_n(u + limbs - y, t + limbs - y, i1 + limbs - y, i2, y);
         u[limbs] = (cy>>1) + i1[limbs];
         t[limbs] = i1[limbs] - (cy&1);
         mpn_addmod_2expp1_1(t + limbs - y, y, cy1 + i2[limbs]);
         mpn_addmod_2expp1_1(u + limbs - y, y, cy2 - i2[limbs]);
      }
   } else if (y == 0) // x != 0
   {
      cy = mpn_sumdiff_n(t, u, i1 + x, i2, limbs - x);
      cy1 = (cy>>1);
      cy2 = -(cy&1);
      cy3 = mpn_neg_n(i1, i1, x);
      cy = mpn_sumdiff_n(t + limbs - x, u + limbs - x, i1, i2 + limbs - x, x);
      u[limbs] = -cy3 - (cy&1) - i2[limbs];
      t[limbs] = -cy3 + i2[limbs] + (cy>>1);
      mpn_addmod_2expp1_1(t + limbs - x, x, cy1 + i1[limbs]);
      mpn_addmod_2expp1_1(u + limbs - x, x, cy2 + i1[limbs]);
   } else if (x == y)
   {
      cy = mpn_sumdiff_n(t, u, i1 + x, i2 + x, limbs - x);
      cy1 = (cy>>1);
      cy2 = -(cy&1);
      cy = mpn_sumdiff_n(t + limbs - x, u + limbs - x, i2, i1, x);
      cy3 = mpn_neg_n(t + limbs - x, t + limbs - x, x);
      u[limbs] = -(cy&1);
      t[limbs] = -(cy>>1) - cy3;
      mpn_addmod_2expp1_1(t + limbs - x, x, cy1 + i1[limbs] + i2[limbs]);
      mpn_addmod_2expp1_1(u + limbs - x, x, cy2 + i1[limbs] - i2[limbs]);
   } else if (x > y)
   {
      cy = mpn_sumdiff_n(t + limbs - y, u + limbs - y, i2, i1 + x - y, y);
      cy3 = mpn_neg_n(t + limbs - y, t + limbs - y, y);
      t[limbs] = -(cy>>1) - cy3;
      u[limbs] = -(cy&1);
      cy3 = mpn_neg_n(i1, i1, x - y);
      cy = mpn_sumdiff_n(t + limbs - x, u + limbs - x, i1, i2 + limbs - x + y, x - y);
      mpn_addmod_2expp1_1(t + limbs - y, y, (cy>>1) + i2[limbs] - cy3);
      mpn_addmod_2expp1_1(u + limbs - y, y, -(cy&1) - i2[limbs] - cy3);
      cy = mpn_sumdiff_n(t, u, i1 + x, i2 + y, limbs - x);
      mpn_addmod_2expp1_1(t + limbs - x, x, (cy>>1) + i1[limbs]);
      mpn_addmod_2expp1_1(u + limbs - x, x, -(cy&1) + i1[limbs]);
   } else //(x < y)
   {
      cy = mpn_sumdiff_n(t + limbs - x, u + limbs - x, i2 + y - x, i1, x);
      cy3 = mpn_neg_n(t + limbs - x, t + limbs - x, x);
      t[limbs] = -(cy>>1) - cy3;
      u[limbs] = -(cy&1);
      cy3 = mpn_neg_n(i2, i2, y - x);
      cy = mpn_sumdiff_n(t + limbs - y, u + limbs - y, i1 + limbs - y + x, i2, y - x);
      mpn_addmod_2expp1_1(t + limbs - x, x, (cy>>1) + i1[limbs] - cy3);
      mpn_addmod_2expp1_1(u + limbs - x, x, -(cy&1) + i1[limbs] + cy3);
      cy = mpn_sumdiff_n(t, u, i1 + x, i2 + y, limbs - y);
      mpn_addmod_2expp1_1(t + limbs - y, y, (cy>>1) + i2[limbs]);
      mpn_addmod_2expp1_1(u + limbs - y, y, -(cy&1) - i2[limbs]);      
   }
}

/* 
   Given an integer i1 modulo 2^wn+1, set t to 2^d*i1 modulo 2^wm+1.
   We must have GMP_LIMB_BITS > d >= 0.
*/
void mpn_mul_2expmod_2expp1(mp_limb_t * t, mp_limb_t * i1, mp_size_t limbs, mp_bitcnt_t d)
{
   mp_limb_signed_t hi, hi2;
   
   if (d == 0)
   {   
      if (t != i1)
         MPN_COPY(t, i1, limbs + 1);
   } else
   {
      hi = i1[limbs]; 
      mpn_lshift(t, i1, limbs + 1, d);
      hi2 = t[limbs];
      t[limbs] = CNST_LIMB(0);
      mpn_sub_1(t, t, limbs + 1, hi2);
      hi >>= (GMP_LIMB_BITS - d);
      mpn_addmod_2expp1_1(t + 1, limbs - 1, -hi);
   }
}

/* 
   Given an integer i1 modulo 2^wn+1, set t to i1/2^d modulo 2^wm+1.
   We must have GMP_LIMB_BITS > d >= 0.
*/
void mpn_div_2expmod_2expp1(mp_limb_t * t, mp_limb_t * i1, mp_size_t limbs, mp_bitcnt_t d)
{
   mp_limb_t lo;
   mp_limb_t * ptr;
   mp_limb_signed_t hi;
   
   if (d == 0)
   {   
      if (t != i1)
         MPN_COPY(t, i1, limbs + 1);
   } else
   {
      hi = i1[limbs];
      lo = mpn_rshift(t, i1, limbs + 1, d);
      t[limbs] = (hi>>d);
      ptr = t + limbs - 1;
      sub_ddmmss(ptr[1], ptr[0], ptr[1], ptr[0], CNST_LIMB(0), lo);
   }
}

/*
   Set u = 2^{ws*tw1}*(s + t), v = 2^{w+ws*tw2}*(s - t)
*/
void FFT_radix2_twiddle_butterfly(mp_limb_t * u, mp_limb_t * v, 
          mp_limb_t * s, mp_limb_t * t, mp_size_t NW, mp_bitcnt_t b1, mp_bitcnt_t b2)
{
   mp_limb_t size = NW/GMP_LIMB_BITS + 1;
   mp_size_t x, y;
   int negate = 0;
   int negate2 = 0;
   
   b1 %= (2*NW);
   if (b1 >= NW) 
   {
      negate2 = 1;
      b1 -= NW;
   }
   x = b1/GMP_LIMB_BITS;
   b1 -= x*GMP_LIMB_BITS;

   b2 %= (2*NW);
   if (b2 >= NW) 
   {
      negate = 1;
      b2 -= NW;
   }
   y = b2/GMP_LIMB_BITS;
   b2 -= y*GMP_LIMB_BITS;
 
   mpn_lshB_sumdiffmod_2expp1(u, v, s, t, size - 1, x, y);
   mpn_mul_2expmod_2expp1(u, u, size - 1, b1);
   if (negate2) mpn_neg_n(u, u, size);
   mpn_mul_2expmod_2expp1(v, v, size - 1, b2);
   if (negate) mpn_neg_n(v, v, size);
}
    
/*
   Set s = i1 + i2, t = z1^i*(i1 - i2) where z1 = exp(2*Pi*I/m) => w bits
*/
void FFT_radix2_butterfly(mp_limb_t * s, mp_limb_t * t, 
                  mp_limb_t * i1, mp_limb_t * i2, mp_size_t i, mp_size_t n, mp_bitcnt_t w)
{
   mp_limb_t size = (w*n)/GMP_LIMB_BITS + 1;
   mp_size_t x, y;
   mp_bitcnt_t b1;
   int negate = 0;

   x = 0;

   b1 = i;
   while (b1 >= n) 
   {
      negate = 1 - negate;
      b1 -= n;
   }
   b1 = b1*w;
   y = b1/GMP_LIMB_BITS;
   b1 -= y*GMP_LIMB_BITS;
 
   mpn_lshB_sumdiffmod_2expp1(s, t, i1, i2, size - 1, x, y);
   mpn_mul_2expmod_2expp1(t, t, size - 1, b1);
   if (negate) mpn_neg_n(t, t, size);
}

/*
   Let w = 2k + 1, i = 2j + 1.
   
   Set s = i1 + i2, t = z1^i*(i1 - i2) where z1 = exp(2*Pi*I/m)

   Here z1 corresponds to multiplication by (2^{3nw/4} - 2^{nw/4})*2^k.

   We see z1^i corresponds to multiplication by
   (2^{3nw/4} - 2^{nw/4})*2^{j+ik}.

   We first multiply by 2^{j + wn/4 + ik}, then again by a further
   2^{wn/2} and subtract.
*/
void FFT_radix2_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t, 
  mp_limb_t * i1, mp_limb_t * i2, mp_size_t i, mp_size_t n, mp_bitcnt_t w, mp_limb_t * temp)
{
   mp_bitcnt_t wn = w*n;
   mp_limb_t size = wn/GMP_LIMB_BITS, cy;
   mp_size_t j = i/2, k = w/2;
   mp_size_t y;
   mp_bitcnt_t b1;
   int negate = 0;

   b1 = j + wn/4 + i*k;
   while (b1 >= wn) 
   {
      negate = 1 - negate;
      b1 -= wn;
   }

   y = b1/GMP_LIMB_BITS;
   b1 -= y*GMP_LIMB_BITS;
 
   /* sumdiff and multiply by 2^{j + wn/4 + i*k} */
   mpn_lshB_sumdiffmod_2expp1(s, t, i1, i2, size, 0, y);
   mpn_mul_2expmod_2expp1(t, t, size, b1);
   if (negate) mpn_neg_n(t, t, size + 1);

   /* multiply by 2^{wn/2} */
   y = size/2;
   
   MPN_COPY(temp + y, t, size - y);
   temp[size] = 0;
   cy = mpn_neg_n(temp, t + size - y, y);
   if ((mp_limb_signed_t) t[size] < 0)
       mpn_add_1(temp + y, temp + y, size - y + 1, -t[size]);
   else
       mpn_sub_1(temp + y, temp + y, size - y + 1, t[size]);
   mpn_sub_1(temp + y, temp + y, size - y + 1, cy); 
   
   /* shift by an additional half limb (rare) */
   if (size & 1) 
       mpn_mul_2expmod_2expp1(temp, temp, size, GMP_LIMB_BITS/2);

   /* subtract */
   mpn_sub_n(t, temp, t, size + 1);
}

/*
   Set s = i1 + z1^i*i2, t = i1 - z1^i*i2 where z1 = exp(-2*Pi*I/m) => w bits
*/
void FFT_radix2_inverse_butterfly(mp_limb_t * s, mp_limb_t * t, 
                  mp_limb_t * i1, mp_limb_t * i2, mp_size_t i, mp_size_t n, mp_bitcnt_t w)
{
   mp_limb_t limbs = (w*n)/GMP_LIMB_BITS;
   mp_size_t y;
   mp_bitcnt_t b1;
   
   b1 = i*w;
   y = b1/GMP_LIMB_BITS;
   b1 -= y*GMP_LIMB_BITS;

   mpn_div_2expmod_2expp1(i2, i2, limbs, b1);
   mpn_sumdiff_rshBmod_2expp1(s, t, i1, i2, limbs, 0, y);
}

/*
   Let w = 2k + 1, i = 2j + 1.
   
   Set s = i1 + z1*i2, t = i1 - z1*i2 where z1 = exp(-2*Pi*I/m)

   Here z1 corresponds to division by 2^k then division by 
   (2^{3nw/4} - 2^{nw/4}).

   We see z1^i corresponds to division by
   (2^{3nw/4} - 2^{nw/4})*2^{j+ik} which is the same as division by
   2^{j+ik + 1} then multiplication by (2^{3nw/4} - 2^{nw/4}).

   Of course, division by 2^{j+ik + 1} is the same as multiplication by
   2^{2*wn - j - ik - 1}. The exponent is positive as i <= 2*n, j < n, k < w/2.

   We first multiply by 2^{2*wn - j - ik - 1 + wn/4} then multiply by an
   additional 2^{nw/2} and subtract.

*/
void FFT_radix2_inverse_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t, 
  mp_limb_t * i1, mp_limb_t * i2, mp_size_t i, mp_size_t n, mp_bitcnt_t w, mp_limb_t * temp)
{
   mp_bitcnt_t wn = w*n;
   mp_limb_t size = wn/GMP_LIMB_BITS, cy;
   mp_size_t j = i/2, k = w/2;
   mp_size_t y2, y;
   mp_size_t b1;
   int negate = 0;

   b1 = 2*wn - j - i*k - 1 + wn/4;
   while (b1 >= wn) 
   {
      negate = 1 - negate;
      b1 -= wn;
   }

   y2 = b1/GMP_LIMB_BITS;
   b1 -= y2*GMP_LIMB_BITS;

   /* multiply by small part of 2^{2*wn - j - ik - 1 + wn/4} */
   if (b1) mpn_mul_2expmod_2expp1(i2, i2, size, b1);
   
   /* multiply by 2^{wn/2} */
   y = size/2;
    
   MPN_COPY(temp + y, i2, size - y);
   temp[size] = 0;
   cy = mpn_neg_n(temp, i2 + size - y, y);
   if ((mp_limb_signed_t) i2[size] < 0)
       mpn_add_1(temp + y, temp + y, size - y + 1, -i2[size]);
   else
       mpn_sub_1(temp + y, temp + y, size - y + 1, i2[size]);
   mpn_sub_1(temp + y, temp + y, size - y + 1, cy); 
   
   /* shift by an additional half limb (rare) */
   if (size & 1) 
      mpn_mul_2expmod_2expp1(temp, temp, size, GMP_LIMB_BITS/2);

   /* subtract and negate... */
   if (negate) mpn_sub_n(i2, temp, i2, size + 1);
   else mpn_sub_n(i2, i2, temp, size + 1);

   /* ...negate and shift **left** by y2 limbs (i.e. shift right by 
   (size - y2) limbs) and sumdiff */
   mpn_sumdiff_rshBmod_2expp1(s, t, i1, i2, size, 0, size - y2);
}

void FFT_radix2_twiddle_inverse_butterfly(mp_limb_t * s, mp_limb_t * t, 
                  mp_limb_t * i1, mp_limb_t * i2, mp_size_t NW, mp_bitcnt_t b1, mp_bitcnt_t b2)
{
   mp_limb_t limbs = NW/GMP_LIMB_BITS;
   mp_size_t x, y;
   int negate = 0;
   int negate2 = 0;
   
   b1 %= (2*NW);
   if (b1 >= NW)
   {
      negate = 1;
      b1 -= NW;
   }
   x = b1/GMP_LIMB_BITS;
   b1 -= x*GMP_LIMB_BITS;

   b2 %= (2*NW);
   if (b2 >= NW)
   {
      negate2 = 1;
      b2 -= NW;
   }
   y = b2/GMP_LIMB_BITS;
   b2 -= y*GMP_LIMB_BITS;

   if (negate) mpn_neg_n(i1, i1, limbs + 1);
   mpn_div_2expmod_2expp1(i1, i1, limbs, b1);
   if (negate2) mpn_neg_n(i2, i2, limbs + 1);
   mpn_div_2expmod_2expp1(i2, i2, limbs, b2);
   mpn_sumdiff_rshBmod_2expp1(s, t, i1, i2, limbs, x, y);
}

/* 
   The radix 2 DIF FFT works as follows:
   Given: inputs [i0, i1, ..., i{m-1}], for m a power of 2

   Output: Fradix2[i0, i1, ..., i{m-1}] = [r0, r1, ..., r{m-1}]
   (Inputs and outputs may be subject to a stride and not in consecutive locations.)

   Algorithm: * Recursively compute [r0, r2, r4, ...., r{m-2}] 
                           = Fradix21[i0+i{m/2}, i1+i{m/2+1}, ..., i{m/2-1}+i{m-1}]
              * Let [t0, t1, ..., t{m/2-1}] 
                    = [i0-i{m/2}, i1-i{m/2+1}, ..., i{m/2-1}-i{m-1}]
              * Let [u0, u1, ..., u{m/2-1}] 
                    = [z1^0*t0, z1^1*t1, ..., z1^{m/2-1}*t{m/2-1}]
                where z1 = exp(2*Pi*I/m)
              * Recursively compute [r1, r3, ..., r{m-1}]
                           = Fradix2[u0, u1, ..., u{m/2-1}]
              
   The parameters are as follows:
              * 2*n is the length of the input and output arrays (m in the above)
              * w is such that 2^w is an 2n-th root of unity in the ring Z/pZ that 
                we are working in, i.e. p = 2^wn + 1 (here n is divisible by 
                GMP_LIMB_BITS)
              * ii is the array of inputs (each input is an array of limbs of length 
                wn/GMP_LIMB_BITS + 1 (the extra limb being a "carry limb")
              * rr is the array of outputs (each being an array of limbs of length
                wn/GMP_LIMB_BITS + 1)
              * rs is the stride with which entries of rr should be written (e.g. 
                rs = 4 means that the output is to be written in every 4th entry of 
                rr)
              * is is the stride with which entries of ii should be read
*/

void FFT_radix2(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (n == 1) 
   {
      FFT_radix2_butterfly(*t1, *t2, ii[0], ii[1], 0, n, w);
      ptr = rr[0];
      rr[0] = *t1;
      *t1 = ptr;
      ptr = rr[rs];
      rr[rs] = *t2;
      *t2 = ptr;
      return;
   }

   // [s0, s1, ..., s{m/2}] = [i0+i{m/2}, i1+i{m/2+1}, ..., i{m/2-1}+i{m-1}]
   // [t0, t1, ..., t{m/2-1}] 
   // = [z1^0*(i0-i{m/2}), z1^1*(i1-i{m/2+1}), ..., z1^{m/2-1}*(i{m/2-1}-i{m-1})]
   // where z1 = exp(2*Pi*I/m), z1 => w bits
   for (i = 0; i < n; i++) 
   {   
      FFT_radix2_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);
   
      ptr = ii[i];
      ii[i] = *t1;
      *t1 = ptr;
      ptr = ii[n+i];
      ii[n+i] = *t2;
      *t2 = ptr;
   }

   // [r0, r2, ..., r{m-2}] = Fradix2[s0, s1, ..., s{m/2-1}]
   FFT_radix2(rr, 1, ii, n/2, 2*w, t1, t2, temp);
   
   // [r1, r3, ..., r{m-1}] = Fradix2[t0, t1, ..., t{m/2-1}]
   FFT_radix2(rr + n, 1, ii+n, n/2, 2*w, t1, t2, temp);
}

/* 
   As for FFT_radix2 except that the length of the input and outputs is 2m = 4n and
   it uses a 2m-th root of unity which is sqrt(2)^w where sqrt(2) is the Schoenhage
   square root of 2 in Z/pZ.

   As p = 2^nw + 1 the square root of 2 is 2^{3nw/4} - 2^{nw/4}. Assuming w = 2k + 1
   the root of unity we are after is (2^{3nw/4} - 2^{nw/4})*2^k.

   Of course if w = 2k then we can just use 2^k as our root of unity.
*/
void FFT_radix2_sqrt2(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp)
{
   mp_limb_t * ptr;
   mp_size_t i;
   
   // [s0, s1, ..., s{m/2}] = [i0+i{m/2}, i1+i{m/2+1}, ..., i{m/2-1}+i{m-1}]
   // [t0, t1, ..., t{m/2-1}] 
   // = [z1^0*(i0-i{m/2}), z1^1*(i1-i{m/2+1}), ..., z1^{m/2-1}*(i{m/2-1}-i{m-1})]
   // where z1 = exp(2*Pi*I/m), z1 => sqrt(2)^w 

   if ((w & 1) == 0)
   {
      FFT_radix2(rr, rs, ii, 2*n, w/2, t1, t2, temp);

      return;
   }
   
   for (i = 0; i < 2*n; i++) 
   {   
      FFT_radix2_butterfly(*t1, *t2, ii[i], ii[2*n+i], i/2, n, w);
   
      ptr = ii[i];
      ii[i] = *t1;
      *t1 = ptr;
      ptr = ii[2*n+i];
      ii[2*n+i] = *t2;
      *t2 = ptr;

      i++;
      
      FFT_radix2_butterfly_sqrt2(*t1, *t2, ii[i], ii[2*n+i], i, n, w, *temp);

      ptr = ii[i];
      ii[i] = *t1;
      *t1 = ptr;
      ptr = ii[2*n+i];
      ii[2*n+i] = *t2;
      *t2 = ptr;
   }

   // [r0, r2, ..., r{m-2}] = Fradix2[s0, s1, ..., s{m/2-1}]
   FFT_radix2(rr, 1, ii, n, w, t1, t2, temp);
   
   // [r1, r3, ..., r{m-1}] = Fradix2[t0, t1, ..., t{m/2-1}]
   FFT_radix2(rr + 2*n, 1, ii+2*n, n, w, t1, t2, temp);
}

int FFT_negacyclic_twiddle(mp_limb_t * r, mp_limb_t * i1, mp_size_t i, mp_size_t n, mp_bitcnt_t w)
{
   mp_limb_t cy;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   int negate = 0;
   while (i >= 2*n)
   {
      negate = 1 - negate;
      i -= 2*n;
   }
   mp_bitcnt_t b1 = (w*i)/2;
   mp_size_t x = b1/GMP_LIMB_BITS;
   b1 -= x*GMP_LIMB_BITS;
   //if ((!x) && (negate)) 
   if (negate) mpn_neg_n(i1, i1, limbs + 1);
   mpn_mul_2expmod_2expp1(i1, i1, limbs, b1);
   if (x)
   {
      /*if (negate)
      {
         r[limbs] = -mpn_neg_n(r + x, i1, limbs - x);
         MPN_COPY(r, i1 + limbs - x, x);
         mpn_addmod_2expp1_1(r + x, limbs - x, i1[limbs]);
      } else
      {*/
         MPN_COPY(r + x, i1, limbs - x);
         r[limbs] = CNST_LIMB(0);
         cy = mpn_neg_n(r, i1 + limbs - x, x);
         mpn_addmod_2expp1_1(r + x, limbs - x, -i1[limbs]);
         mpn_sub_1(r + x, r + x, limbs - x + 1, cy); 
      //}
      return 1;
   }
   return 0;
}

/* 
   Set r to i1*z1^i for i < 2n, where z1 corresponds to shifting by w bits
*/
void FFT_twiddle(mp_limb_t * r, mp_limb_t * i1, mp_size_t i, mp_size_t n, mp_bitcnt_t w)
{
   mp_limb_t cy;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   int negate = 0;
   while (i >= n)
   {
      negate = 1 - negate;
      i -= n;
   }
   mp_bitcnt_t b1 = w*i;
   mp_size_t x = b1/GMP_LIMB_BITS;
   b1 -= x*GMP_LIMB_BITS;
   if (x)
   {
      MPN_COPY(r + x, i1, limbs - x);
      r[limbs] = 0;
      cy = mpn_neg_n(r, i1 + limbs - x, x);
      mpn_addmod_2expp1_1(r + x, limbs - x, - i1[limbs]);
      mpn_sub_1(r + x, r + x, limbs - x + 1, cy); 
      if (negate) mpn_neg_n(r, r, limbs + 1);
      mpn_mul_2expmod_2expp1(r, r, limbs, b1);
   } else
   {
      if (negate) 
      {
         mpn_neg_n(r, i1, limbs + 1);
         mpn_mul_2expmod_2expp1(r, r, limbs, b1);
      } else
         mpn_mul_2expmod_2expp1(r, i1, limbs, b1);
   }
}

/*
   Let w = 2k + 1, i = 2j + 1 for i < 4n.
   
   Set r = z1^i*i1 where z1 = exp(2*Pi*I/m)

   Here z1 corresponds to multiplication by (2^{3nw/4} - 2^{nw/4})*2^k.

   We see z1^i corresponds to multiplication by
   (2^{3nw/4} - 2^{nw/4})*2^{j+ik}.

   We first multiply by 2^{j + wn/4 + ik}, then again by a further
   2^{wn/2} and subtract.
*/
void FFT_twiddle_sqrt2(mp_limb_t * r, 
       mp_limb_t * i1, mp_size_t i, mp_size_t n, mp_bitcnt_t w, mp_limb_t * temp)
{
   mp_bitcnt_t wn = w*n;
   mp_limb_t size = wn/GMP_LIMB_BITS, cy;
   mp_size_t j = i/2, k = w/2;
   mp_size_t y;
   mp_bitcnt_t b1;
   int negate = 0;

   b1 = j + wn/4 + i*k;
   while (b1 >= wn) 
   {
      negate = 1 - negate;
      b1 -= wn;
   }

   y = b1/GMP_LIMB_BITS;
   b1 -= y*GMP_LIMB_BITS;
 
   /* multiply by 2^{j + wn/4 + i*k} */
   if (y)
   {
      mpn_copyi(temp + y, i1, size - y);
      cy = mpn_neg_n(temp, i1 + size - y, y);
      temp[size] = 0;
      mpn_addmod_2expp1_1(temp + y, size - y, -i1[size]);
      mpn_sub_1(temp + y, temp + y, size - y + 1, cy); 
      mpn_mul_2expmod_2expp1(r, temp, size, b1);
      if (negate) mpn_neg_n(r, r, size + 1);
   } else
   {
      mpn_mul_2expmod_2expp1(r, i1, size, b1);
      if (negate) mpn_neg_n(r, r, size + 1);
   }
   /* multiply by 2^{wn/2} */
   y = size/2;
   
   MPN_COPY(temp + y, r, size - y);
   temp[size] = 0;
   cy = mpn_neg_n(temp, r + size - y, y);
   mpn_addmod_2expp1_1(temp + y, size - y, -r[size]);
   mpn_sub_1(temp + y, temp + y, size - y + 1, cy); 
   
   /* shift by an additional half limb (rare) */
   if (size & 1) 
       mpn_mul_2expmod_2expp1(temp, temp, size, GMP_LIMB_BITS/2);

   /* subtract */
   mpn_sub_n(r, temp, r, size + 1);
}

/* 
   Truncate FFT to any length given by trunc, so long as trunc is 
   divisible by 8.
*/
void FFT_radix2_truncate1(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, 
      mp_size_t trunc)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (trunc == 2*n)
   {
      FFT_radix2(rr, rs, ii, n, w, t1, t2, temp);
      return;
   }

   if (trunc <= n)
   {
      for (i = 0; i < n; i++)
         mpn_add_n(ii[i], ii[i], ii[i+n], size);
      
      FFT_radix2_truncate1(rr, rs, ii, n/2, 2*w, t1, t2, temp, trunc);
   } else
   {

      // [s0, s1, ..., s{m/2}] = [i0+i{m/2}, i1+i{m/2+1}, ..., i{m/2-1}+i{m-1}]
      // [t0, t1, ..., t{m/2-1}] 
      // = [z1^0*(i0-i{m/2}), z1^1*(i1-i{m/2+1}), ..., z1^{m/2-1}*(i{m/2-1}-i{m-1})]
      // where z1 = exp(2*Pi*I/m), z1 => w bits
      for (i = 0; i < n; i++) 
      {   
         FFT_radix2_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);
   
         ptr = ii[i];
         ii[i] = *t1;
         *t1 = ptr;
         ptr = ii[n+i];
         ii[n+i] = *t2;
         *t2 = ptr;
      }

      // [r0, r2, ..., r{m-2}] = Fradix2[s0, s1, ..., s{m/2-1}]
      FFT_radix2(rr, 1, ii, n/2, 2*w, t1, t2, temp);
   
      // [r1, r3, ..., r{m-1}] = Fradix2[t0, t1, ..., t{m/2-1}]
      FFT_radix2_truncate1(rr + n, 1, ii+n, n/2, 2*w, t1, t2, temp, trunc - n);
   }
}

void FFT_radix2_truncate1_twiddle(mp_limb_t ** ii, mp_size_t is,
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
      mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs, mp_size_t trunc)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (trunc == 2*n)
   {
      FFT_radix2_twiddle(ii, is, n, w, t1, t2, temp, ws, r, c, rs);
      return;
   }

   if (trunc <= n)
   {
      for (i = 0; i < n; i++)
         mpn_add_n(ii[i*is], ii[i*is], ii[(i+n)*is], size);
      
      FFT_radix2_truncate1_twiddle(ii, is, n/2, 2*w, t1, t2, temp, ws, r, c, 2*rs, trunc);
   } else
   {

      // [s0, s1, ..., s{m/2}] = [i0+i{m/2}, i1+i{m/2+1}, ..., i{m/2-1}+i{m-1}]
      // [t0, t1, ..., t{m/2-1}] 
      // = [z1^0*(i0-i{m/2}), z1^1*(i1-i{m/2+1}), ..., z1^{m/2-1}*(i{m/2-1}-i{m-1})]
      // where z1 = exp(2*Pi*I/m), z1 => w bits
      for (i = 0; i < n; i++) 
      {   
         FFT_radix2_butterfly(*t1, *t2, ii[i*is], ii[(n+i)*is], i, n, w);
   
         ptr = ii[i*is];
         ii[i*is] = *t1;
         *t1 = ptr;
         ptr = ii[(n+i)*is];
         ii[(n+i)*is] = *t2;
         *t2 = ptr;
      }

      // [r0, r2, ..., r{m-2}] = Fradix2[s0, s1, ..., s{m/2-1}]
      FFT_radix2_twiddle(ii, is, n/2, 2*w, t1, t2, temp, ws, r, c, 2*rs);
   
      // [r1, r3, ..., r{m-1}] = Fradix2[t0, t1, ..., t{m/2-1}]
      FFT_radix2_truncate1_twiddle(ii + n*is, is, n/2, 2*w, t1, t2, temp, ws, r + rs, c, 2*rs, trunc - n);
   }
}

/* 
   Truncate FFT to any length given by trunc, so long as trunc is 
   divisible by 8. Assumes zeros from trunc to n.
*/
void FFT_radix2_truncate(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, 
      mp_size_t trunc)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (trunc == 2*n)
   {
      FFT_radix2(rr, rs, ii, n, w, t1, t2, temp);
      return;
   }

   if (trunc <= n)
   {
      // PASS
      FFT_radix2_truncate(rr, rs, ii, n/2, 2*w, t1, t2, temp, trunc);
   } else
   {
      // PASS
      // [s0, s1, ..., s{m/2}] = [i0+i{m/2}, i1+i{m/2+1}, ..., i{m/2-1}+i{m-1}]
      // [t0, t1, ..., t{m/2-1}] 
      // = [z1^0*(i0-i{m/2}), z1^1*(i1-i{m/2+1}), ..., z1^{m/2-1}*(i{m/2-1}-i{m-1})]
      // where z1 = exp(2*Pi*I/m), z1 => w bits
      for (i = 0; i < trunc - n; i++) 
      {   
         FFT_radix2_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);
   
         ptr = ii[i];
         ii[i] = *t1;
         *t1 = ptr;
         ptr = ii[n+i];
         ii[n+i] = *t2;
         *t2 = ptr;
      }

      for (i = trunc; i < 2*n; i++)
      {
         FFT_twiddle(ii[i], ii[i-n], i - n, n, w); 
      }
   
      // [r0, r2, ..., r{m-2}] = Fradix2[s0, s1, ..., s{m/2-1}]
      FFT_radix2(rr, 1, ii, n/2, 2*w, t1, t2, temp);

      // [r1, r3, ..., r{m-1}] = Fradix2[t0, t1, ..., t{m/2-1}]
      FFT_radix2_truncate1(rr + n, 1, ii+n, n/2, 2*w, t1, t2, temp, trunc - n);
   }
}

void FFT_radix2_truncate_twiddle(mp_limb_t ** ii, mp_size_t is,
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
      mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs, mp_size_t trunc)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (trunc == 2*n)
   {
      FFT_radix2_twiddle(ii, is, n, w, t1, t2, temp, ws, r, c, rs);
      return;
   }

   if (trunc <= n)
   {
      // PASS
      FFT_radix2_truncate_twiddle(ii, is, n/2, 2*w, t1, t2, temp, ws, r, c, 2*rs, trunc);
   } else
   {
      // PASS
      // [s0, s1, ..., s{m/2}] = [i0+i{m/2}, i1+i{m/2+1}, ..., i{m/2-1}+i{m-1}]
      // [t0, t1, ..., t{m/2-1}] 
      // = [z1^0*(i0-i{m/2}), z1^1*(i1-i{m/2+1}), ..., z1^{m/2-1}*(i{m/2-1}-i{m-1})]
      // where z1 = exp(2*Pi*I/m), z1 => w bits
      for (i = 0; i < trunc - n; i++) 
      {   
         FFT_radix2_butterfly(*t1, *t2, ii[i*is], ii[(n+i)*is], i, n, w);
   
         ptr = ii[i*is];
         ii[i*is] = *t1;
         *t1 = ptr;
         ptr = ii[(n+i)*is];
         ii[(n+i)*is] = *t2;
         *t2 = ptr;
      }

      for (i = trunc; i < 2*n; i++)
      {
         FFT_twiddle(ii[i*is], ii[(i-n)*is], i - n, n, w); 
      }
   
      // [r0, r2, ..., r{m-2}] = Fradix2[s0, s1, ..., s{m/2-1}]
      FFT_radix2_twiddle(ii, is, n/2, 2*w, t1, t2, temp, ws, r, c, 2*rs);

      // [r1, r3, ..., r{m-1}] = Fradix2[t0, t1, ..., t{m/2-1}]
      FFT_radix2_truncate1_twiddle(ii + n*is, is, n/2, 2*w, t1, t2, temp, ws, r + rs, c, 2*rs, trunc - n);
   }
}

void FFT_radix2_truncate_sqrt2(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, 
      mp_size_t trunc)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (trunc == 4*n)
   {
      FFT_radix2_sqrt2(rr, rs, ii, n, w, t1, t2, temp);
      return;
   }

   if ((w & 1) == 0)
   {
      FFT_radix2_truncate(rr, rs, ii, 2*n, w/2, t1, t2, temp, trunc);
      return;
   }
   
   for (i = 0; i < trunc - 2*n; i++) 
   {   
      FFT_radix2_butterfly(*t1, *t2, ii[i], ii[2*n+i], i/2, n, w);
   
      ptr = ii[i];
      ii[i] = *t1;
      *t1 = ptr;
      ptr = ii[2*n+i];
      ii[2*n+i] = *t2;
      *t2 = ptr;

      i++;
      
      FFT_radix2_butterfly_sqrt2(*t1, *t2, ii[i], ii[2*n+i], i, n, w, *temp);

      ptr = ii[i];
      ii[i] = *t1;
      *t1 = ptr;
      ptr = ii[2*n+i];
      ii[2*n+i] = *t2;
      *t2 = ptr;
   }

   for (i = trunc; i < 4*n; i++)
   {
      FFT_twiddle(ii[i], ii[i - 2*n], i/2 - n, n, w); 
         
      i++;

      FFT_twiddle_sqrt2(ii[i], ii[i - 2*n], i - 2*n, n, w, *temp); 
   }
   
   // [r0, r2, ..., r{m-2}] = Fradix2[s0, s1, ..., s{m/2-1}]
   FFT_radix2(rr, 1, ii, n, w, t1, t2, temp);

   // [r1, r3, ..., r{m-1}] = Fradix2[t0, t1, ..., t{m/2-1}]
   FFT_radix2_truncate1(rr + 2*n, 1, ii + 2*n, n, w, t1, t2, temp, trunc - 2*n);
}

void FFT_radix2_negacyclic(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp)
{
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t limbs = (w*n)/GMP_LIMB_BITS;
   
   /*
      we first apply twiddle factors corresponding to shifts of w*i/2 bits
   */
   
   if (w & 1)
   {
      for (i = 0; i < n; i++) 
      {   
          FFT_twiddle(*t1, ii[i], i/2, n, w);
          ptr = ii[i];
          ii[i] = *t1;
          *t1 = ptr;
            
          FFT_twiddle(*t1, ii[n+i], (n+i)/2, n, w);
          ptr = ii[n+i];
          ii[n+i] = *t1;
          *t1 = ptr;

          FFT_radix2_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);
   
          ptr = ii[i];
          ii[i] = *t1;
          *t1 = ptr;
          ptr = ii[n+i];
          ii[n+i] = *t2;
          *t2 = ptr;

          i++;
          
          FFT_twiddle_sqrt2(*t1, ii[i], i, n, w, *temp);
          ptr = ii[i];
          ii[i] = *t1;
          *t1 = ptr;
          
          FFT_twiddle_sqrt2(*t1, ii[n+i], n+i, n, w, *temp);
          ptr = ii[n+i];
          ii[n+i] = *t1;
          *t1 = ptr;
      
          FFT_radix2_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);
   
          ptr = ii[i];
          ii[i] = *t1;
          *t1 = ptr;
          ptr = ii[n+i];
          ii[n+i] = *t2;
          *t2 = ptr;
       }
   } else
   {
       for (i = 0; i < n; i++) 
       {   
          FFT_twiddle(*t1, ii[i], i, 2*n, w/2);
          ptr = ii[i];
          ii[i] = *t1;
          *t1 = ptr;
            
          FFT_twiddle(*t1, ii[n+i], n+i, 2*n, w/2);
          ptr = ii[n+i];
          ii[n+i] = *t1;
          *t1 = ptr;

      
          FFT_radix2_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);
   
          ptr = ii[i];
          ii[i] = *t1;
          *t1 = ptr;
          ptr = ii[n+i];
          ii[n+i] = *t2;
          *t2 = ptr;
       }
   }

   FFT_radix2(rr, 1, ii, n/2, 2*w, t1, t2, temp);
   
   FFT_radix2(rr + n, 1, ii+n, n/2, 2*w, t1, t2, temp);

/*   j = 0; k = 1; l = n+1;
   temp[0] = ii[1];
   ii[1] = ii[n];
   for (i = 2; i < n; i += 2)
   {
      temp[k++] = ii[i];
      ii[i] = temp[j++];
      temp[k++] = ii[i+1];
      ii[i+1] = ii[l++];
   }
   for ( ; i < 2*n; i+=2)
   {
      ii[i] = temp[j++];
      ii[i+1] = ii[l++];
   }*/
}

/*
   FFT of length 2*n on entries of ii with stride is.
   Apply additional twists by z^{c*i} where i starts at r and increases by rs
   for each coeff. Note z => ws bits.
*/
void FFT_radix2_twiddle(mp_limb_t ** ii, mp_size_t is,
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
      mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (n == 1) 
   {
      mp_size_t tw1, tw2;
      tw1 = r*c;
      tw2 = tw1 + rs*c;
      FFT_radix2_twiddle_butterfly(*t1, *t2, ii[0], ii[is], n*w, tw1*ws, tw2*ws);
      ptr = ii[0];
      ii[0] = *t1;
      *t1 = ptr;
      ptr = ii[is];
      ii[is] = *t2;
      *t2 = ptr;
      return;
   }

   // [s0, s1, ..., s{m/2}] = [i0+i{m/2}, i1+i{m/2+1}, ..., i{m/2-1}+i{m-1}]
   // [t0, t1, ..., t{m/2-1}] 
   // = [z1^0*(i0-i{m/2}), z1^1*(i1-i{m/2+1}), ..., z1^{m/2-1}*(i{m/2-1}-i{m-1})]
   // where z1 = exp(2*Pi*I/m), z1 => w bits
   for (i = 0; i < n; i++) 
   {   
      FFT_radix2_butterfly(*t1, *t2, ii[i*is], ii[(n+i)*is], i, n, w);
   
      ptr = ii[i*is];
      ii[i*is] = *t1;
      *t1 = ptr;
      ptr = ii[(n+i)*is];
      ii[(n+i)*is] = *t2;
      *t2 = ptr;
   }

   // [r0, r2, ..., r{m-2}] = Fradix2[s0, s1, ..., s{m/2-1}]
   FFT_radix2_twiddle(ii, is, n/2, 2*w, t1, t2, temp, ws, r, c, 2*rs);
   
   // [r1, r3, ..., r{m-1}] = Fradix2[t0, t1, ..., t{m/2-1}]
   FFT_radix2_twiddle(ii+n*is, is, n/2, 2*w, t1, t2, temp, ws, r + rs, c, 2*rs);
}

void IFFT_radix2(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (n == 1) 
   {
      FFT_radix2_inverse_butterfly(*t1, *t2, ii[0], ii[1], 0, n, w);
      ptr = rr[0];
      rr[0] = *t1;
      *t1 = ptr;
      ptr = rr[rs];
      rr[rs] = *t2;
      *t2 = ptr;
      return;
   }

   // [s0, s1, ..., s{m/2-1}] = Fradix2_inverse[i0, i2, ..., i{m-2}]
   IFFT_radix2(ii, 1, ii, n/2, 2*w, t1, t2, temp);
   
   // [t{m/2}, t{m/2+1}, ..., t{m-1}] = Fradix2_inverse[i1, i3, ..., i{m-1}]
   IFFT_radix2(ii+n, 1, ii+n, n/2, 2*w, t1, t2, temp);

   // [r0, r1, ..., r{m/2-1}] 
   // = [s0+z1^0*t0, s1+z1^1*t1, ..., s{m/2-1}+z1^{m/2-1}*t{m-1}]
   // [r{m/2}, r{m/2+1}, ..., r{m-1}] 
   // = [s0-z1^0*t{m/2}, s1-z1^1*t{m/2+1}, ..., s{m/2-1}-z1^{m/2-1}*t{m/2-1}]
   // where z1 = exp(-2*Pi*I/m), z1 => w bits
   for (i = 0; i < n; i++) 
   {   
      FFT_radix2_inverse_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);
   
      ptr = rr[i];
      rr[i] = *t1;
      *t1 = ptr;
      ptr = rr[n+i];
      rr[n+i] = *t2;
      *t2 = ptr;
   }
}

void IFFT_radix2_sqrt2(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if ((w & 1) == 0)
   {
      IFFT_radix2(rr, rs, ii, 2*n, w/2, t1, t2, temp);

      return;
   }

   // [s0, s1, ..., s{m/2-1}] = Fradix2_inverse[i0, i2, ..., i{m-2}]
   IFFT_radix2(ii, 1, ii, n, w, t1, t2, temp);
   
   // [t{m/2}, t{m/2+1}, ..., t{m-1}] = Fradix2_inverse[i1, i3, ..., i{m-1}]
   IFFT_radix2(ii+2*n, 1, ii+2*n, n, w, t1, t2, temp);

   // [r0, r1, ..., r{m/2-1}] 
   // = [s0+z1^0*t0, s1+z1^1*t1, ..., s{m/2-1}+z1^{m/2-1}*t{m-1}]
   // [r{m/2}, t{m/2+1}, ..., r{m-1}] 
   // = [s0-z1^0*t{m/2}, s1-z1^1*t{m/2+1}, ..., s{m/2-1}-z1^{m/2-1}*t{m-1}]
   // where z1 = exp(-2*Pi*I/m), z1 => w bits
   for (i = 0; i < 2*n; i++) 
   {   
      FFT_radix2_inverse_butterfly(*t1, *t2, ii[i], ii[2*n+i], i/2, n, w);
   
      ptr = rr[i];
      rr[i] = *t1;
      *t1 = ptr;
      ptr = rr[2*n+i];
      rr[2*n+i] = *t2;
      *t2 = ptr;

      i++;

      FFT_radix2_inverse_butterfly_sqrt2(*t1, *t2, ii[i], ii[2*n+i], i, n, w, *temp);
   
      ptr = rr[i];
      rr[i] = *t1;
      *t1 = ptr;
      ptr = rr[2*n+i];
      rr[2*n+i] = *t2;
      *t2 = ptr;
   }
}

void IFFT_radix2_truncate1(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
      mp_size_t trunc)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (trunc == 2*n)
   {
      IFFT_radix2(rr, rs, ii, n, w, t1, t2, temp);
      return;
   }

   if (trunc <= n)
   {
      // PASS
      for (i = trunc; i < n; i++)
      {
         mpn_add_n(ii[i], ii[i], ii[i+n], size);
         mpn_div_2expmod_2expp1(ii[i], ii[i], size - 1, 1);
      }
      
      IFFT_radix2_truncate1(rr, rs, ii, n/2, 2*w, t1, t2, temp, trunc);

      for (i = 0; i < trunc; i++)
         mpn_addsub_n(ii[i], ii[i], ii[i], ii[n+i], size);

      return;
   }

   // [s0, s1, ..., s{m/2-1}] = Fradix2_inverse[i0, i2, ..., i{m-2}]
   IFFT_radix2(ii, 1, ii, n/2, 2*w, t1, t2, temp);

   for (i = trunc - n; i < n; i++)
   {
      mpn_sub_n(ii[i+n], ii[i], ii[i+n], size);
      FFT_twiddle(*t1, ii[i+n], i, n, w);
      mpn_add_n(ii[i], ii[i], ii[i+n], size);
      ptr = ii[i+n];
      ii[i+n] = *t1;
      *t1 = ptr;
   }
   
    // [t{m/2}, t{m/2+1}, ..., t{m-1}] = Fradix2_inverse[i1, i3, ..., i{m-1}]
   IFFT_radix2_truncate1(ii+n, 1, ii+n, n/2, 2*w, t1, t2, temp, trunc - n);

   // [r0, r1, ..., r{m/2-1}] 
   // = [s0+z1^0*t0, s1+z1^1*t1, ..., s{m/2-1}+z1^{m/2-1}*t{m-1}]
   // [r{m/2}, t{m/2+1}, ..., r{m-1}] 
   // = [s0-z1^0*t{m/2}, s1-z1^1*t{m/2+1}, ..., s{m/2-1}-z1^{m/2-1}*t{m-1}]
   // where z1 = exp(-2*Pi*I/m), z1 => w bits
   for (i = 0; i < trunc - n; i++) 
   {   
      FFT_radix2_inverse_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);
   
      ptr = rr[i];
      rr[i] = *t1;
      *t1 = ptr;
      ptr = rr[n+i];
      rr[n+i] = *t2;
      *t2 = ptr;
   }
}

void IFFT_radix2_truncate1_twiddle(mp_limb_t ** ii, mp_size_t is,
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
      mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs, mp_size_t trunc)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (trunc == 2*n)
   {
      IFFT_radix2_twiddle(ii, is, n, w, t1, t2, temp, ws, r, c, rs);
      return;
   }

   if (trunc <= n)
   {
      // PASS
      for (i = trunc; i < n; i++)
      {
         mpn_add_n(ii[i*is], ii[i*is], ii[(i+n)*is], size);
         mpn_div_2expmod_2expp1(ii[i*is], ii[i*is], size - 1, 1);
      }
      
      IFFT_radix2_truncate1_twiddle(ii, is, n/2, 2*w, t1, t2, temp, ws, r, c, 2*rs, trunc);

      for (i = 0; i < trunc; i++)
         mpn_addsub_n(ii[i*is], ii[i*is], ii[i*is], ii[(n+i)*is], size);

      return;
   }

   // [s0, s1, ..., s{m/2-1}] = Fradix2_inverse[i0, i2, ..., i{m-2}]
   IFFT_radix2_twiddle(ii, is, n/2, 2*w, t1, t2, temp, ws, r, c, 2*rs);

   for (i = trunc - n; i < n; i++)
   {
      mpn_sub_n(ii[(i+n)*is], ii[i*is], ii[(i+n)*is], size);
      FFT_twiddle(*t1, ii[(i+n)*is], i, n, w);
      mpn_add_n(ii[i*is], ii[i*is], ii[(i+n)*is], size);
      ptr = ii[(i+n)*is];
      ii[(i+n)*is] = *t1;
      *t1 = ptr;
   }

   // [t{m/2}, t{m/2+1}, ..., t{m-1}] = Fradix2_inverse[i1, i3, ..., i{m-1}]
   IFFT_radix2_truncate1_twiddle(ii + n*is, is, n/2, 2*w, t1, t2, temp, ws, r + rs, c, 2*rs, trunc - n);

   // [r0, r1, ..., r{m/2-1}] 
   // = [s0+z1^0*t0, s1+z1^1*t1, ..., s{m/2-1}+z1^{m/2-1}*t{m-1}]
   // [r{m/2}, t{m/2+1}, ..., r{m-1}] 
   // = [s0-z1^0*t{m/2}, s1-z1^1*t{m/2+1}, ..., s{m/2-1}-z1^{m/2-1}*t{m-1}]
   // where z1 = exp(-2*Pi*I/m), z1 => w bits
   for (i = 0; i < trunc - n; i++) 
   {   
      FFT_radix2_inverse_butterfly(*t1, *t2, ii[i*is], ii[(n+i)*is], i, n, w);
   
      ptr = ii[i*is];
      ii[i*is] = *t1;
      *t1 = ptr;
      ptr = ii[(n+i)*is];
      ii[(n+i)*is] = *t2;
      *t2 = ptr;
   }
}

/* 
   Truncate IFFT to given length. Requires trunc a multiple of 8.
   Assumes (conceptually) zeroes from trunc to n.
*/
void IFFT_radix2_truncate(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
      mp_size_t trunc)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (trunc == 2*n)
   {
      IFFT_radix2(rr, rs, ii, n, w, t1, t2, temp);
      return;
   }

   if (trunc <= n)
   {
      // PASS
      IFFT_radix2_truncate(rr, rs, ii, n/2, 2*w, t1, t2, temp, trunc);

      for (i = 0; i < trunc; i++)
         mpn_add_n(ii[i], ii[i], ii[i], size);

      return;
   }

   //PASS
   // [s0, s1, ..., s{m/2-1}] = Fradix2_inverse[i0, i2, ..., i{m-2}]
   IFFT_radix2(ii, 1, ii, n/2, 2*w, t1, t2, temp);

   for (i = trunc; i < 2*n; i++)
   {
      FFT_twiddle(ii[i], ii[i-n], i - n, n, w);
   }
   
   // [t{m/2}, t{m/2+1}, ..., t{m-1}] = Fradix2_inverse[i1, i3, ..., i{m-1}]
   IFFT_radix2_truncate1(ii+n, 1, ii+n, n/2, 2*w, t1, t2, temp, trunc - n);

   // [r0, r1, ..., r{m/2-1}] 
   // = [s0+z1^0*t0, s1+z1^1*t1, ..., s{m/2-1}+z1^{m/2-1}*t{m-1}]
   // [r{m/2}, t{m/2+1}, ..., r{m-1}] 
   // = [s0-z1^0*t{m/2}, s1-z1^1*t{m/2+1}, ..., s{m/2-1}-z1^{m/2-1}*t{m-1}]
   // where z1 = exp(-2*Pi*I/m), z1 => w bits
   for (i = 0; i < trunc - n; i++) 
   {   
      FFT_radix2_inverse_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);
   
      ptr = rr[i];
      rr[i] = *t1;
      *t1 = ptr;
      ptr = rr[n+i];
      rr[n+i] = *t2;
      *t2 = ptr;
   }

   for (i = trunc - n; i < n; i++)
         mpn_add_n(ii[i], ii[i], ii[i], size);
}

void IFFT_radix2_truncate_twiddle(mp_limb_t ** ii, mp_size_t is,
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
      mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs, mp_size_t trunc)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (trunc == 2*n)
   {
      IFFT_radix2_twiddle(ii, is, n, w, t1, t2, temp, ws, r, c, rs);
      return;
   }

   if (trunc <= n)
   {
      // PASS
      IFFT_radix2_truncate_twiddle(ii, is, n/2, 2*w, t1, t2, temp, ws, r, c, 2*rs, trunc);

      for (i = 0; i < trunc; i++)
         mpn_add_n(ii[i*is], ii[i*is], ii[i*is], size);

      return;
   }

   //PASS
   // [s0, s1, ..., s{m/2-1}] = Fradix2_inverse[i0, i2, ..., i{m-2}]
   IFFT_radix2_twiddle(ii, is, n/2, 2*w, t1, t2, temp, ws, r, c, 2*rs);

   for (i = trunc; i < 2*n; i++)
   {
      FFT_twiddle(ii[i*is], ii[(i-n)*is], i - n, n, w);
   }
   
   // [t{m/2}, t{m/2+1}, ..., t{m-1}] = Fradix2_inverse[i1, i3, ..., i{m-1}]
   IFFT_radix2_truncate1_twiddle(ii+n*is, is, n/2, 2*w, t1, t2, temp, ws, r + rs, c, 2*rs, trunc - n);

   // [r0, r1, ..., r{m/2-1}] 
   // = [s0+z1^0*t0, s1+z1^1*t1, ..., s{m/2-1}+z1^{m/2-1}*t{m-1}]
   // [r{m/2}, t{m/2+1}, ..., r{m-1}] 
   // = [s0-z1^0*t{m/2}, s1-z1^1*t{m/2+1}, ..., s{m/2-1}-z1^{m/2-1}*t{m-1}]
   // where z1 = exp(-2*Pi*I/m), z1 => w bits
   for (i = 0; i < trunc - n; i++) 
   {   
      FFT_radix2_inverse_butterfly(*t1, *t2, ii[i*is], ii[(n+i)*is], i, n, w);
   
      ptr = ii[i*is];
      ii[i*is] = *t1;
      *t1 = ptr;
      ptr = ii[(n+i)*is];
      ii[(n+i)*is] = *t2;
      *t2 = ptr;
   }

   for (i = trunc - n; i < n; i++)
         mpn_add_n(ii[i*is], ii[i*is], ii[i*is], size);
}

void IFFT_radix2_truncate_sqrt2(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
      mp_size_t trunc)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (trunc == 4*n)
   {
      IFFT_radix2_sqrt2(rr, rs, ii, n, w, t1, t2, temp);
      return;
   }

   if ((w & 1) == 0)
   {
      IFFT_radix2_truncate(rr, rs, ii, 2*n, w/2, t1, t2, temp, trunc);
      return;
   }

   //PASS
   // [s0, s1, ..., s{m/2-1}] = Fradix2_inverse[i0, i2, ..., i{m-2}]
   IFFT_radix2(ii, 1, ii, n, w, t1, t2, temp);

   for (i = trunc; i < 4*n; i++)
   {
      FFT_twiddle(ii[i], ii[i - 2*n], i/2 - n, n, w);

      i++;

      FFT_twiddle_sqrt2(ii[i], ii[i - 2*n], i - 2*n, n, w, *temp);
   }
   
   // [t{m/2}, t{m/2+1}, ..., t{m-1}] = Fradix2_inverse[i1, i3, ..., i{m-1}]
   IFFT_radix2_truncate1(ii + 2*n, 1, ii + 2*n, n, w, t1, t2, temp, trunc - 2*n);

   // [r0, r1, ..., r{m/2-1}] 
   // = [s0+z1^0*t0, s1+z1^1*t1, ..., s{m/2-1}+z1^{m/2-1}*t{m-1}]
   // [r{m/2}, t{m/2+1}, ..., r{m-1}] 
   // = [s0-z1^0*t{m/2}, s1-z1^1*t{m/2+1}, ..., s{m/2-1}-z1^{m/2-1}*t{m-1}]
   // where z1 = exp(-2*Pi*I/m), z1 => w bits
   for (i = 0; i < trunc - 2*n; i++) 
   {   
      FFT_radix2_inverse_butterfly(*t1, *t2, ii[i], ii[2*n+i], i/2, n, w);
   
      ptr = rr[i];
      rr[i] = *t1;
      *t1 = ptr;
      ptr = rr[2*n+i];
      rr[2*n+i] = *t2;
      *t2 = ptr;

      i++;

      FFT_radix2_inverse_butterfly_sqrt2(*t1, *t2, ii[i], ii[2*n+i], i, n, w, *temp);
   
      ptr = rr[i];
      rr[i] = *t1;
      *t1 = ptr;
      ptr = rr[2*n+i];
      rr[2*n+i] = *t2;
      *t2 = ptr;
   }

  for (i = trunc - 2*n; i < 2*n; i++)
     mpn_add_n(ii[i], ii[i], ii[i], size);
}

void IFFT_radix2_negacyclic(mp_limb_t ** rr, mp_size_t rs, mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp)
{
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;

   /*j = 0; k = 0;
   for (i = 0; i < n; i += 2)
   {
      ii[i] = ii[2*i];
      temp[k++] = ii[i+1];
      ii[i+1] = ii[2*(i+1)];
   }
   for ( ; i < 2*n; i+=2)
   {
      ii[i] = temp[j++];
      temp[k++] = ii[i+1];
      ii[i+1] = temp[j++];
   }*/

   // [s0, s1, ..., s{m/2-1}] = Fradix2_inverse[i0, i2, ..., i{m-2}]
   IFFT_radix2(ii, 1, ii, n/2, 2*w, t1, t2, temp);
   
   // [t{m/2}, t{m/2+1}, ..., t{m-1}] = Fradix2_inverse[i1, i3, ..., i{m-1}]
   IFFT_radix2(ii+n, 1, ii+n, n/2, 2*w, t1, t2, temp);

   // [r0, r1, ..., r{m/2-1}] 
   // = [s0+z1^0*t0, s1+z1^1*t1, ..., s{m/2-1}+z1^{m/2-1}*t{m-1}]
   // [r{m/2}, t{m/2+1}, ..., r{m-1}] 
   // = [s0-z1^0*t{m/2}, s1-z1^1*t{m/2+1}, ..., s{m/2-1}-z1^{m/2-1}*t{m-1}]
   // where z1 = exp(-2*Pi*I/m), z1 => w bits
   if (w & 1)
   {
      for (i = 0; i < n; i++) 
      {   
          FFT_radix2_inverse_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);
   
          ptr = rr[i];
          rr[i] = *t1;
          *t1 = ptr;
          ptr = rr[n+i];
          rr[n+i] = *t2;
          *t2 = ptr;

          FFT_twiddle(*t1, ii[i], 2*n - i/2, n, w);
          ptr = ii[i];
          ii[i] = *t1;
          *t1 = ptr;
            
          FFT_twiddle(*t1, ii[n+i], 2*n - (n+i)/2, n, w);
          ptr = ii[n+i];
          ii[n+i] = *t1;
          *t1 = ptr;
            
          i++;

          FFT_radix2_inverse_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);

          ptr = rr[i];
          rr[i] = *t1;
          *t1 = ptr;
          ptr = rr[n+i];
          rr[n+i] = *t2;
          *t2 = ptr;

          FFT_twiddle_sqrt2(*t1, ii[i], 4*n-i, n, w, *temp);
          ptr = ii[i];
          ii[i] = *t1;
          *t1 = ptr;
          
          FFT_twiddle_sqrt2(*t1, ii[n+i], 3*n-i, n, w, *temp);
          ptr = ii[n+i];
          ii[n+i] = *t1;
          *t1 = ptr;
      
       }
   } else
   {
       for (i = 0; i < n; i++) 
       {   
          FFT_radix2_inverse_butterfly(*t1, *t2, ii[i], ii[n+i], i, n, w);
   
          ptr = rr[i];
          rr[i] = *t1;
          *t1 = ptr;
          ptr = rr[n+i];
          rr[n+i] = *t2;
          *t2 = ptr;

          FFT_twiddle(*t1, ii[i], 4*n-i, 2*n, w/2);
          ptr = ii[i];
          ii[i] = *t1;
          *t1 = ptr;
            
          FFT_twiddle(*t1, ii[n+i], 3*n-i, 2*n, w/2);
          ptr = ii[n+i];
          ii[n+i] = *t1;
          *t1 = ptr; 
       }
   }
}

void IFFT_radix2_twiddle(mp_limb_t ** ii, mp_size_t is,
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
              mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs)
{
   mp_limb_t ** ss, ** tt;
   mp_limb_t * ptr;
   mp_size_t i, j, k, l;
   mp_size_t size = (w*n)/GMP_LIMB_BITS + 1;
   
   if (n == 1) 
   {
      mp_size_t tw1, tw2;
      tw1 = r*c;
      tw2 = tw1 + rs*c;
      FFT_radix2_twiddle_inverse_butterfly(*t1, *t2, ii[0], ii[is], n*w, tw1*ws, tw2*ws);
      ptr = ii[0];
      ii[0] = *t1;
      *t1 = ptr;
      ptr = ii[is];
      ii[is] = *t2;
      *t2 = ptr;
      return;
   }

   // [s0, s1, ..., s{m/2-1}] = Fradix2_inverse[i0, i2, ..., i{m-2}]
   IFFT_radix2_twiddle(ii, is, n/2, 2*w, t1, t2, temp, ws, r, c, 2*rs);
   
   // [t{m/2}, t{m/2+1}, ..., t{m-1}] = Fradix2_inverse[i1, i3, ..., i{m-1}]
   IFFT_radix2_twiddle(ii+n*is, is, n/2, 2*w, t1, t2, temp, ws, r + rs, c, 2*rs);

   // [r0, r1, ..., r{m/2-1}] 
   // = [s0+z1^0*t0, s1+z1^1*t1, ..., s{m/2-1}+z1^{m/2-1}*t{m-1}]
   // [r{m/2}, t{m/2+1}, ..., r{m-1}] 
   // = [s0-z1^0*t{m/2}, s1-z1^1*t{m/2+1}, ..., s{m/2-1}-z1^{m/2-1}*t{m-1}]
   // where z1 = exp(-2*Pi*I/m), z1 => w bits
   for (i = 0; i < n; i++) 
   {   
      FFT_radix2_inverse_butterfly(*t1, *t2, ii[i*is], ii[(n+i)*is], i, n, w);
   
      ptr = ii[i*is];
      ii[i*is] = *t1;
      *t1 = ptr;
      ptr = ii[(n+i)*is];
      ii[(n+i)*is] = *t2;
      *t2 = ptr;
   }
}

/*
   The matrix Fourier algorithm for a 1D Fourier transform of length m = 2*n,
   works as follows:

   * Split the coefficients into R rows of C columns (here C = n1, R = m/n1 = n2)
   * Perform a length R FFT on each column, i.e. with an input stride of n1
   * Multiply each coefficient by z^{r*c} where z = exp(2*Pi*I/m), note z=>w bits
   * Perform a length C FFT on each row, i.e. with an input stride of 1
*/
void FFT_radix2_mfa(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w, 
                    mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1)
{
   mp_size_t i, j;
   mp_size_t n2 = (2*n)/n1;
   mp_bitcnt_t depth = 0;
   mp_bitcnt_t depth2 = 0;
   mp_limb_t * ptr;

   while ((1UL<<depth) < n2) depth++;
   while ((1UL<<depth2) < n1) depth2++;
   
   // n2 rows, n1 cols

   for (i = 0; i < n1; i++)
   {   
      // FFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      
      FFT_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1);
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
   }

   for (i = 0; i < n2; i++)
   {
      FFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
      
      for (j = 0; j < n1; j++)
      {
         mp_size_t s = mpir_revbin(j, depth2);
         if (j < s)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + s];
            ii[i*n1 + s] = ptr;
         }
      }
   }
}

/*
   The matrix Fourier algorithm for a 1D Fourier transform of length m = 4*n,
   works as follows:

   * Perform 1 row of the FFT
   
   * Perform ordinary MFA on both halves
*/
void FFT_radix2_mfa_sqrt2(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w, 
                    mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1)
{
   mp_size_t i, j;
   mp_size_t n2 = (2*n)/n1;
   mp_bitcnt_t depth = 0;
   mp_bitcnt_t depth2 = 0;
   mp_limb_t * ptr;

   while ((1UL<<depth) < n2) depth++;
   while ((1UL<<depth2) < n1) depth2++;

   /* first half FFT */
   // n2 rows, n1 cols

   for (i = 0; i < n1; i++)
   {   
      /* first row of FFT */
      if ((w & 1) == 1)
      {
         for (j = i; j < i + n1*n2; j+=n1) 
         {   
            if ((j & 1) == 0)
            {
               FFT_radix2_butterfly(*t1, *t2, ii[j], ii[2*n+j], j/2, n, w);
   
               ptr = ii[j];
               ii[j] = *t1;
               *t1 = ptr;
               ptr = ii[2*n+j];
               ii[2*n+j] = *t2;
               *t2 = ptr;
            } else
            {       
               FFT_radix2_butterfly_sqrt2(*t1, *t2, ii[j], ii[2*n+j], j, n, w, *temp);

               ptr = ii[j];
               ii[j] = *t1;
               *t1 = ptr;
               ptr = ii[2*n+j];
               ii[2*n+j] = *t2;
               *t2 = ptr;
            }
         }
      } else
      {
         for (j = i; j < n1*n2; j+=n1) 
         {   
            FFT_radix2_butterfly(*t1, *t2, ii[j], ii[2*n+j], j, 2*n, w/2);
   
            ptr = ii[j];
            ii[j] = *t1;
            *t1 = ptr;
            ptr = ii[2*n+j];
            ii[2*n+j] = *t2;
            *t2 = ptr;
         }
      }
   
      // FFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      
      FFT_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1);
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
   }

   for (i = 0; i < n2; i++)
   {
      FFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
      
      for (j = 0; j < n1; j++)
      {
         mp_size_t s = mpir_revbin(j, depth2);
         if (j < s)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + s];
            ii[i*n1 + s] = ptr;
         }
      }
   }

   ii += 2*n;

   /* second half FFT */
   // n2 rows, n1 cols

   for (i = 0; i < n1; i++)
   {   
      // FFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      
      FFT_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1);
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
   }

   for (i = 0; i < n2; i++)
   {
      FFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
      
      for (j = 0; j < n1; j++)
      {
         mp_size_t s = mpir_revbin(j, depth2);
         if (j < s)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + s];
            ii[i*n1 + s] = ptr;
         }
      }
   }
}

/* 
    trunc must be a multiple of 2*n1
*/
void FFT_radix2_mfa_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w, 
      mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)
{
   mp_size_t i, j;
   mp_size_t n2 = (2*n)/n1;
   mp_size_t trunc2 = (trunc - 2*n)/n1;
   mp_size_t s;
   mp_bitcnt_t depth = 0;
   mp_bitcnt_t depth2 = 0;
   mp_limb_t * ptr;

   while ((1UL<<depth) < n2) depth++;
   while ((1UL<<depth2) < n1) depth2++;

   /* first half FFT */
   // n2 rows, n1 cols
   
   for (i = 0; i < n1; i++)
   {   
      /* first row of FFT */
      if ((w & 1) == 1)
      {
         for (j = i; j < trunc - 2*n; j+=n1) 
         {   
            if ((j & 1) == 0)
            {
               FFT_radix2_butterfly(*t1, *t2, ii[j], ii[2*n+j], j/2, n, w);
   
               ptr = ii[j];
               ii[j] = *t1;
               *t1 = ptr;
               ptr = ii[2*n+j];
               ii[2*n+j] = *t2;
               *t2 = ptr;
            } else
            {       
               FFT_radix2_butterfly_sqrt2(*t1, *t2, ii[j], ii[2*n+j], j, n, w, *temp);

               ptr = ii[j];
               ii[j] = *t1;
               *t1 = ptr;
               ptr = ii[2*n+j];
               ii[2*n+j] = *t2;
               *t2 = ptr;
            }
         }

         for ( ; j < 2*n; j+=n1)
         {
             if ((i & 1) == 0)
                FFT_twiddle(ii[j + 2*n], ii[j], j/2, n, w); 
             else
                FFT_twiddle_sqrt2(ii[j + 2*n], ii[j], j, n, w, *temp); 
         }
      } else
      {
         for (j = i; j < trunc - 2*n; j+=n1) 
         {   
            FFT_radix2_butterfly(*t1, *t2, ii[j], ii[2*n+j], j, 2*n, w/2);
   
            ptr = ii[j];
            ii[j] = *t1;
            *t1 = ptr;
            ptr = ii[2*n+j];
            ii[2*n+j] = *t2;
            *t2 = ptr;
         }

         for ( ; j < 2*n; j+=n1)
            FFT_twiddle(ii[j + 2*n], ii[j], j, 2*n, w/2);
      }
   
      // FFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      
      FFT_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1);
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
   }
   
   for (i = 0; i < n2; i++)
   {
      FFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
      
      for (j = 0; j < n1; j++)
      {
         mp_size_t t = mpir_revbin(j, depth2);
         if (j < t)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + t];
            ii[i*n1 + t] = ptr;
         }
      }
   }
   
   ii += 2*n;

   /* second half FFT */
   // n2 rows, n1 cols

   for (i = 0; i < n1; i++)
   {   
      // FFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      
      FFT_radix2_truncate1_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1, trunc2);
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
   }

   for (s = 0; s < trunc2; s++)
   {
      i = mpir_revbin(s, depth);
      FFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
      
      for (j = 0; j < n1; j++)
      {
         mp_size_t t = mpir_revbin(j, depth2);
         if (j < t)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + t];
            ii[i*n1 + t] = ptr;
         }
      }
   }
}

void FFT_radix2_mfa_truncate(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w, 
        mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)
{
   mp_size_t i, j, s;
   mp_size_t n2 = (2*n)/n1;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_bitcnt_t depth = 0;
   mp_bitcnt_t depth2 = 0;
   mp_limb_t * ptr;

   while ((1UL<<depth) < n2) depth++;
   while ((1UL<<depth2) < n1) depth2++;
   
   trunc /= n1;

   // n2 rows, n1 cols

   for (i = 0; i < n1; i++)
   {   
      // FFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      
      FFT_radix2_truncate_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1, trunc);
      for (j = 0; j < n2; j++)
      {
         s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
   }
   
   for (s = 0; s < trunc; s++)
   {
      i = mpir_revbin(s, depth);
      FFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
      
      for (j = 0; j < n1; j++)
      {
         mp_size_t t = mpir_revbin(j, depth2);
         if (j < t)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + t];
            ii[i*n1 + t] = ptr;
         }
         mpn_normmod_2expp1(ii[i*n1 + j], limbs);
       }
   }
}

void IFFT_radix2_mfa(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w, 
                    mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1)
{
   mp_size_t i, j;
   mp_size_t n2 = (2*n)/n1;
   mp_bitcnt_t depth = 0;
   mp_bitcnt_t depth2 = 0;
   mp_limb_t * ptr;

   while ((1UL<<depth) < n2) depth++;
   while ((1UL<<depth2) < n1) depth2++;

   // n2 rows, n1 cols

   for (i = 0; i < n2; i++)
   {
      for (j = 0; j < n1; j++)
      {
         mp_size_t s = mpir_revbin(j, depth2);
         if (j < s)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + s];
            ii[i*n1 + s] = ptr;
         }
      }      
      
      IFFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
   }

   for (i = 0; i < n1; i++)
   {   
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
      
      // IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      IFFT_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1);
   }

}

void IFFT_radix2_mfa_sqrt2(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w, 
                    mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1)
{
   mp_size_t i, j;
   mp_size_t n2 = (2*n)/n1;
   mp_bitcnt_t depth = 0;
   mp_bitcnt_t depth2 = 0;
   mp_limb_t * ptr;

   while ((1UL<<depth) < n2) depth++;
   while ((1UL<<depth2) < n1) depth2++;

   /* first half IFFT */
   // n2 rows, n1 cols

   for (i = 0; i < n2; i++)
   {
      for (j = 0; j < n1; j++)
      {
         mp_size_t s = mpir_revbin(j, depth2);
         if (j < s)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + s];
            ii[i*n1 + s] = ptr;
         }
      }      
      
      IFFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
   }

   for (i = 0; i < n1; i++)
   {   
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
      
      // IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      IFFT_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1);
   }

   ii += 2*n;

   /* second half IFFT */
   // n2 rows, n1 cols

   for (i = 0; i < n2; i++)
   {
      for (j = 0; j < n1; j++)
      {
         mp_size_t s = mpir_revbin(j, depth2);
         if (j < s)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + s];
            ii[i*n1 + s] = ptr;
         }
      }      
      
      IFFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
   }

   for (i = 0; i < n1; i++)
   {   
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
      
      // IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      IFFT_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1);

      /* final row of IFFT */

      if ((w & 1) == 1)
      {
         for (j = i; j < n1*n2; j+=n1) 
         {   
            if ((j & 1) == 0)
            {
               FFT_radix2_inverse_butterfly(*t1, *t2, ii[j - 2*n], ii[j], j/2, n, w);
   
               ptr = ii[j - 2*n];
               ii[j - 2*n] = *t1;
               *t1 = ptr;
               ptr = ii[j];
               ii[j] = *t2;
               *t2 = ptr;
            } else
            {
               FFT_radix2_inverse_butterfly_sqrt2(*t1, *t2, ii[j - 2*n], ii[j], j, n, w, *temp);
   
               ptr = ii[j - 2*n];
               ii[j - 2*n] = *t1;
               *t1 = ptr;
               ptr = ii[j];
               ii[j] = *t2;
               *t2 = ptr;
            }
         }
      } else
      {
         for (j = i; j < n1*n2; j+=n1) 
         {   
            FFT_radix2_inverse_butterfly(*t1, *t2, ii[j - 2*n], ii[j], j, 2*n, w/2);
   
            ptr = ii[j - 2*n];
            ii[j - 2*n] = *t1;
            *t1 = ptr;
            ptr = ii[j];
            ii[j] = *t2;
            *t2 = ptr;
         }
      }
   }
}

void IFFT_radix2_mfa_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w, 
      mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)
{
   mp_size_t i, j;
   mp_size_t n2 = (2*n)/n1;
   mp_size_t trunc2 = (trunc - 2*n)/n1;
   mp_size_t s;
   mp_bitcnt_t depth = 0;
   mp_bitcnt_t depth2 = 0;
   mp_bitcnt_t size = (w*n)/GMP_LIMB_BITS + 1;
   mp_limb_t * ptr;

   while ((1UL<<depth) < n2) depth++;
   while ((1UL<<depth2) < n1) depth2++;

   /* first half IFFT */
   // n2 rows, n1 cols

   for (i = 0; i < n2; i++)
   {
      for (j = 0; j < n1; j++)
      {
         mp_size_t s = mpir_revbin(j, depth2);
         if (j < s)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + s];
            ii[i*n1 + s] = ptr;
         }
      }      
      
      IFFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
   }
   
   for (i = 0; i < n1; i++)
   {   
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
      
      // IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      IFFT_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1);
   }
   
   ii += 2*n;

   /* second half IFFT */
   // n2 rows, n1 cols

   for (s = 0; s < trunc2; s++)
   {
      i = mpir_revbin(s, depth);
      for (j = 0; j < n1; j++)
      {
         mp_size_t t = mpir_revbin(j, depth2);
         if (j < t)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + t];
            ii[i*n1 + t] = ptr;
         }
      }      
      
      IFFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
   }

   for (i = 0; i < n1; i++)
   {   
      for (j = 0; j < trunc2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }

      for ( ; j < n2; j++)
      {
         mp_size_t u = i + j*n1;
         if ((w & 1) == 1)
         {
            if ((i & 1) == 0)
               FFT_twiddle(ii[i + j*n1], ii[u - 2*n], u/2, n, w); 
            else
               FFT_twiddle_sqrt2(ii[i + j*n1], ii[u - 2*n], u, n, w, *temp); 
         } else
            FFT_twiddle(ii[i + j*n1], ii[u - 2*n], u, 2*n, w/2);
      }

      // IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      IFFT_radix2_truncate1_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1, trunc2);
      
      /* final row of IFFT */

      if ((w & 1) == 1)
      {
         for (j = i; j < trunc - 2*n; j+=n1) 
         {   
            if ((j & 1) == 0)
            {
               FFT_radix2_inverse_butterfly(*t1, *t2, ii[j - 2*n], ii[j], j/2, n, w);
   
               ptr = ii[j - 2*n];
               ii[j - 2*n] = *t1;
               *t1 = ptr;
               ptr = ii[j];
               ii[j] = *t2;
               *t2 = ptr;
            } else
            {
               FFT_radix2_inverse_butterfly_sqrt2(*t1, *t2, ii[j - 2*n], ii[j], j, n, w, *temp);
   
               ptr = ii[j - 2*n];
               ii[j - 2*n] = *t1;
               *t1 = ptr;
               ptr = ii[j];
               ii[j] = *t2;
               *t2 = ptr;
            }
         }
      } else
      {
         for (j = i; j < trunc - 2*n; j+=n1) 
         {   
            FFT_radix2_inverse_butterfly(*t1, *t2, ii[j - 2*n], ii[j], j, 2*n, w/2);
   
            ptr = ii[j - 2*n];
            ii[j - 2*n] = *t1;
            *t1 = ptr;
            ptr = ii[j];
            ii[j] = *t2;
            *t2 = ptr;
         }
      }

      for (j = trunc + i - 2*n; j < 2*n; j+=n1)
           mpn_add_n(ii[j - 2*n], ii[j - 2*n], ii[j - 2*n], size);
   }
}

void IFFT_radix2_mfa_truncate_sqrt2_combined(mp_limb_t ** ii, mp_limb_t ** jj, 
                  mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, 
                mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc, mp_limb_t * tt)
{
   mp_size_t i, j;
   mp_size_t n2 = (2*n)/n1;
   mp_size_t trunc2 = (trunc - 2*n)/n1;
   mp_size_t s;
   mp_bitcnt_t depth = 0;
   mp_bitcnt_t depth2 = 0;
   mp_size_t limbs = (w*n)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_limb_t * ptr;
   int k;

   while ((1UL<<depth) < n2) depth++;
   while ((1UL<<depth2) < n1) depth2++;

   /* first half IFFT */
   // n2 rows, n1 cols

   k = mpn_fft_best_k(limbs, 0);
       
   for (i = 0; i < n2; i++)
   {
      for (j = 0; j < n1; j++)
      {
          mp_limb_t c;
          mpn_normmod_2expp1(jj[i*n1 + j], limbs);
          mpn_normmod_2expp1(ii[i*n1 + j], limbs);
          c = ii[i*n1 + j][limbs] + 2*jj[i*n1 + j][limbs];
          ii[i*n1 + j][limbs] = new_mpn_mulmod_2expp1(ii[i*n1 + j], ii[i*n1 + j], jj[i*n1 + j], c, n*w, tt);
          //ii[i*n1 + j][limbs] = mpn_mul_fft_aux(ii[i*n1 + j], limbs, ii[i*n1 + j], limbs, jj[i*n1 + j], limbs, k, 1);
      }
      
      for (j = 0; j < n1; j++)
      {
         mp_size_t s = mpir_revbin(j, depth2);
 
         if (j < s)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + s];
            ii[i*n1 + s] = ptr;
         }
      }      
      IFFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
   }
   
   for (i = 0; i < n1; i++)
   {   
      for (j = 0; j < n2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
      
      // IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      IFFT_radix2_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1);
   }
   
   ii += 2*n;
   jj += 2*n;

   /* second half IFFT */
   // n2 rows, n1 cols

   for (s = 0; s < trunc2; s++)
   {
      i = mpir_revbin(s, depth);
     
      for (j = 0; j < n1; j++)
      {
          mp_limb_t c;
          mpn_normmod_2expp1(jj[i*n1 + j], limbs);
          mpn_normmod_2expp1(ii[i*n1 + j], limbs);
          c = ii[i*n1 + j][limbs] + 2*jj[i*n1 + j][limbs];
          ii[i*n1 + j][limbs] = new_mpn_mulmod_2expp1(ii[i*n1 + j], ii[i*n1 + j], jj[i*n1 + j], c, n*w, tt);
          //ii[i*n1 + j][limbs] = mpn_mul_fft_aux(ii[i*n1 + j], limbs, ii[i*n1 + j], limbs, jj[i*n1 + j], limbs, k, 1);
      }
      
      for (j = 0; j < n1; j++)
      {
         mp_size_t t = mpir_revbin(j, depth2);
         
         if (j < t)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + t];
            ii[i*n1 + t] = ptr;
         }
      }      
      
      IFFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
   }

   for (i = 0; i < n1; i++)
   {   
      for (j = 0; j < trunc2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }

      for ( ; j < n2; j++)
      {
         mp_size_t u = i + j*n1;
         if ((w & 1) == 1)
         {
            if ((i & 1) == 0)
               FFT_twiddle(ii[i + j*n1], ii[u - 2*n], u/2, n, w); 
            else
               FFT_twiddle_sqrt2(ii[i + j*n1], ii[u - 2*n], u, n, w, *temp); 
         } else
            FFT_twiddle(ii[i + j*n1], ii[u - 2*n], u, 2*n, w/2);
      }

      // IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      IFFT_radix2_truncate1_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1, trunc2);
      
      /* final row of IFFT */

      if ((w & 1) == 1)
      {
         for (j = i; j < trunc - 2*n; j+=n1) 
         {   
            if ((j & 1) == 0)
            {
               FFT_radix2_inverse_butterfly(*t1, *t2, ii[j - 2*n], ii[j], j/2, n, w);
   
               ptr = ii[j - 2*n];
               ii[j - 2*n] = *t1;
               *t1 = ptr;
               ptr = ii[j];
               ii[j] = *t2;
               *t2 = ptr;
            } else
            {
               FFT_radix2_inverse_butterfly_sqrt2(*t1, *t2, ii[j - 2*n], ii[j], j, n, w, *temp);
   
               ptr = ii[j - 2*n];
               ii[j - 2*n] = *t1;
               *t1 = ptr;
               ptr = ii[j];
               ii[j] = *t2;
               *t2 = ptr;
            }
         }
      } else
      {
         for (j = i; j < trunc - 2*n; j+=n1) 
         {   
            FFT_radix2_inverse_butterfly(*t1, *t2, ii[j - 2*n], ii[j], j, 2*n, w/2);
   
            ptr = ii[j - 2*n];
            ii[j - 2*n] = *t1;
            *t1 = ptr;
            ptr = ii[j];
            ii[j] = *t2;
            *t2 = ptr;
         }
      }

      for (j = trunc + i - 2*n; j < 2*n; j+=n1)
           mpn_add_n(ii[j - 2*n], ii[j - 2*n], ii[j - 2*n], size);
   }
}

void IFFT_radix2_mfa_truncate(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w, 
     mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)
{
   mp_size_t i, j, s;
   mp_size_t n2 = (2*n)/n1;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_bitcnt_t depth = 0;
   mp_bitcnt_t depth2 = 0;
   mp_limb_t * ptr;

   while ((1UL<<depth) < n2) depth++;
   while ((1UL<<depth2) < n1) depth2++;

   trunc /= n1;

   // n2 rows, n1 cols

   for (s = 0; s < trunc; s++)
   {
      i = mpir_revbin(s, depth);
      for (j = 0; j < n1; j++)
      {
         mp_size_t t = mpir_revbin(j, depth2);
         if (j < t)
         {
            ptr = ii[i*n1 + j];
            ii[i*n1 + j] = ii[i*n1 + t];
            ii[i*n1 + t] = ptr;
         }
      }
      
      IFFT_radix2(ii + i*n1, 1, ii + i*n1, n1/2, w*n2, t1, t2, temp);
   }

   for (i = 0; i < n1; i++)
   {   
      for (j = 0; j < n2; j++)
      {
         s = mpir_revbin(j, depth);
         if (j < s)
         {
            ptr = ii[i + j*n1];
            ii[i + j*n1] = ii[i + s*n1];
            ii[i + s*n1] = ptr;
         }
      }
      
      // IFFT of length n2 on column i, applying z^{r*i} for rows going up in steps 
      // of 1 starting at row 0, where z => w bits
      IFFT_radix2_truncate_twiddle(ii + i, n1, n2/2, w*n1, t1, t2, temp, w, 0, i, 1, trunc);
      for (j = 0; j < trunc; j++)
         mpn_normmod_2expp1(ii[i + j], limbs);
   }

}

void fft_naive_convolution_1(mp_limb_t * r, mp_limb_t * ii, mp_limb_t * jj, mp_size_t m)
{
   mp_size_t i, j;

   for (i = 0; i < m; i++)
      r[i] = ii[0]*jj[i];

   for (i = 1; i < m; i++)
   {
      for (j = 0; j < m - i; j++)
         r[i+j] += ii[i]*jj[j];

      for ( ; j < m; j++)
         r[i+j-m] -=ii[i]*jj[j];
   }
}

void FFT_mulmod_2expp1(mp_limb_t * r1, mp_limb_t * i1, mp_limb_t * i2, 
                 mp_size_t r_limbs, mp_bitcnt_t depth, mp_bitcnt_t w)
{
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t bits1 = (r_limbs*GMP_LIMB_BITS)/(2*n);
   
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j;

   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *u1, *u2, *s1, *r, *ii0, *jj0;
   mp_limb_t c;
   
   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 4*n + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
   }
   ii0 = ptr;
   t1 = ii0 + 2*n;
   t2 = t1 + size;
   s1 = t2 + size;
   r = s1 + size;
   
   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size) 
   {
      jj[i] = ptr;
   }
   jj0 = ptr;
   u1 = jj0 + 2*n;
   u2 = u1 + size;
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);

   j = FFT_split_bits(ii, i1, r_limbs, bits1, limbs);
   for ( ; j < 2*n; j++)
      MPN_ZERO(ii[j], limbs + 1);

   for (i = 0; i < 2*n; i++)
      ii0[i] = ii[i][0];
 
   FFT_radix2_negacyclic(ii, 1, ii, n, w, &t1, &t2, &s1);
   for (j = 0; j < 2*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);

   j = FFT_split_bits(jj, i2, r_limbs, bits1, limbs);
   for ( ; j < 2*n; j++)
      MPN_ZERO(jj[j], limbs + 1);
   
   for (i = 0; i < 2*n; i++)
      jj0[i] = jj[i][0];
   
   FFT_radix2_negacyclic(jj, 1, jj, n, w, &u1, &u2, &s1);
      
   for (j = 0; j < 2*n; j++)
   {
      mpn_normmod_2expp1(jj[j], limbs);
      c = ii[j][limbs] + 2*jj[j][limbs];
      ii[j][limbs] = mpn_mulmod_2expp1(ii[j], ii[j], jj[j], c, n*w, tt);
   }
   
   IFFT_radix2_negacyclic(ii, 1, ii, n, w, &t1, &t2, &s1);
   
   fft_naive_convolution_1(r, ii0, jj0, 2*n);

   for (j = 0; j < 2*n; j++)
   {
      mp_limb_t t, cy2;
      
      mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 1);
      mpn_normmod_2expp1(ii[j], limbs);

      t = ii[j][limbs];
      ii[j][limbs] = r[j] - ii[j][0];
      cy2 = mpn_add_1(ii[j], ii[j], limbs + 1, ii[j][limbs]);
      add_ssaaaa(r[j], ii[j][limbs], 0, ii[j][limbs], 0, t);
      if (cy2) r[j]++;
   }
   
   MPN_ZERO(r1, r_limbs + 1);
   FFT_combine_bits(r1, ii, 2*n - 1, bits1, limbs + 1, r_limbs + 1);
   
   /* 
      as the negacyclic convolution has effectively done subtractions
      some of the coefficients will be negative, so need to subtract p
   */
   mp_bitcnt_t b = 0;
   mp_size_t ll = 0;
   mp_size_t limb_add = bits1/GMP_LIMB_BITS;
   
   for (j = 0; j < 2*n - 2; j++)
   {   
      if (r[j]) 
         mpn_sub_1(r1 + ll + 1, r1 + ll + 1, r_limbs - ll, 1);
      else if ((mp_limb_signed_t) ii[j][limbs] < 0) // coefficient was -ve
      {
         mpn_sub_1(r1 + ll + 1, r1 + ll + 1, r_limbs - ll, 1);
         mpn_sub_1(r1 + ll + limbs + 1, r1 + ll + limbs + 1, r_limbs - limbs - ll, 1);
      }

      ll += limb_add;
   }
   /* penultimate coefficient, top bit was already ignored */
   if (r[j] || (mp_limb_signed_t) ii[j][limbs] < 0) // coefficient was -ve
      mpn_sub_1(r1 + ll + 1, r1 + ll + 1, r_limbs - ll, 1);
   
   /* final coefficient wraps around */
   r1[r_limbs] += mpn_add_n(r1 + r_limbs - limb_add, r1 + r_limbs - limb_add, ii[2*n - 1], limb_add);
   c = mpn_sub_n(r1, r1, ii[2*n - 1] + limb_add, limbs + 1 - limb_add);
   mpn_addmod_2expp1_1(r1 + limbs + 1 - limb_add, r_limbs - limbs - 1 + limb_add, -c);
   mpn_normmod_2expp1(r1, r_limbs);
   
   TMP_FREE;
}

mp_limb_t new_mpn_mulmod_2expp1(mp_limb_t * r, mp_limb_t * i1, mp_limb_t * i2, 
                           mp_limb_t c, mp_limb_t bits, mp_limb_t * tt)
{
   return mpn_mulmod_2expp1(r, i1, i2, c, bits, tt);
}

mp_limb_t fft_mulmod_2expp1(mp_limb_t * r, mp_limb_t * i1, mp_limb_t * i2, 
                           mp_size_t n, mp_size_t w, mp_limb_t * tt)
{
   mp_size_t bits = n*w;
   mp_size_t limbs = bits/GMP_LIMB_BITS;
   mp_bitcnt_t depth = 1;

   mp_size_t n1, w1;
   mp_bitcnt_t bits1;

   if (limbs < 250) 
   {
      mp_limb_t c = i1[limbs] + 2*i2[limbs];
      return mpn_mulmod_2expp1(r, i1, i2, c, bits, tt);
   }
   
   while ((1UL<<(2*depth)) < bits) depth++;
   depth--;

   n1 = (1UL<<depth);
   w1 = bits/(1UL<<(2*depth));

   bits1 = (n1*w1)/2;
   
   depth -= 3;
   w1 *= 64;

   if (n > (1UL<<15) || (n == (1UL<<15) && w == 2))
   {
      depth++;
      w1/=4;
   }

   if (n > (1UL<<17))
   {
      depth++;
      w1/=4;
   }

   FFT_mulmod_2expp1(r, i1, i2, bits/GMP_LIMB_BITS, depth, w1);

   return 0;
}

/*
   The main integer multiplication routine. Multiplies i1 of n1 limbs by i2 of
   n2 limbs and puts the result in r1, which must have space for n1 + n2 limbs.
   
   Currently there is no sqrt2 trick, thus a convolution of depth d gives a 
   convolution length 2*n where n = 2^d with coefficients of size n*w, where 
   w = w2^2. 

   FIXME: there is no reason for w to be a square, nor for us to pass in w2. 
   Just pass in w.

   Note the *output* polynomial must fit in that convolution size. 
   
   The computation, bits1 = (n*w - depth)/2, gives the maximum size of input 
   coefficients for a given n and w if you want the output to be meaningful.

   The function automatically truncates the FFT, pointwise mults and IFFT, 
   thus the inputs can be different lengths. The function will just segfault
   if n and w2 are not sufficiently large. Except for the smallest multiplications
   (where w2 should be 2), the value of w2 should probably always just be 1.
*/
void new_mpn_mul(mp_limb_t * r1, mp_limb_t * i1, mp_size_t n1, mp_limb_t * i2, mp_size_t n2,
                 mp_bitcnt_t depth, mp_bitcnt_t w)
{
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t bits1 = (n*w - depth)/2; 
   mp_size_t sqrt = (1UL<<(depth/2));

   mp_size_t r_limbs = n1 + n2;
   mp_size_t j1 = (n1*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t j2 = (n2*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t trunc = ((j1 + j2 - 2 + 2*sqrt)/(2*sqrt)) * 2*sqrt;
   
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s, t, u;

   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *u1, *u2, **s1, **s2;
   mp_limb_t c;
   
   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = (mp_limb_t **) t2 + size;
   
   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size) 
   {
      jj[i] = ptr;
   }
   u1 = ptr;
   u2 = u1 + size;
   s2 = (mp_limb_t **) u2 + size;
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   j = FFT_split_bits(ii, i1, n1, bits1, limbs);
   for ( ; j < trunc; j++)
      MPN_ZERO(ii[j], limbs + 1);
   FFT_radix2_mfa_truncate(ii, n, w, &t1, &t2, s1, sqrt, trunc);
            
   j = FFT_split_bits(jj, i2, n2, bits1, limbs);
   for ( ; j < trunc; j++)
      MPN_ZERO(jj[j], limbs + 1);
   FFT_radix2_mfa_truncate(jj, n, w, &u1, &u2, s2, sqrt, trunc);
    
   for (s = 0; s < trunc/sqrt; s++)
   {
      u = mpir_revbin(s, (depth + 1)/2)*sqrt;
      for (t = 0; t < sqrt; t++)
      {
         j = u + t; 
         c = ii[j][limbs] + 2*jj[j][limbs];
         ii[j][limbs] = new_mpn_mulmod_2expp1(ii[j], ii[j], jj[j], c, n*w, tt);
      }
   }

   IFFT_radix2_mfa_truncate(ii, n, w, &t1, &t2, s1, sqrt, trunc);
   for (j = 0; j < trunc; j++)
   {
      mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 1);
      mpn_normmod_2expp1(ii[j], limbs);
   }
   MPN_ZERO(r1, r_limbs);
   FFT_combine_bits(r1, ii, j1 + j2 - 1, bits1, limbs, r_limbs);
      
   TMP_FREE;
}

void new_mpn_mul2(mp_limb_t * r1, mp_limb_t * i1, mp_size_t n1, mp_limb_t * i2, mp_size_t n2,
                 mp_bitcnt_t depth, mp_bitcnt_t w)
{
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t bits1 = (n*w - (depth+1))/2; 
   
   mp_size_t r_limbs = n1 + n2;
   mp_size_t j1 = (n1*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t j2 = (n2*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s, t, u;

   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *u1, *u2, *s1, *s2;
   mp_limb_t c;
   
   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size) + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
   {
      ii[i] = ptr;
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = t2 + size;
   
   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size) + 3*size);
   for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
   {
      jj[i] = ptr;
   }
   u1 = ptr;
   u2 = u1 + size;
   s2 = u2 + size;
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   j1 = FFT_split_bits(ii, i1, n1, bits1, limbs);
   for (j = j1 ; j < 4*n; j++)
      MPN_ZERO(ii[j], limbs + 1);
   FFT_radix2_sqrt2(ii, 1, ii, n, w, &t1, &t2, &s1);
    
   j2 = FFT_split_bits(jj, i2, n2, bits1, limbs);
   for (j = j2 ; j < 4*n; j++)
      MPN_ZERO(jj[j], limbs + 1);
   FFT_radix2_sqrt2(jj, 1, jj, n, w, &u1, &u2, &s2);      

   for (j = 0; j < 4*n; j++)
   {
      mpn_normmod_2expp1(ii[j], limbs);
      mpn_normmod_2expp1(jj[j], limbs);
      c = ii[j][limbs] + 2*jj[j][limbs];
      ii[j][limbs] = new_mpn_mulmod_2expp1(ii[j], ii[j], jj[j], c, n*w, tt);
   }

   IFFT_radix2_sqrt2(ii, 1, ii, n, w, &t1, &t2, &s2);
   for (j = 0; j < 4*n; j++)
   {
      mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
      mpn_normmod_2expp1(ii[j], limbs);
   }
   MPN_ZERO(r1, r_limbs);
   FFT_combine_bits(r1, ii, j1 + j2 - 1, bits1, limbs, r_limbs);
     
   TMP_FREE;
}

void new_mpn_mul3(mp_limb_t * r1, mp_limb_t * i1, mp_size_t n1, mp_limb_t * i2, mp_size_t n2,
                 mp_bitcnt_t depth, mp_bitcnt_t w, mp_size_t sqrt)
{
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t bits1 = (n*w - (depth+1))/2; 
   
   mp_size_t r_limbs = n1 + n2;
   mp_size_t j1 = (n1*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t j2 = (n2*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s, t, u;

   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *s1;
   mp_limb_t c;
   
   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size) + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
   {
      ii[i] = ptr;
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = t2 + size;
   
   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size));
   for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
   {
      jj[i] = ptr;
   }
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   j1 = FFT_split_bits(ii, i1, n1, bits1, limbs);
   for (j = j1 ; j < 4*n; j++)
      MPN_ZERO(ii[j], limbs + 1);
   FFT_radix2_mfa_sqrt2(ii, n, w, &t1, &t2, &s1, sqrt);
    
   j2 = FFT_split_bits(jj, i2, n2, bits1, limbs);
   for (j = j2 ; j < 4*n; j++)
      MPN_ZERO(jj[j], limbs + 1);
   FFT_radix2_mfa_sqrt2(jj, n, w, &t1, &t2, &s1, sqrt);      

   {
      int k;
      for (k = 0; k < 20; k++)
         if (mpn_fft_next_size (limbs, k) == limbs) break;
      for (j = 0; j < 4*n; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         mpn_normmod_2expp1(jj[j], limbs);
         c = ii[j][limbs] + 2*jj[j][limbs];
         //ii[j][limbs] = new_mpn_mulmod_2expp1(ii[j], ii[j], jj[j], c, n*w, tt);
         ii[j][limbs] = mpn_mul_fft_aux (ii[j], limbs, ii[j], limbs, jj[j], limbs, k + 6, 1);
      }
   }

   IFFT_radix2_mfa_sqrt2(ii, n, w, &t1, &t2, &s1, sqrt);
   for (j = 0; j < 4*n; j++)
   {
      mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
      mpn_normmod_2expp1(ii[j], limbs);
   }
   
   MPN_ZERO(r1, r_limbs);
   FFT_combine_bits(r1, ii, j1 + j2 - 1, bits1, limbs, r_limbs);
     
   TMP_FREE;
}

void new_mpn_mul4(mp_limb_t * r1, mp_limb_t * i1, mp_size_t n1, mp_limb_t * i2, mp_size_t n2,
                 mp_bitcnt_t depth, mp_bitcnt_t w)
{
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t bits1 = (n*w - (depth+1))/2; 
   
   mp_size_t r_limbs = n1 + n2;
   mp_size_t j1 = (n1*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t j2 = (n2*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s, t, u, trunc;

   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *s1;
   mp_limb_t c;
   
   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size) + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
   {
      ii[i] = ptr;
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = t2 + size;
   
   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size));
   for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
   {
      jj[i] = ptr;
   }
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   trunc = 2*((j1 + j2)/2); /* trunc must be divisible by 8 */

   j1 = FFT_split_bits(ii, i1, n1, bits1, limbs);
   for (j = j1 ; j < 4*n; j++)
      MPN_ZERO(ii[j], limbs + 1);
   
   FFT_radix2_truncate_sqrt2(ii, 1, ii, n, w, &t1, &t2, &s1, trunc);
    
   j2 = FFT_split_bits(jj, i2, n2, bits1, limbs);
   for (j = j2 ; j < 4*n; j++)
      MPN_ZERO(jj[j], limbs + 1);
   FFT_radix2_truncate_sqrt2(jj, 1, jj, n, w, &t1, &t2, &s1, trunc);      

   {
      int k;
      for (k = 0; k < 20; k++)
         if (mpn_fft_next_size (limbs, k) == limbs) break;
      for (j = 0; j < trunc; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         mpn_normmod_2expp1(jj[j], limbs);
         c = ii[j][limbs] + 2*jj[j][limbs];
         ii[j][limbs] = new_mpn_mulmod_2expp1(ii[j], ii[j], jj[j], c, n*w, tt);
         //ii[j][limbs] = mpn_mul_fft_aux (ii[j], limbs, ii[j], limbs, jj[j], limbs, k + 4, 1);
      }
   }

   IFFT_radix2_truncate_sqrt2(ii, 1, ii, n, w, &t1, &t2, &s1, trunc);
   for (j = 0; j < trunc; j++)
   {
      mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
      mpn_normmod_2expp1(ii[j], limbs);
   }
   
   MPN_ZERO(r1, r_limbs);
   FFT_combine_bits(r1, ii, j1 + j2 - 1, bits1, limbs, r_limbs);
     
   TMP_FREE;
}

void new_mpn_mul5(mp_limb_t * r1, mp_limb_t * i1, mp_size_t n1, mp_limb_t * i2, mp_size_t n2,
                 mp_bitcnt_t depth, mp_bitcnt_t w)
{
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t bits1 = (n*w - depth)/2; 
   
   mp_size_t r_limbs = n1 + n2;
   mp_size_t j1 = (n1*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t j2 = (n2*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s, t, u, trunc;

   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *s1;
   mp_limb_t c;
   
   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = t2 + size;
   
   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size));
   for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size) 
   {
      jj[i] = ptr;
   }
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   trunc = 2*((j1 + j2)/2); /* trunc must be divisible by 8 */

   j1 = FFT_split_bits(ii, i1, n1, bits1, limbs);
   for (j = j1 ; j < 2*n; j++)
      MPN_ZERO(ii[j], limbs + 1);
   
   FFT_radix2_truncate(ii, 1, ii, n, w, &t1, &t2, &s1, trunc);
    
   j2 = FFT_split_bits(jj, i2, n2, bits1, limbs);
   for (j = j2 ; j < 2*n; j++)
      MPN_ZERO(jj[j], limbs + 1);
   FFT_radix2_truncate(jj, 1, jj, n, w, &t1, &t2, &s1, trunc);      

   {
      int k;
      for (k = 0; k < 20; k++)
         if (mpn_fft_next_size (limbs, k) == limbs) break;
      for (j = 0; j < trunc; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         mpn_normmod_2expp1(jj[j], limbs);
         c = ii[j][limbs] + 2*jj[j][limbs];
         ii[j][limbs] = new_mpn_mulmod_2expp1(ii[j], ii[j], jj[j], c, n*w, tt);
         //ii[j][limbs] = mpn_mul_fft_aux (ii[j], limbs, ii[j], limbs, jj[j], limbs, k + 4, 1);
      }
   }

   IFFT_radix2_truncate(ii, 1, ii, n, w, &t1, &t2, &s1, trunc);
   for (j = 0; j < trunc; j++)
   {
      mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 1);
      mpn_normmod_2expp1(ii[j], limbs);
   }
   
   MPN_ZERO(r1, r_limbs);
   FFT_combine_bits(r1, ii, j1 + j2 - 1, bits1, limbs, r_limbs);
     
   TMP_FREE;
}

void new_mpn_mul6(mp_limb_t * r1, mp_limb_t * i1, mp_size_t n1, mp_limb_t * i2, mp_size_t n2,
                 mp_bitcnt_t depth, mp_bitcnt_t w)
{
   mp_size_t n = (1UL<<depth);
   mp_size_t sqrt = (1UL<<(depth/2));
   mp_bitcnt_t bits1 = (n*w - (depth+1))/2; 
   
   mp_size_t r_limbs = n1 + n2;
   mp_size_t j1 = (n1*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t j2 = (n2*GMP_LIMB_BITS - 1)/bits1 + 1;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s, t, u, trunc;

   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *s1;
   mp_limb_t c;
   
   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size) + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
   {
      ii[i] = ptr;
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = t2 + size;
   
   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size));
   for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
   {
      jj[i] = ptr;
   }
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   trunc = 2*sqrt*((j1 + j2 + 2*sqrt - 2)/(2*sqrt)); /* trunc must be divisible by sqrt */

   j1 = FFT_split_bits(ii, i1, n1, bits1, limbs);
   for (j = j1; j < 4*n; j++)
      MPN_ZERO(ii[j], limbs + 1);
   
   FFT_radix2_mfa_truncate_sqrt2(ii, n, w, &t1, &t2, &s1, sqrt, trunc);
    
   j2 = FFT_split_bits(jj, i2, n2, bits1, limbs);
   for (j = j2; j < 4*n; j++)
      MPN_ZERO(jj[j], limbs + 1);
   FFT_radix2_mfa_truncate_sqrt2(jj, n, w, &t1, &t2, &s1, sqrt, trunc);      

   {
      int k = mpn_fft_best_k(limbs, 0);
      mp_size_t trunc2 = (trunc - 2*n)/sqrt;
      mp_size_t depth2 = depth - (depth/2);
      mp_size_t t, u;
      for (j = 0; j < 2*n; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         mpn_normmod_2expp1(jj[j], limbs);
         //c = ii[j][limbs] + 2*jj[j][limbs];
         //ii[j][limbs] = new_mpn_mulmod_2expp1(ii[j], ii[j], jj[j], c, n*w, tt);
         //ii[j][limbs] = mpn_mul_fft_aux(ii[j], limbs, ii[j], limbs, jj[j], limbs, k, 1);
         fft_mulmod_2expp1(ii[j], ii[j], jj[j], n, w, tt);
      }
      for (j = 0; j < trunc2; j++)
      {
         mp_size_t s = mpir_revbin(j, depth2 +1);
         for (t = 0; t < sqrt; t++)
         {
            u = 2*n + s*sqrt+t;
            mpn_normmod_2expp1(ii[u], limbs);
            mpn_normmod_2expp1(jj[u], limbs);
            //c = ii[s][limbs] + 2*jj[s][limbs];
            //ii[s][limbs] = new_mpn_mulmod_2expp1(ii[s], ii[s], jj[s], c, n*w, tt);
            //ii[u][limbs] = mpn_mul_fft_aux(ii[u], limbs, ii[u], limbs, jj[u], limbs, k, 1);
            fft_mulmod_2expp1(ii[u], ii[u], jj[u], n, w, tt);
         }
      }
   }
   IFFT_radix2_mfa_truncate_sqrt2(ii, n, w, &t1, &t2, &s1, sqrt, trunc);
   
   //IFFT_radix2_mfa_truncate_sqrt2_combined(ii, jj, n, w, &t1, &t2, &s1, sqrt, trunc, tt);
   for (j = 0; j < trunc; j++)
   {
      mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
      mpn_normmod_2expp1(ii[j], limbs);
   }
   
   MPN_ZERO(r1, r_limbs);
   FFT_combine_bits(r1, ii, j1 + j2 - 1, bits1, limbs, r_limbs);
     
   TMP_FREE;
}


/************************************************************************************

   Test code

************************************************************************************/

void mpn_to_mpz(mpz_t m, mp_limb_t * i, mp_size_t limbs)
{
   mp_limb_signed_t hi;
   
   mpz_realloc(m, limbs + 1);
   MPN_COPY(m->_mp_d, i, limbs + 1);
   hi = i[limbs];
   if (hi < 0L)
   {
      mpn_neg_n(m->_mp_d, m->_mp_d, limbs + 1);
      m->_mp_size = limbs + 1;
      while ((m->_mp_size) && (!m->_mp_d[m->_mp_size - 1])) 
         m->_mp_size--;
      m->_mp_size = -m->_mp_size;
   } else
   {
      m->_mp_size = limbs + 1;
      while ((m->_mp_size) && (!m->_mp_d[m->_mp_size - 1])) 
         m->_mp_size--;
   }
}

void ref_norm(mpz_t m, mpz_t p)
{
   mpz_mod(m, m, p);
}

void ref_submod_i(mpz_t m, mpz_t i1, mpz_t i2, mpz_t p, mp_size_t n, mp_bitcnt_t w)
{
   mpz_sub(m, i1, i2);
   mpz_mul_2exp(m, m, (n*w)/2);
   mpz_mod(m, m, p);
}

void ref_mul_2expmod(mpz_t m, mpz_t i2, mpz_t p, mp_size_t n, mp_bitcnt_t w, mp_bitcnt_t d)
{
   mpz_mul_2exp(m, i2, d);
   mpz_mod(m, m, p);
}

void ref_div_2expmod(mpz_t m, mpz_t i2, mpz_t p, mp_size_t n, mp_bitcnt_t w, mp_bitcnt_t d)
{
   mpz_t temp;
   mpz_init(temp);
   mpz_set_ui(temp, 1);
   mpz_mul_2exp(temp, temp, d);
   mpz_invert(temp, temp, p);
   mpz_mul(m, i2, temp);
   mpz_mod(m, m, p);
   mpz_clear(temp);
}

void ref_lshB_sumdiffmod(mpz_t t, mpz_t u, mpz_t i1, 
                      mpz_t i2, mpz_t p, mp_size_t n, mp_bitcnt_t w, mp_bitcnt_t x, mp_bitcnt_t y)
{
   mpz_add(t, i1, i2);
   mpz_sub(u, i1, i2);
   mpz_mul_2exp(t, t, x*GMP_LIMB_BITS);
   mpz_mul_2exp(u, u, y*GMP_LIMB_BITS);
   mpz_mod(t, t, p);
   mpz_mod(u, u, p);
}

void ref_sumdiff_rshBmod(mpz_t t, mpz_t u, mpz_t i1, 
                      mpz_t i2, mpz_t p, mp_size_t n, mp_bitcnt_t w, mp_bitcnt_t x, mp_bitcnt_t y)
{
   mpz_t mult1, mult2;
   mpz_init(mult1);
   mpz_init(mult2);
   mpz_set_ui(mult1, 1);
   mpz_mul_2exp(mult1, mult1, x*GMP_LIMB_BITS);
   mpz_invert(mult1, mult1, p);
   mpz_set_ui(mult2, 1);
   mpz_mul_2exp(mult2, mult2, y*GMP_LIMB_BITS);
   mpz_invert(mult2, mult2, p);
   mpz_mul(mult1, mult1, i1);
   mpz_mul(mult2, mult2, i2);
   mpz_add(t, mult1, mult2);
   mpz_sub(u, mult1, mult2);
   mpz_mod(t, t, p);
   mpz_mod(u, u, p);
   mpz_clear(mult1);
   mpz_clear(mult2);
}

/* set p = 2^wn + 1 */
void set_p(mpz_t p, mp_size_t n, mp_bitcnt_t w)
{
   mpz_set_ui(p, 1);
   mpz_mul_2exp(p, p, n*w);
   mpz_add_ui(p, p, 1);
}

void rand_n(mp_limb_t * n, gmp_randstate_t state, mp_size_t limbs)
{
   mpn_rrandom(n, state, limbs);
   n[limbs] = gmp_urandomm_ui(state, 10);
   if (gmp_urandomm_ui(state, 2)) n[limbs] = -n[limbs];
}

void test_norm()
{
   mp_size_t i, j, k, l, n, w, limbs;
   mpz_t p, m, m2;
   mpz_init(p);
   mpz_init(m);
   mpz_init(m2);
   mp_limb_t * nn;
   gmp_randstate_t state;
   gmp_randinit_default(state);
   TMP_DECL;

   for (i = GMP_LIMB_BITS; i < 32*GMP_LIMB_BITS; i += GMP_LIMB_BITS)
   {
      for (j = 1; j < 32; j++)
      {
         for (k = 1; k <= GMP_NUMB_BITS; k <<= 1)
         {
            n = i/k;
            w = j*k;
            limbs = (n*w)/GMP_LIMB_BITS;
            TMP_MARK;
            nn = TMP_BALLOC_LIMBS(limbs + 1);
            mpn_rrandom(nn, state, limbs + 1);
            mpn_to_mpz(m, nn, limbs);
            set_p(p, n, w);
            
            mpn_normmod_2expp1(nn, limbs);
            mpn_to_mpz(m2, nn, limbs);
            ref_norm(m, p);

            if (mpz_cmp(m, m2) != 0)
            {
               printf("mpn_normmod_2expp1 error\n");
               gmp_printf("want %Zx\n\n", m);
               gmp_printf("got  %Zx\n", m2);
               abort();
            }
            TMP_FREE;
         }
      }
   }
   mpz_clear(p);
   mpz_clear(m);
   mpz_clear(m2);
   gmp_randclear(state);
}

void test_mul_2expmod()
{
   mp_size_t i, j, k, l, n, w, limbs, d;
   mpz_t p, m, m2, mn1, mn2;
   mpz_init(p);
   mpz_init(m);
   mpz_init(m2);
   mpz_init(mn1);
   mp_limb_t * nn1, * r;
   gmp_randstate_t state;
   gmp_randinit_default(state);
   TMP_DECL;

   for (i = 2*GMP_LIMB_BITS; i < 64*GMP_LIMB_BITS; i += 2*GMP_LIMB_BITS)
   {
      for (j = 1; j < 32; j++)
      {
         for (k = 1; k <= 2*GMP_NUMB_BITS; k <<= 1)
         {
            for (d = 0; d < GMP_LIMB_BITS; d++)
            {
               n = i/k;
               w = j*k;
               limbs = (n*w)/GMP_LIMB_BITS;
               TMP_MARK;
               nn1 = TMP_BALLOC_LIMBS(limbs + 1);
               r = TMP_BALLOC_LIMBS(limbs + 1);
               rand_n(nn1, state, limbs);
               mpn_to_mpz(mn1, nn1, limbs);
               set_p(p, n, w);
               
               mpn_mul_2expmod_2expp1(r, nn1, limbs, d);
               mpn_to_mpz(m2, r, limbs);
               ref_norm(m2, p);
               ref_mul_2expmod(m, mn1, p, n, w, d);
               
               if (mpz_cmp(m, m2) != 0)
               {
                  printf("mpn_mul_2expmod_2expp1 error\n");
                  gmp_printf("want %Zx\n\n", m);
                  gmp_printf("got  %Zx\n", m2);
                  abort();
               }
               TMP_FREE;
            }
         }
      }
   }
   mpz_clear(p);
   mpz_clear(m);
   mpz_clear(m2);
   mpz_clear(mn1);
   gmp_randclear(state);
}

void test_FFT_negacyclic_twiddle()
{
   mp_size_t i, j, k, l, n, w, limbs, d;
   mpz_t p, m, m2, mn1, mn2;
   mpz_init(p);
   mpz_init(m);
   mpz_init(m2);
   mpz_init(mn1);
   mp_limb_t * nn1, * r;
   gmp_randstate_t state;
   gmp_randinit_default(state);
   TMP_DECL;

   for (i = 2*GMP_LIMB_BITS; i < 20*GMP_LIMB_BITS; i += 2*GMP_LIMB_BITS)
   {
      for (j = 1; j < 10; j++)
      {
         for (k = 1; k <= 2*GMP_NUMB_BITS; k <<= 1)
         {
            n = i/k;
            w = 2*j*k;
            for (d = 0; d < 2*n; d++)
            {
               limbs = (n*w)/GMP_LIMB_BITS;
               TMP_MARK;
               nn1 = TMP_BALLOC_LIMBS(limbs + 1);
               r = TMP_BALLOC_LIMBS(limbs + 1);
               rand_n(nn1, state, limbs);
               mpn_to_mpz(mn1, nn1, limbs);
               set_p(p, n, w);
               
               if (!FFT_negacyclic_twiddle(r, nn1, d, n, w))
                  MPN_COPY(r, nn1, limbs + 1);
               mpn_to_mpz(m2, r, limbs);
               ref_norm(m2, p);
               ref_mul_2expmod(m, mn1, p, n, w, d*w/2);
               
               if (mpz_cmp(m, m2) != 0)
               {
                  printf("FFT_negacyclic_twiddle error\n");
                  gmp_printf("want %Zx\n\n", m);
                  gmp_printf("got  %Zx\n", m2);
                  abort();
               }
               TMP_FREE;
            }
         }
      }
   }

   for (i = 2*GMP_LIMB_BITS; i < 20*GMP_LIMB_BITS; i += 2*GMP_LIMB_BITS)
   {
      for (j = 1; j < 10; j++)
      {
         for (k = 1; k <= 2*GMP_NUMB_BITS; k <<= 1)
         {
            n = i/k;
            w = 2*j*k;
            for (d = 0; d < 2*n; d++)
            {
               limbs = (n*w)/GMP_LIMB_BITS;
               TMP_MARK;
               nn1 = TMP_BALLOC_LIMBS(limbs + 1);
               r = TMP_BALLOC_LIMBS(limbs + 1);
               rand_n(nn1, state, limbs);
               mpn_to_mpz(mn1, nn1, limbs);
               set_p(p, n, w);
               
               if (!FFT_negacyclic_twiddle(r, nn1, 4*n - d, n, w))
                  MPN_COPY(r, nn1, limbs + 1);
               mpn_to_mpz(m2, r, limbs);
               ref_norm(m2, p);
               ref_div_2expmod(m, mn1, p, n, w, d*w/2);
               
               if (mpz_cmp(m, m2) != 0)
               {
                  printf("FFT_negacyclic_twiddle error\n");
                  gmp_printf("want %Zx\n\n", m);
                  gmp_printf("got  %Zx\n", m2);
                  abort();
               }
               TMP_FREE;
            }
         }
      }
   }
   mpz_clear(p);
   mpz_clear(m);
   mpz_clear(m2);
   mpz_clear(mn1);
   gmp_randclear(state);
}

void test_div_2expmod()
{
   mp_size_t i, j, k, l, n, w, limbs, d;
   mpz_t p, m, m2, mn1, mn2;
   mpz_init(p);
   mpz_init(m);
   mpz_init(m2);
   mpz_init(mn1);
   mp_limb_t * nn1, * r;
   gmp_randstate_t state;
   gmp_randinit_default(state);
   TMP_DECL;

   for (i = 2*GMP_LIMB_BITS; i < 64*GMP_LIMB_BITS; i += 2*GMP_LIMB_BITS)
   {
      for (j = 1; j < 32; j++)
      {
         for (k = 1; k <= 2*GMP_NUMB_BITS; k <<= 1)
         {
            for (d = 0; d < GMP_LIMB_BITS; d++)
            {
               n = i/k;
               w = j*k;
               limbs = (n*w)/GMP_LIMB_BITS;
               TMP_MARK;
               nn1 = TMP_BALLOC_LIMBS(limbs + 1);
               r = TMP_BALLOC_LIMBS(limbs + 1);
               rand_n(nn1, state, limbs);
            
               mpn_to_mpz(mn1, nn1, limbs);
               set_p(p, n, w);
            
               mpn_div_2expmod_2expp1(r, nn1, limbs, d);
               mpn_to_mpz(m2, r, limbs);
               ref_norm(m2, p);
               ref_norm(mn1, p);
               ref_mul_2expmod(m, m2, p, n, w, d);

               if (mpz_cmp(m, mn1) != 0)
               {
                  printf("mpn_div_2expmod_2expp1 error\n");
                  gmp_printf("want %Zx\n\n", mn1);
                  gmp_printf("got  %Zx\n", m);
                  abort();
               }
               TMP_FREE;
            }
         }
      }
   }
   mpz_clear(p);
   mpz_clear(m);
   mpz_clear(m2);
   mpz_clear(mn1);
   gmp_randclear(state);
}

void test_lshB_sumdiffmod()
{
   mp_size_t c, i, j, k, l, x, y, n, w, limbs;
   mpz_t p, ma, mb, m2a, m2b, mn1, mn2;
   mpz_init(p);
   mpz_init(ma);
   mpz_init(mb);
   mpz_init(m2a);
   mpz_init(m2b);
   mpz_init(mn1);
   mpz_init(mn2);
   mp_limb_t * nn1, * nn2, * r1, * r2;
   gmp_randstate_t state;
   gmp_randinit_default(state);
   TMP_DECL;

   for (i = 2*GMP_LIMB_BITS; i < 20*GMP_LIMB_BITS; i += 2*GMP_LIMB_BITS)
   {
      for (j = 1; j < 10; j++)
      {
         for (k = 1; k <= 2*GMP_NUMB_BITS; k <<= 1)
         {
            n = i/k;
            w = j*k;
            limbs = (n*w)/GMP_LIMB_BITS;
            for (c = 0; c < limbs; c++)
            {
               x = gmp_urandomm_ui(state, limbs + 1);
               y = gmp_urandomm_ui(state, limbs + 1);
               TMP_MARK;
               nn1 = TMP_BALLOC_LIMBS(limbs + 1);
               nn2 = TMP_BALLOC_LIMBS(limbs + 1);
               r1 = TMP_BALLOC_LIMBS(limbs + 1);
               r2 = TMP_BALLOC_LIMBS(limbs + 1);
               rand_n(nn1, state, limbs);
               rand_n(nn2, state, limbs);
            
               mpn_to_mpz(mn1, nn1, limbs);
               mpn_to_mpz(mn2, nn2, limbs);
               set_p(p, n, w);
            
               mpn_lshB_sumdiffmod_2expp1(r1, r2, nn1, nn2, limbs, x, y);
               mpn_to_mpz(m2a, r1, limbs);
               mpn_to_mpz(m2b, r2, limbs);
               ref_norm(m2a, p);
               ref_norm(m2b, p);
               ref_lshB_sumdiffmod(ma, mb, mn1, mn2, p, n, w, x, y);

               if (mpz_cmp(ma, m2a) != 0)
               {
                  printf("mpn_lshB_sumdiffmod_2expp1 error a\n");
                  printf("x = %ld, y = %ld\n", x, y);
                  gmp_printf("want %Zx\n\n", ma);
                  gmp_printf("got  %Zx\n", m2a);
                  abort();
               }
               if (mpz_cmp(mb, m2b) != 0)
               {
                  printf("mpn_lshB_sumdiffmod_2expp1 error b\n");
                  printf("x = %ld, y = %ld\n", x, y);
                  gmp_printf("want %Zx\n\n", mb);
                  gmp_printf("got  %Zx\n", m2b);
                  abort();
               }
               TMP_FREE;
            }
         }
      }
   }
   mpz_clear(p);
   mpz_clear(ma);
   mpz_clear(mb);
   mpz_clear(m2a);
   mpz_clear(m2b);
   mpz_clear(mn1);
   mpz_clear(mn2);
   gmp_randclear(state);
}

void test_sumdiff_rshBmod()
{
   mp_size_t c, i, j, k, l, x, y, n, w, limbs;
   mpz_t p, ma, mb, m2a, m2b, mn1, mn2;
   mpz_init(p);
   mpz_init(ma);
   mpz_init(mb);
   mpz_init(m2a);
   mpz_init(m2b);
   mpz_init(mn1);
   mpz_init(mn2);
   mp_limb_t * nn1, * nn2, * r1, * r2;
   gmp_randstate_t state;
   gmp_randinit_default(state);
   TMP_DECL;

   for (i = 2*GMP_LIMB_BITS; i < 20*GMP_LIMB_BITS; i += 2*GMP_LIMB_BITS)
   {
      for (j = 1; j < 10; j++)
      {
         for (k = 1; k <= 2*GMP_NUMB_BITS; k <<= 1)
         {
            n = i/k;
            w = j*k;
            limbs = (n*w)/GMP_LIMB_BITS;
            for (c = 0; c < limbs; c++)
            {
               x = gmp_urandomm_ui(state, limbs);
               y = gmp_urandomm_ui(state, limbs);
               TMP_MARK;
               nn1 = TMP_BALLOC_LIMBS(limbs + 1);
               nn2 = TMP_BALLOC_LIMBS(limbs + 1);
               r1 = TMP_BALLOC_LIMBS(limbs + 1);
               r2 = TMP_BALLOC_LIMBS(limbs + 1);
               rand_n(nn1, state, limbs);
               rand_n(nn2, state, limbs);
            
               mpn_to_mpz(mn1, nn1, limbs);
               mpn_to_mpz(mn2, nn2, limbs);
               set_p(p, n, w);
            
               mpn_sumdiff_rshBmod_2expp1(r1, r2, nn1, nn2, limbs, x, y);
               mpn_to_mpz(m2a, r1, limbs);
               mpn_to_mpz(m2b, r2, limbs);
               ref_norm(m2a, p);
               ref_norm(m2b, p);
               ref_sumdiff_rshBmod(ma, mb, mn1, mn2, p, n, w, x, y);

               if (mpz_cmp(ma, m2a) != 0)
               {
                  printf("mpn_sumdiff_rshBmod_2expp1 error a\n");
                  printf("x = %ld, y = %ld, limbs = %ld\n", x, y, limbs);
                  gmp_printf("want %Zx\n\n", ma);
                  gmp_printf("got  %Zx\n", m2a);
                  abort();
               }
               if (mpz_cmp(mb, m2b) != 0)
               {
                  printf("mpn_sumdiff_rshBmod_2expp1 error b\n");
                  printf("x = %ld, y = %ld, limbs = %ld\n", x, y, limbs);
                  gmp_printf("want %Zx\n\n", mb);
                  gmp_printf("got  %Zx\n", m2b);
                  abort();
               }
               TMP_FREE;
            }
         }
      }
   }
   mpz_clear(p);
   mpz_clear(ma);
   mpz_clear(mb);
   mpz_clear(m2a);
   mpz_clear(m2b);
   mpz_clear(mn1);
   mpz_clear(mn2);
   gmp_randclear(state);
}

void time_mul_with_negacyclic()
{
   mp_bitcnt_t depth = 17UL;
   mp_bitcnt_t w = 1UL;
   mp_size_t iters = 1;

   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t bits1 = (n*w - depth)/2; 
   mp_bitcnt_t bits = n*bits1;
   mp_size_t int_limbs = (bits - 1UL)/GMP_LIMB_BITS + 1;
   
   mp_size_t i, j;
   mp_limb_t *i1, *i2, *r1;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   i1 = TMP_BALLOC_LIMBS(4*int_limbs);
   i2 = i1 + int_limbs;
   r1 = i2 + int_limbs;
   
   mpn_urandomb(i1, state, bits);
   mpn_urandomb(i2, state, bits);
  
   for (i = 0; i < iters; i++)
   {
      new_mpn_mul(r1, i1, int_limbs, i2, int_limbs, depth, w);
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void test_mulmod()
{
   mp_bitcnt_t depth = 15UL; //should be at least 6 
   mp_bitcnt_t w = 1UL; /* should be 1 or a power of 2 */
   mp_size_t iters = 10000;

   mp_size_t n = (1UL<<depth);
   
   mp_bitcnt_t bits = n*w;
   mp_size_t int_limbs = bits/GMP_LIMB_BITS;
   
   mp_size_t i, j;
   mp_limb_t *i1, *i2, *r1, *r2, *tt;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   i1 = TMP_BALLOC_LIMBS(6*(int_limbs+1));
   i2 = i1 + int_limbs + 1;
   r1 = i2 + int_limbs + 1;
   r2 = r1 + int_limbs + 1;
   tt = r2 + int_limbs + 1;
   
   for (i = 0; i < iters; i++)
   {
      mpn_rrandom(i1, state, int_limbs);
      i1[int_limbs] = CNST_LIMB(0);
      mpn_rrandom(i2, state, int_limbs);
      i2[int_limbs] = CNST_LIMB(0);
      mpn_mulmod_2expp1(r2, i1, i2, 0, bits, tt);
      fft_mulmod_2expp1(r1, i1, i2, n, w, tt);
      
      mp_size_t wrong = 0;
      for (j = 0; j < int_limbs; j++)
      {
         if (r1[j] != r2[j]) 
         {
            if (wrong < 10) 
               printf("error in limb %ld, %lx != %lx\n", j, r1[j], r2[j]);
            wrong++;
         } 
      }
      if (wrong) printf("%ld limbs wrong\n", wrong);
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void test_fft_ifft()
{
   mp_bitcnt_t depth = 10UL;
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t w = 1;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, ** kk, *tt, *t1, *t2, *u1, *u2, *s1, *s2;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      rand_n(ii[i], state, limbs);
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = t2 + size;

   for (j = 0; j < 2*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);

   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size));
   for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size) 
   {
      jj[i] = ptr;
      MPN_COPY(jj[i], ii[i], limbs + 1);
   }
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (i = 0; i < 1; i++)
   {
      FFT_radix2(ii, 1, ii, n, w, &t1, &t2, &s1);
      
      IFFT_radix2(ii, 1, ii, n, w, &t1, &t2, &s1);
      for (j = 0; j < 2*n; j++)
      {
         mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 1);
         mpn_normmod_2expp1(ii[j], limbs);
      }

      for (j = 0; j < 2*n; j++)
      {
         if (mpn_cmp(ii[j], jj[j], limbs + 1) != 0)
         {
            printf("Error in entry %ld\n", j);
            abort();
         }
      }
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void test_fft_ifft_negacyclic()
{
   mp_bitcnt_t depth = 11UL;
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t w = 1;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, ** kk, *tt, *t1, *t2, *u1, *u2, *s1, *s2;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      rand_n(ii[i], state, limbs);
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = t2 + size;

   for (j = 0; j < 2*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);

   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size));
   for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size) 
   {
      jj[i] = ptr;
      MPN_COPY(jj[i], ii[i], limbs + 1);
   }
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (i = 0; i < 1; i++)
   {
      FFT_radix2_negacyclic(ii, 1, ii, n, w, &t1, &t2, &s1);
      
      IFFT_radix2_negacyclic(ii, 1, ii, n, w, &t1, &t2, &s1);
      for (j = 0; j < 2*n; j++)
      {
         mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth+1);
         mpn_normmod_2expp1(ii[j], limbs);
      }

      for (j = 0; j < 2*n; j++)
      {
         if (mpn_cmp(ii[j], jj[j], limbs + 1) != 0)
         {
            printf("Error in entry %ld\n", j);
            abort();
         }
      }
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void test_fft_ifft_sqrt2()
{
   mp_bitcnt_t depth = 6UL;
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t w = 1;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *t1, *t2, *s1;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size) + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      rand_n(ii[i], state, limbs);
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = t2 + size;

   for (j = 0; j < 4*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);
  
   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size));
   for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
   {
      jj[i] = ptr;
      MPN_COPY(jj[i], ii[i], size);
   }
   
   FFT_radix2_sqrt2(ii, 1, ii, n, w, &t1, &t2, &s1);
   
   for (j = 0; j < 4*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);

   IFFT_radix2_sqrt2(ii, 1, ii, n, w, &t1, &t2, &s1);
   
   for (j = 0; j < 4*n; j++)
   {
      mpn_mul_2expmod_2expp1(jj[j], jj[j], limbs, depth + 2);
      mpn_normmod_2expp1(jj[j], limbs);
      mpn_normmod_2expp1(ii[j], limbs);
   }

   for (j = 0; j < 4*n; j++)
   {
      if (mpn_cmp(ii[j], jj[j], size) != 0)
      {
         printf("Error in entry %ld\n", j);
         gmp_printf("%Nx != \n%Nx\n", ii[j], size, jj[j], size);
         abort();
      }
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void test_fft_ifft_truncate()
{
   mp_bitcnt_t depth = 10UL;
   mp_size_t n = (1UL<<depth);            
   mp_bitcnt_t w = 1;
   mp_size_t iter = 1000;

   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s, count;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, ** kk, *tt, *t1, *t2, *u1, *u2, *v1, *v2, **s1, **s2, **s3;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;
 
   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*n + 2*n*size + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      rand_n(ii[i], state, limbs);
   }
   t1 = (mp_limb_t *) ii + 2*n + 2*n*size;
   t2 = t1 + size;
   s1 = (mp_limb_t **) t2 + size;

   for (j = 0; j < 2*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);

   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(2*n + 2*n*size + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size) 
   {
      jj[i] = ptr;
      MPN_COPY(jj[i], ii[i], limbs + 1);
   }
   u1 = (mp_limb_t *) jj + 2*n + 2*n*size;
   u2 = u1 + size;
   s2 = (mp_limb_t **) u2 + size;
   
   kk = (mp_limb_t **) TMP_BALLOC_LIMBS(2*n + 2*n*size + 2*n + 2*size);
   for (i = 0, j = 0, ptr = (mp_limb_t *) kk + 2*n; i < 2*n; i++, ptr += size) 
   {
      kk[i] = ptr;
   }
   v1 = (mp_limb_t *) kk + 2*n + 2*n*size;
   v2 = v1 + size;
   s3 = (mp_limb_t **) v2 + size;
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (count = 0; count < iter; count++)
   {
      for (i = 0; i < 2*n; i++) 
      {
         rand_n(ii[i], state, limbs);
      }
   
      for (j = 0; j < 2*n; j++)
         mpn_normmod_2expp1(ii[j], limbs);

      for (i = 0; i < 2*n; i++) 
      {
         MPN_COPY(jj[i], ii[i], limbs + 1);
      }
      
      mp_size_t trunc = gmp_urandomm_ui(state, 2*n) + 1;
      trunc = ((trunc + 1)/2)*2;
   
      FFT_radix2_truncate(ii, 1, ii, n, w, &t1, &t2, s1, trunc);
      for (j = 0; j < trunc; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         MPN_COPY(kk[j], ii[j], limbs + 1);
      }
      
      IFFT_radix2_truncate(kk, 1, kk, n, w, &v1, &v2, s3, trunc);
      for (j = 0; j < trunc; j++)
      {
         mpn_mul_2expmod_2expp1(jj[j], jj[j], limbs, depth + 1);
         mpn_normmod_2expp1(jj[j], limbs);
         mpn_normmod_2expp1(kk[j], limbs);
         if (mpn_cmp(kk[j], jj[j], limbs + 1) != 0)
         {
            gmp_printf("Error in entry %ld, %Nx != %Nx\n", j, kk[j], limbs + 1, jj[j], limbs + 1);
            abort();
         }
      }
   }

   TMP_FREE;

   gmp_randclear(state);
}

void test_fft_ifft_truncate_sqrt2()
{
   mp_bitcnt_t depth = 15UL;
   mp_size_t n = (1UL<<depth);            
   mp_bitcnt_t w = 1;
   mp_size_t iter = 1;

   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s, count;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, ** kk, *tt, *t1, *t2, *u1, *u2, *v1, *v2, *s1, *s2, *s3;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;
 
   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(4*n + 4*n*size + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      rand_n(ii[i], state, limbs);
   }
   t1 = (mp_limb_t *) ii + 4*n + 4*n*size;
   t2 = t1 + size;
   s1 = t2 + size;

   for (j = 0; j < 4*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);

   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(4*n + 4*n*size + 3*size);
   for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
   {
      jj[i] = ptr;
      MPN_COPY(jj[i], ii[i], limbs + 1);
   }
   u1 = (mp_limb_t *) jj + 4*n + 4*n*size;
   u2 = u1 + size;
   s2 = u2 + size;
   
   kk = (mp_limb_t **) TMP_BALLOC_LIMBS(4*n + 4*n*size + 3*size);
   for (i = 0, j = 0, ptr = (mp_limb_t *) kk + 4*n; i < 4*n; i++, ptr += size) 
   {
      kk[i] = ptr;
   }
   v1 = (mp_limb_t *) kk + 4*n + 4*n*size;
   v2 = v1 + size;
   s3 = v2 + size;
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (count = 0; count < iter; count++)
   {
      for (i = 0; i < 4*n; i++) 
      {
         rand_n(ii[i], state, limbs);
      }
   
      for (j = 0; j < 4*n; j++)
         mpn_normmod_2expp1(ii[j], limbs);

      for (i = 0; i < 4*n; i++) 
      {
         MPN_COPY(jj[i], ii[i], limbs + 1);
      }
      
      mp_size_t trunc = gmp_urandomm_ui(state, 2*n) + 2*n + 1;
      trunc = ((trunc + 7)/8)*8;
   
      FFT_radix2_truncate_sqrt2(ii, 1, ii, n, w, &t1, &t2, &s1, trunc);
      for (j = 0; j < trunc; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         MPN_COPY(kk[j], ii[j], limbs + 1);
      }
      
      IFFT_radix2_truncate_sqrt2(kk, 1, kk, n, w, &v1, &v2, &s3, trunc);
      for (j = 0; j < trunc; j++)
      {
         mpn_mul_2expmod_2expp1(jj[j], jj[j], limbs, depth + 2);
         mpn_normmod_2expp1(jj[j], limbs);
         mpn_normmod_2expp1(kk[j], limbs);
         if (mpn_cmp(kk[j], jj[j], limbs + 1) != 0)
         {
            gmp_printf("Error in entry %ld, %Nx != %Nx\n", j, kk[j], limbs + 1, jj[j], limbs + 1);
            abort();
         }
      }
   }

   TMP_FREE;

   gmp_randclear(state);
}

void test_fft_ifft_mfa_truncate_sqrt2()
{
   mp_bitcnt_t depth = 15UL;
   mp_size_t n = (1UL<<depth);            
   mp_size_t sqrt = (1UL<<(depth/2));
   mp_bitcnt_t w = 1;
   mp_size_t iter = 1;

   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s, count;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, ** kk, *tt, *t1, *t2, *u1, *u2, *v1, *v2, *s1, *s2, *s3;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;
 
   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(4*n + 4*n*size + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      rand_n(ii[i], state, limbs);
   }
   t1 = (mp_limb_t *) ii + 4*n + 4*n*size;
   t2 = t1 + size;
   s1 = t2 + size;

   for (j = 0; j < 4*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);

   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(4*n + 4*n*size + 3*size);
   for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
   {
      jj[i] = ptr;
      MPN_COPY(jj[i], ii[i], limbs + 1);
   }
   u1 = (mp_limb_t *) jj + 4*n + 4*n*size;
   u2 = u1 + size;
   s2 = u2 + size;
   
   kk = (mp_limb_t **) TMP_BALLOC_LIMBS(4*n + 4*n*size + 3*size);
   for (i = 0, j = 0, ptr = (mp_limb_t *) kk + 4*n; i < 4*n; i++, ptr += size) 
   {
      kk[i] = ptr;
   }
   v1 = (mp_limb_t *) kk + 4*n + 4*n*size;
   v2 = v1 + size;
   s3 = v2 + size;
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (count = 0; count < iter; count++)
   {
      for (i = 0; i < 4*n; i++) 
      {
         rand_n(ii[i], state, limbs);
      }
   
      for (j = 0; j < 4*n; j++)
         mpn_normmod_2expp1(ii[j], limbs);

      for (i = 0; i < 4*n; i++) 
      {
         MPN_COPY(jj[i], ii[i], limbs + 1);
      }
      
      mp_size_t trunc = gmp_urandomm_ui(state, 2*n) + 2*n + 1;
      trunc = ((trunc + sqrt - 1)/sqrt)*sqrt;
   
      FFT_radix2_mfa_truncate_sqrt2(ii, n, w, &t1, &t2, &s1, sqrt, trunc);
      for (j = 0; j < 4*n; j++)
      {
         mpn_normmod_2expp1(ii[j], limbs);
         MPN_COPY(kk[j], ii[j], limbs + 1);
      }
      
      IFFT_radix2_mfa_truncate_sqrt2(kk, n, w, &v1, &v2, &s3, sqrt, trunc);
      for (j = 0; j < trunc; j++)
      {
         mpn_mul_2expmod_2expp1(jj[j], jj[j], limbs, depth + 2);
         mpn_normmod_2expp1(jj[j], limbs);
         mpn_normmod_2expp1(kk[j], limbs);
         if (mpn_cmp(kk[j], jj[j], limbs + 1) != 0)
         {
            gmp_printf("Error in entry %ld, %Nx != %Nx\n", j, kk[j], limbs + 1, jj[j], limbs + 1);
            abort();
         }
      }
   }

   TMP_FREE;

   gmp_randclear(state);
}

void test_fft_ifft_mfa()
{
   mp_bitcnt_t depth = 12UL;
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t w = 1;
   mp_size_t sqrt = (1UL<<(depth/2));
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, ** kk, *tt, *t1, *t2, *u1, *u2, *v1, *v2, **s1, **s2, **s3;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      rand_n(ii[i], state, limbs);
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = (mp_limb_t **) t2 + size;

   for (j = 0; j < 2*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);

   
   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size) 
   {
      jj[i] = ptr;
      MPN_COPY(jj[i], ii[i], limbs + 1);
   }
   u1 = ptr;
   u2 = u1 + size;
   s2 = (mp_limb_t **) u2 + size;
   
   kk = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) kk + 2*n; i < 2*n; i++, ptr += size) 
   {
      kk[i] = ptr;
   }
   v1 = ptr;
   v2 = v1 + size;
   s3 = (mp_limb_t **) v2 + size;
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (i = 0; i < 10; i++)
   {
      for (j = 0; j < 2*n; j++) 
      {
          rand_n(ii[j], state, limbs);
          mpn_normmod_2expp1(ii[j], limbs);
          MPN_COPY(jj[j], ii[j], limbs + 1);
      }

      FFT_radix2_mfa(ii, n, w, &t1, &t2, s1, sqrt);
      for (j = 0; j < 2*n; j++)
         mpn_normmod_2expp1(ii[j], limbs);
      for (j = 0; j < 2*n; j++)
      {
         MPN_COPY(kk[j], ii[j], limbs + 1);
      }
      IFFT_radix2_mfa(kk, n, w, &v1, &v2, s3, sqrt);
      for (j = 0; j < 2*n; j++)
      {
         mpn_mul_2expmod_2expp1(jj[j], jj[j], limbs, depth + 1);
         mpn_normmod_2expp1(jj[j], limbs);
         mpn_normmod_2expp1(kk[j], limbs);
      }

      for (j = 0; j < 2*n; j++)
      {
         if (mpn_cmp(kk[j], jj[j], limbs + 1) != 0)
         {
            gmp_printf("Error in entry %ld, %Nx != %Nx\n", j, kk[j], limbs + 1, jj[j], limbs + 1);
            abort();
         }
      }
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void test_fft_ifft_mfa_sqrt2()
{
   mp_bitcnt_t depth = 13UL;
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t w = 4;
   mp_size_t sqrt = (1UL<<(depth/2));
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, ** kk, *tt, *t1, *t2, *s1;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size) + 3*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      rand_n(ii[i], state, limbs);
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = t2 + size;
   
   for (j = 0; j < 4*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);

   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size));
   for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
   {
      jj[i] = ptr;
      MPN_COPY(jj[i], ii[i], limbs + 1);
   }
   
   kk = (mp_limb_t **) TMP_BALLOC_LIMBS(4*(n + n*size));
   for (i = 0, ptr = (mp_limb_t *) kk + 4*n; i < 4*n; i++, ptr += size) 
   {
      kk[i] = ptr;
   }
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (i = 0; i < 1; i++)
   {
      FFT_radix2_mfa_sqrt2(ii, n, w, &t1, &t2, &s1, sqrt);

      for (j = 0; j < 4*n; j++)
         mpn_normmod_2expp1(ii[j], limbs);
      
      for (j = 0; j < 4*n; j++)
         MPN_COPY(kk[j], ii[j], limbs + 1);

      IFFT_radix2_mfa_sqrt2(kk, n, w, &t1, &t2, &s1, sqrt);
      for (j = 0; j < 4*n; j++)
      {
         mpn_mul_2expmod_2expp1(jj[j], jj[j], limbs, depth + 2);
         mpn_normmod_2expp1(jj[j], limbs);
         mpn_normmod_2expp1(kk[j], limbs);
      }

      for (j = 0; j < 4*n; j++)
      {
         if (mpn_cmp(kk[j], jj[j], limbs + 1) != 0)
         {
            gmp_printf("Error in entry %ld, %Nx != %Nx\n", j, kk[j], limbs + 1, jj[j], limbs + 1);
            abort();
         }
      }
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void test_fft_ifft_mfa_truncate()
{
   mp_bitcnt_t depth = 12UL;
   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t w = 1;
   mp_size_t iters = 100;
   mp_size_t sqrt = (1UL<<(depth/2));
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s, count, trunc;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, ** kk, *tt, *t1, *t2, *u1, *u2, *v1, *v2, **s1, **s2, **s3;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      rand_n(ii[i], state, limbs);
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = (mp_limb_t **) t2 + size;

   for (j = 0; j < 2*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);

   
   jj = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size) 
   {
      jj[i] = ptr;
      MPN_COPY(jj[i], ii[i], limbs + 1);
   }
   u1 = ptr;
   u2 = u1 + size;
   s2 = (mp_limb_t **) u2 + size;
   
   kk = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) kk + 2*n; i < 2*n; i++, ptr += size) 
   {
      kk[i] = ptr;
   }
   v1 = ptr;
   v2 = v1 + size;
   s3 = (mp_limb_t **) v2 + size;
   
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (count = 0; count < iters; count++)
   {
      mp_size_t trunc = (gmp_urandomm_ui(state, n/sqrt) + 1)*sqrt*2;
      for (i = 0; i < 2*n; i++) 
      {
         rand_n(ii[i], state, limbs);
         mpn_normmod_2expp1(ii[i], limbs);
         MPN_COPY(jj[i], ii[i], limbs + 1);
      }

      FFT_radix2_mfa_truncate(ii, n, w, &t1, &t2, s1, sqrt, trunc);
      
      for (j = 0; j < 2*n; j++)
         mpn_normmod_2expp1(ii[j], limbs);
      
      IFFT_radix2_mfa_truncate(ii, n, w, &v1, &v2, s3, sqrt, trunc);
      
      for (j = 0; j < 2*n; j++)
      {
         mpn_mul_2expmod_2expp1(jj[j], jj[j], limbs, depth + 1);
         mpn_normmod_2expp1(jj[j], limbs);
         mpn_normmod_2expp1(ii[j], limbs);
      }

      for (j = 0; j < trunc; j++)
      {
         if (mpn_cmp(ii[j], jj[j], limbs + 1) != 0)
         {
            gmp_printf("Error in entry %ld, %Nx != %Nx\n", j, kk[j], limbs + 1, jj[j], limbs + 1);
            abort();
         }
      }
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void test_fft_truncate()
{
   mp_bitcnt_t depth = 10UL;
   mp_size_t n = (1UL<<depth);            
   mp_bitcnt_t w = 1;
   mp_size_t iter = 1000;

   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s, count;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, ** kk, *tt, *t1, *t2, *u1, *u2, *v1, *v2, **s1, **s2, **s3;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;
 
   for (count = 0; count < iter; count++)
   {
      mp_size_t trunc = gmp_urandomm_ui(state, 2*n) + 1;
      trunc = ((trunc + 7)/8)*8;
   
      TMP_MARK;

      ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*n + 2*n*size + 2*n + 2*size);
      for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
      {
         ii[i] = ptr;
         if (i < trunc) rand_n(ii[i], state, limbs);
         else MPN_ZERO(ii[i], limbs + 1);
      }
      t1 = (mp_limb_t *) ii + 2*n + 2*n*size;
      t2 = t1 + size;
      s1 = (mp_limb_t **) t2 + size;

      for (j = 0; j < 2*n; j++)
         mpn_normmod_2expp1(ii[j], limbs);

      jj = (mp_limb_t **) TMP_BALLOC_LIMBS(2*n + 2*n*size + 2*n + 2*size);
      for (i = 0, ptr = (mp_limb_t *) jj + 2*n; i < 2*n; i++, ptr += size) 
      {
         jj[i] = ptr;
         MPN_COPY(jj[i], ii[i], limbs + 1);
      }
      u1 = (mp_limb_t *) jj + 2*n + 2*n*size;
      u2 = u1 + size;
      s2 = (mp_limb_t **) u2 + size;
   
      tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
      for (i = 0; i < 1; i++)
      {
         FFT_radix2_truncate(ii, 1, ii, n, w, &t1, &t2, s1, trunc);
         FFT_radix2(jj, 1, jj, n, w, &u1, &u2, s2);
         
         for (j = 0; j < trunc; j++)
         {
            mpn_normmod_2expp1(jj[j], limbs);
            mpn_normmod_2expp1(ii[j], limbs);
            if (mpn_cmp(ii[j], jj[j], limbs + 1) != 0)
            {
               gmp_printf("Error in entry %ld, %Nx != %Nx\n", j, ii[j], limbs + 1, jj[j], limbs + 1);
               abort();
            }
         }
      }
      
      TMP_FREE;
   }

   gmp_randclear(state);
}

void time_mfa()
{
   mp_bitcnt_t depth = 12L;
   mp_size_t iters = 1000;
   mp_bitcnt_t w2 = 1;

   mp_size_t n = (1UL<<depth)/w2;
   mp_bitcnt_t w = w2*w2;
   
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *u1, *u2, **s1, **s2;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      if (i < n) rand_n(ii[i], state, limbs);
      else MPN_ZERO(ii[i], limbs + 1);
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = (mp_limb_t **) t2 + size;
   
   for (j = 0; j < 2*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);
     
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (i = 0; i < iters; i++)
   {
      FFT_radix2_mfa(ii, n, w, &t1, &t2, s1, (1UL<<(depth/2))/w2);
   }
     
   TMP_FREE;
   gmp_randclear(state);
}

void time_ifft()
{
   mp_bitcnt_t depth = 16L;
   mp_size_t iters = 1;
   mp_bitcnt_t w = 1;
   
   mp_size_t n = (1UL<<depth);
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *u1, *u2, **s1, **s2;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      if (i < n) rand_n(ii[i], state, limbs);
      else MPN_ZERO(ii[i], limbs + 1);
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = (mp_limb_t **) t2 + size;
   
   for (j = 0; j < 2*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);
     
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (i = 0; i < iters; i++)
   {
      IFFT_radix2(ii, 1, ii, n, w, &t1, &t2, s1);
   }
     
   TMP_FREE;
   gmp_randclear(state);
}

void time_negacyclic_fft()
{
   mp_bitcnt_t depth = 10L;
   mp_size_t iters = 10000;
   mp_bitcnt_t w = 4;
   
   mp_size_t n = 512;//(1UL<<depth);
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *u1, *u2, **s1, **s2;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      if (i < n) rand_n(ii[i], state, limbs);
      else MPN_ZERO(ii[i], limbs + 1);
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = (mp_limb_t **) t2 + size;

   for (j = 0; j < 2*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);
     
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (i = 0; i < iters; i++)
   {
      FFT_radix2_negacyclic(ii, 1, ii, n, w, &t1, &t2, s1);
   }
     
   TMP_FREE;
   gmp_randclear(state);
}

void time_imfa()
{
   mp_bitcnt_t depth = 16L;
   mp_size_t iters = 1;
   mp_bitcnt_t w2 = 1;

   mp_bitcnt_t w = w2*w2;
   mp_size_t n = (1UL<<depth)/w2;
   mp_size_t limbs = (n*w)/GMP_LIMB_BITS;
   mp_size_t size = limbs + 1;
   mp_size_t i, j, s;
   mp_limb_t * ptr;
   mp_limb_t ** ii, ** jj, *tt, *t1, *t2, *u1, *u2, **s1, **s2;
   mp_limb_t c;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   ii = (mp_limb_t **) TMP_BALLOC_LIMBS(2*(n + n*size) + 2*n + 2*size);
   for (i = 0, ptr = (mp_limb_t *) ii + 2*n; i < 2*n; i++, ptr += size) 
   {
      ii[i] = ptr;
      if (i < n) rand_n(ii[i], state, limbs);
      else MPN_ZERO(ii[i], limbs + 1);
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = (mp_limb_t **) t2 + size;

   for (j = 0; j < 2*n; j++)
      mpn_normmod_2expp1(ii[j], limbs);
     
   tt = (mp_limb_t *) TMP_BALLOC_LIMBS(2*size);
   
   for (i = 0; i < iters; i++)
   {
      IFFT_radix2_mfa(ii, n, w, &t1, &t2, s1, (1UL<<(depth/2))/w2);
   }
     
   TMP_FREE;
   gmp_randclear(state);
}

void time_mul()
{
   mp_bitcnt_t depth = 10UL;
   mp_bitcnt_t w = 3UL;
   mp_size_t iters = 100;

   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t bits1 = (n*w - depth)/2; 
   mp_bitcnt_t bits = (8364032*8)/8; //n*bits1;
   printf("bits = %ld\n", bits);
   mp_size_t int_limbs = (bits - 1UL)/GMP_LIMB_BITS + 1;
   
   mp_size_t i, j;
   mp_limb_t *i1, *i2, *r1, *r2;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   i1 = TMP_BALLOC_LIMBS(6*int_limbs);
   i2 = i1 + int_limbs;
   r1 = i2 + int_limbs;
   r2 = r1 + 2*int_limbs;
   
   mpn_urandomb(i1, state, bits);
   mpn_urandomb(i2, state, bits);
  
   for (i = 0; i < iters; i++)
   {
      new_mpn_mul(r1, i1, int_limbs, i2, int_limbs, depth, w);
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void time_mul2()
{
   mp_bitcnt_t depth = 17UL;
   mp_bitcnt_t w = 1UL;
   mp_size_t iters = 1;

   mp_size_t n = (1UL<<depth);
   mp_size_t sqrt = (1UL<<(depth/2));
   mp_bitcnt_t bits1 = (n*w - (depth + 1))/2; 
   mp_bitcnt_t bits = 2*n*bits1;
   printf("bits = %ld\n", bits);
   mp_size_t int_limbs = bits/GMP_LIMB_BITS;
   
   mp_size_t i, j;
   mp_limb_t *i1, *i2, *r1, *r2;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   i1 = TMP_BALLOC_LIMBS(6*int_limbs);
   i2 = i1 + int_limbs;
   r1 = i2 + int_limbs;
   r2 = r1 + 2*int_limbs;
   
   mpn_urandomb(i1, state, bits);
   mpn_urandomb(i2, state, bits);
  
   for (i = 0; i < iters; i++)
   {
      new_mpn_mul3(r1, i1, int_limbs, i2, int_limbs, depth, w, sqrt);
      //mpn_mul(r1, i1, int_limbs, i2, int_limbs);
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void time_mul4()
{
   mp_bitcnt_t depth = 13UL;
   mp_bitcnt_t w = 1UL;
   mp_size_t iters = 1;

   mp_size_t n = (1UL<<depth);
   mp_size_t sqrt = (1UL<<(depth/2));
   mp_bitcnt_t bits1 = (n*w - (depth + 1))/2; 
   mp_bitcnt_t bits = 2*n*bits1;
   mp_size_t int_limbs = bits/GMP_LIMB_BITS;
   //mp_size_t n1 = (3*int_limbs)/4;
   //mp_size_t n2 = (3*int_limbs)/4;
   mp_size_t n1 = int_limbs;
   mp_size_t n2 = int_limbs;
   mp_bitcnt_t b1 = n1 * GMP_LIMB_BITS;
   mp_bitcnt_t b2 = n2 * GMP_LIMB_BITS;
   
   mp_size_t i, j;
   mp_limb_t *i1, *i2, *r1, *r2;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   i1 = TMP_BALLOC_LIMBS(6*int_limbs);
   i2 = i1 + int_limbs;
   r1 = i2 + int_limbs;
   r2 = r1 + 2*int_limbs;
   
   mpn_urandomb(i1, state, b1);
   mpn_urandomb(i2, state, b2);
  
   printf("b1 = %ld, b2 = %ld\n", b1, b2);

   for (i = 0; i < iters; i++)
   {
      new_mpn_mul4(r1, i1, n1, i2, n2, depth, w);
      //mpn_mul(r1, i1, n1, i2, n2);
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void time_mul6()
{
   mp_bitcnt_t depth = 13UL;
   mp_bitcnt_t w = 2UL;
   mp_size_t iters = 1;

   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t bits1 = (n*w - (depth + 1))/2; 
   mp_bitcnt_t bits = 2*n*bits1;
   mp_size_t int_limbs = bits/GMP_LIMB_BITS;
   mp_size_t n1 = (3*int_limbs)/4;
   mp_size_t n2 = (3*int_limbs)/4;
   //mp_size_t n1 = int_limbs;
   //mp_size_t n2 = int_limbs;
   mp_bitcnt_t b1 = n1 * GMP_LIMB_BITS;
   mp_bitcnt_t b2 = n2 * GMP_LIMB_BITS;
   
   mp_size_t i, j;
   mp_limb_t *i1, *i2, *r1, *r2;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   i1 = TMP_BALLOC_LIMBS(6*int_limbs);
   i2 = i1 + int_limbs;
   r1 = i2 + int_limbs;
   r2 = r1 + 2*int_limbs;
   
   mpn_urandomb(i1, state, b1);
   mpn_urandomb(i2, state, b2);
  
   printf("b1 = %ld, b2 = %ld\n", b1, b2);

   for (i = 0; i < iters; i++)
   {
      new_mpn_mul6(r1, i1, n1, i2, n2, depth, w);
      //mpn_mul(r1, i1, n1, i2, n2);
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void test_mul()
{
   mp_bitcnt_t depth = 15UL;
   mp_bitcnt_t w = 2UL;
   mp_size_t iters = 1;

   mp_size_t n = (1UL<<depth);
   mp_size_t sqrt = (1UL<<(depth/2));
   mp_bitcnt_t bits1 = (n*w - (depth + 1))/2; 
   mp_bitcnt_t bits = 2*n*bits1;
   mp_size_t int_limbs = bits/GMP_LIMB_BITS;
   
   mp_size_t i, j;
   mp_limb_t *i1, *i2, *r1, *r2;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   i1 = TMP_BALLOC_LIMBS(6*int_limbs);
   i2 = i1 + int_limbs;
   r1 = i2 + int_limbs;
   r2 = r1 + 2*int_limbs;
   
   for (i = 0; i < iters; i++)
   {
      mpn_urandomb(i1, state, bits);
      mpn_urandomb(i2, state, bits);
  
      mpn_mul_n(r2, i1, i2, int_limbs);
      new_mpn_mul3(r1, i1, int_limbs, i2, int_limbs, depth, w, sqrt);
      
      for (j = 0; j < 2*int_limbs; j++)
      {
         if (r1[j] != r2[j]) 
         {
            printf("error in limb %ld, %lx != %lx\n", j, r1[j], r2[j]);
            abort();
         } 
      }
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void test_mul5()
{
   mp_bitcnt_t depth = 14UL;
   mp_bitcnt_t w = 1UL;
   mp_size_t iters = 1;

   mp_size_t n = (1UL<<depth);
   mp_size_t sqrt = (1UL<<(depth/2));
   mp_bitcnt_t bits1 = (n*w - (depth + 1))/2; 
   mp_bitcnt_t bits = n*bits1;
   mp_size_t int_limbs = bits/GMP_LIMB_BITS;
   mp_size_t n1 = (3*int_limbs)/4;
   mp_size_t n2 = (3*int_limbs)/4;
   mp_bitcnt_t b1 = n1 * GMP_LIMB_BITS;
   mp_bitcnt_t b2 = n2 * GMP_LIMB_BITS;
   
   mp_size_t i, j;
   mp_limb_t *i1, *i2, *r1, *r2;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   i1 = TMP_BALLOC_LIMBS(6*int_limbs);
   i2 = i1 + int_limbs;
   r1 = i2 + int_limbs;
   r2 = r1 + 2*int_limbs;
   
   for (i = 0; i < iters; i++)
   {
      mpn_urandomb(i1, state, b1);
      mpn_urandomb(i2, state, b2);
  
      mpn_mul(r2, i1, n1, i2, n2);
      new_mpn_mul5(r1, i1, n1, i2, n2, depth, w);
      
      for (j = 0; j < n1+n2; j++)
      {
         if (r1[j] != r2[j]) 
         {
            printf("error in limb %ld, %lx != %lx\n", j, r1[j], r2[j]);
            abort();
         } 
      }
   }
      
   TMP_FREE;
   gmp_randclear(state);
}

void test_mul4()
{
   mp_bitcnt_t depth = 14UL;
   mp_bitcnt_t w = 1UL;
   mp_size_t iters = 1;

   mp_size_t n = (1UL<<depth);
   mp_bitcnt_t bits1 = (n*w - (depth + 1))/2; 
   mp_bitcnt_t bits = 2*n*bits1;
   mp_size_t int_limbs = bits/GMP_LIMB_BITS;
   mp_size_t n1 = (3*int_limbs)/4;
   mp_size_t n2 = (3*int_limbs)/4;
   mp_bitcnt_t b1 = n1 * GMP_LIMB_BITS;
   mp_bitcnt_t b2 = n2 * GMP_LIMB_BITS;
   
   mp_size_t i, j;
   mp_limb_t *i1, *i2, *r1, *r2;
   gmp_randstate_t state;
   gmp_randinit_default(state);

   TMP_DECL;

   TMP_MARK;

   i1 = TMP_BALLOC_LIMBS(6*int_limbs);
   i2 = i1 + int_limbs;
   r1 = i2 + int_limbs;
   r2 = r1 + 2*int_limbs;
   
   for (i = 0; i < iters; i++)
   {
      mpn_urandomb(i1, state, b1);
      mpn_urandomb(i2, state, b2);
  
      mpn_mul(r2, i1, n1, i2, n2);
      new_mpn_mul6(r1, i1, n1, i2, n2, depth, w);
      
      for (j = 0; j < n1+n2; j++)
      {
         if (r1[j] != r2[j]) 
         {
            printf("error in limb %ld, %lx != %lx\n", j, r1[j], r2[j]);
            abort();
         } 
      }
   }
      
   TMP_FREE;
   gmp_randclear(state);
} 

int main(void)
{
#if TEST
   test_mulmod(); printf("MULMOD....PASS\n");
   test_fft_ifft_negacyclic(); printf("FFT_IFFT_NEGACYCLIC...PASS\n");
   test_mul4(); printf("MUL4...PASS\n");
   test_fft_ifft_mfa_truncate(); printf("FFT_IFFT_MFA_TRUNCATE...PASS\n");
   test_fft_ifft_mfa(); printf("FFT_IFFT_MFA...PASS\n");
   test_fft_ifft_mfa_sqrt2(); printf("FFT_IFFT_MFA_SQRT2...PASS\n");
   test_fft_ifft_mfa_truncate_sqrt2(); printf("FFT_IFFT_MFA_TRUNCATE_SQRT2...PASS\n");
   test_mul5(); printf("MUL5...PASS\n");
   
   test_fft_ifft_sqrt2(); printf("FFT_IFFT_SQRT...PASS\n");
   test_norm(); printf("mpn_normmod_2expp1...PASS\n");
   test_mul_2expmod(); printf("mpn_mul_2expmod_2expp1...PASS\n");
   test_div_2expmod(); printf("mpn_div_2expmod_2expp1...PASS\n");
   test_lshB_sumdiffmod(); printf("mpn_lshB_sumdiffmod_2expp1...PASS\n");
   test_sumdiff_rshBmod(); printf("mpn_sumdiff_rshBmod_2expp1...PASS\n");
   
   test_FFT_negacyclic_twiddle(); printf("FFT_negacyclic_twiddle...PASS\n");
   test_fft_ifft(); printf("FFT_IFFT...PASS\n");
   test_fft_truncate(); printf("FFT_TRUNCATE...PASS\n");
   test_fft_ifft_truncate(); printf("FFT_IFFT_TRUNCATE...PASS\n");
   test_fft_ifft_truncate_sqrt2(); printf("FFT_IFFT_TRUNCATE_SQRT2...PASS\n");
   
#endif

#if TIME
   //time_ifft();
   //time_mfa();
   //time_imfa();
   //time_mul_with_negacyclic(); // negacyclic is currently *disabled*
   //time_negacyclic_fft();
   time_mul6();
#endif

   return 0;
}
