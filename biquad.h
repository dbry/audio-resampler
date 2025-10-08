////////////////////////////////////////////////////////////////////////////
//                           **** BIQUAD ****                             //
//                     Simple Biquad Filter Library                       //
//                Copyright (c) 2021 - 2022 David Bryant.                 //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// biquad.h

#ifndef BIQUAD_H
#define BIQUAD_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

typedef struct {
    float a0, a1, a2, a3, a4, b1, b2, b3, b4;
} BiquadCoefficients;

typedef struct {
    float a[5], b[5];               // coefficients
    float x[4], y[4];               // delayed input/output
    int order, index;
} Biquad;

#ifdef __cplusplus
extern "C" {
#endif

void biquad_init (Biquad *f, const BiquadCoefficients *coeffs, float gain);

void biquad_lowpass (BiquadCoefficients *filter, double frequency);
void biquad_highpass (BiquadCoefficients *filter, double frequency);

void biquad_apply_buffer (Biquad *f, float *buffer, int num_samples, int stride);
float biquad_apply_sample (Biquad *f, float input);

#ifdef __cplusplus
}
#endif

#endif

