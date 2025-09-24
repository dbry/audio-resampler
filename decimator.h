////////////////////////////////////////////////////////////////////////////
//                           **** DECIMATOR ****                          //
//                    Float to Integer Audio Decimation                   //
//                 Copyright (c) 2006 - 2024 David Bryant                 //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// decimator.h

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "biquad.h"

#ifdef ENABLE_THREADS
#include "workers.h"
#endif

#define DITHER_HIGHPASS     0x1
#define DITHER_FLAT         0x2
#define DITHER_LOWPASS      0x4
#define DITHER_ENABLED      (DITHER_HIGHPASS | DITHER_FLAT | DITHER_LOWPASS)

#define SHAPING_1ST_ORDER   0x100
#define SHAPING_2ND_ORDER   0x200
#define SHAPING_3RD_ORDER   0x400
#define SHAPING_ATH_CURVE   0x800
#define SHAPING_ENABLED     (SHAPING_1ST_ORDER | SHAPING_2ND_ORDER | SHAPING_3RD_ORDER | SHAPING_ATH_CURVE)

#define DECIMATE_MULTITHREADED  0x1000

typedef struct {
    int numChannels, outputBits, outputBytes, dither_type, flags;
    float outputGain, *feedback;
    uint32_t *tpdf_generators;
    Biquad *noise_shapers;

#ifdef ENABLE_THREADS
    Workers *workers;
    const float *input;
    int numInputFrames;
    unsigned char *output;
    int stride, clips;

    uint32_t tpdf_generator;
    Biquad noise_shaper;
    float feedback_val;
#endif
} Decimate;

#ifdef __cplusplus
extern "C" {
#endif

void floatIntegersLE (unsigned char *input, float inputGain, int inputBits, int inputBytes, int inputStride, float *output, int numSamples);

Decimate *decimateInit (int numChannels, int outputBits, int outputBytes, float outputGain, int sampleRate, int flags);
int decimateProcessLE (Decimate *cxt, const float *const *input, int numInputFrames, unsigned char *const *output);
int decimateProcessInterleavedLE (Decimate *cxt, const float *input, int numInputFrames, unsigned char *output);
void decimateFree (Decimate *cxt);

#ifdef __cplusplus
}
#endif
