////////////////////////////////////////////////////////////////////////////
//                         **** RESAMPLER ****                            //
//                     Sinc-based Audio Resampling                        //
//                Copyright (c) 2006 - 2023 David Bryant.                 //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// resampler.h

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#ifdef ENABLE_THREADS
#include "workers.h"
#endif

#define SUBSAMPLE_INTERPOLATE   0x1
#define BLACKMAN_HARRIS         0x2
#define INCLUDE_LOWPASS         0x4
#define RESAMPLE_MULTITHREADED  0x8
#define NO_FILTER_REDUCTION     0x10
#define RESAMPLE_FIXED_RATIO    0x20    // internal use only, do not set

typedef struct {
    unsigned int input_used, output_generated;
} ResampleResult;

typedef struct {
    int numChannels, numSamples, numFilters, numTaps, inputIndex, flags;
    double *tempFilter, outputOffset, fixedRatio, lowpassRatio;
    float **buffers, **filters;

#ifdef ENABLE_THREADS
    Workers *workers;
    const float *input;
    float *cbuffer, *output;
    int numInputFrames, numOutputFrames, stride;
    ResampleResult res;
    double ratio;
#endif
} Resample;

#ifdef __cplusplus
extern "C" {
#endif

Resample *resampleInit (int numChannels, int numTaps, int numFilters, double lowpassRatio, int flags);
Resample *resampleFixedRatioInit (int numChannels, int numTaps, int maxFilters, double sourceRate, double destinRate, int lowpassFreq, int flags);
ResampleResult resampleProcess (Resample *cxt, const float *const *input, int numInputFrames, float *const *output, int numOutputFrames, double ratio);
ResampleResult resampleProcessInterleaved (Resample *cxt, const float *input, int numInputFrames, float *output, int numOutputFrames, double ratio);
unsigned int resampleGetRequiredSamples (Resample *cxt, int numOutputFrames, double ratio);
unsigned int resampleGetExpectedOutput (Resample *cxt, int numInputFrames, double ratio);
void resampleAdvancePosition (Resample *cxt, double delta);
double resampleGetLowpassRatio (Resample *cxt);
double resampleGetPosition (Resample *cxt);
int resampleGetNumFilters (Resample *cxt);
int resampleInterpolationUsed (Resample *cxt);
void resampleReset (Resample *cxt);
void resampleFree (Resample *cxt);

#ifdef __cplusplus
}
#endif
