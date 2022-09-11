////////////////////////////////////////////////////////////////////////////
//                         **** RESAMPLER ****                            //
//                     Sinc-based Audio Resampling                        //
//                Copyright (c) 2006 - 2022 David Bryant.                 //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// resampler.c

#include "resampler.h"

#ifndef M_PI
#define M_PI 3.14159265358979324
#endif

// This is the basic convolution operation that is the core of the resampler and utilizes the
// bulk of the CPU load (assuming reasonably long filters). The first version is the canonical
// form, followed by three variations that may or may not be faster depending on your compiler,
// options, and system. Try 'em and use the fastest, or rewrite them using SIMD. Note that on
// gcc and clang, -Ofast can make a huge difference.

#if 1   // Version 1 (canonical)
static double apply_filter (float *A, float *B, int num_taps)
{
    float sum = 0.0;

    do sum += *A++ * *B++;
    while (--num_taps);

    return sum;
}
#endif

#if 0   // Version 2 (2x unrolled loop)
static double apply_filter (float *A, float *B, int num_taps)
{
    int num_loops = num_taps >> 1;
    float sum = 0.0;

    do {
        sum += (A[0] * B[0]) + (A[1] * B[1]);
        A += 2; B += 2;
    } while (--num_loops);

    return sum;
}
#endif

#if 0   // Version 3 (4x unrolled loop)
static double apply_filter (float *A, float *B, int num_taps)
{
    int num_loops = num_taps >> 2;
    float sum = 0.0;

    do {
        sum += (A[0] * B[0]) + (A[1] * B[1]) + (A[2] * B[2]) + (A[3] * B[3]);
        A += 4; B += 4;
    } while (--num_loops);

    return sum;
}
#endif

#if 0   // Version 4 (outside-in order, may be more accurate)
static double apply_filter (float *A, float *B, int num_taps)
{
    int i = num_taps - 1;
    float sum = 0.0;

    do {
        sum += (A[0] * B[0]) + (A[i] * B[i]);
        A++; B++;
    } while ((i -= 2) > 0);

    return sum;
}
#endif

static void init_filter (Resample *cxt, float *filter, double fraction, double lowpass_ratio)
{
    const double a0 = 0.35875;
    const double a1 = 0.48829;
    const double a2 = 0.14128;
    const double a3 = 0.01168;
    double filter_sum = 0.0;
    int i;

    // "dist" is the absolute distance from the sinc maximum to the filter tap to be calculated, in radians
    // "ratio" is that distance divided by half the tap count such that it reaches π at the window extremes

    // Note that with this scaling, the odd terms of the Blackman-Harris calculation appear to be negated
    // with respect to the reference formula version.

    for (i = 0; i < cxt->numTaps; ++i) {
        double dist = fabs ((cxt->numTaps / 2 - 1) + fraction - i) * M_PI;
        double ratio = dist / (cxt->numTaps / 2);
        double value;

        if (dist != 0.0) {
            value = sin (dist * lowpass_ratio) / (dist * lowpass_ratio);

            if (cxt->flags & BLACKMAN_HARRIS)
                value *= a0 + a1 * cos (ratio) + a2 * cos (2 * ratio) + a3 * cos (3 * ratio);
            else
                value *= 0.5 * (1.0 + cos (ratio));     // Hann window
        }
        else
            value = 1.0;

        filter_sum += cxt->tempFilter [i] = value;
    }

    // filter should have unity DC gain

    double scaler = 1.0 / filter_sum, error = 0.0;

    for (i = cxt->numTaps / 2; i < cxt->numTaps; i = cxt->numTaps - i - (i >= cxt->numTaps / 2)) {
        filter [i] = (cxt->tempFilter [i] *= scaler) - error;
        error += filter [i] - cxt->tempFilter [i];
    }
}

static double subsample_no_interpolate (Resample *cxt, float *source, double offset)
{
    source += (int) floor (offset);
    offset -= floor (offset);

    if (offset == 0.0 && !(cxt->flags & INCLUDE_LOWPASS))
        return *source;

    return apply_filter (
        cxt->filters [(int) floor (offset * cxt->numFilters + 0.5)],
        source - cxt->numTaps / 2 + 1, cxt->numTaps);
}

static double subsample_interpolate (Resample *cxt, float *source, double offset)
{
    double sum1, sum2;
    int i;

    source += (int) floor (offset);
    offset -= floor (offset);

     if (offset == 0.0 && !(cxt->flags & INCLUDE_LOWPASS))
        return *source;

    i = (int) floor (offset *= cxt->numFilters);
    sum1 = apply_filter (cxt->filters [i], source - cxt->numTaps / 2 + 1, cxt->numTaps);

    if ((offset -= i) == 0.0 && !(cxt->flags & INCLUDE_LOWPASS))
        return sum1;

    sum2 = apply_filter (cxt->filters [i+1], source - cxt->numTaps / 2 + 1, cxt->numTaps);

    return sum2 * offset + sum1 * (1.0 - offset);
}

double subsample (Resample *cxt, float *source, double offset)
{
    if (cxt->flags & SUBSAMPLE_INTERPOLATE)
        return subsample_interpolate (cxt, source, offset);
    else
        return subsample_no_interpolate (cxt, source, offset);
}

Resample *resampleInit (int numChannels, int numTaps, int numFilters, double lowpassRatio, int flags)
{
    Resample *cxt = calloc (1, sizeof (Resample));
    int i;

    if (lowpassRatio > 0.0 && lowpassRatio < 1.0)
        flags |= INCLUDE_LOWPASS;
    else {
        flags &= ~INCLUDE_LOWPASS;
        lowpassRatio = 1.0;
    }

    if ((numTaps & 3) || numTaps <= 0 || numTaps > 1024) {
        fprintf (stderr, "must 4-1024 filter taps, and a multiple of 4!\n");
        return NULL;
    }

    if (numFilters < 2 || numFilters > 1024) {
        fprintf (stderr, "must be 2-1024 filters!\n");
        return NULL;
    }

    cxt->numChannels = numChannels;
    cxt->numSamples = numTaps * 16;
    cxt->numFilters = numFilters;
    cxt->numTaps = numTaps;
    cxt->flags = flags;

    // note that we actually have one more than the specified number of filters

    cxt->filters = calloc (cxt->numFilters + 1, sizeof (float*));
    cxt->tempFilter = malloc (numTaps * sizeof (double));

    for (i = 0; i <= cxt->numFilters; ++i) {
        cxt->filters [i] = calloc (cxt->numTaps, sizeof (float));
        init_filter (cxt, cxt->filters [i], (double) i / cxt->numFilters, lowpassRatio);
    }

    free (cxt->tempFilter); cxt->tempFilter = NULL;
    cxt->buffers = calloc (numChannels, sizeof (float*));

    for (i = 0; i < numChannels; ++i)
        cxt->buffers [i] = calloc (cxt->numSamples, sizeof (float));

    cxt->outputOffset = numTaps / 2;
    cxt->inputIndex = numTaps;

    return cxt;
}

void resampleReset (Resample *cxt)
{
    int i;

    for (i = 0; i < cxt->numChannels; ++i)
        memset (cxt->buffers [i], 0,  cxt->numSamples * sizeof (float));

    cxt->outputOffset = cxt->numTaps / 2;
    cxt->inputIndex = cxt->numTaps;
}

ResampleResult resampleProcess (Resample *cxt, const float *const *input, int numInputFrames, float *const *output, int numOutputFrames, double ratio)
{
    int half_taps = cxt->numTaps / 2, i;
    ResampleResult res = { 0, 0 };

    while (numOutputFrames > 0) {
        if (cxt->outputOffset >= cxt->inputIndex - half_taps) {
            if (numInputFrames > 0) {
                if (cxt->inputIndex == cxt->numSamples) {
                    for (i = 0; i < cxt->numChannels; ++i)
                        memmove (cxt->buffers [i], cxt->buffers [i] + cxt->numSamples - cxt->numTaps, cxt->numTaps * sizeof (float));

                    cxt->outputOffset -= cxt->numSamples - cxt->numTaps;
                    cxt->inputIndex -= cxt->numSamples - cxt->numTaps;
                }

                for (i = 0; i < cxt->numChannels; ++i)
                    cxt->buffers [i] [cxt->inputIndex] = input [i] [res.input_used];

                cxt->inputIndex++;
                res.input_used++;
                numInputFrames--;
            }
            else
                break;
        }
        else {
            for (i = 0; i < cxt->numChannels; ++i)
                output [i] [res.output_generated] = subsample (cxt, cxt->buffers [i], cxt->outputOffset);

            cxt->outputOffset += (1.0 / ratio);
            res.output_generated++;
            numOutputFrames--;
        }
    }

    return res;
}

ResampleResult resampleProcessInterleaved (Resample *cxt, const float *input, int numInputFrames, float *output, int numOutputFrames, double ratio)
{
    int half_taps = cxt->numTaps / 2, i;
    ResampleResult res = { 0, 0 };

    while (numOutputFrames > 0) {
        if (cxt->outputOffset >= cxt->inputIndex - half_taps) {
            if (numInputFrames > 0) {
                if (cxt->inputIndex == cxt->numSamples) {
                    for (i = 0; i < cxt->numChannels; ++i)
                        memmove (cxt->buffers [i], cxt->buffers [i] + cxt->numSamples - cxt->numTaps, cxt->numTaps * sizeof (float));

                    cxt->outputOffset -= cxt->numSamples - cxt->numTaps;
                    cxt->inputIndex -= cxt->numSamples - cxt->numTaps;
                }

                for (i = 0; i < cxt->numChannels; ++i)
                    cxt->buffers [i] [cxt->inputIndex] = *input++;

                cxt->inputIndex++;
                res.input_used++;
                numInputFrames--;
            }
            else
                break;
        }
        else {
            for (i = 0; i < cxt->numChannels; ++i)
                *output++ = subsample (cxt, cxt->buffers [i], cxt->outputOffset);

            cxt->outputOffset += (1.0 / ratio);
            res.output_generated++;
            numOutputFrames--;
        }
    }

    return res;
}

unsigned int resampleGetRequiredSamples (Resample *cxt, int numOutputFrames, double ratio)
{
    int half_taps = cxt->numTaps / 2;
    int input_index = cxt->inputIndex;
    double offset = cxt->outputOffset;
    ResampleResult res = { 0, 0 };

    while (numOutputFrames > 0) {
        if (offset >= input_index - half_taps) {
            if (input_index == cxt->numSamples) {
                offset -= cxt->numSamples - cxt->numTaps;
                input_index -= cxt->numSamples - cxt->numTaps;
            }

            input_index++;
            res.input_used++;
        }
        else {
            offset += (1.0 / ratio);
            numOutputFrames--;
        }
    }

    return res.input_used;
}

unsigned int resampleGetExpectedOutput (Resample *cxt, int numInputFrames, double ratio)
{
    int half_taps = cxt->numTaps / 2;
    int input_index = cxt->inputIndex;
    double offset = cxt->outputOffset;
    ResampleResult res = { 0, 0 };

    while (1) {
        if (offset >= input_index - half_taps) {
            if (numInputFrames > 0) {
                if (input_index == cxt->numSamples) {
                    offset -= cxt->numSamples - cxt->numTaps;
                    input_index -= cxt->numSamples - cxt->numTaps;
                }

                input_index++;
                numInputFrames--;
            }
            else
                break;
        }
        else {
            offset += (1.0 / ratio);
            res.output_generated++;
        }
    }

    return res.output_generated;
}

void resampleAdvancePosition (Resample *cxt, double delta)
{
    if (delta < 0.0)
        fprintf (stderr, "resampleAdvancePosition() can only advance forward!\n");
    else
        cxt->outputOffset += delta;
}

double resampleGetPosition (Resample *cxt)
{
    return cxt->outputOffset + (cxt->numTaps / 2.0) - cxt->inputIndex;
}

void resampleFree (Resample *cxt)
{
    int i;

    for (i = 0; i <= cxt->numFilters; ++i)
        free (cxt->filters [i]);

    free (cxt->filters);

    for (i = 0; i < cxt->numChannels; ++i)
        free (cxt->buffers [i]);

    free (cxt->buffers);
    free (cxt);
}
