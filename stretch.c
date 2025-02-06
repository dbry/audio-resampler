////////////////////////////////////////////////////////////////////////////
//                        **** AUDIO-STRETCH ****                         //
//                      Time Domain Harmonic Scaler                       //
//                    Copyright (c) 2022 David Bryant                     //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// stretch.c

// Time Domain Harmonic Compression and Expansion
//
// This library performs time domain harmonic scaling with pitch detection
// to stretch the timing of a 32-bit PCM signal (either mono or stereo) from
// 1/4 to 4 times its original length. This is done without altering any of
// the tonal characteristics.
//
// Use stereo (num_chans = 2), when both channels are from same source
// and should contain approximately similar content. For independent channels,
// prefer using multiple instances.
// see https://github.com/dbry/audio-stretch/issues/6

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "stretch.h"

static void merge_blocks (float *output, float *input1, float *input2, int samples);
static int find_period_fast (const Stretch *cxt, const float *samples);
static int find_period (const Stretch *cxt, float *samples);

/*
 * Initialize a context of the time stretching code. The shortest and longest periods
 * are specified here. The longest period determines the lowest fundamental frequency
 * that can be handled correctly. Note that higher frequencies can be handled than the
 * shortest period would suggest because multiple periods can be combined, and the
 * worst-case performance will suffer if too short a period is selected. The flags are:
 *
 * STRETCH_FAST_FLAG    0x1     Use the "fast" version of the period calculation
 *
 * STRETCH_DUAL_FLAG    0x2     Cascade two instances of the stretcher to expand
 *                              available ratios to 0.25X to 4.00X
 */

Stretch *stretchInit (int shortest_period, int longest_period, int num_channels, int flags)
{
    int max_periods = 3;
    Stretch *cxt;

    if (flags & STRETCH_FAST_FLAG) {
        longest_period = (longest_period + 1) & ~1;
        shortest_period &= ~1;
        max_periods = 4;
    }

    if (longest_period <= shortest_period || shortest_period < MIN_PERIOD || longest_period > MAX_PERIOD) {
        fprintf (stderr, "stretchInit(): invalid periods!\n");
        return NULL;
    }

    cxt = (Stretch *) calloc (1, sizeof (Stretch));

    if (cxt) {
        cxt->inbuff_samples = longest_period * num_channels * max_periods;
        cxt->inbuff = calloc (cxt->inbuff_samples, sizeof (*cxt->inbuff));

        if (num_channels == 2 || (flags & STRETCH_FAST_FLAG))
            cxt->calcbuff = calloc (longest_period * num_channels, sizeof (*cxt->calcbuff));

        if ((flags & STRETCH_FAST_FLAG))
            cxt->results = calloc (longest_period, sizeof (*cxt->results));
    }

    if (!cxt || !cxt->inbuff || (num_channels == 2 && (flags & STRETCH_FAST_FLAG) && !cxt->calcbuff) || ((flags & STRETCH_FAST_FLAG) && !cxt->results)) {
        fprintf (stderr, "stretchInit(): out of memory!\n");
        return NULL;
    }

    cxt->head = cxt->tail = cxt->longest = longest_period * num_channels;
    cxt->fast_mode = (flags & STRETCH_FAST_FLAG) ? 1 : 0;
    cxt->shortest = shortest_period * num_channels;
    cxt->num_chans = num_channels;

    if (flags & STRETCH_DUAL_FLAG) {
        cxt->next = stretchInit (shortest_period, longest_period, num_channels, flags & ~STRETCH_DUAL_FLAG);
        cxt->intermediate = calloc (longest_period * num_channels * max_periods, sizeof (*cxt->intermediate));
    }

    return cxt;
}

/*
 * Re-Initialize a context of the time stretching code - as if freshly created
 * with stretchInit(). This drops all internal state.
 */

void stretchReset (Stretch *cxt)
{
    cxt->head = cxt->tail = cxt->longest;
    memset (cxt->inbuff, 0, cxt->tail * sizeof (*cxt->inbuff));

    if (cxt->next)
        stretchReset (cxt->next);
}

/*
 * Determine how many samples (per channel) should be reserved in 'output'-array
 * for stretchProcess() and stretchFlush(). max_num_samples and max_ratio are the
 * maximum values that will be passed to stretchProcess().
 */

int stretchGetOutputCapacity (Stretch *cxt, int max_num_samples, double max_ratio)
{
    int max_period = cxt->longest / cxt->num_chans;
    int max_expected_samples;
    double next_ratio = 1.0;

    if (cxt->next) {
        if (max_ratio < 0.5) {
            next_ratio = max_ratio / 0.5;
            max_ratio = 0.5;
        }
        else if (max_ratio > 2.0) {
            next_ratio = max_ratio / 2.0;
            max_ratio = 2.0;
        }
        else
            next_ratio = 1.0;
    }

    max_expected_samples = (int) ceil (max_num_samples * ceil (max_ratio * 2.0) / 2.0) +
        max_period * (cxt->fast_mode ? 4 : 3);

    if (cxt->next)
        max_expected_samples = stretchGetOutputCapacity (cxt->next, max_expected_samples, next_ratio);

    return max_expected_samples;
}

/*
 * Process the specified samples with the given ratio (which is normally clipped to
 * the range 0.5 to 2.0, or 0.25 to 4.00 for the "dual" mode). Note that in stereo
 * the number of samples refers to the samples for one channel (i.e., not the total
 * number of values passed) and can be as large as desired (samples are buffered here).
 * The ratio may change between calls, but there is some latency to consider because
 * audio is buffered here and a new ratio may be applied to previously sent samples.
 *
 * The exact number of samples output is not easy to determine in advance, so a function
 * is provided (stretch_output_capacity()) that calculates the maximum number of samples
 * that can be generated from a single call to this function (or stretchFlush()) given
 * a number of samples and maximum ratio. It is reccomended that that function be used
 * after initialization to allocate in advance the buffer size required. Be sure to
 * multiply the return value by the number channels!
 */

int stretchProcess (Stretch *cxt, const float *samples, int num_samples, float *output, double ratio)
{
    int out_samples = 0, next_samples = 0;
    float *outbuf = output;
    double next_ratio = 1.0;

    /* if there's a cascaded instance after this one, try to do as much of the ratio here and the rest in "next" */

    if (cxt->next) {
        outbuf = cxt->intermediate;

        if (ratio < 0.5) {
            next_ratio = ratio / 0.5;
            ratio = 0.5;
        }
        else if (ratio > 2.0) {
            next_ratio = ratio / 2.0;
            ratio = 2.0;
        }
        else
            next_ratio = 1.0;
    }

    num_samples *= cxt->num_chans;

    /* this really should not happen, but a good idea to clamp in case */

    if (ratio < 0.5)
        ratio = 0.5;
    else if (ratio > 2.0)
        ratio = 2.0;

    /* while we have pending samples to read into our buffer */

    while (num_samples) {

        /* copy in as many samples as we have room for */

        int samples_to_copy = num_samples;

        if (samples_to_copy > cxt->inbuff_samples - cxt->head)
            samples_to_copy = cxt->inbuff_samples - cxt->head;

        memcpy (cxt->inbuff + cxt->head, samples, samples_to_copy * sizeof (cxt->inbuff [0]));
        num_samples -= samples_to_copy;
        samples += samples_to_copy;
        cxt->head += samples_to_copy;

        /* while there are enough samples to process (3 or 4 times the longest period), do so */

        while (cxt->tail >= cxt->longest && cxt->head - cxt->tail >= cxt->longest * (cxt->fast_mode ? 3 : 2)) {
            double process_ratio;
            int period;

            if (ratio != 1.0 || cxt->outsamples_error)
                period = cxt->fast_mode ? find_period_fast (cxt, cxt->inbuff + cxt->tail) :
                    find_period (cxt, cxt->inbuff + cxt->tail);
            else
                period = cxt->longest;

            /*
             * Once we have calculated the best-match period, there are 4 possible transformations
             * available to convert the input samples to output samples. Obviously we can simply
             * copy the samples verbatim (1:1). Standard TDHS provides algorithms for 2:1 and
             * 1:2 scaling, and I have created an obvious extension for 2:3 scaling. To achieve
             * intermediate ratios we maintain a "error" term (in samples) and use that here to
             * calculate the actual transformation to apply.
             */

            if (cxt->outsamples_error == 0.0)
                process_ratio = floor (ratio * 2.0 + 0.5) / 2.0;
            else if (cxt->outsamples_error > 0.0)
                process_ratio = floor (ratio * 2.0) / 2.0;
            else
                process_ratio = ceil (ratio * 2.0) / 2.0;

            if (process_ratio == 0.5) {
                merge_blocks (outbuf + out_samples, cxt->inbuff + cxt->tail,
                    cxt->inbuff + cxt->tail + period, period);
                cxt->outsamples_error += period - (period * 2.0 * ratio);
                out_samples += period;
                cxt->tail += period * 2;
            }
            else if (process_ratio == 1.0) {
                memcpy (outbuf + out_samples, cxt->inbuff + cxt->tail, period * 2 * sizeof (cxt->inbuff [0]));

                if (ratio != 1.0)
                    cxt->outsamples_error += (period * 2.0) - (period * 2.0 * ratio);
                else
                    cxt->outsamples_error = 0; /* if the ratio is 1.0, we can never cancel the error, so just do it now */

                out_samples += period * 2;
                cxt->tail += period * 2;
            }
            else if (process_ratio == 1.5) {
                memcpy (outbuf + out_samples, cxt->inbuff + cxt->tail, period * sizeof (cxt->inbuff [0]));
                merge_blocks (outbuf + out_samples + period, cxt->inbuff + cxt->tail + period,
                    cxt->inbuff + cxt->tail, period);
                memcpy (outbuf + out_samples + period * 2, cxt->inbuff + cxt->tail + period, period * sizeof (cxt->inbuff [0]));
                cxt->outsamples_error += (period * 3.0) - (period * 2.0 * ratio);
                out_samples += period * 3;
                cxt->tail += period * 2;
            }
            else if (process_ratio == 2.0) {
                merge_blocks (outbuf + out_samples, cxt->inbuff + cxt->tail,
                    cxt->inbuff + cxt->tail - period, period * 2);

                cxt->outsamples_error += (period * 2.0) - (period * ratio);
                out_samples += period * 2;
                cxt->tail += period;

                if (cxt->fast_mode) {
                    merge_blocks (outbuf + out_samples, cxt->inbuff + cxt->tail,
                        cxt->inbuff + cxt->tail - period, period * 2);

                    cxt->outsamples_error += (period * 2.0) - (period * ratio);
                    out_samples += period * 2;
                    cxt->tail += period;
                }
            }
            else
                fprintf (stderr, "stretchProcess: fatal programming error: process_ratio == %g\n", process_ratio);

            /* if there's another cascaded instance after this, pass the just stretched samples into that */

            if (cxt->next) {
                next_samples += stretchProcess (cxt->next, outbuf, out_samples / cxt->num_chans, output + next_samples * cxt->num_chans, next_ratio);
                out_samples = 0;
            }

            /* finally, left-justify the samples in the buffer leaving one longest period of history */

            int samples_to_move = cxt->inbuff_samples - cxt->tail + cxt->longest;

            memmove (cxt->inbuff, cxt->inbuff + cxt->tail - cxt->longest,
                samples_to_move * sizeof (cxt->inbuff [0]));

            cxt->head -= cxt->tail - cxt->longest;
            cxt->tail = cxt->longest;
        }
    }

    /*
     * This code is not strictly required, but will reduce latency, especially in the dual-instance case, by
     * always flushing all pending samples if no actual stretching is desired (i.e., ratio is 1.0 and there's
     * no error to compensate for). This case is more common now than previously because of the gap detection
     * and cascaded instances.
     */

    if (ratio == 1.0 && !cxt->outsamples_error && cxt->head != cxt->tail) {
        int samples_leftover = cxt->head - cxt->tail;

        if (cxt->next)
            next_samples += stretchProcess (cxt->next, cxt->inbuff + cxt->tail, samples_leftover / cxt->num_chans,
                output + next_samples * cxt->num_chans, next_ratio);
        else {
            memcpy (outbuf + out_samples, cxt->inbuff + cxt->tail, samples_leftover * sizeof (*output));
            out_samples += samples_leftover;
        }

        memmove (cxt->inbuff, cxt->inbuff + cxt->head - cxt->longest, cxt->longest * sizeof (cxt->inbuff [0]));
        cxt->head = cxt->tail = cxt->longest;
    }

    return cxt->next ? next_samples : out_samples / cxt->num_chans;
}  

/*
 * Flush any leftover samples out at normal speed. For cascaded dual instances this must be called
 * twice to completely flush, or simply call it until it returns zero samples. The maximum number
 * of samples that can be returned from each call of this function can be determined in advance with
 * stretch_output_capacity().
 */

int stretchFlush (Stretch *cxt, float *output)
{
    int samples_leftover = cxt->head - cxt->tail;
    int samples_flushed = 0;

    if (cxt->next) {
        if (samples_leftover)
            samples_flushed = stretchProcess (cxt->next, cxt->inbuff + cxt->tail, samples_leftover / cxt->num_chans, output, 1.0);

        if (!samples_flushed)
            samples_flushed = stretchFlush (cxt->next, output);
    }
    else {
        memcpy (output, cxt->inbuff + cxt->tail, samples_leftover * sizeof (*output));
        samples_flushed = samples_leftover / cxt->num_chans;
    }

    cxt->tail = cxt->head;
    memset (cxt->inbuff, 0, cxt->tail * sizeof (*cxt->inbuff));

    return samples_flushed;
}

/* free instance */

void stretchFree (Stretch *cxt)
{
    if (cxt) {
        free (cxt->calcbuff);
        free (cxt->results);
        free (cxt->inbuff);

        if (cxt->next) {
            stretchFree (cxt->next);
            free (cxt->intermediate);
        }

        free (cxt);
    }
}

/*
 * The pitch detection is done by finding the period that produces the
 * maximum value for the following correlation formula applied to two
 * consecutive blocks of the given period length:
 *
 *         sum of the absolute values of each sample in both blocks
 *   ---------------------------------------------------------------------
 *   sum of the absolute differences of each corresponding pair of samples
 *
 * This formula was chosen for two reasons.  First, it produces output values
 * that can directly compared regardless of the pitch period.  Second, the
 * numerator can be accumulated for successive periods, and only the
 * denominator need be completely recalculated.
 */

static int find_period (const Stretch *cxt, float *samples)
{
    float sum, factor, best_factor = 0;
    float *calcbuff = samples;
    int period, best_period;
    int i, j;

    period = best_period = cxt->shortest / cxt->num_chans;

    // convert stereo to mono, and accumulate sum for longest period

    if (cxt->num_chans == 2) {
        calcbuff = cxt->calcbuff;

        for (sum = i = j = 0; i < cxt->longest * 2; i += 2)
            sum += fabs (calcbuff [j++] = (samples [i] + samples [i+1]) / 2.0);
    }
    else
        for (sum = i = 0; i < cxt->longest; ++i)
            sum += fabs (calcbuff [i]) + fabs (calcbuff [i+cxt->longest]);

    // if silence return longest period, else calculate scaler based on largest sum

    if (!sum)
        return cxt->longest;

    /* accumulate sum for shortest period size */

    for (sum = i = 0; i < period; ++i)
        sum += fabs (calcbuff [i]) + fabs (calcbuff [i+period]);

    /* this loop actually cycles through all period lengths */

    while (1) {
        float *comp = calcbuff + period * 2;
        float *ref = calcbuff + period;
        float diff = 0;

        /* compute sum of absolute differences */

        while (ref != calcbuff)
            diff += fabs (*--ref - *--comp);

        /*
         * Here we calculate and store the resulting correlation
         * factor.  Note that we must watch for a difference of
         * zero, meaning a perfect match.  Also, for increased
         * precision using integer math, we scale the sum.
         */

        factor = (diff == 0.0) ? FLT_MAX : sum / diff;

        if (factor >= best_factor) {
            best_factor = factor;
            best_period = period;
        }

        /* see if we're done */

        if (period * cxt->num_chans == cxt->longest)
            break;

        /* update accumulating sum and current period */

        sum += fabs (calcbuff [period * 2]) + fabs (calcbuff [period * 2 + 1]);
        period++;
    }

    return best_period * cxt->num_chans;
}

/*
 * This pitch detection function is similar to find_period() above, except that it
 * is optimized for speed. The audio data corresponding to two maximum periods is
 * averaged 2:1 into the calculation buffer, and then the calulations are done
 * for every other period length. Because the time is essentially proportional to
 * both the number of samples and the number of period lengths to try, this scheme
 * can reduce the time by a factor approaching 4x. The correlation results on either
 * side of the peak are compared to calculate a more accurate center of the period.
 */

#ifndef M_E
#define M_E		2.7182818284590452354	/* e */
#endif

static int find_period_fast (const Stretch *cxt, const float *samples)
{
    float sum, best_factor = 0;
    int period, best_period;
    int i, j;

    best_period = period = cxt->shortest / (cxt->num_chans * 2);

    /* first step is compressing data 2:1 into calcbuff, and calculating maximum sum */

    if (cxt->num_chans == 2)
        for (sum = i = j = 0; i < cxt->longest * 2; i += 4)
            sum += fabs (cxt->calcbuff [j++] = (samples [i] + samples [i+1] + samples [i+2] + samples [i+3]) / 2.0);
    else
        for (sum = i = j = 0; i < cxt->longest * 2; i += 2)
            sum += fabs (cxt->calcbuff [j++] = (samples [i] + samples [i+1]) / 2.0);

    // if silence return longest period, else calculate scaler based on largest sum

    if (!sum)
        return cxt->longest;

    /* accumulate sum for shortest period */

    for (sum = i = 0; i < period; ++i)
        sum += fabs (cxt->calcbuff [i]) + fabs (cxt->calcbuff [i+period]);

    /* this loop actually cycles through all period lengths */

    while (1) {
        float *comp = cxt->calcbuff + period * 2;
        float *ref = cxt->calcbuff + period;
        float diff = 0.0;

        /* compute sum of absolute differences */

        while (ref != cxt->calcbuff)
            diff += fabs (*--ref - *--comp);

        /*
         * Here we calculate and store the resulting correlation
         * factor.  Note that we must watch for a difference of
         * zero, meaning a perfect match.  Also, for increased
         * precision using integer math, we scale the sum.
         */

        cxt->results [period] = diff == 0.0 ? FLT_MAX : sum / diff;

        if (cxt->results [period] >= best_factor) {    /* check if best yet */
            best_factor = cxt->results [period];
            best_period = period;
        }

        /* see if we're done */

        if (period * cxt->num_chans * 2 == cxt->longest)
            break;

        /* update accumulating sum and current period */

        sum += fabs (cxt->calcbuff [period * 2]) + fabs (cxt->calcbuff [period * 2 + 1]);
        period++;
    }

    if (best_period * cxt->num_chans * 2 != cxt->shortest && best_period * cxt->num_chans * 2 != cxt->longest) {
        float high_side_diff = cxt->results [best_period] - cxt->results [best_period+1];
        float low_side_diff = cxt->results [best_period] - cxt->results [best_period-1];

        if (low_side_diff > high_side_diff * M_E)
            best_period = best_period * 2 + 1;
        else if (high_side_diff > low_side_diff * M_E)
            best_period = best_period * 2 - 1;
        else
            best_period *= 2;
    }
    else
        best_period *= 2;           /* shortest or longest use as is */

    return best_period * cxt->num_chans;
}

/*
 * To combine the two periods into one, each corresponding pair of samples
 * are averaged with a linearly sliding scale.  At the beginning of the period
 * the first sample dominates, and at the end the second sample dominates.  In
 * this way the resulting block blends with the previous and next blocks.
 */

static void merge_blocks (float *output, float *input1, float *input2, int samples)
{
    int i;

    for (i = 0; i < samples; ++i)
        output [i] = ((input1 [i]) * (samples - i) + input2 [i] * i) / samples;
}
