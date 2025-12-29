////////////////////////////////////////////////////////////////////////////
//                        **** AUDIO-STRETCH ****                         //
//                      Time Domain Harmonic Scaler                       //
//                    Copyright (c) 2022 David Bryant                     //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// stretch.h

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

#ifndef STRETCH_H
#define STRETCH_H

#include <stdint.h>

#if defined(PATH_WIDTH) && (PATH_WIDTH==64)
typedef double artsample_t;
#else
typedef float artsample_t;
#endif

#define MIN_PERIOD  24          /* minimum allowable pitch period */
#define MAX_PERIOD  2400        /* maximum allowable pitch period */

#define STRETCH_FAST_FLAG    0x1    // use "fast" version of period determination code
#define STRETCH_DUAL_FLAG    0x2    // cascade two instances (doubles usable ratio range)

typedef struct stretch {
    int num_chans, inbuff_samples, shortest, longest, tail, head, fast_mode;
    artsample_t *inbuff, *calcbuff, *results;
    double outsamples_error;

    struct stretch *next;
    artsample_t *intermediate;
} Stretch;

#ifdef __cplusplus
extern "C" {
#endif

Stretch *stretchInit (int shortest_period, int longest_period, int num_channels, int flags);
int stretchGetOutputCapacity (Stretch *cxt, int max_num_samples, double max_ratio);
int stretchProcess (Stretch *cxt, const artsample_t *samples, int num_samples, artsample_t *output, double ratio);
int stretchFlush (Stretch *cxt, artsample_t *output);
void stretchReset (Stretch *cxt);
void stretchFree (Stretch *cxt);

#ifdef __cplusplus
}
#endif

#endif

