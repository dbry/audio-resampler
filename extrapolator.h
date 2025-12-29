////////////////////////////////////////////////////////////////////////////
//                        **** EXTRAPOLATOR ****                          //
//                        LPC-Based Extrapolator                          //
//                   Copyright (c) 2025 David Bryant                      //
//                         All Rights Reserved                            //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// extrapolator.h

#ifndef EXTRAPOLATOR_H
#define EXTRAPOLATOR_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#if defined(PATH_WIDTH) && (PATH_WIDTH==64)
typedef double artsample_t;
#else
typedef float artsample_t;
#endif

#ifndef NCOEFFS
#define NCOEFFS 4
#endif

#ifndef MAXLOOPS
#define MAXLOOPS 100000
#endif

#ifdef __cplusplus
extern "C" {
#endif

double extrapolate_forward (artsample_t *values, int nvalues, int num_to_extrapolate);
double extrapolate_reverse (artsample_t *values, int nvalues, int num_to_extrapolate);

#ifdef __cplusplus
}
#endif

#endif

