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
#include <math.h>

#ifndef NCOEFFS
#define NCOEFFS 4
#endif

#ifndef MAXLOOPS
#define MAXLOOPS 100000
#endif

#ifdef __cplusplus
extern "C" {
#endif

double extrapolate_forward (float *values, int nvalues, int num_to_extrapolate);
double extrapolate_reverse (float *values, int nvalues, int num_to_extrapolate);

#ifdef __cplusplus
}
#endif

#endif

