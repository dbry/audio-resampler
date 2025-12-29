////////////////////////////////////////////////////////////////////////////
//                        **** EXTRAPOLATOR ****                          //
//                        LPC-Based Extrapolator                          //
//                   Copyright (c) 2025 David Bryant                      //
//                         All Rights Reserved                            //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// extrapolator.c

#include "extrapolator.h"

// Generalized extrapolation using LPC on a single array of floats. A set of LPC coefficients are generated from the
// the first 'nvalues' and then that filter is used to generate 'num_to_extrapolate' samples that are written just
// after the initial values in the array. See description of calc_lpc_coeffs() for a desciption of the return
// value and an explanation of "maxloops".
//
// Also see extrapolate_reverse() for a version with the data going backward.

static double calc_lpc_coeffs (const artsample_t *values, int nvalues, float *coeffs, int maxloops);

double extrapolate_forward (artsample_t *values, int nvalues, int num_to_extrapolate)
{
    artsample_t *src = values + nvalues - NCOEFFS;
    artsample_t *dst = values + nvalues;
    float coeffs [NCOEFFS];
    double quality;

    memset (values + nvalues, 0, num_to_extrapolate * sizeof (artsample_t));
    quality = calc_lpc_coeffs (values, nvalues, coeffs, MAXLOOPS);

    for (int i = 0; i < num_to_extrapolate; ++i) {
        double sum = 0.0;

        for (int c = 0; c < NCOEFFS; ++c)
            sum += src [c] * coeffs [NCOEFFS - c - 1];

        *dst++ = -sum;
        src++;
    }

    return quality;
}

// Similar to extraplate_forward() above except backwards. The 'values' parameter points to the
// first value PAST the actual samples to be extrapolated (backwards). In other words, it is valid
// for the pointer to reference an invalid address because only prior values will be accessed.

double extrapolate_reverse (artsample_t *values, int nvalues, int num_to_extrapolate)
{
    artsample_t *rbuffer = calloc (sizeof (artsample_t), nvalues + num_to_extrapolate);
    double quality;
    int i;

    for (i = 0; i < nvalues; ++i)
        rbuffer [i] = values [-1 - i];

    quality = extrapolate_forward (rbuffer, nvalues, num_to_extrapolate);

    for (i = nvalues; i < nvalues + num_to_extrapolate; ++i)
        values [-1 - i] = rbuffer [i];

    free (rbuffer);
    return quality;
}

// Calculate optimum LPC coefficients for predicting the specified values. Note the the first 'n' values are not
// predicted, where 'n' is the number of coefficients. This function operates iteratively with a decreasing
// step size and can become very slow in some situations. For this reason a parameter (maxloops) is included
// to limit the maximum delay. A good starting point is 100,000 loops.
//
// The initial set of calculated coefficients may be wildly unstable and therefore unsuitable for extrapolation
// (while being perfectly fine when not presented with novel input). Therefore an additional step is performed
// converting the LPC coefficients to PARCOR coefficients (the reflection coefficients) and these are checked
// for instabilities (i.e., coefficients outside +/- 1.0). If these are found then they are clipped and the
// coefficients are converted back to LPC coefficients.
//
// The return value is an evaluation of the "quality" of the coefficients. Because this was originally
// targeting lossless compression, it represents the approximate number of bits per sample that the
// prediction could remove from the data (i.e., used as a decorrelator).

static void lpc_to_parcor (const double *lpc, double *parcor, int ncoeffs);
static void parcor_to_lpc (const double *parcor, double *lpc, int ncoeffs);

static double calc_lpc_coeffs (const artsample_t *values, int nvalues, float *coeffs, int maxloops)
{
    double values_rms = 0.0, deltas_rms = 0.0, filter_rms_error = 0.0, step = 3.0 / (1 << 4), quality_factor = 20.0;
    double *sums = malloc (sizeof (double) * (nvalues - NCOEFFS));
    int nevals = nvalues - NCOEFFS, loops = 0, changes = 0;

    memset (coeffs, 0, sizeof (coeffs [0]) * NCOEFFS);      // start with a nothing filter

    // calculate the RMS level of the values and the value deltas

    for (int i = 0; i < nevals; ++i) {
        deltas_rms += (values [i + NCOEFFS] - values [i + NCOEFFS - 1]) *
                      (values [i + NCOEFFS] - values [i + NCOEFFS - 1]);

        values_rms += values [i + NCOEFFS] *
                      values [i + NCOEFFS];
    }

    if (values_rms == 0.0) {
        free (sums);
        return quality_factor;
    }

    filter_rms_error = values_rms;

    // loop until either the error is zero, or not being reduced, or we exceed a given max loop count

    while (filter_rms_error > 0.0 && (!maxloops || loops < maxloops)) {
        int tcoeff;

        for (int k = 0; k < nevals; ++k) {
            const artsample_t *kvalues = values + k;
            double zero_sum = 0.0;

            for (int c = 0; c < NCOEFFS; ++c)
                zero_sum += coeffs [NCOEFFS - c - 1] * kvalues [c];

            sums [k] = zero_sum + kvalues [NCOEFFS];
        }

        for (tcoeff = 0; loops++, tcoeff < NCOEFFS; tcoeff++) {
            double low_rms_error = 0.0, hi_rms_error = 0.0;

            for (int k = 0; k < nevals; ++k) {
                double delta = values [k + NCOEFFS - tcoeff - 1] * step;

                low_rms_error += (sums [k] - delta) * (sums [k] - delta);
                hi_rms_error += (sums [k] + delta) * (sums [k] + delta);
            }

            if (low_rms_error < filter_rms_error || hi_rms_error < filter_rms_error) {
                if (low_rms_error < hi_rms_error) {
                    filter_rms_error = low_rms_error;
                    coeffs [tcoeff] -= step;
                }
                else {
                    filter_rms_error = hi_rms_error;
                    coeffs [tcoeff] += step;
                }

                changes++;
                break;
            }
        }

        // if there's no improvement for this step size, reduce it or exit if we're done

        if (tcoeff == NCOEFFS) {
            if (step > 3.0 / (1 << 22))
                step *= 0.5;
            else
                break;
        }
    }

    free (sums);

    // assuming we actually generated something, check for and correct instabilities

    if (changes) {
        double dcoeffs [NCOEFFS], parcor [NCOEFFS];
        int outliers = 0;

        for (int i = 0; i < NCOEFFS; ++i)
            dcoeffs [i] = coeffs [i];

        lpc_to_parcor (dcoeffs, parcor, NCOEFFS);

        for (int i = 0; i < NCOEFFS; ++i)
            if (fabs (parcor [i]) > 0.9999) {
                parcor [i] = parcor [i] < 0.0 ? -0.9999 : 0.9999;
                outliers++;
            }

        if (outliers) {
            parcor_to_lpc (parcor, dcoeffs, NCOEFFS);

            for (int i = 0; i < NCOEFFS; ++i)
                coeffs [i] = dcoeffs [i];
        }
    }

    // check the filter again for effectiveness at predicting the given values

    filter_rms_error = 0.0;

    for (int k = 0; k < nevals; ++k) {
        const artsample_t *kvalues = values + k;
        double sum = 0.0;

        for (int c = 0; c < NCOEFFS; ++c)
            sum += coeffs [NCOEFFS - c - 1] * kvalues [c];

        filter_rms_error += (sum + kvalues [NCOEFFS]) * (sum + kvalues [NCOEFFS]);
    }

    // we might just decide to use a very simple filter or no filter if we get bad results

    if (deltas_rms < filter_rms_error && deltas_rms < values_rms) {
        memset (coeffs, 0, sizeof (coeffs [0]) * NCOEFFS);
        filter_rms_error = deltas_rms;
        coeffs [0] = -1.0;
    }
    else if (values_rms <= filter_rms_error) {
        memset (coeffs, 0, sizeof (coeffs [0]) * NCOEFFS);
        filter_rms_error = values_rms;
    }

    // calculate the quality and clip to 20.0 (which would be amazing)
    // note that it should not be <= 0.0 because we caught that case just above

    if (filter_rms_error != 0.0)
        quality_factor = (log (values_rms / filter_rms_error) * 0.5) / log (2.0);

    if (quality_factor > 20.0)
        quality_factor = 20.0;

    // sanitize the quality factor just in case, including for NaN

    if (quality_factor < 0.0 || quality_factor > 20.0 || quality_factor != quality_factor) {
        fprintf (stderr, "quality factor = %g\n", quality_factor);
        exit (1);
    }

    return quality_factor;
}

// Function to convert LPC coefficients to PARCOR coefficients

static void lpc_to_parcor (const double *lpc, double *parcor, int ncoeffs)
{
    double *temp_lpc = calloc (sizeof (double), ncoeffs);
    double *next_lpc = calloc (sizeof (double), ncoeffs);

    for (int i = 0; i < ncoeffs; ++i)
        temp_lpc [i] = lpc [i];

    for (int m = ncoeffs - 1; m >= 0; --m) {
        parcor [m] = temp_lpc [m];

        // if parcor[m] is close to 1 or -1, numerical instability can occur
        double denominator = 1.0 - (parcor [m] * parcor [m]);

        if (fabs (denominator) < 1e-6) {
            parcor [m] = parcor [m] < 0.0 ? -0.9999995 : 0.9999995;
            denominator = 1.0 - (parcor [m] * parcor [m]);
        }

        if (m > 0) {
            for (int i = 0; i < m; ++i)
                next_lpc [i] = (temp_lpc [i] - parcor [m] * temp_lpc [m - i - 1]) / denominator;

            for (int i = 0; i < m; ++i)
                temp_lpc [i] = next_lpc [i];
        }
    }

    free (temp_lpc);
    free (next_lpc);
}

// Function to convert PARCOR coefficients to LPC coefficients

static void parcor_to_lpc (const double *parcor, double *lpc, int ncoeffs)
{
    for (int i = 0; i < ncoeffs; i++) {
        lpc [i] = parcor [i];

        for (int j = 0; j < i / 2; j++) {
            double tmp = lpc [j];

            lpc [j] += parcor [i] * lpc [i - 1 - j];
            lpc [i - 1 - j] += parcor [i] * tmp;
        }

        if (i & 1)
            lpc [i >> 1] += lpc [i >> 1] * parcor [i];
    }
}
