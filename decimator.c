////////////////////////////////////////////////////////////////////////////
//                           **** DECIMATOR ****                          //
//                    Float to Integer Audio Decimation                   //
//                 Copyright (c) 2006 - 2025 David Bryant                 //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// decimator.c

#include "decimator.h"

// Initialize a decimator context with the specified characteristics. The returned context
// pointer is used for all subsequent calls to the decimator (and should not be dereferenced).
// A NULL return indicates an error. A gain parameter is provided as a convenience to avoid
// another pass through the data. The output bitcount parameter is obvious because it
// determines the magnitude of the output, but the output bytecountis less obvious (it allows
// for scenarios like 24-bit audio in 32-bit containers). Only little-endian support for now.
//
// Uses the decoupled form ("H(z)" in the referenced Wiki article) of the noise-shaping filter,
// which is refactored here from the direct form N(z). The Lame ATH curves for 44.1 Khz and 48
// KHz by Sebastian Gesemann are included, along with some others.
//
// see https://wiki.hydrogenaud.io/index.php?title=Noise_shaping

static void shaper_init (Biquad *f, float a0, float a1, float a2, float a3, float a4, float b1, float b2, float b3, float b4);

Decimate *decimateInit (int numChannels, int outputBits, int outputBytes, float outputGain, int sampleRate, int flags)
{
    Decimate *cxt = calloc (1, sizeof (Decimate));

    cxt->numChannels = numChannels;
    cxt->outputBytes = outputBytes;
    cxt->outputBits = outputBits;
    cxt->outputGain = outputGain;
    cxt->flags = flags;

    cxt->feedback = calloc (numChannels, sizeof (float));

    if (flags & DITHER_ENABLED) {
        int generator_bytes = numChannels * sizeof (uint32_t);
        unsigned char *seed = malloc (generator_bytes);
        uint32_t random = 0x31415926;

        cxt->tpdf_generators = (uint32_t *) seed;

        while (generator_bytes--) {
            *seed++ = random >> 24;
            random = ((random << 4) - random) ^ 1;
            random = ((random << 4) - random) ^ 1;
            random = ((random << 4) - random) ^ 1;
        }

        if (flags & DITHER_HIGHPASS)
            cxt->dither_type = -1;
        else if (flags & DITHER_LOWPASS)
            cxt->dither_type = 1;
        else if (flags & DITHER_FLAT)
            cxt->dither_type = 0;
    }

    if (flags & SHAPING_ENABLED) {
        int ch;

        cxt->noise_shapers = calloc (numChannels, sizeof (Biquad));

        for (ch = 0; ch < numChannels; ++ch) {
            if (flags & SHAPING_ATH_CURVE) {
                if (sampleRate == 32000)
                    shaper_init (cxt->noise_shapers + ch, 1.0, -0.780459, +0.569358, -0.348221, +0.466316, +0.950797, +0.282052, +0.004337, +1.76209e-5);
                else if (sampleRate == 44100)
                    shaper_init (cxt->noise_shapers + ch, 1.0, -1.1474, 0.5383, -0.3530, 0.3475, 1.0587, 0.0676, -0.6054, -0.2738);
                else if (sampleRate == 48000)
                    shaper_init (cxt->noise_shapers + ch, 1.0, -1.3344, 0.7455, -0.4602, 0.4363, 0.9030, 0.0116, -0.5853, -0.2571);
                else if (sampleRate == 88200)
                    shaper_init (cxt->noise_shapers + ch, 1.0, -2.150679, +2.1402057, -1.042712, +0.206838, +0.67433, +1.017047, +0.4028633, +0.098656);
                else if (sampleRate == 96000)
                    shaper_init (cxt->noise_shapers + ch, 1.0, -2.16994, +2.01986, -0.894857, +0.1557738, +0.517789, +1.1062189, +0.4825786, +0.244994);
                else
                    shaper_init (cxt->noise_shapers + ch, +1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);   // default to first-order for non-standard rates
            }
            else if (flags & SHAPING_1ST_ORDER)
                shaper_init (cxt->noise_shapers + ch, +1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            else if (flags & SHAPING_2ND_ORDER)
                shaper_init (cxt->noise_shapers + ch, +1.0, -2.0, +1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            else if (flags & SHAPING_3RD_ORDER)
                shaper_init (cxt->noise_shapers + ch, +1.0, -3.0, +3.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
    }

#ifdef ENABLE_THREADS
    if (numChannels > 1 && (flags & DECIMATE_MULTITHREADED))
        cxt->workers = workersInit (numChannels - 1);
#endif

    return cxt;
}

// Run the decimation process for the specified number of frames.
//
// This is the "non-interleaved" version of the decimator where the audio sample buffers for
// different channels are passed in as an array of float pointers and similarly the audio is
// output using an array of unsigned character pointers, one for each channel. There is also
// an "interleaved" version (see below).

#ifdef ENABLE_THREADS
static int decimateProcessSingleChanLE (void *ptr, void *sync_not_used);
#endif

static inline double tpdf_dither (uint32_t *generator, int type);

int decimateProcessLE (Decimate *cxt, const float *const *input, int numInputFrames, unsigned char *const *output)
{
#ifdef ENABLE_THREADS
    if (cxt->workers) {
        Decimate *worker_contexts = calloc (cxt->numChannels, sizeof (Decimate));
        int clipped_samples = 0, ch;

        for (ch = 0; ch < cxt->numChannels; ++ch) {
            Decimate *wcxt = worker_contexts + ch;

            *wcxt = *cxt;
            wcxt->input = input [ch] - 1;
            wcxt->numInputFrames = numInputFrames;
            wcxt->output = output [ch];
            wcxt->stride = 1;
            wcxt->clips = 0;
            wcxt->feedback_val = cxt->feedback [ch];
            wcxt->tpdf_generator = cxt->tpdf_generators [ch];
            memcpy (&wcxt->noise_shaper, cxt->noise_shapers + ch, sizeof (Biquad));

            workersEnqueueJob (cxt->workers, decimateProcessSingleChanLE, wcxt,
                ch < cxt->numChannels - 1 ? WaitForAvailableWorkerThread : DontUseWorkerThread);
        }

        workersWaitAllJobs (cxt->workers);

        for (ch = 0; ch < cxt->numChannels; ++ch) {
            Decimate *wcxt = worker_contexts + ch;

            cxt->feedback [ch] = wcxt->feedback_val;
            cxt->tpdf_generators [ch] = wcxt->tpdf_generator;
            memcpy (cxt->noise_shapers + ch, &wcxt->noise_shaper, sizeof (Biquad));
            clipped_samples += wcxt->clips;
        }

        free (worker_contexts);
        return clipped_samples;
    }
    else {
#endif
    float scaler = (1 << cxt->outputBits) / 2.0 * cxt->outputGain, codevalue;
    int pre_zeros = cxt->outputBytes - ((cxt->outputBits + 7) / 8);
    int32_t offset = (cxt->outputBits <= 8) * 128;
    int32_t highclip = (1 << (cxt->outputBits - 1)) - 1;
    int32_t lowclip = ~highclip;
    int leftshift = (24 - cxt->outputBits) % 8;
    int clipped_samples = 0, ch, i, j;

    for (i = 0; i < numInputFrames; ++i)
        for (ch = 0; ch < cxt->numChannels; ++ch) {
            float dither_value = (cxt->flags & DITHER_ENABLED) ? tpdf_dither (cxt->tpdf_generators + ch, cxt->dither_type) : 0.0;
            unsigned char *outp = output [ch] + i * cxt->outputBytes;
            int32_t outvalue;

            for (j = 0; j < pre_zeros; ++j)
                *outp++ = 0;

            codevalue = (input [ch] [i] * scaler) - cxt->feedback [ch];
            outvalue = floor (codevalue + dither_value + 0.5);

            if (cxt->flags & SHAPING_ENABLED)
                cxt->feedback [ch] = biquad_apply_sample (cxt->noise_shapers + ch, outvalue - codevalue);

            if (outvalue > highclip) {
                outvalue = highclip;
                clipped_samples++;
            }
            else if (outvalue < lowclip) {
                outvalue = lowclip;
                clipped_samples++;
            }

            *outp++ = outvalue = ((uint32_t) outvalue << leftshift) + offset;

            if (cxt->outputBits > 8) {
                *outp++ = outvalue >> 8;

                if (cxt->outputBits > 16)
                    *outp++ = outvalue >> 16;
            }
        }

    return clipped_samples;

#ifdef ENABLE_THREADS
    }
#endif
}

// This is the "interleaved" version of the decimator where the audio samples for different
// channels are passed in sequence in a single buffer. There is also a "non-interleaved"
// version for independent buffers, which is otherwise identical (see above).

int decimateProcessInterleavedLE (Decimate *cxt, const float *input, int numInputFrames, unsigned char *output)
{
#ifdef ENABLE_THREADS
    if (cxt->workers) {
        Decimate *worker_contexts = calloc (cxt->numChannels, sizeof (Decimate));
        int clipped_samples = 0, ch;

        for (ch = 0; ch < cxt->numChannels; ++ch) {
            Decimate *wcxt = worker_contexts + ch;

            *wcxt = *cxt;
            wcxt->input = input + ch - cxt->numChannels;
            wcxt->numInputFrames = numInputFrames;
            wcxt->output = output + (ch * cxt->outputBytes);
            wcxt->stride = cxt->numChannels;
            wcxt->clips = 0;
            wcxt->tpdf_generator = cxt->tpdf_generators [ch];
            wcxt->feedback_val = cxt->feedback [ch];
            memcpy (&wcxt->noise_shaper, cxt->noise_shapers + ch, sizeof (Biquad));

            workersEnqueueJob (cxt->workers, decimateProcessSingleChanLE, wcxt,
                ch < cxt->numChannels - 1 ? WaitForAvailableWorkerThread : DontUseWorkerThread);
        }

        workersWaitAllJobs (cxt->workers);

        for (ch = 0; ch < cxt->numChannels; ++ch) {
            Decimate *wcxt = worker_contexts + ch;

            cxt->tpdf_generators [ch] = wcxt->tpdf_generator;
            cxt->feedback [ch] = wcxt->feedback_val;
            memcpy (cxt->noise_shapers + ch, &wcxt->noise_shaper, sizeof (Biquad));
            clipped_samples += wcxt->clips;
        }

        free (worker_contexts);
        return clipped_samples;
    }
    else {
#endif
    float scaler = (1 << cxt->outputBits) / 2.0 * cxt->outputGain, codevalue;
    int pre_zeros = cxt->outputBytes - ((cxt->outputBits + 7) / 8);
    int32_t offset = (cxt->outputBits <= 8) * 128;
    int32_t highclip = (1 << (cxt->outputBits - 1)) - 1;
    int32_t lowclip = ~highclip;
    int leftshift = (24 - cxt->outputBits) % 8;
    int clipped_samples = 0, ch, i, j;

    for (i = 0; i < numInputFrames; ++i)
        for (ch = 0; ch < cxt->numChannels; ++ch) {
            float dither_value = (cxt->flags & DITHER_ENABLED) ? tpdf_dither (cxt->tpdf_generators + ch, cxt->dither_type) : 0.0;
            int32_t outvalue;

            for (j = 0; j < pre_zeros; ++j)
                *output++ = 0;

            codevalue = (*input++ * scaler) - cxt->feedback [ch];
            outvalue = floor (codevalue + dither_value + 0.5);

            if (cxt->flags & SHAPING_ENABLED)
                cxt->feedback [ch] = biquad_apply_sample (cxt->noise_shapers + ch, outvalue - codevalue);

            if (outvalue > highclip) {
                outvalue = highclip;
                clipped_samples++;
            }
            else if (outvalue < lowclip) {
                outvalue = lowclip;
                clipped_samples++;
            }

            *output++ = outvalue = ((uint32_t) outvalue << leftshift) + offset;

            if (cxt->outputBits > 8) {
                *output++ = outvalue >> 8;

                if (cxt->outputBits > 16)
                    *output++ = outvalue >> 16;
            }
        }

    return clipped_samples;

#ifdef ENABLE_THREADS
    }
#endif
}

#ifdef ENABLE_THREADS

static int decimateProcessSingleChanLE (void *ptr, void *sync_not_used)
{
    Decimate *cxt = ptr;
    float scaler = (1 << cxt->outputBits) / 2.0 * cxt->outputGain, codevalue;
    int pre_zeros = cxt->outputBytes - ((cxt->outputBits + 7) / 8);
    int stride_bytes = (cxt->stride - 1) * cxt->outputBytes;
    int32_t offset = (cxt->outputBits <= 8) * 128;
    int32_t highclip = (1 << (cxt->outputBits - 1)) - 1;
    int32_t lowclip = ~highclip;
    int leftshift = (24 - cxt->outputBits) % 8;
    int i, j;

    for (i = 0; i < cxt->numInputFrames; ++i) {
        float dither_value = (cxt->flags & DITHER_ENABLED) ? tpdf_dither (&cxt->tpdf_generator, cxt->dither_type) : 0.0;
        int32_t outvalue;

        for (j = 0; j < pre_zeros; ++j)
            *cxt->output++ = 0;

        codevalue = (*(cxt->input += cxt->stride) * scaler) - cxt->feedback_val;
        outvalue = floor (codevalue + dither_value + 0.5);

        if (cxt->flags & SHAPING_ENABLED)
            cxt->feedback_val = biquad_apply_sample (&cxt->noise_shaper, outvalue - codevalue);

        if (outvalue > highclip) {
            outvalue = highclip;
            cxt->clips++;
        }
        else if (outvalue < lowclip) {
            outvalue = lowclip;
            cxt->clips++;
        }

        *cxt->output++ = outvalue = ((uint32_t) outvalue << leftshift) + offset;

        if (cxt->outputBits > 8) {
            *cxt->output++ = outvalue >> 8;

            if (cxt->outputBits > 16)
                *cxt->output++ = outvalue >> 16;
        }

        cxt->output += stride_bytes;
    }

    return 0;
}

#endif

// Free all resources associated with the decimator context, including the context pointer
// itself. Do not use the context after this call.

void decimateFree (Decimate *cxt)
{
    if (cxt) {
        free (cxt->tpdf_generators);
        free (cxt->noise_shapers);
        free (cxt->feedback);

#ifdef ENABLE_THREADS
        if (cxt->workers)
            workersDeinit (cxt->workers);
#endif

        free (cxt);
    }
}

// Return a tpdf random value in the range: -1.0 <= n < 1.0
// type: -1: negative intersample correlation (HF boost)
//        0: no correlation (independent samples, flat spectrum)
//        1: positive intersample correlation (LF boost)

static inline double tpdf_dither (uint32_t *generator, int type)
{
    uint32_t random = *generator, first;

    random = ((random << 4) - random) ^ 1;
    random = ((random << 4) - random) ^ 1;
    first = type ? *generator ^ ((int32_t) type >> 31) : ~random;
    random = ((random << 4) - random) ^ 1;
    random = ((random << 4) - random) ^ 1;
    *generator = random = ((random << 4) - random) ^ 1;

    return (((first >> 1) + (random >> 1)) / 2147483648.0) - 1.0;
}

// Specify H(z) filter indirectly with N(z). Note, this is passed the actual noise-shaping
// transfer function, and so a0 needs to be 1.0. This function translates the filter to the
// H(z) form. The input to the resulting filter is the unfiltered quantization noise and the
// output, delayed by one sample, is subtracted from the [next] value to be quantized.

static void shaper_init (Biquad *f, float a0, float a1, float a2, float a3, float a4, float b1, float b2, float b3, float b4)
{
    BiquadCoefficients coeffs = { 0 };

    if (a0 != 1.0) {
        fprintf (stderr, "shaper_init() error: a0 = %g, should be one!\n", a0);
        exit (1);
    }

    coeffs.a0 = b1 - a1;
    coeffs.a1 = b2 - a2;
    coeffs.a2 = b3 - a3;
    coeffs.a3 = b4 - a4;

    coeffs.b1 = b1;
    coeffs.b2 = b2;
    coeffs.b3 = b3;
    coeffs.b4 = b4;

    biquad_init (f, &coeffs, 1.0);
}

// This function is provided as a complmentary function to the decimation, but is much simpler
// because the operation is naturally lossless and therefore does not require any state or
// context or channel awareness (i.e., samples are samples). A gain parameter is provided
// as a convenience.

void floatIntegersLE (unsigned char *input, float inputGain, int inputBits, int inputBytes, int inputStride, float *output, int numSamples)
{
    int post_skip = inputStride * inputBytes - ((inputBits + 7) / 8);

    input += inputBytes - ((inputBits + 7) / 8);

    if (inputBits <= 8) {
        float gain_factor = inputGain / 128.0;
        int i;

        for (i = 0; i < numSamples; ++i, input += post_skip)
            *output++ = ((int) *input++ - 128) * gain_factor;
    }
    else if (inputBits <= 16) {
        float gain_factor = inputGain / 32768.0;
        int i;

        for (i = 0; i < numSamples; ++i, input += post_skip) {
            int16_t value = *input++;
            value += *input++ << 8;
            *output++ = value * gain_factor;
        }
    }
    else if (inputBits <= 24) {
        float gain_factor = inputGain / 8388608.0;
        int i;

        for (i = 0; i < numSamples; ++i, input += post_skip) {
            int32_t value = *input++;
            value += *input++ << 8;
            value += (uint32_t) (signed char) *input++ << 16;
            *output++ = value * gain_factor;
        }
    }
}
