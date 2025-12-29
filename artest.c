////////////////////////////////////////////////////////////////////////////
//                           **** ARTEST ****                             //
//                        Audio Resampling Tester                         //
//                 Copyright (c) 2006-2025 David Bryant                   //
//                         All Rights Reserved                            //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#define _USE_MATH_DEFINES
#include <math.h>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include "resampler.h"
#include "decimator.h"
#include "stretch.h"
#include "biquad.h"

// This program is used to benchmark and test the audio resampler library
// in two major ways:
//
// 1. Reverse resample and inverse paste for measuring resampling fidelity
// 2. Testing of interleaved and non-interleaved versions, and decimation
// 3. Benchmarking speed without concern for I/O latency

#define FIXED_RATIO     // define if resampler module has resampleFixedRatioInit()

static ResampleResult resampleProcessInterleavedSimulator (Resample *cxt, const artsample_t *input, int numInputFrames, artsample_t *output, int numOutputFrames, double ratio);
static ResampleResult resampleProcessAndFlushInterleavedSimulator (Resample *cxt, const artsample_t *input, int numInputFrames, artsample_t *output, int numOutputFrames, double ratio);
static int decimateProcessInterleavedSimulatorLE (Decimate *cxt, const artsample_t *input, int numInputFrames, unsigned char *output);
static void fill_buffer_with_noise (artsample_t *data, int count), fill_buffer_with_tone (artsample_t *data, int count, int chans, double freq);
static void fade_in (artsample_t *data, int count), fade_out (artsample_t *data, int count);

static const char *usage =
" Usage:    ARTEST [-options] [< infile.raw] [> outfile.raw]\n\n"
" Version:  %d-bit\n\n"
" Options:  -1|2|3|4    = quality presets, default = 3\n"
"           -b<num>     = inbuffer samples (default 4096)\n"
"           -c<num>     = number of channels (1-256, default 2)\n"
"           -n<num>     = number of seconds (1-3600, default 60)\n"
"           -h[<Hz>]    = use 1kc freq tone instead of white noise\n"
"           -s<Hz>      = source sample rate in Hz\n"
"                          (follow rate with 'k' for kHz)\n"
"           -d<Hz>      = destin sample rate in Hz\n"
"                          (follow rate with 'k' for kHz)\n"
"           -l<Hz>      = lowpass frequency in Hz\n"
"                          (follow rate with 'k' for kHz)\n"
"           -f<num>     = number of sinc filters (1-1024, default 380)\n"
"           -t<num>     = number of sinc taps (4-1024, default 380)\n"
"           -o<bits>    = change output file bitdepth (4-24)\n"
#ifdef FIXED_RATIO
"           -e          = calc exact filters / no interpolation\n"
#endif
"           -r          = read input from stdin\n"
"           -w<num>     = write specified raw stream to stdout\n"
"             <num>=1   = source stream\n"
"             <num>=2   = destination stream\n"
"             <num>=3   = decimated stream (from -o)\n"
"             <num>=4   = inverse resampled stream (from -i)\n"
"             <num>=5   = subtracted (error) stream (from -i)\n"
#ifdef ENABLE_THREADS
"           -m          = run resampler and decimator multi-threaded\n"
#endif
"           -i          = inverse-resample and compare to source\n"
"           -a          = do not fade-in and fade-out audio endpoints\n"
#ifdef ENABLE_EXTRAPOLATION
"           -x          = extrapolate audio endpoints to reduce glitches\n"
#endif
#if !defined(PATH_WIDTH) || (PATH_WIDTH==32)
"           -p          = precise, use doubles not floats for convolution\n"
#endif
"           -v          = test non-interleaved versions (if they exist)\n"
"           -o<bits>    = change output file bitdepth (4-24 or 32)\n\n";

typedef struct {
    uint64_t count, checksum;
    artsample_t min, max;
    double rms;
    int chans;
} Stats;

void update_stats (artsample_t *data, int samples, int chans, Stats *stats)
{
    int tsamples = samples * chans;

    stats->count += tsamples;
    stats->chans = chans;

    while (tsamples--) {
        stats->checksum = (stats->checksum * 3) + *(uint32_t*)(data);
        if (*data > stats->max) stats->max = *data;
        if (*data < stats->min) stats->min = *data;
        stats->rms += *data * *data;
        data++;
    }
}

char *display_stats (Stats *stats)
{
    static char string [128];

    sprintf (string, "count = %9llu, checksum = %016llx, range = %.7f to %.7f, RMS = %.2f dB", (unsigned long long) stats->count / stats->chans,
        (unsigned long long) stats->checksum, stats->min, stats->max, log10 (stats->rms / stats->count * 2.0) * 10.0);

    return string;
}

int main (int argc, char **argv)
{
    int inbuffer_samples = 4096, outbuffer_samples = 0, invbuffer_samples = 0, rembuffer_samples = 0, max_rembuffer_samples = 0;
    int dither = DITHER_HIGHPASS, noise_shaping = SHAPING_ATH_CURVE, multithreading = 0, fades = 1;
    int read_stdin = 0, write_stdout = 0, exact = 0, non_interleaved = 0, inv_resample = 0, buffers;
    int chans = 2, taps = 380, filters = 380, seconds = 60, outbits = 32, outbytes = 4;
    artsample_t *inbuffer = NULL, *outbuffer = NULL, *invbuffer = NULL, *rembuffer = NULL;
    int source_rate = 0, destin_rate = 0, lowpass_freq = 0;
    int flags = BLACKMAN_HARRIS | SUBSAMPLE_INTERPOLATE;
    Resample *resampler = NULL, *inv_resampler = NULL;
    Stats out_stats = { 0, 0, 1e20, -1e20, 0, chans };
    Stats in_stats = { 0, 0, 1e20, -1e20, 0, chans };
    Stats inv_stats = { 0, 0, 1e20, -1e20, 0, chans };
    Stats diff_stats = { 0, 0, 1e20, -1e20, 0, chans };
    unsigned char *decimate_buffer;
    Decimate *decimator = NULL;
    uint64_t dec_checksum = 0;
    int clipped_samples = 0;
    double ratio, inv_ratio;
    uint64_t out_bytes = 0;
    double tone_freq = 0.0;

    if (argc < 3) {
        fprintf (stderr, usage, sizeof (artsample_t) * 8);
        return 0;
    }

    // loop through command-line arguments

    while (--argc) {
#if defined (_WIN32)
        if ((**++argv == '-' || **argv == '/') && (*argv)[1])
#else
        if ((**++argv == '-') && (*argv)[1])
#endif
            while (*++*argv)
                switch (**argv) {

                    case '1':
                        filters = taps = 48;
                        break;

                    case '2':
                        filters = 320;
                        taps = 156;
                        break;

                    case '3':
                        filters = taps = 380;
                        break;

                    case '4':
                        filters = taps = 988;
                        break;

                    case 'a':
                        fades = 0;
                        break;
#ifdef FIXED_RATIO
                    case 'e':
                        exact = 1;
                        break;
#endif
                    case 'r':
                        read_stdin = 1;
                        break;

                    case 'w':
                        write_stdout = strtod (++*argv, argv);

                        if (write_stdout < 0 || write_stdout > 5) {
                            fprintf (stderr, "\nwritten stream must be 0 - 5!\n");
                            return 1;
                        }

                        --*argv;
                        break;

                    case 'i':
                        inv_resample = 1;
                        break;

                    case 'v':
                        non_interleaved = 1;
                        break;
#ifdef ENABLE_EXTRAPOLATION
                    case 'x':
                        flags |= EXTRAPOLATE_ENDPOINTS;
                        break;
#endif
                    case 'p':
                        flags |= EXTEND_CONVOLUTION_MATH;
                        break;
#ifdef ENABLE_THREADS
                    case 'm':
                        multithreading = DECIMATE_MULTITHREADED;
                        flags |= RESAMPLE_MULTITHREADED;
                        break;
#endif
                    case 'H': case 'h':
                        {
                            double freq = strtod (++*argv, argv);

                            if ((**argv & 0xdf) == 'K')
                                freq *= 1000.0;
                            else
                                --*argv;

                            tone_freq = freq == 0.0 ? 1000.0 : freq;
                        }

                        break;

                    case 'S': case 's':
                        {
                            double rate = strtod (++*argv, argv);

                            if ((**argv & 0xdf) == 'K')
                                rate *= 1000.0;
                            else
                                --*argv;

                            source_rate = rate;
                        }

                        break;

                    case 'D': case 'd':
                        {
                            double freq = strtod (++*argv, argv);

                            if ((**argv & 0xdf) == 'K')
                                freq *= 1000.0;
                            else
                                --*argv;

                            destin_rate = freq;
                        }

                        break;

                    case 'L': case 'l':
                        {
                            double freq = strtod (++*argv, argv);

                            if ((**argv & 0xdf) == 'K')
                                freq *= 1000.0;
                            else
                                --*argv;

                            lowpass_freq = freq;
                            flags |= INCLUDE_LOWPASS;
                        }

                        break;

                    case 'B': case 'b':
                        inbuffer_samples = strtod (++*argv, argv);

                        if (inbuffer_samples < 256 || inbuffer_samples > 65536) {
                            fprintf (stderr, "\ninbuffer samples must be 256 - 65536!\n");
                            return 1;
                        }

                        --*argv;
                        break;

                    case 'C': case 'c':
                        chans = strtod (++*argv, argv);

                        if (chans < 1 || chans > 256) {
                            fprintf (stderr, "\nnum of chans must be 1 - 256!\n");
                            return 1;
                        }

                        --*argv;
                        break;

                    case 'F': case 'f':
                        filters = strtod (++*argv, argv);

                        if (filters < 1 || filters > 1024) {
                            fprintf (stderr, "\nnum of filters must be 1 - 1024!\n");
                            return 1;
                        }

                        --*argv;
                        break;

                    case 'N': case 'n':
                        seconds = strtod (++*argv, argv);

                        if (seconds < 1 || seconds > 36000) {
                            fprintf (stderr, "\nnumber of seconds must be 1 - 36000!\n");
                            return 1;
                        }

                        --*argv;
                        break;

                    case 'O': case 'o':
                        outbits = strtod (++*argv, argv);

                        if (outbits != 32 && (outbits < 4 || outbits > 24)) {
                            fprintf (stderr, "\noutbits must be 4 - 24 (for integer) or 32 (for float)!\n");
                            return 1;
                        }

                        outbytes = (outbits + 7) / 8;
                        --*argv;
                        break;

                    case 'T': case 't':
                        taps = strtod (++*argv, argv);

                        if ((taps & 3) || taps < 4 || taps > 1024) {
                            fprintf (stderr, "\nnum of taps must be 4 - 1024 and a multiple of 4!\n");
                            return 1;
                        }

                        --*argv;
                        break;

                    default:
                        fprintf (stderr, "\nillegal option: %c !\n", **argv);
                        return 1;
                }
        else {
            fprintf (stderr, "\nextra unknown argument: %s !\n", *argv);
            return 1;
        }
    }

    if (!(destin_rate && source_rate) || !filters || !taps || !chans) {
        fprintf (stderr, "\nsomething is missing!\n\n");
        exit (1);
    }

    if ((flags & INCLUDE_LOWPASS) && !lowpass_freq && !exact) {
        fprintf (stderr, "\nspecify lowpass frequency, auto lowpass can only be used with exact resampling (-e)!\n\n");
        exit (1);
    }

#if defined(_WIN32)
        _setmode (_fileno (stdout), O_BINARY);
        _setmode (_fileno (stdin), O_BINARY);
#endif

    ratio = (double) destin_rate / source_rate;
    outbuffer_samples = floor ((inbuffer_samples + taps / 2) * ratio + 10);
    inbuffer = malloc (inbuffer_samples * chans * sizeof (artsample_t));
    outbuffer = malloc (outbuffer_samples * chans * sizeof (artsample_t));
    buffers = ceil ((double) seconds * source_rate / inbuffer_samples);

    if (inv_resample) {
        invbuffer_samples = floor ((outbuffer_samples + taps / 2) / ratio + 10);
        invbuffer = malloc (invbuffer_samples * chans * sizeof (artsample_t));
        inv_ratio = (double) source_rate / destin_rate;
    }

    if (ratio != 1.0 || lowpass_freq) {
#ifdef FIXED_RATIO
        if (exact) {
            int actual_filters;

            resampler = resampleFixedRatioInit (chans, taps, filters, source_rate, destin_rate, lowpass_freq, flags);
            actual_filters = resampleGetNumFilters (resampler);

            if (resampleGetLowpassRatio (resampler) == 1.0)
                fprintf (stderr, "w1 --> w2: %d %d-tap fixed-ratio sinc resampler%s, no lowpass, %s interpolation\n", actual_filters, taps,
                    actual_filters > 1 ? "s" : "", resampleInterpolationUsed (resampler) ? "with" : "no");
            else
                fprintf (stderr, "w1 --> w2: %d %d-tap fixed-rate sinc resampler%s with lowpass at %lu Hz, %s interpolation\n", actual_filters, taps,
                    actual_filters > 1 ? "s" : "", (unsigned long)(resampleGetLowpassRatio (resampler) * source_rate / 2.0), resampleInterpolationUsed (resampler) ? "with" : "no");

            if (inv_resample) {
                inv_resampler = resampleFixedRatioInit (chans, taps, filters, destin_rate, source_rate, lowpass_freq, flags);
                actual_filters = resampleGetNumFilters (inv_resampler);

                if (resampleGetLowpassRatio (inv_resampler) == 1.0)
                    fprintf (stderr, "w2 --> w4: %d %d-tap fixed-ratio sinc resampler%s, no lowpass, %s interpolation\n", actual_filters, taps,
                        actual_filters > 1 ? "s" : "", resampleInterpolationUsed (inv_resampler) ? "with" : "no");
                else
                    fprintf (stderr, "w2 --> w4: %d %d-tap fixed-rate sinc resampler%s with lowpass at %lu Hz, %s interpolation\n", actual_filters, taps,
                        actual_filters > 1 ? "s" : "", (unsigned long)(resampleGetLowpassRatio (inv_resampler) * destin_rate / 2.0), resampleInterpolationUsed (inv_resampler) ? "with" : "no");
            }

            inv_ratio = ratio = 0.0;
        }
        else
#endif
        {
            resampler = resampleInit (chans, taps, filters, lowpass_freq * 2.0 / source_rate, flags);

            if (resampleGetLowpassRatio (resampler) == 1.0)
                fprintf (stderr, "w1 --> w2: %d %d-tap fixed-ratio sinc resampler%s, no lowpass, %s interpolation\n", filters, taps,
                    filters > 1 ? "s" : "", resampleInterpolationUsed (resampler) ? "with" : "no");
            else
                fprintf (stderr, "w1 --> w2: %d %d-tap fixed-rate sinc resampler%s with lowpass at %lu Hz, %s interpolation\n", filters, taps,
                    filters > 1 ? "s" : "", (unsigned long)(resampleGetLowpassRatio (resampler) * source_rate / 2.0), resampleInterpolationUsed (resampler) ? "with" : "no");

            if (inv_resample) {
                inv_resampler = resampleInit (chans, taps, filters, lowpass_freq * 2.0 / destin_rate, flags);

                if (resampleGetLowpassRatio (inv_resampler) == 1.0)
                    fprintf (stderr, "w2 --> w4: %d %d-tap fixed-ratio sinc resampler%s, no lowpass, %s interpolation\n", filters, taps,
                        filters > 1 ? "s" : "", resampleInterpolationUsed (inv_resampler) ? "with" : "no");
                else
                    fprintf (stderr, "w2 --> w4: %d %d-tap fixed-rate sinc resampler%s with lowpass at %lu Hz, %s interpolation\n", filters, taps,
                        filters > 1 ? "s" : "", (unsigned long)(resampleGetLowpassRatio (inv_resampler) * destin_rate / 2.0), resampleInterpolationUsed (inv_resampler) ? "with" : "no");
            }
        }

        resampleAdvancePosition (resampler, taps / 2.0);

        if (inv_resample)
            resampleAdvancePosition (inv_resampler, taps / 2.0);
    }

    if (outbits != 32) {
        decimator = decimateInit (chans, outbits, outbytes, 1.0, destin_rate, dither | noise_shaping | multithreading);
        decimate_buffer = malloc (outbuffer_samples * chans * outbytes);
    }

    // everything is set up, so generate (or read) the requested data and process it

    for (int bi = 0; (bi < buffers || read_stdin) && inbuffer_samples; ++bi) {
        ResampleResult res, inv_res;

        if (read_stdin)
            inbuffer_samples = fread (inbuffer, sizeof (artsample_t) * chans, inbuffer_samples, stdin);
        else {
            if (tone_freq)
                fill_buffer_with_tone (inbuffer, inbuffer_samples, chans, tone_freq / source_rate);
            else
                fill_buffer_with_noise (inbuffer, inbuffer_samples * chans);

            if (fades) {
                if (bi == 0)
                    fade_in (inbuffer, inbuffer_samples * chans);
                else if (bi == buffers - 1)
                    fade_out (inbuffer, inbuffer_samples * chans);
            }
        }

        if (!inbuffer_samples)
            break;

        update_stats (inbuffer, inbuffer_samples, chans, &in_stats);

        if (write_stdout == 1)
            fwrite (inbuffer, sizeof (artsample_t) * chans, inbuffer_samples, stdout);

        if (!resampler) {
            memcpy (outbuffer, inbuffer, inbuffer_samples * chans * sizeof (artsample_t));
            res.output_generated = res.input_used = inbuffer_samples;
        }
        else if (bi < buffers - 1)
            res = non_interleaved ?
                resampleProcessInterleavedSimulator (resampler, inbuffer, inbuffer_samples, outbuffer, outbuffer_samples, ratio) :
                resampleProcessInterleaved (resampler, inbuffer, inbuffer_samples, outbuffer, outbuffer_samples, ratio);
        else
            res = non_interleaved ?
                resampleProcessAndFlushInterleavedSimulator (resampler, inbuffer, inbuffer_samples, outbuffer, outbuffer_samples, ratio) :
                resampleProcessAndFlushInterleaved (resampler, inbuffer, inbuffer_samples, outbuffer, outbuffer_samples, ratio);

        if (res.input_used != inbuffer_samples || res.output_generated == outbuffer_samples) {
            fprintf (stderr, "fatal error in resample results!\n");
            exit (1);
        }

        update_stats (outbuffer, res.output_generated, chans, &out_stats);

        if (write_stdout == 2)
            fwrite (outbuffer, sizeof (artsample_t) * chans, res.output_generated, stdout);

        if (inv_resample) {
            if (!inv_resampler) {
                memcpy (invbuffer, outbuffer, res.output_generated * chans * sizeof (artsample_t));
                inv_res.output_generated = inv_res.input_used = res.output_generated;
            }
            else if (bi < buffers - 1) {
                inv_res = non_interleaved ?
                    resampleProcessInterleavedSimulator (inv_resampler, outbuffer, res.output_generated, invbuffer, invbuffer_samples, inv_ratio) :
                    resampleProcessInterleaved (inv_resampler, outbuffer, res.output_generated, invbuffer, invbuffer_samples, inv_ratio);
            }
            else
                inv_res = non_interleaved ?
                    resampleProcessAndFlushInterleavedSimulator (inv_resampler, outbuffer, res.output_generated, invbuffer, invbuffer_samples, inv_ratio) :
                    resampleProcessAndFlushInterleaved (inv_resampler, outbuffer, res.output_generated, invbuffer, invbuffer_samples, inv_ratio);

            // this handles the case where, at the end, we've generated one or two more samples than we started with because of rounding
            if (inv_res.output_generated > rembuffer_samples + inbuffer_samples) {
                fprintf (stderr, "info: we generated %d extra sample(s) on round-trip resample\n", inv_res.output_generated - (rembuffer_samples + inbuffer_samples));
                inv_res.output_generated = rembuffer_samples + inbuffer_samples;
            }
            else if (bi == buffers - 1 && inv_res.output_generated < rembuffer_samples + inbuffer_samples)
                fprintf (stderr, "info: we generated %d fewer sample(s) on round-trip resample\n", rembuffer_samples + inbuffer_samples - inv_res.output_generated);

            if (inv_res.input_used != res.output_generated || inv_res.output_generated == invbuffer_samples) {
                fprintf (stderr, "fatal error in inverse resample results!\n");
                exit (1);
            }

            update_stats (invbuffer, inv_res.output_generated, chans, &inv_stats);

            if (write_stdout == 4)
                fwrite (invbuffer, sizeof (artsample_t) * chans, inv_res.output_generated, stdout);

            int next_rembuffer_samples = rembuffer_samples + inbuffer_samples - inv_res.output_generated;

            if (next_rembuffer_samples > max_rembuffer_samples) {
                rembuffer = realloc (rembuffer, next_rembuffer_samples * chans * sizeof (artsample_t));
                max_rembuffer_samples = next_rembuffer_samples;
            }

            int rembuffer_diff_count = 0, inbuffer_diff_count = 0, rembuffer_move_count = 0, inbuffer_move_count = 0;
            artsample_t *rembuffer_src_ptr = rembuffer;
            artsample_t *rembuffer_dst_ptr = rembuffer;
            artsample_t *invbuffer_ptr = invbuffer;
            artsample_t *inbuffer_ptr = inbuffer;

            if (inv_res.output_generated >= rembuffer_samples) {
                rembuffer_diff_count = rembuffer_samples * chans;
                inbuffer_diff_count = (inv_res.output_generated - rembuffer_samples) * chans;
                inbuffer_move_count = next_rembuffer_samples * chans;
            }
            else {
                rembuffer_diff_count = inv_res.output_generated * chans;
                rembuffer_move_count = (rembuffer_samples - inv_res.output_generated) * chans;
                inbuffer_move_count = inbuffer_samples * chans;
            }

            while (rembuffer_diff_count--)
                *invbuffer_ptr++ -= *rembuffer_src_ptr++;

            while (inbuffer_diff_count--)
                *invbuffer_ptr++ -= *inbuffer_ptr++;

            while (rembuffer_move_count--)
                *rembuffer_dst_ptr++ = *rembuffer_src_ptr++;

            while (inbuffer_move_count--)
                *rembuffer_dst_ptr++ = *inbuffer_ptr++;

            rembuffer_samples = next_rembuffer_samples;
            update_stats (invbuffer, inv_res.output_generated, chans, &diff_stats);

            if (write_stdout == 5)
                fwrite (invbuffer, sizeof (artsample_t) * chans, inv_res.output_generated, stdout);
        }

        if (outbits != 32) {
            int bytes_to_check;

            if (non_interleaved)
                clipped_samples += decimateProcessInterleavedSimulatorLE (decimator, outbuffer, res.output_generated, decimate_buffer);
            else
                clipped_samples += decimateProcessInterleavedLE (decimator, outbuffer, res.output_generated, decimate_buffer);

            out_bytes += res.output_generated * chans * outbytes;

            if (write_stdout == 3)
                fwrite (decimate_buffer, chans * outbytes, res.output_generated, stdout);

            bytes_to_check = chans * outbytes * res.output_generated;

            for (int si = 0; si < bytes_to_check; ++si)
                dec_checksum = (dec_checksum * 3) + decimate_buffer [si];
        }
    }

    resampleFree (inv_resampler);
    resampleFree (resampler);
    decimateFree (decimator);
    free (rembuffer);
    free (outbuffer);
    free (invbuffer);
    free (inbuffer);

    fprintf (stderr, "\n");

    fprintf (stderr, "   input (-w1): %s\n", display_stats (&in_stats));
    fprintf (stderr, "  output (-w2): %s\n", display_stats (&out_stats));

    if (inv_resample) {
        fprintf (stderr, " inverse (-w4): %s\n", display_stats (&inv_stats));
        fprintf (stderr, "    diff (-w5): %s\n", display_stats (&diff_stats));
    }

    if (out_bytes)
        fprintf (stderr, "decimate (-w3): count = %9llu, checksum = %016llx, clipped samples = %d\n",
            (unsigned long long) out_bytes, (unsigned long long) dec_checksum, clipped_samples);

    fprintf (stderr, "\n");
    return 0;
}

// simulate the interleaved decimator API but use the non-interleaved version

static int decimateProcessInterleavedSimulatorLE (Decimate *cxt, const artsample_t *input, int numInputFrames, unsigned char *output)
{
    artsample_t **input_array = calloc (sizeof (artsample_t*), cxt->numChannels);
    unsigned char **output_array = calloc (sizeof (unsigned char*), cxt->numChannels);

    for (int c = 0; c < cxt->numChannels; ++c) {
        input_array [c] = malloc (numInputFrames * sizeof (artsample_t));
        output_array [c] = malloc (numInputFrames * cxt->outputBytes);
    }

    for (int i = 0; i < numInputFrames; ++i)
        for (int c = 0; c < cxt->numChannels; ++c)
            input_array [c] [i] = *input++;

    int res = decimateProcessLE (cxt, (const artsample_t * const *) input_array, numInputFrames, output_array);

    for (int i = 0; i < numInputFrames; ++i)
        for (int c = 0; c < cxt->numChannels; ++c) {
            unsigned char *cp = output_array [c] + i * cxt->outputBytes;

            for (int b = 0; b < cxt->outputBytes; ++b)
                *output++ = *cp++;
        }

    for (int c = 0; c < cxt->numChannels; ++c) {
        free (input_array [c]);
        free (output_array [c]);
    }

    free (input_array);
    free (output_array);

    return res;
}

#if 1   // this version of the simulator is for determining whether the two versions generate the same values,
        // but hurts performance so should not be used to evaluate the relative performance of the versions.

static ResampleResult resampleProcessInterleavedSimulator (Resample *cxt, const artsample_t *input, int numInputFrames, artsample_t *output, int numOutputFrames, double ratio)
{
    artsample_t **input_array = calloc (sizeof (artsample_t*), cxt->numChannels);
    artsample_t **output_array = calloc (sizeof (artsample_t*), cxt->numChannels);

    for (int c = 0; c < cxt->numChannels; ++c) {
        input_array [c] = input ? malloc (numInputFrames * sizeof (artsample_t)) : NULL;
        output_array [c] = malloc (numOutputFrames * sizeof (artsample_t));
    }

    for (int i = 0; i < numInputFrames; ++i)
        for (int c = 0; c < cxt->numChannels; ++c)
            input_array [c] [i] = *input++;

    ResampleResult res = resampleProcess (cxt, (const artsample_t * const *) input_array, numInputFrames, output_array, numOutputFrames, ratio);

    for (int i = 0; i < res.output_generated; ++i)
        for (int c = 0; c < cxt->numChannels; ++c)
            *output++ = output_array [c] [i];

    for (int c = 0; c < cxt->numChannels; ++c) {
        free (input_array [c]);
        free (output_array [c]);
    }

    free (input_array);
    free (output_array);

    return res;
}

#else   // This version of the simulator is for evaluating the relative performance of the two versions because it does almost
        // nothing, but the resulting audio will have discontinuities and not match, but from run to run it should be consistent.

static ResampleResult resampleProcessInterleavedSimulator (Resample *cxt, const artsample_t *input, int numInputFrames, artsample_t *output, int numOutputFrames, double ratio)
{
    artsample_t *input_array [cxt->numChannels];
    artsample_t *output_array [cxt->numChannels];

    for (int c = 0; c < cxt->numChannels; ++c) {
        input_array [c] = (artsample_t* const) (input + numInputFrames * c);
        output_array [c] = (artsample_t* const) (output + numOutputFrames * c);
    }

    for (int i = 0; i < numInputFrames; ++i)
        for (int c = 0; c < cxt->numChannels; ++c)
            input_array [c] [i] = *input++;

    return resampleProcess (cxt, (const artsample_t * const *) input_array, numInputFrames, output_array, numOutputFrames, ratio);
}

#endif

static ResampleResult resampleProcessAndFlushInterleavedSimulator (Resample *cxt, const artsample_t *input, int numInputFrames, artsample_t *output, int numOutputFrames, double ratio)
{
    artsample_t **input_array = calloc (sizeof (artsample_t*), cxt->numChannels);
    artsample_t **output_array = calloc (sizeof (artsample_t*), cxt->numChannels);

    for (int c = 0; c < cxt->numChannels; ++c) {
        input_array [c] = malloc (numInputFrames * sizeof (artsample_t));
        output_array [c] = malloc (numOutputFrames * sizeof (artsample_t));
    }

    for (int i = 0; i < numInputFrames; ++i)
        for (int c = 0; c < cxt->numChannels; ++c)
            input_array [c] [i] = *input++;

    ResampleResult res = resampleProcessAndFlush (cxt, input ? (const artsample_t * const *) input_array : NULL, numInputFrames, output_array, numOutputFrames, ratio);

    for (int i = 0; i < res.output_generated; ++i)
        for (int c = 0; c < cxt->numChannels; ++c)
            *output++ = output_array [c] [i];

    for (int c = 0; c < cxt->numChannels; ++c) {
        free (input_array [c]);
        free (output_array [c]);
    }

    free (input_array);
    free (output_array);

    return res;
}

// fill buffer with +/-0.5 white noise

static void fill_buffer_with_noise (artsample_t *data, int count)
{
    static uint64_t random = 0x3141592653589793;

    while (count--) {
        random = ((random << 4) - random) ^ 1;
        random = ((random << 4) - random) ^ 1;
        random = ((random << 4) - random) ^ 1;
        *data++ = (int32_t)(random >> 32) / 4294967296.0;
    }
}

// fill buffer with +/-0.5 tone at specified frequency

static void fill_buffer_with_tone (artsample_t *data, int count, int chans, double freq)
{
    static double phase_angle;
    double chan_offset;

    if (chans > 2)
        chan_offset = 2.0 * M_PI / chans;
    else
        chan_offset = M_PI / 2.0;

    while (count--) {
        *data++ = sin (phase_angle += 2 * M_PI * freq) * 0.5;

        for (int c = 1; c < chans; ++c)
            *data++ = sin (phase_angle + chan_offset * c) * 0.5;
    }
}

static void fade_in (artsample_t *data, int count)
{
    int zcount = count / 4;
    int fcount = count - zcount;

    for (int i = 0; i < zcount; ++i)
        *data++ = 0.0;

    for (int i = 0; i < fcount; ++i)
        *data++ *= (cos ((fcount - i) * M_PI / fcount) + 1.0) / 2.0;
}

static void fade_out (artsample_t *data, int count)
{
    int zcount = count / 4;
    int fcount = count - zcount;

    for (int i = 0; i < fcount; ++i)
        *data++ *= (cos (i * M_PI / fcount) + 1.0) / 2.0;

    for (int i = 0; i < zcount; ++i)
        *data++ = 0.0;
}
