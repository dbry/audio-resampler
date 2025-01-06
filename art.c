////////////////////////////////////////////////////////////////////////////
//                            **** ART ****                               //
//                        Audio Resampling Tool                           //
//                 Copyright (c) 2006-2023 David Bryant                   //
//                         All Rights Reserved                            //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <math.h>

#include "resampler.h"
#include "decimator.h"
#include "biquad.h"

#define VERSION         0.3

#define IS_BIG_ENDIAN (*(uint16_t *)"\0\xff" < 0x0100)

static const char *sign_on = "\n"
" ART  Audio Resampling Tool  Version %.1f\n"
" Copyright (c) 2006 - 2025 David Bryant.\n\n";

static const char *usage =
" Usage:     ART [-options] infile.wav outfile.wav\n\n"
" Options:  -1|2|3|4    = quality presets, default = 3\n"
"           -r<Hz>      = resample to specified rate\n"
"           -g<dB>      = apply gain (default = 0 dB)\n"
"           -s<degrees> = add specified phase shift (+/-360 degrees)\n"
"           -l<Hz>      = specify alternate lowpass frequency\n"
"           -f<num>     = number of sinc filters (2-1024)\n"
"           -t<num>     = number of sinc taps (4-1024, multiples of 4)\n"
"           -o<bits>    = change output file bitdepth (4-24 or 32)\n"
"           -d<sel>     = override default dither (which is HP tpdf):\n"
"                           sel = 0 for no dither\n"
"                           sel = 1 for flat tpdf dither\n"
"                           sel = 2 for LP tpdf dither\n"
"           -n<sel>     = override default noise-shaping (which is ATH)\n"
"                           sel = 0 for no noise-shaping\n"
"                           sel = 1 for 1st-order shaping\n"
"                           sel = 2 for 2nd-order shaping\n"
"                           sel = 3 for 3rd-order shaping\n"
"           -i<bytes>   = test non-interleaved decimator (debug only)\n"
"           -0          = disable sinc resampler completely (debug only)\n"
"           -b          = Blackman-Harris windowing (best stopband)\n"
"           -h          = Hann windowing (fastest transition)\n"
"           -p          = pre/post filtering (cascaded biquads)\n"
"           -q          = quiet mode (display errors only)\n"
"           -v          = verbose (display lots of info)\n"
"           -y          = overwrite outfile if it exists\n\n"
" Web:       Visit www.github.com/dbry/audio-resampler for latest version and info\n\n";

static int wav_process (char *infilename, char *outfilename);
static int bh4_window, hann_window, num_taps = 256, num_filters = 256, outbits, verbosity, pre_post_filter;
static int dither = DITHER_HIGHPASS, noise_shaping = SHAPING_ATH_CURVE;
static int test_non_interleaved, non_interleaved_bytes;
static unsigned long resample_rate, lowpass_freq;
static double phase_shift, gain = 1.0;

int main (argc, argv) int argc; char **argv;
{
    int overwrite = 0;
    char *infilename = NULL, *outfilename = NULL;
    FILE *outfile;

    // loop through command-line arguments

    while (--argc) {
#if defined (_WIN32)
        if ((**++argv == '-' || **argv == '/') && (*argv)[1])
#else
        if ((**++argv == '-') && (*argv)[1])
#endif
            while (*++*argv)
                switch (**argv) {

		    case '0':                           // for debug only (forces resampler off)
			num_filters = num_taps = 0;
			break;

		    case '1':
			num_filters = num_taps = 16;
			break;

		    case '2':
			num_filters = num_taps = 64;
			break;

		    case '3':
			num_filters = num_taps = 256;
			break;

		    case '4':
			num_filters = num_taps = 1024;
			break;

                    case 'P': case 'p':
                        pre_post_filter = 1;
                        break;

                    case 'Q': case 'q':
                        verbosity = -1;
                        break;

                    case 'V': case 'v':
                        verbosity = 1;
                        break;

                    case 'Y': case 'y':
                        overwrite = 1;
                        break;

		    case 'R': case 'r':
			resample_rate = strtod (++*argv, argv);
			--*argv;
			break;

		    case 'D': case 'd':
                        {
                            int dither_select = strtod (++*argv, argv);

                            switch (dither_select) {
                                case 0:
                                    dither = 0;
                                    break;
                                case 1:
                                    dither = DITHER_FLAT;
                                    break;
                                case 2:
                                    dither = DITHER_LOWPASS;
                                    break;
                                default:
                                    fprintf (stderr, "\ndither override must be 0, 1, or 2!\n");
                                    return 1;
                            }
                        }

			--*argv;
			break;

		    case 'N': case 'n':
                        {
			    int noise_shaping_select = strtod (++*argv, argv);

                            switch (noise_shaping_select) {
                                case 0:
                                    noise_shaping = 0;
                                    break;
                                case 1:
                                    noise_shaping = SHAPING_1ST_ORDER;
                                    break;
                                case 2:
                                    noise_shaping = SHAPING_2ND_ORDER;
                                    break;
                                case 3:
                                    noise_shaping = SHAPING_3RD_ORDER;
                                    break;
                                default:
                                    fprintf (stderr, "\nnoise-shaping override must be 0, 1, 2, or 3!\n");
                                    return 1;
                            }
                        }

			--*argv;
			break;

		    case 'S': case 's':
			phase_shift = strtod (++*argv, argv) / 360.0;

                        if (phase_shift <= -1.0 || phase_shift >= 1.0) {
                            fprintf (stderr, "\nphase shift must be less than +/- 1 sample!\n");
                            return 1;
                        }

			--*argv;
			break;

		    case 'G': case 'g':
			gain = pow (10.0, strtod (++*argv, argv) / 20.0);
			--*argv;
			break;

		    case 'L': case 'l':
			lowpass_freq = strtod (++*argv, argv);
			--*argv;
			break;

		    case 'F': case 'f':
			num_filters = strtod (++*argv, argv);

                        if (num_filters < 2 || num_filters > 1024) {
                            fprintf (stderr, "\nnum of filters must be 2 - 1024!\n");
                            return 1;
                        }

			--*argv;
			break;

		    case 'I': case 'i':
			non_interleaved_bytes = strtod (++*argv, argv);
                        test_non_interleaved = 1;
			--*argv;
			break;

		    case 'O': case 'o':
			outbits = strtod (++*argv, argv);

                        if (outbits != 32 && (outbits < 4 || outbits > 24)) {
                            fprintf (stderr, "\noutbits must be 4 - 24 (for integer) or 32 (for float)!\n");
                            return 1;
                        }

			--*argv;
			break;

		    case 'T': case 't':
			num_taps = strtod (++*argv, argv);

                        if ((num_taps & 3) || num_taps < 4 || num_taps > 1024) {
                            fprintf (stderr, "\nnum of taps must be 4 - 1024 and a multiple of 4!\n");
                            return 1;
                        }

			--*argv;
			break;

		    case 'B': case 'b':
			bh4_window = 1;
			break;

		    case 'H': case 'h':
			hann_window = 1;
			break;

                    default:
                        fprintf (stderr, "\nillegal option: %c !\n", **argv);
                        return 1;
                }
        else if (!infilename) {
            infilename = malloc (strlen (*argv) + 10);
            strcpy (infilename, *argv);
        }
        else if (!outfilename) {
            outfilename = malloc (strlen (*argv) + 10);
            strcpy (outfilename, *argv);
        }
        else {
            fprintf (stderr, "\nextra unknown argument: %s !\n", *argv);
            return 1;
        }
    }

    if (verbosity >= 0)
        fprintf (stderr, sign_on, VERSION);

    if (!outfilename) {
        printf ("%s", usage);
        return 0;
    }

    if (!strcmp (infilename, outfilename)) {
        fprintf (stderr, "can't overwrite input file (specify different/new output file name)\n");
        return -1;
    }

    if (!overwrite && (outfile = fopen (outfilename, "r"))) {
        fclose (outfile);
        fprintf (stderr, "output file \"%s\" exists (use -y to overwrite)\n", outfilename);
        return -1;
    }

    int res = wav_process (infilename, outfilename);

    free (infilename);
    free (outfilename);

    return res;
}

typedef struct {
    char ckID [4];
    uint32_t ckSize;
    char formType [4];
} RiffChunkHeader;

typedef struct {
    char ckID [4];
    uint32_t ckSize;
} ChunkHeader;

#define ChunkHeaderFormat "4L"

typedef struct {
    uint16_t FormatTag, NumChannels;
    uint32_t SampleRate, BytesPerSecond;
    uint16_t BlockAlign, BitsPerSample;
    uint16_t cbSize;
    union {
        uint16_t ValidBitsPerSample;
        uint16_t SamplesPerBlock;
        uint16_t Reserved;
    } Samples;
    int32_t ChannelMask;
    uint16_t SubFormat;
    char GUID [14];
} WaveHeader;

#define WaveHeaderFormat "SSLLSSSSLS"

#define WAVE_FORMAT_PCM         0x1
#define WAVE_FORMAT_IEEE_FLOAT  0x3
#define WAVE_FORMAT_EXTENSIBLE  0xfffe

static unsigned int process_audio (FILE *infile, FILE *outfile, unsigned long sample_rate,
    unsigned long num_samples, int num_channels, int inbits);

static int write_pcm_wav_header (FILE *outfile, int bps, int num_channels, unsigned long num_samples, unsigned long sample_rate, uint32_t channel_mask);
static void little_endian_to_native (void *data, char *format);
static void native_to_little_endian (void *data, char *format);

static int wav_process (char *infilename, char *outfilename)
{
    int format = 0, res = 0, inbits = 0, num_channels = 0;
    unsigned long num_samples = 0, sample_rate = 0;
    uint32_t channel_mask = 0;
    FILE *infile, *outfile;
    RiffChunkHeader riff_chunk_header;
    ChunkHeader chunk_header;
    WaveHeader WaveHeader;

    // open both input and output files

    if (!(infile = fopen (infilename, "rb"))) {
        fprintf (stderr, "can't open file \"%s\" for reading!\n", infilename);
        return -1;
    }

    if (!(outfile = fopen (outfilename, "wb"))) {
        fprintf (stderr, "can't open file \"%s\" for writing!\n", outfilename);
        fclose (infile);
        return -1;
    }

    // read (and write) initial RIFF form header

    if (!fread (&riff_chunk_header, sizeof (RiffChunkHeader), 1, infile) ||
        strncmp (riff_chunk_header.ckID, "RIFF", 4) ||
        strncmp (riff_chunk_header.formType, "WAVE", 4)) {
            fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
            fclose (outfile);
            fclose (infile);
            return -1;
    }

    // loop through all elements of the RIFF wav header (until the data chuck)

    while (1) {

        if (!fread (&chunk_header, sizeof (ChunkHeader), 1, infile)) {
            fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
            fclose (outfile);
            fclose (infile);
            return -1;
        }

        little_endian_to_native (&chunk_header, ChunkHeaderFormat);

        // if it's the format chunk, we want to get some info out of there and
        // make sure it's a .wav file we can handle

        if (!strncmp (chunk_header.ckID, "fmt ", 4)) {
            int supported = 1;

            if (chunk_header.ckSize < 16 || chunk_header.ckSize > sizeof (WaveHeader) ||
                !fread (&WaveHeader, chunk_header.ckSize, 1, infile)) {
                    fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
                    fclose (outfile);
                    fclose (infile);
                    return -1;
            }

            little_endian_to_native (&WaveHeader, WaveHeaderFormat);

            format = (WaveHeader.FormatTag == WAVE_FORMAT_EXTENSIBLE && chunk_header.ckSize == 40) ?
                WaveHeader.SubFormat : WaveHeader.FormatTag;

            if (WaveHeader.FormatTag == WAVE_FORMAT_EXTENSIBLE && chunk_header.ckSize == 40)
                channel_mask = WaveHeader.ChannelMask;
            else if (WaveHeader.NumChannels <= 2)
                channel_mask = 0x5 - WaveHeader.NumChannels;
            else if (WaveHeader.NumChannels < 32)
                channel_mask = (1U << WaveHeader.NumChannels) - 1;
            else
                channel_mask = 0xffffffff;

            inbits = (chunk_header.ckSize == 40 && WaveHeader.Samples.ValidBitsPerSample) ?
                WaveHeader.Samples.ValidBitsPerSample : WaveHeader.BitsPerSample;

            if (WaveHeader.NumChannels < 1 || WaveHeader.NumChannels > 32)
                supported = 0;
            else if (format == WAVE_FORMAT_PCM) {

                if (inbits < 4 || inbits > 24)
                    supported = 0;

                if (WaveHeader.BlockAlign != WaveHeader.NumChannels * ((inbits + 7) / 8))
                    supported = 0;
            }
            else if (format == WAVE_FORMAT_IEEE_FLOAT) {

                if (inbits != 32)
                    supported = 0;

                if (WaveHeader.BlockAlign != WaveHeader.NumChannels * 4)
                    supported = 0;
            }
            else
                supported = 0;

            if (!supported) {
                fprintf (stderr, "\"%s\" is an unsupported .WAV format!\n", infilename);
                fclose (outfile);
                fclose (infile);
                return -1;
            }

            if (verbosity > 0) {
                fprintf (stderr, "format tag size = %d\n", chunk_header.ckSize);
                fprintf (stderr, "FormatTag = 0x%x, NumChannels = %u, BitsPerSample = %u\n",
                    WaveHeader.FormatTag, WaveHeader.NumChannels, WaveHeader.BitsPerSample);
                fprintf (stderr, "BlockAlign = %u, SampleRate = %lu, BytesPerSecond = %lu\n",
                    WaveHeader.BlockAlign, (unsigned long) WaveHeader.SampleRate, (unsigned long) WaveHeader.BytesPerSecond);

                if (chunk_header.ckSize > 16)
                    fprintf (stderr, "cbSize = %d, ValidBitsPerSample = %d\n", WaveHeader.cbSize,
                        WaveHeader.Samples.ValidBitsPerSample);

                if (chunk_header.ckSize > 20)
                    fprintf (stderr, "ChannelMask = %x, SubFormat = %d\n",
                        WaveHeader.ChannelMask, WaveHeader.SubFormat);
            }
        }
        else if (!strncmp (chunk_header.ckID, "data", 4)) {

            // on the data chunk, get size and exit parsing loop

            if (!WaveHeader.NumChannels) {      // make sure we saw a "fmt" chunk...
                fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
                fclose (outfile);
                fclose (infile);
                return -1;
            }

            if (!chunk_header.ckSize) {
                fprintf (stderr, "this .WAV file has no audio samples, probably is corrupt!\n");
                fclose (outfile);
                fclose (infile);
                return -1;
            }

            if (chunk_header.ckSize % WaveHeader.BlockAlign) {
                fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
                fclose (outfile);
                fclose (infile);
                return -1;
            }

            num_samples = chunk_header.ckSize / WaveHeader.BlockAlign;

            if (!num_samples) {
                fprintf (stderr, "this .WAV file has no audio samples, probably is corrupt!\n");
                fclose (outfile);
                fclose (infile);
                return -1;
            }

            if (verbosity > 0)
                fprintf (stderr, "num samples = %lu\n", num_samples);

            num_channels = WaveHeader.NumChannels;
            sample_rate = WaveHeader.SampleRate;
            break;
        }
        else {          // just ignore/copy unknown chunks
            unsigned int bytes_to_copy = (chunk_header.ckSize + 1) & ~1L;

            if (verbosity > 0)
                fprintf (stderr, "extra unknown chunk \"%c%c%c%c\" of %u bytes\n",
                    chunk_header.ckID [0], chunk_header.ckID [1], chunk_header.ckID [2],
                    chunk_header.ckID [3], bytes_to_copy);

            while (bytes_to_copy) {
                unsigned int bytes_to_read = bytes_to_copy, bytes_read;
                char temp_buffer [256];

                if (bytes_to_read > sizeof (temp_buffer))
                    bytes_to_read = sizeof (temp_buffer);

                bytes_read = fread (temp_buffer, 1, bytes_to_read, infile);

                if (bytes_read != bytes_to_read) {
                    fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
                    fclose (outfile);
                    fclose (infile);
                    return -1;
                }

                bytes_to_copy -= bytes_read;
            }
        }
    }

    if (!num_channels || !sample_rate || !inbits || !num_samples) {
        fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
        fclose (outfile);
        fclose (infile);
        return -1;
    }

    // if not specified, preserve sample rate and bitdepth of input

    if (!resample_rate) resample_rate = sample_rate;
    if (!outbits) outbits = inbits;

    if (verbosity >= 0)
        fprintf (stderr, "resampling %d-channel file \"%s\" (%db/%dk) to \"%s\" (%db/%dk)...\n",
            num_channels, infilename, inbits, (int)((sample_rate + 500) / 1000),
            outfilename, outbits, (int)((resample_rate + 500) / 1000));

    if (!write_pcm_wav_header (outfile, outbits, num_channels, num_samples, resample_rate, channel_mask)) {
        fprintf (stderr, "can't write to file \"%s\"!\n", outfilename);
        fclose (outfile);
        fclose (infile);
        return -1;
    }

    unsigned int output_samples = process_audio (infile, outfile, sample_rate, num_samples, num_channels, inbits);

    // write an extra padding zero byte if the data chunk is not an even size

    if ((output_samples * num_channels * ((outbits + 7) / 8)) & 1)
        fwrite ("", 1, 1, outfile);

    rewind (outfile);

    if (!write_pcm_wav_header (outfile, outbits, num_channels, output_samples, resample_rate, channel_mask)) {
        fprintf (stderr, "can't write to file \"%s\"!\n", outfilename);
        fclose (outfile);
        fclose (infile);
        return -1;
    }

    fclose (outfile);
    fclose (infile);
    return res;
}

#define BUFFER_SAMPLES          4096

static unsigned int process_audio (FILE *infile, FILE *outfile, unsigned long sample_rate,
    unsigned long num_samples, int num_channels, int inbits)
{
    double sample_ratio = (double) resample_rate / sample_rate, lowpass_ratio = 1.0;
    unsigned int outbuffer_samples = (int) floor (BUFFER_SAMPLES * sample_ratio * 1.1 + 100.0);
    unsigned long remaining_samples = num_samples, output_samples = 0, clipped_samples = 0;
    float *const outbuffer = malloc (outbuffer_samples * num_channels * sizeof (float));
    float *inbuffer = malloc (BUFFER_SAMPLES * num_channels * sizeof (float));
    int samples_to_append = 0, pre_filter = 0, post_filter = 0;
    Biquad *lowpass1 = NULL, *lowpass2 = NULL;
    int flags = SUBSAMPLE_INTERPOLATE;
    BiquadCoefficients lowpass_coeff;
    unsigned char *tmpbuffer = NULL;
    void *readbuffer = inbuffer;
    Decimate *decimator = NULL;
    Resample *resampler = NULL;

    // when downsampling, calculate the optimum lowpass based on resample filter
    // length (i.e., more taps allow us to lowpass closer to Nyquist)

    if (sample_ratio < 1.0) {
        lowpass_ratio -= (10.24 / num_taps);

        if (lowpass_ratio < 0.84)           // limit the lowpass for very short filters
            lowpass_ratio = 0.84;

        if (lowpass_ratio < sample_ratio)   // avoid discontinuities near unity sample ratios
            lowpass_ratio = sample_ratio;
    }

    if (lowpass_freq) {
        double user_lowpass_ratio;

        if (sample_ratio < 1.0)
            user_lowpass_ratio = lowpass_freq / (resample_rate / 2.0);
        else
            user_lowpass_ratio = lowpass_freq / (sample_rate / 2.0);

        if (user_lowpass_ratio >= 1.0)
            fprintf (stderr, "warning: ignoring invalid lowpass frequency specification (at or over Nyquist)\n");
        else
            lowpass_ratio = user_lowpass_ratio;
    }

    if (bh4_window || !hann_window)
        flags |= BLACKMAN_HARRIS;

    if (lowpass_ratio * sample_ratio < 0.98 && pre_post_filter) {
        double cutoff = lowpass_ratio * sample_ratio / 2.0;
        biquad_lowpass (&lowpass_coeff, cutoff);
        pre_filter = 1;

        if (verbosity > 0)
            fprintf (stderr, "cascaded biquad pre-filter at %g Hz\n", sample_rate * cutoff);
    }

    if (num_filters && (sample_ratio != 1.0 || lowpass_ratio != 1.0 || phase_shift != 0.0)) {
        if (sample_ratio < 1.0) {
            resampler = resampleInit (num_channels, num_taps, num_filters, sample_ratio * lowpass_ratio, flags | INCLUDE_LOWPASS);

            if (verbosity > 0)
                fprintf (stderr, "%d-tap sinc downsampler with lowpass at %g Hz\n", num_taps, sample_ratio * lowpass_ratio * sample_rate / 2.0);
        }
        else if (lowpass_ratio < 1.0) {
            resampler = resampleInit (num_channels, num_taps, num_filters, lowpass_ratio, flags | INCLUDE_LOWPASS);

            if (verbosity > 0)
                fprintf (stderr, "%d-tap sinc resampler with lowpass at %g Hz\n", num_taps, lowpass_ratio * sample_rate / 2.0);
        }
        else {
            resampler = resampleInit (num_channels, num_taps, num_filters, 1.0, flags);

            if (verbosity > 0)
                fprintf (stderr, "%d-tap pure sinc resampler (no lowpass), %g Hz Nyquist\n", num_taps, sample_rate / 2.0);
        }

        samples_to_append = num_taps / 2;
    }

    if (lowpass_ratio / sample_ratio < 0.98 && pre_post_filter && !pre_filter) {
        double cutoff = lowpass_ratio / sample_ratio / 2.0;
        biquad_lowpass (&lowpass_coeff, cutoff);
        post_filter = 1;

        if (verbosity > 0)
            fprintf (stderr, "cascaded biquad post-filter at %g Hz\n", resample_rate * cutoff);
    }

    if (pre_filter || post_filter) {
        lowpass1 = calloc (num_channels, sizeof (Biquad));
        lowpass2 = calloc (num_channels, sizeof (Biquad));

        for (int i = 0; i < num_channels; ++i) {
            biquad_init (lowpass1 + i, &lowpass_coeff, 1.0);
            biquad_init (lowpass2 + i, &lowpass_coeff, 1.0);
        }
    }

    if (outbits != 32) {
        if (test_non_interleaved && non_interleaved_bytes)
            decimator = decimateInit (num_channels, outbits, non_interleaved_bytes, 1.0, resample_rate, dither | noise_shaping);
        else
            decimator = decimateInit (num_channels, outbits, (outbits + 7) / 8, 1.0, resample_rate, dither | noise_shaping);
    }

    if (inbits != 32 || outbits != 32) {
        int max_samples = BUFFER_SAMPLES, max_bytes = 2;

        if (outbuffer_samples > BUFFER_SAMPLES)
            max_samples = outbuffer_samples;

        if (inbits > 16 || outbits > 16)
            max_bytes = 3;

        tmpbuffer = malloc (max_samples * num_channels * max_bytes);

        if (inbits != 32)
            readbuffer = tmpbuffer;
    }

    // this takes care of the filter delay and any user-specified phase shift
    if (resampler)
        resampleAdvancePosition (resampler, num_taps / 2.0 + phase_shift);

    uint32_t progress_divider = 0, percent;

    if (verbosity >= 0 && remaining_samples > 1000) {
        progress_divider = (remaining_samples + 50) / 100;
        fprintf (stderr, "\rprogress: %d%% ", percent = 0); fflush (stderr);
    }

    while (remaining_samples + samples_to_append) {

        // first we read the audio data, converting to 32-bit float (if not already) and applying gain

        unsigned long samples_to_read = remaining_samples, samples_read, samples_generated;
        ResampleResult res;

        if (samples_to_read > BUFFER_SAMPLES)
            samples_to_read = BUFFER_SAMPLES;

        samples_read = fread (readbuffer, num_channels * ((inbits + 7) / 8), samples_to_read, infile);
        remaining_samples -= samples_read;

        if (!samples_read) {
            int samples_to_append_now = samples_to_append;

            if (!samples_to_append_now)
                break;

            if (samples_to_append_now > BUFFER_SAMPLES)
                samples_to_append_now = BUFFER_SAMPLES;

            memset (readbuffer, (inbits <= 8) * 128, samples_to_append_now * num_channels * ((inbits + 7) / 8));
            samples_read = samples_to_append_now;
            samples_to_append -= samples_to_append_now;
        }

        if (inbits > 24) {
            if (IS_BIG_ENDIAN) {
                unsigned char *bptr = (unsigned char *) inbuffer, word [4];
                int wcount = samples_read * num_channels;

                while (wcount--) {
                    memcpy (word, bptr, 4);
                    *bptr++ = word [3];
                    *bptr++ = word [2];
                    *bptr++ = word [1];
                    *bptr++ = word [0];
                }
            }

            if (gain != 1.0)
                for (int i = 0; i < samples_read * num_channels; ++i)
                    inbuffer [i] *= gain;
        }
        else
            floatIntegersLE (tmpbuffer, gain, inbits, (inbits + 7) / 8, 1, inbuffer, samples_read * num_channels);

        // common code to process the audio in 32-bit floats

        if (pre_filter)
            for (int i = 0; i < num_channels; ++i) {
                biquad_apply_buffer (lowpass1 + i, inbuffer + i, samples_read, num_channels);
                biquad_apply_buffer (lowpass2 + i, inbuffer + i, samples_read, num_channels);
            }

        if (resampler) {
            res = resampleProcessInterleaved (resampler, inbuffer, samples_read, outbuffer, outbuffer_samples, sample_ratio);
            samples_generated = res.output_generated;
        }
        else {
            memcpy (outbuffer, inbuffer, samples_read * num_channels * sizeof (float));
            samples_generated = samples_read;
        }

        if (post_filter)
            for (int i = 0; i < num_channels; ++i) {
                biquad_apply_buffer (lowpass1 + i, outbuffer + i, samples_generated, num_channels);
                biquad_apply_buffer (lowpass2 + i, outbuffer + i, samples_generated, num_channels);
            }

        // finally write the audio, converting to appropriate integer format if requested

        if (outbits != 32) {
            if (test_non_interleaved) {
                unsigned char *const output_array [num_channels];
                const float* const input_array [num_channels];
                unsigned char *destinptr = tmpbuffer;
                int outbytes = (outbits + 7) / 8;
                float *sourceptr = outbuffer;

                for (int c = 0; c < num_channels; ++c) {
                    if (non_interleaved_bytes)
                        ((unsigned char**) output_array) [c] = malloc (samples_generated * non_interleaved_bytes);
                    else
                        ((unsigned char**) output_array) [c] = malloc (samples_generated * sizeof (unsigned char) * outbytes);

                    ((float**) input_array) [c] = malloc (samples_generated * sizeof (float));
                }

                for (int i = 0; i < samples_generated; ++i)
                    for (int c = 0; c < num_channels; ++c)
                        ((float**) input_array) [c] [i] = *sourceptr++;

                clipped_samples += decimateProcessLE (decimator, input_array, samples_generated, output_array);

                for (int i = 0; i < samples_generated; ++i)
                    for (int c = 0; c < num_channels; ++c) {
                        unsigned char *src = ((unsigned char *) output_array [c]) + i * outbytes;

                        if (non_interleaved_bytes)
                            src = ((unsigned char *) output_array [c]) + i * non_interleaved_bytes + (non_interleaved_bytes - outbytes);

                        for (int j = 0; j < outbytes; ++j)
                            *destinptr++ = *src++;
                   }

                for (int c = 0; c < num_channels; ++c) {
                    free ((void *) output_array [c]);
                    free ((void *) input_array [c]);
                }
            }
            else
                clipped_samples += decimateProcessInterleavedLE (decimator, outbuffer, samples_generated, tmpbuffer);

            fwrite (tmpbuffer, num_channels * ((outbits + 7) / 8), samples_generated, outfile);
        }
        else {
            if (IS_BIG_ENDIAN) {
                unsigned char *bptr = (unsigned char *) outbuffer, word [4];
                int wcount = samples_generated * num_channels;

                while (wcount--) {
                    memcpy (word, bptr, 4);
                    *bptr++ = word [3];
                    *bptr++ = word [2];
                    *bptr++ = word [1];
                    *bptr++ = word [0];
                }
            }

            fwrite (outbuffer, num_channels * sizeof (float), samples_generated, outfile);
        }

        output_samples += samples_generated;

        if (progress_divider) {
            int new_percent = 100 - remaining_samples / progress_divider;

            if (new_percent != percent) {
                fprintf (stderr, "\rprogress: %d%% ", percent = new_percent);
                fflush (stderr);
            }
        }
    }

    if (verbosity >= 0)
        fprintf (stderr, "\r...completed successfully\n");

    decimateFree (decimator);
    resampleFree (resampler);
    free (inbuffer);
    free (outbuffer);
    free (tmpbuffer);
    free (lowpass1);
    free (lowpass2);

    if (clipped_samples)
        fprintf (stderr, "warning: %lu samples were clipped, suggest reducing gain!\n", clipped_samples);

    if (remaining_samples)
        fprintf (stderr, "warning: file terminated early!\n");

    return output_samples;
}

static int write_pcm_wav_header (FILE *outfile, int bps, int num_channels, unsigned long num_samples, unsigned long sample_rate, uint32_t channel_mask)
{
    RiffChunkHeader riffhdr;
    ChunkHeader datahdr, fmthdr;
    WaveHeader wavhdr;

    int wavhdrsize = 16;
    int bytes_per_sample = (bps + 7) / 8;
    int format = (bps == 32) ? WAVE_FORMAT_IEEE_FLOAT : WAVE_FORMAT_PCM;
    uint32_t total_data_bytes = num_samples * bytes_per_sample * num_channels;

    memset (&wavhdr, 0, sizeof (wavhdr));

    wavhdr.FormatTag = format;
    wavhdr.NumChannels = num_channels;
    wavhdr.SampleRate = sample_rate;
    wavhdr.BytesPerSecond = sample_rate * num_channels * bytes_per_sample;
    wavhdr.BlockAlign = bytes_per_sample * num_channels;
    wavhdr.BitsPerSample = bps;

    // write an extended header if more than 2 channels or a non-standard channel-mask

    if (num_channels > 2 || channel_mask != 0x5 - num_channels) {
        wavhdrsize = sizeof (wavhdr);
        wavhdr.cbSize = 22;
        wavhdr.Samples.ValidBitsPerSample = bps;
        wavhdr.SubFormat = format;
        wavhdr.ChannelMask = channel_mask;
        wavhdr.FormatTag = WAVE_FORMAT_EXTENSIBLE;
        wavhdr.BitsPerSample = bps;
        wavhdr.GUID [4] = 0x10;
        wavhdr.GUID [6] = 0x80;
        wavhdr.GUID [9] = 0xaa;
        wavhdr.GUID [11] = 0x38;
        wavhdr.GUID [12] = 0x9b;
        wavhdr.GUID [13] = 0x71;
    }

    memcpy (riffhdr.ckID, "RIFF", sizeof (riffhdr.ckID));
    memcpy (riffhdr.formType, "WAVE", sizeof (riffhdr.formType));
    riffhdr.ckSize = (sizeof (riffhdr) + wavhdrsize + sizeof (datahdr) + total_data_bytes + 1) & ~1U;
    memcpy (fmthdr.ckID, "fmt ", sizeof (fmthdr.ckID));
    fmthdr.ckSize = wavhdrsize;

    memcpy (datahdr.ckID, "data", sizeof (datahdr.ckID));
    datahdr.ckSize = total_data_bytes;

    // write the RIFF chunks up to just before the data starts

    native_to_little_endian (&riffhdr, ChunkHeaderFormat);
    native_to_little_endian (&fmthdr, ChunkHeaderFormat);
    native_to_little_endian (&wavhdr, WaveHeaderFormat);
    native_to_little_endian (&datahdr, ChunkHeaderFormat);

    return fwrite (&riffhdr, sizeof (riffhdr), 1, outfile) &&
        fwrite (&fmthdr, sizeof (fmthdr), 1, outfile) &&
        fwrite (&wavhdr, wavhdrsize, 1, outfile) &&
        fwrite (&datahdr, sizeof (datahdr), 1, outfile);
}

static void little_endian_to_native (void *data, char *format)
{
    unsigned char *cp = (unsigned char *) data;
    int32_t temp;

    while (*format) {
        switch (*format) {
            case 'L':
                temp = cp [0] + ((int32_t) cp [1] << 8) + ((int32_t) cp [2] << 16) + ((int32_t) cp [3] << 24);
                * (int32_t *) cp = temp;
                cp += 4;
                break;

            case 'S':
                temp = cp [0] + (cp [1] << 8);
                * (short *) cp = (short) temp;
                cp += 2;
                break;

            default:
                if (isdigit ((unsigned char) *format))
                    cp += *format - '0';

                break;
        }

        format++;
    }
}

static void native_to_little_endian (void *data, char *format)
{
    unsigned char *cp = (unsigned char *) data;
    int32_t temp;

    while (*format) {
        switch (*format) {
            case 'L':
                temp = * (int32_t *) cp;
                *cp++ = (unsigned char) temp;
                *cp++ = (unsigned char) (temp >> 8);
                *cp++ = (unsigned char) (temp >> 16);
                *cp++ = (unsigned char) (temp >> 24);
                break;

            case 'S':
                temp = * (short *) cp;
                *cp++ = (unsigned char) temp;
                *cp++ = (unsigned char) (temp >> 8);
                break;

            default:
                if (isdigit ((unsigned char) *format))
                    cp += *format - '0';

                break;
        }

        format++;
    }
}

