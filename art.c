////////////////////////////////////////////////////////////////////////////
//                            **** ART ****                               //
//                        Audio Resampling Tool                           //
//                 Copyright (c) 2006-2025 David Bryant                   //
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
#include "stretch.h"
#include "biquad.h"

#define VERSION         0.5

#define IS_BIG_ENDIAN (*(uint16_t *)"\0\xff" < 0x0100)

static const char *sign_on = "\n"
" ART  Audio Resampling Tool  Version %.1f\n"
" Copyright (c) 2006 - 2025 David Bryant.\n\n";

static const char *usage =
" Usage:     ART [-options] infile.wav outfile.wav\n\n"
" Options:  -1|2|3|4    = quality presets, default = 3\n"
"           -r<Hz>      = resample to specified rate in Hz\n"
"                           (follow rate with 'k' for kHz)\n"
"           -g<dB>      = apply gain (default = 0 dB)\n"
"           -s<degrees> = add specified phase shift (+/-360 degrees)\n"
"           -l<Hz>      = specify alternate lowpass frequency in Hz\n"
"                           (follow freq with 'k' for kHz)\n"
"           -f<num>     = number of sinc filters (1-1024)\n"
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
"           -a          = allpass sinc (no lowpass, even downsampling)\n"
"           -b          = Blackman-Harris windowing (best stopband)\n"
"           -h          = Hann windowing (fastest transition)\n"
#ifdef ENABLE_THREADS
"           -m          = use multithreading on stereo & multichannel files\n"
#endif
"           -p          = pre/post filtering (cascaded biquads)\n"
"           -q          = quiet mode (display errors only)\n"
"           -v          = verbose (display lots of info)\n"
"           -y          = overwrite outfile if it exists\n\n"
"           Using any of the following options will invoke the audio-stretch\n"
"           functionality that, unlike regular resampling, can result in very\n"
"           audible and annoying artifacts, especially when used with large\n"
"           stretch ratios or with highly polyphonic source material. Also,\n"
"           these only work with mono or stereo audio (not multichannel).\n\n"
"           --pitch=<cents>   = set pitch shift in cents (+/-2400)\n"
"           --tempo=<ratio>   = set tempo ratio (0.25x - 4.0x)\n"
"           --duration=<[+|-][[hh:]mm:]ss.ss> = set a target duration\n"
"                                 (absolute or +/-relative to source)\n\n"
" Web:       Visit www.github.com/dbry/audio-resampler for latest version and info\n\n";

static int wav_process (char *infilename, char *outfilename);

static int bh4_window, hann_window, num_taps = 256, num_filters = 320, outbits, verbosity, pre_post_filter, allpass, enable_threads;
static int dither = DITHER_HIGHPASS, noise_shaping = SHAPING_ATH_CURVE;
static double pitch_ratio = 1.0, tempo_ratio = 1.0;
static unsigned long resample_rate, lowpass_freq;
static double phase_shift, gain = 1.0;

static struct time_spec {
    int value_is_relative, value_is_valid;
    double value;
} duration;

static void parse_time_spec (struct time_spec *dst, char *src);

int main (int argc, char **argv)
{
    int overwrite = 0, res;
    char *infilename = NULL, *outfilename = NULL;
    FILE *outfile;

    // loop through command-line arguments

    while (--argc) {
        if (**++argv == '-' && (*argv)[1] == '-' && (*argv)[2]) {
            char *long_option = *argv + 2, *long_param = long_option;

            while (*long_param)
                if (*long_param++ == '=')
                    break;

            if (!strncmp (long_option, "pitch", 5)) {                   // --pitch
                double pitch_cents = strtod (long_param, NULL);

                if (pitch_cents < -2400 || pitch_cents > 2400) {
                    fprintf (stderr, "invalid pitch shift, must be +/- 2400 cents (2 octaves)!\n");
                    return 1;
                }

                pitch_ratio = pow (2.0, pitch_cents / 1200.0);
            }
            else if (!strncmp (long_option, "tempo", 5)) {              // --tempo
                tempo_ratio = strtod (long_param, NULL);

                if (tempo_ratio < 0.25 || tempo_ratio > 4.0) {
                    fprintf (stderr, "invalid tempo, must be 0.25 to 4.0!\n");
                    return 1;
                }
            }
            else if (!strncmp (long_option, "duration", 5)) {           // --duration
                parse_time_spec (&duration, long_param);

                if (!duration.value_is_valid) {
                    fprintf (stderr, "invalid --duration parameter!\n");
                    return 1;
                }
            }
            else {
                fprintf (stderr, "unknown option: %s !\n", long_option);
                return 1;
            }
        }
#if defined (_WIN32)
        else if ((**argv == '-' || **argv == '/') && (*argv)[1])
#else
        else if ((**argv == '-') && (*argv)[1])
#endif
            while (*++*argv)
                switch (**argv) {

                    case '1':
                        num_filters = num_taps = 16;
                        break;

                    case '2':
                        num_filters = num_taps = 64;
                        break;

                    case '3':
                        num_filters = 320;      // hack to allow optimized 44.1k --> 96k at default quality
                        num_taps = 256;
                        break;

                    case '4':
                        num_filters = num_taps = 1024;
                        break;

                    case 'A': case 'a':
                        allpass = 1;
                        break;

                    case 'M': case 'm':
                        enable_threads = 1;
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
                        {
                            double rate = strtod (++*argv, argv);

                            if ((**argv & 0xdf) == 'K')
                                rate *= 1000.0;
                            else
                                --*argv;

                            resample_rate = rate;
                        }

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
                        {
                            double freq = strtod (++*argv, argv);

                            if ((**argv & 0xdf) == 'K')
                                freq *= 1000.0;
                            else
                                --*argv;

                            lowpass_freq = freq;
                        }

                        break;

                    case 'F': case 'f':
                        num_filters = strtod (++*argv, argv);

                        if (num_filters < 1 || num_filters > 1024) {
                            fprintf (stderr, "\nnum of filters must be 1 - 1024!\n");
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

    if (lowpass_freq && allpass) {
        fprintf (stderr, "error: can't specify BOTH the allpass option and a lowpass frequency!\n");
        return 1;
    }

    if (duration.value_is_valid && tempo_ratio != 1.0) {
        fprintf (stderr, "error: can't specify BOTH a tempo change and a target duration!\n");
        return 1;
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

    res = wav_process (infilename, outfilename);

    free (infilename);
    free (outfilename);

    return res;
}

// Parse the parameter of the --duration option, which is of the form:
//   [+|-][[hh:]mm:]ss.ss
// The value is returned in a double (in the "dst" struct) as absolute seconds,
// and the sign and abs/rel selection is in the "value_is_relative" field.

static void parse_time_spec (struct time_spec *dst, char *src)
{
    int colons = 0;

    memset (dst, 0, sizeof (*dst));

    if (*src == '+' || *src == '-')
        dst->value_is_relative = (*src++ == '+') ? 1 : -1;

    while (*src)
        if (*src == ':') {
            if (++colons == 3 || dst->value != floor (dst->value))
                return;

            src++;
            dst->value *= 60.0;
            continue;
        }
        else if (*src == '.' || isdigit (*src)) {
            double temp = strtod (src, &src);

            if (temp < 0.0 || (colons && temp >= 60.0))
                return;

            dst->value += temp;
        }
        else
            return;

    dst->value_is_valid = 1;
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
    unsigned char GUID [14];
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
    unsigned int output_samples;
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

    output_samples = process_audio (infile, outfile, sample_rate, num_samples, num_channels, inbits);

    if (output_samples == (unsigned int) -1) {
        fclose (outfile);
        fclose (infile);
        return -1;
    }

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

#define BUFFER_SAMPLES          16384

static unsigned int process_audio (FILE *infile, FILE *outfile, unsigned long sample_rate,
    unsigned long num_samples, int num_channels, int inbits)
{
    unsigned long remaining_samples = num_samples, output_samples = 0, clipped_samples = 0, target_output_samples;
    float *inbuffer = malloc (BUFFER_SAMPLES * num_channels * sizeof (float)), *stretch_buffer = NULL, *outbuffer;
    double sample_ratio = (double) resample_rate / sample_rate, stretch_ratio = 1.0;
    int upper_frequency = 350, lower_frequency = 50, pre_filter = 0, post_filter = 0;
    Biquad *lowpass1 = NULL, *lowpass2 = NULL;
    uint32_t progress_divider = 0, percent;
    float *resample_buffer = inbuffer;
    BiquadCoefficients lowpass_coeff;
    unsigned char *tmpbuffer = NULL;
    unsigned int outbuffer_samples;
    void *readbuffer = inbuffer;
    Decimate *decimator = NULL;
    Resample *resampler = NULL;
    Stretch *stretcher = NULL;

    // if the user specified a target duration (absolute or relative), we calculate a tempo ratio here

    if (duration.value_is_valid) {
        double source_seconds = (double) num_samples / sample_rate, target_seconds;

        switch (duration.value_is_relative) {
            case -1:
                target_seconds = source_seconds - duration.value;
                break;

            case 1:
                target_seconds = source_seconds + duration.value;
                break;

            default:
                target_seconds = duration.value;
                break;
        }

        if (target_seconds <= 0.0) {
            fprintf (stderr, "error: invalid relative duration specified!\n");
            return -1;
        }

        tempo_ratio = source_seconds / target_seconds;
    }

    // if a stretch is probably required, initialize that here

    if (pitch_ratio != 1.0 || tempo_ratio != 1.0) {
        stretch_ratio = pitch_ratio / tempo_ratio;
        sample_ratio /= pitch_ratio;

        if (stretch_ratio != 1.0) {
            int stretch_flags = (stretch_ratio < 0.5 || stretch_ratio > 2.0) ? STRETCH_DUAL_FLAG : 0;
            int stretch_samples;

            if (num_channels > 2) {
                fprintf (stderr, "error: audio stretch only works with mono or stereo, not %d-channel\n", num_channels);
                return -1;
            }

            if (stretch_ratio < 0.25 || stretch_ratio > 4.0) {
                fprintf (stderr, "error: audio stretch requires excessive ratio %g\n", stretch_ratio);
                return -1;
            }

            stretcher = stretchInit (sample_rate / upper_frequency, sample_rate / lower_frequency, num_channels, stretch_flags);
            stretch_samples = stretchGetOutputCapacity (stretcher, BUFFER_SAMPLES, stretch_ratio);
            resample_buffer = stretch_buffer = malloc (stretch_samples * num_channels * sizeof (float));
            outbuffer_samples = (int) floor (stretch_samples * sample_ratio * 1.1 + 100.0);

            if (verbosity > 0)
                fprintf (stderr, "audio stretch initialized with ratio %g\n", stretch_ratio);
        }
        else 
            outbuffer_samples = (int) floor (BUFFER_SAMPLES * sample_ratio * 1.1 + 100.0);
    }
    else 
        outbuffer_samples = (int) floor (BUFFER_SAMPLES * sample_ratio * 1.1 + 100.0);

    outbuffer = malloc (outbuffer_samples * num_channels * sizeof (float));
    target_output_samples = (unsigned long) floor ((double) num_samples * stretch_ratio * sample_ratio + 0.5);

    if (num_filters && (sample_ratio != 1.0 || lowpass_freq || phase_shift != 0.0)) {
        int flags = SUBSAMPLE_INTERPOLATE | INCLUDE_LOWPASS;

#ifdef ENABLE_THREADS
        if (enable_threads)
            flags |= RESAMPLE_MULTITHREADED;
#endif

        if (bh4_window || !hann_window)
            flags |= BLACKMAN_HARRIS;

        if (phase_shift != 0.0)
            flags |= NO_FILTER_REDUCTION;

        if (allpass)
            flags &= ~INCLUDE_LOWPASS;

        resampler = resampleFixedRatioInit (num_channels, num_taps, num_filters, sample_rate * pitch_ratio, resample_rate, lowpass_freq, flags);

        if (!resampler) {
            fprintf (stderr, "error: resampler initialization failed!\n");
            return -1;
        }

        lowpass_freq = resampleGetLowpassRatio (resampler) * (sample_rate * pitch_ratio / 2.0);
        num_filters = resampleGetNumFilters (resampler);

        if (verbosity > 0) {
            if (resampleGetLowpassRatio (resampler) == 1.0)
                fprintf (stderr, "%d %d-tap fixed-ratio sinc resampler%s, no lowpass, %s interpolation\n", num_filters, num_taps,
                    num_filters > 1 ? "s" : "", resampleInterpolationUsed (resampler) ? "with" : "no");
            else
                fprintf (stderr, "%d %d-tap fixed-rate sinc resampler%s with lowpass at %lu Hz, %s interpolation\n", num_filters, num_taps,
                    num_filters > 1 ? "s" : "", lowpass_freq, resampleInterpolationUsed (resampler) ? "with" : "no");
        }
    }

    if (pre_post_filter) {
        if (resample_rate <= sample_rate) {
            double cutoff = resample_rate * 0.45 / sample_rate;
            biquad_lowpass (&lowpass_coeff, cutoff);
            pre_filter = 1;

            if (verbosity > 0)
                fprintf (stderr, "cutoff = %g, cascaded biquad pre-filter at %g Hz\n", cutoff, sample_rate * cutoff);
        }
        else {
            double cutoff = (double) sample_rate * 0.45 / resample_rate;
            biquad_lowpass (&lowpass_coeff, cutoff);
            post_filter = 1;

            if (verbosity > 0)
                fprintf (stderr, "cascaded biquad post-filter at %g Hz\n", resample_rate * cutoff);
        }
    }

    if (pre_filter || post_filter) {
        int  i;

        lowpass1 = calloc (num_channels, sizeof (Biquad));
        lowpass2 = calloc (num_channels, sizeof (Biquad));

        for (i = 0; i < num_channels; ++i) {
            biquad_init (lowpass1 + i, &lowpass_coeff, 1.0);
            biquad_init (lowpass2 + i, &lowpass_coeff, 1.0);
        }
    }

    if (outbits != 32) {
        int decimate_flags = dither | noise_shaping;

#ifdef ENABLE_THREADS
        if (enable_threads)
            decimate_flags |= DECIMATE_MULTITHREADED;
#endif

        decimator = decimateInit (num_channels, outbits, (outbits + 7) / 8, 1.0, resample_rate, decimate_flags);
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

    if (verbosity >= 0 && remaining_samples > 1000) {
        progress_divider = (remaining_samples + 50) / 100;
        fprintf (stderr, "\rprogress: %d%% ", percent = 0); fflush (stderr);
    }

    // loop until we have generated the target number of samples

    while (output_samples < target_output_samples) {

        // first we read the audio data, converting to 32-bit float (if not already) and applying gain

        unsigned long samples_to_read = remaining_samples, samples_read, samples_generated;
        ResampleResult res;

        if (samples_to_read > BUFFER_SAMPLES)
            samples_to_read = BUFFER_SAMPLES;

        samples_read = fread (readbuffer, num_channels * ((inbits + 7) / 8), samples_to_read, infile);
        remaining_samples -= samples_read;

        if (!samples_read) {
            memset (readbuffer, (inbits <= 8) * 128, BUFFER_SAMPLES * num_channels * ((inbits + 7) / 8));
            samples_read = BUFFER_SAMPLES;
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

            if (gain != 1.0) {
                int  i;
                for (i = 0; i < samples_read * num_channels; ++i)
                    inbuffer [i] *= gain;
            }
        }
        else
            floatIntegersLE (tmpbuffer, gain, inbits, (inbits + 7) / 8, 1, inbuffer, samples_read * num_channels);

        // common code to process the audio in 32-bit floats
        // first step is any audio stretching, which is done into stretch_buffer

        if (stretcher)
            samples_read = stretchProcess (stretcher, inbuffer, samples_read, stretch_buffer, stretch_ratio);

        // then any pre-filtering, which is done inline

        if (pre_filter) {
            int  i;
            for (i = 0; i < num_channels; ++i) {
                biquad_apply_buffer (lowpass1 + i, inbuffer + i, samples_read, num_channels);
                biquad_apply_buffer (lowpass2 + i, inbuffer + i, samples_read, num_channels);
            }
        }

        // next is the actual re-sampling, which is done into the output buffer

        if (resampler) {
            res = resampleProcessInterleaved (resampler, resample_buffer, samples_read, outbuffer, outbuffer_samples, sample_ratio);
            samples_generated = res.output_generated;
        }
        else {
            memcpy (outbuffer, resample_buffer, samples_read * num_channels * sizeof (float));
            samples_generated = samples_read;
        }

        // final processing is any post-filtering, which is done inline

        if (post_filter) {
            int  i;
            for (i = 0; i < num_channels; ++i) {
                biquad_apply_buffer (lowpass1 + i, outbuffer + i, samples_generated, num_channels);
                biquad_apply_buffer (lowpass2 + i, outbuffer + i, samples_generated, num_channels);
            }
        }

        // finally write the audio, converting to appropriate integer format if requested

        if (output_samples + samples_generated > target_output_samples)
            samples_generated = target_output_samples - output_samples;

        if (outbits != 32) {
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
    stretchFree (stretcher);
    free (stretch_buffer);
    free (outbuffer);
    free (tmpbuffer);
    free (inbuffer);
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
