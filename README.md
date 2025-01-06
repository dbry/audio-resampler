## AUDIO-RESAMPLER

Audio Resampling Engine & Command-Line Tool

Copyright (c) 2025 David Bryant.

All Rights Reserved.

Distributed under the [BSD Software License](https://github.com/dbry/audio-resampler/blob/master/license.txt).

## What is this?

This is a simple audio resampler, written entirely in C and specifically targeting embedded systems. It
provides fine control over both the CPU load and memory footprint so it can be easily adapted to a wide
range of hardware (e.g., ESP32 to high-end ARM). It is also well suited for ASRC (asynchronous sample rate
converter) applications because it provides a function to query the exact phase position of the resampler
(which is required in the feedback loop of an ASRC).

The package includes a command-line program (**ART**) to experiment with the resampler and serve as example
code for the engine API. The resampling and filtering code works with only 32-bit float audio data, however
there is also code provided to convert to and from integer audio samples. The decimation includes several
configurable options for dither and noise-shaping, including strong ATH filters tailored for popular
sampling rates.

**ART** works with Microsoft WAV files and includes four quality presets that set the number and size of
the sinc filters and serve as a starting point for experimentation:

Preset|Number of sinc filters|Number of taps per filter|RAM use (stereo)
------|----------------------|-------------------------|----------------
-1    |       16             |            16           | 3.4 Kbytes
-2    |       64             |            64           | 25.3 Kbytes
-3    |      256             |           256           | 293 Kbytes
-4    |     1024             |          1024           | 4244 Kbytes

Preset **-3** is the default and is a reasonable compromise for high-quality resampling on a PC. Presets **-1**
and **-2** are more suited for realtime use on embedded systems and **-4** represents the highest quality
available for this tool.

[Infinite Wave](https://infinitewave.ca/) has a very useful site comparing the audio performance of various sample rate converters and now includes this resampler utilizing presets **-2**, **-3**, and **-4**. [Here's the comparison of **ART** preset **-2** and preset **-4**](https://src.infinitewave.ca/?Top=ART_Preset2&Bot=ART_Preset4&Spec=0100), and you can scroll through to many other converters to see that **ART** is competetive in quality with virtually anything out there.

**ART** supports integer samples from 4-bits to 24-bits, as well as 32-bit floating-point samples. Any
number of channels are supported. Normally the output bitdepth is set to the same as the input file, however
this can be forced to one of the other supported bitdepths with the **-o** option.

Both the resampling and filter engines are endian-safe and the command-line program is endian-aware
(although, of course, WAV files are always little-endian).

## Technical Description

The resampler uses windowed sinc interpolation filters. A configurable number of filters are generated on
initialization representing subdivisions of the unit circle. Normally, the output sample value is calculated
by convolving the input samples (the required history is stored in the resampler) with the two sinc filters
on either side of the desired phase angle, and then linear interpolating to get the precise result. If it is
known that there will always be an exact sinc filter (or there's enough room for many filters) then the
interpolation can be skipped and only the nearest sinc filter is used.

The sinc filters are generated with either Hann or Blackman-Harris (4 term) windowing functions. The
Blackman-Harris is usually the best choice (and the default in the CLI) because it has very good stopband
(side-lobe) rejection. However, in some situations (e.g., short filters) the Hann window might be better
because it has a sharper transition (this is controlled in the CLI with the **-b** and **-h** options).

For upsampling or simple resampling operations the sinc filters preserve the original waveform (i.e., they
are purely interpolative). However, for downsampling applications this is not sufficient because aliasing
occurs if content at the old sampling rate moves above the Nyquist frequency of the new sampling rate. For
these cases the sinc filters are constructed with an implicit lowpass, and this lowpass frequency is
optimized for the length of the interpolation filters (i.e., longer filters allow the lowpass to be closer
to the Nyquist frequency). This is also available for upsampling if desired to eliminate frequencies close
to the Nyquist frequency of the input and reduce aliasing further at the expense of some HF loss (e.g., it
might be desirable to set a lowpass of 20 kHz when resampling up from 44.1 kHz even though it's not
strictly required). The lowpass option is enabled with the **-l** option in the CLI).

It is sometimes desirable to reduce aliasing further with lowpass filters either before downsampling
or after upsampling as this can be more efficient than increasing the length of the sinc filters. This is
enabled with the **-p** option in the CLI and implements a cascaded pair of 2nd-order biquads. Note that
unlike the sinc filters, these filters are not linear-phase and will introduce group delay.

For version 0.3 the decimation code has been moved from the command-line program into its own module
and header file (decimator.[ch]) so that it can be utilized by other applications. Also, the dither and
noise-shaping are now configurable and the noise-shaping defaults to strong ATH filters for popular
sampling rates based on [this excellent guide](https://wiki.hydrogenaud.io/index.php?title=Noise_shaping).

## Building

To build an optimized version of the command-line tool (**ART**) on Linux or OS-X:

> $ gcc -O3 -mavx2 -fno-signed-zeros -fno-trapping-math -fassociative-math art.c resampler.c decimator.c biquad.c -lm -o art

The "help" display from the command-line app:

```
 Usage:     ART [-options] infile.wav outfile.wav

 Options:  -1|2|3|4    = quality presets, default = 3
           -r<Hz>      = resample to specified rate
           -g<dB>      = apply gain (default = 0 dB)
           -s<degrees> = add specified phase shift (+/-360 degrees)
           -l<Hz>      = specify alternate lowpass frequency
           -f<num>     = number of sinc filters (2-1024)
           -t<num>     = number of sinc taps (4-1024, multiples of 4)
           -o<bits>    = change output file bitdepth (4-24 or 32)
           -d<sel>     = override default dither (which is HP tpdf):
                           sel = 0 for no dither
                           sel = 1 for flat tpdf dither
                           sel = 2 for LP tpdf dither
           -n<sel>     = override default noise-shaping (which is ATH)
                           sel = 0 for no noise-shaping
                           sel = 1 for 1st-order shaping
                           sel = 2 for 2nd-order shaping
                           sel = 3 for 3rd-order shaping
           -b          = Blackman-Harris windowing (best stopband)
           -h          = Hann windowing (fastest transition)
           -p          = pre/post filtering (cascaded biquads)
           -q          = quiet mode (display errors only)
           -v          = verbose (display lots of info)
           -y          = overwrite outfile if it exists

 Web:       Visit www.github.com/dbry/audio-resampler for latest version and info
```

## Caveats

- The resampling engine is a single C file, with another C file for the biquad filters. Don't expect
the quality and performance of more advanced libraries, but also don't expect much difficulty integrating
it. The simplicity and flexibility of this code might make it appealing for many applications, especially
on limited-resource systems.
- In the command-line program, unknown RIFF chunk types are correctly parsed on input files, but are
*not* passed to the output file, and pipes are not supported.
- The command-line program is not very restrictive about the option parameters, so it's very easy to
get bad results or even crashes with crazy input.
