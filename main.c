////////////////////////////////////////////////////////////////////////////
//                        **** AUDIO-STRETCH ****                         //
//                      Time Domain Harmonic Scaler                       //
//                    Copyright (c) 2022 David Bryant                     //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// main.c

// This module provides a demo for the TDHS library using WAV files.

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "stretch.h"

static const char *sign_on = "\n"
" AUDIO-STRETCH  Time Domain Harmonic Scaling Demo  Version 0.3\n"
" Copyright (c) 2022 David Bryant. All Rights Reserved.\n\n";

static const char *usage =
" Usage:     AUDIO-STRETCH [-options] infile.wav outfile.wav\n\n"
" Options:  -r<n.n> = stretch ratio (0.5 to 2.0, default = 1.0)\n"
"           -u<n>   = upper freq period limit (default = 333 Hz)\n"
"           -l<n>   = lower freq period limit (default = 55 Hz)\n"
"           -c      = cycle through all ratios, starting higher\n"
"           -cc     = cycle through all ratios, starting lower\n"
"           -s      = scale rate to preserve duration (not pitch)\n"
"           -f      = fast pitch detection (default >= 32 kHz)\n"
"           -n      = normal pitch detection (default < 32 kHz)\n"
"           -q      = quiet mode (display errors only)\n"
"           -v      = verbose (display lots of info)\n"
"           -y      = overwrite outfile if it exists\n\n"
" Web:      Visit www.github.com/dbry/audio-stretch for latest version\n\n";

typedef struct {
    char ckID [4];
    uint32_t ckSize;
    char formType [4];
} RiffChunkHeader;

typedef struct {
    char ckID [4];
    uint32_t ckSize;
} ChunkHeader;

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

#define WAVE_FORMAT_PCM         0x1
#define WAVE_FORMAT_EXTENSIBLE  0xfffe

static int write_pcm_wav_header (FILE *outfile, uint32_t num_samples, int num_channels, int bytes_per_sample, uint32_t sample_rate);

#define BUFFER_SAMPLES 1024

static int verbose_mode, quiet_mode;

int main (argc, argv) int argc; char **argv;
{
    int asked_help = 0, overwrite = 0, scale_rate = 0, force_fast = 0, force_normal = 0, cycle_ratio = 0;
    int upper_frequency = 333, lower_frequency = 55, min_period, max_period;
    uint32_t samples_to_process, insamples = 0, outsamples = 0;
    char *infilename = NULL, *outfilename = NULL;
    RiffChunkHeader riff_chunk_header;
    WaveHeader WaveHeader = { 0 };
    ChunkHeader chunk_header;
    StretchHandle stretcher;
    FILE *infile, *outfile;
    float ratio = 1.0;

    // loop through command-line arguments

    while (--argc) {
#ifdef _WIN32
        if ((**++argv == '-' || **argv == '/') && (*argv)[1])
#else
        if ((**++argv == '-') && (*argv)[1])
#endif
            while (*++*argv)
                switch (**argv) {

                    case 'U': case 'u':
                        upper_frequency = strtol (++*argv, argv, 10);

                        if (upper_frequency <= 40) {
                            fprintf (stderr, "\nupper frequency must be at least 40 Hz!\n");
                            return -1;
                        }

                        --*argv;
                        break;

                    case 'L': case 'l':
                        lower_frequency = strtol (++*argv, argv, 10);

                        if (lower_frequency < 20) {
                            fprintf (stderr, "\nlower frequency must be at least 20 Hz!\n");
                            return -1;
                        }

                        --*argv;
                        break;

                    case 'R': case 'r':
                        ratio = strtod (++*argv, argv);

                        if (ratio < 0.25 || ratio > 4.0) {
                            fprintf (stderr, "\nratio must be from 0.25 to 4.0!\n");
                            return -1;
                        }

                        --*argv;
                        break;

                    case 'S': case 's':
                        scale_rate = 1;
                        break;

                    case 'C': case 'c':
                        cycle_ratio++;
                        break;

                    case 'F': case 'f':
                        force_fast = 1;
                        break;

                    case 'N': case 'n':
                        force_normal = 1;
                        break;

                    case 'H': case 'h':
                        asked_help = 1;
                        break;

                    case 'V': case 'v':
                        verbose_mode = 1;
                        break;

                    case 'Q': case 'q':
                        quiet_mode = 1;
                        break;

                    case 'Y': case 'y':
                        overwrite = 1;
                        break;

                    default:
                        fprintf (stderr, "\nillegal option: %c !\n", **argv);
                        return -1;
                }
        else if (!infilename)
            infilename = *argv;
        else if (!outfilename)
            outfilename = *argv;
        else {
            fprintf (stderr, "\nextra unknown argument: %s !\n", *argv);
            return -1;
        }
    }

    if (!quiet_mode)
        fprintf (stderr, "%s", sign_on);

    if (!outfilename || asked_help) {
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

    if (!(infile = fopen (infilename, "rb"))) {
        fprintf (stderr, "can't open file \"%s\" for reading!\n", infilename);
        return 1;
    }

    // read initial RIFF form header

    if (!fread (&riff_chunk_header, sizeof (RiffChunkHeader), 1, infile) ||
        strncmp (riff_chunk_header.ckID, "RIFF", 4) ||
        strncmp (riff_chunk_header.formType, "WAVE", 4)) {
            fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
            return 1;
    }

    // loop through all elements of the RIFF wav header (until the data chuck)

    while (1) {
        if (!fread (&chunk_header, sizeof (ChunkHeader), 1, infile)) {
            fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
            return 1;
        }

        // if it's the format chunk, we want to get some info out of there and
        // make sure it's a .wav file we can handle

        if (!strncmp (chunk_header.ckID, "fmt ", 4)) {
            int format, bits_per_sample;

            if (chunk_header.ckSize < 16 || chunk_header.ckSize > sizeof (WaveHeader) ||
                !fread (&WaveHeader, chunk_header.ckSize, 1, infile)) {
                    fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
                    return 1;
            }

            format = (WaveHeader.FormatTag == WAVE_FORMAT_EXTENSIBLE && chunk_header.ckSize == 40) ?
                WaveHeader.SubFormat : WaveHeader.FormatTag;

            bits_per_sample = (chunk_header.ckSize == 40 && WaveHeader.Samples.ValidBitsPerSample) ?
                WaveHeader.Samples.ValidBitsPerSample : WaveHeader.BitsPerSample;

            if (bits_per_sample != 16) {
                fprintf (stderr, "\"%s\" is not a 16-bit .WAV file!\n", infilename);
                return 1;
            }

            if (WaveHeader.NumChannels < 1 || WaveHeader.NumChannels > 2) {
                fprintf (stderr, "\"%s\" is not a mono or stereo .WAV file!\n", infilename);
                return 1;
            }

            if (WaveHeader.BlockAlign != WaveHeader.NumChannels * 2) {
                fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
                return 1;
            }

            if (format == WAVE_FORMAT_PCM) {
                if (WaveHeader.SampleRate < 8000 || WaveHeader.SampleRate > 48000) {
                    fprintf (stderr, "\"%s\" sample rate is %lu, must be 8000 to 48000!\n", infilename, (unsigned long) WaveHeader.SampleRate);
                    return 1;
                }
            }
            else {
                fprintf (stderr, "\"%s\" is not a PCM .WAV file!\n", infilename);
                return 1;
            }
        }
        else if (!strncmp (chunk_header.ckID, "data", 4)) {

            // on the data chunk, get size and exit parsing loop

            if (!WaveHeader.SampleRate) {      // make sure we saw a "fmt" chunk...
                fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
                return 1;
            }

            if (!chunk_header.ckSize) {
                fprintf (stderr, "this .WAV file has no audio samples, probably is corrupt!\n");
                return 1;
            }

            if (chunk_header.ckSize % WaveHeader.BlockAlign) {
                fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
                return 1;
            }

            samples_to_process = chunk_header.ckSize / WaveHeader.BlockAlign;

            if (!samples_to_process) {
                fprintf (stderr, "this .WAV file has no audio samples, probably is corrupt!\n");
                return 1;
            }

            break;
        }
        else {          // just ignore unknown chunks
            uint32_t bytes_to_eat = (chunk_header.ckSize + 1) & ~1L;
            char dummy;

            while (bytes_to_eat--)
                if (!fread (&dummy, 1, 1, infile)) {
                    fprintf (stderr, "\"%s\" is not a valid .WAV file!\n", infilename);
                    return 1;
                }
        }
    }

    if (upper_frequency < lower_frequency * 2 || upper_frequency >= WaveHeader.SampleRate / 2) {
        fprintf (stderr, "invalid frequencies specified!\n");
        fclose (infile);
        return 1;
    }

    min_period = WaveHeader.SampleRate / upper_frequency;
    max_period = WaveHeader.SampleRate / lower_frequency;
    int flags = 0;

    if (ratio < 0.5 || ratio > 2.0)
        flags |= STRETCH_DUAL_FLAG;

    if ((force_fast || WaveHeader.SampleRate >= 32000) && !force_normal)
        flags |= STRETCH_FAST_FLAG;

    if (verbose_mode)
        fprintf (stderr, "initializing stretch library with period range = %d to %d, %d channels, %s\n",
            min_period, max_period, WaveHeader.NumChannels, (flags & STRETCH_FAST_FLAG) ? "fast mode" : "normal mode");

    if (!quiet_mode && ratio == 1.0 && !cycle_ratio)
        fprintf (stderr, "warning: a ratio of 1.0 will do nothing but copy the WAV file!\n");

    if (!quiet_mode && ratio != 1.0 && cycle_ratio && !scale_rate)
        fprintf (stderr, "warning: specifying ratio with cycling doesn't do anything (unless scaling rate)\n");

    stretcher = stretch_init (min_period, max_period, WaveHeader.NumChannels, flags);

    if (!stretcher) {
        fprintf (stderr, "can't initialize stretcher\n");
        fclose (infile);
        return 1;
    }

    if (!(outfile = fopen (outfilename, "wb"))) {
        fprintf (stderr, "can't open file \"%s\" for writing!\n", outfilename);
        fclose (infile);
        return 1;
    }

    uint32_t scaled_rate = scale_rate ? (uint32_t)(WaveHeader.SampleRate * ratio + 0.5) : WaveHeader.SampleRate;
    write_pcm_wav_header (outfile, 0, WaveHeader.NumChannels, 2, scaled_rate);

    int16_t *inbuffer = malloc (BUFFER_SAMPLES * WaveHeader.BlockAlign);
    int16_t *outbuffer = malloc ( stretch_output_capacity (stretcher, BUFFER_SAMPLES) * WaveHeader.BlockAlign);

    if (!inbuffer || !outbuffer) {
        fprintf (stderr, "can't allocate required memory!\n");
        fclose (infile);
        return 1;
    }

    while (1) {
        int samples_read = fread (inbuffer, WaveHeader.BlockAlign,
            samples_to_process >= BUFFER_SAMPLES ? BUFFER_SAMPLES : samples_to_process, infile);
        int samples_generated;

        insamples += samples_read;
        samples_to_process -= samples_read;

        if (cycle_ratio)
            ratio = (sin ((double) outsamples / WaveHeader.SampleRate) * (cycle_ratio & 1 ? 0.75 : -0.75)) + 1.25;

        if (samples_read)
            samples_generated = stretch_samples (stretcher, inbuffer, samples_read, outbuffer, ratio);
        else
            samples_generated = stretch_flush (stretcher, outbuffer);

        if (samples_generated) {
            fwrite (outbuffer, WaveHeader.BlockAlign, samples_generated, outfile);
            outsamples += samples_generated;
        }

        if (!samples_read && !samples_generated)
            break;
    }

    free (inbuffer);
    free (outbuffer);
    stretch_deinit (stretcher);

    fclose (infile);

    rewind (outfile);
    write_pcm_wav_header (outfile, outsamples, WaveHeader.NumChannels, 2, scaled_rate);
    fclose (outfile);

    if (insamples && verbose_mode) {
        fprintf (stderr, "done, %lu samples --> %lu samples (ratio = %.3f)\n",
            (unsigned long) insamples, (unsigned long) outsamples, (float) outsamples / insamples);
        if (scale_rate)
            fprintf (stderr, "sample rate changed from %lu Hz to %lu Hz\n",
                (unsigned long) WaveHeader.SampleRate, (unsigned long) scaled_rate);
    }

    return 0;
}

static int write_pcm_wav_header (FILE *outfile, uint32_t num_samples, int num_channels, int bytes_per_sample, uint32_t sample_rate)
{
    RiffChunkHeader riffhdr;
    ChunkHeader datahdr, fmthdr;
    WaveHeader wavhdr;

    int wavhdrsize = 16;
    uint32_t total_data_bytes = num_samples * bytes_per_sample * num_channels;

    memset (&wavhdr, 0, sizeof (wavhdr));

    wavhdr.FormatTag = WAVE_FORMAT_PCM;
    wavhdr.NumChannels = num_channels;
    wavhdr.SampleRate = sample_rate;
    wavhdr.BytesPerSecond = sample_rate * num_channels * bytes_per_sample;
    wavhdr.BlockAlign = bytes_per_sample * num_channels;
    wavhdr.BitsPerSample = bytes_per_sample * 8;

    memcpy (riffhdr.ckID, "RIFF", sizeof (riffhdr.ckID));
    memcpy (riffhdr.formType, "WAVE", sizeof (riffhdr.formType));
    riffhdr.ckSize = sizeof (riffhdr) + wavhdrsize + sizeof (datahdr) + total_data_bytes;
    memcpy (fmthdr.ckID, "fmt ", sizeof (fmthdr.ckID));
    fmthdr.ckSize = wavhdrsize;

    memcpy (datahdr.ckID, "data", sizeof (datahdr.ckID));
    datahdr.ckSize = total_data_bytes;

    return fwrite (&riffhdr, sizeof (riffhdr), 1, outfile) &&
        fwrite (&fmthdr, sizeof (fmthdr), 1, outfile) &&
        fwrite (&wavhdr, wavhdrsize, 1, outfile) &&
        fwrite (&datahdr, sizeof (datahdr), 1, outfile);
}
