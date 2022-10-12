////////////////////////////////////////////////////////////////////////////
//                        **** AUDIO-STRETCH ****                         //
//                      Time Domain Harmonic Scaler                       //
//                    Copyright (c) 2022 David Bryant                     //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// stretch.c

// Time Domain Harmonic Compression and Expansion
//
// This library performs time domain harmonic scaling with pitch detection
// to stretch the timing of a 16-bit PCM signal (either mono or stereo) from
// 1/2 to 2 times its original length. This is done without altering any of
// the tonal characteristics.
//
// Use stereo (num_chans = 2), when both channels are from same source
// and should contain approximately similar content.
// For independent channels, prefer using multiple StretchHandle-instances.
// see https://github.com/dbry/audio-stretch/issues/6
// Multiple instances, of course, will consume more CPU load.
// In addition, different output amounts need to be handled.


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "stretch.h"

#define MIN_PERIOD  24          /* minimum allowable pitch period */
#define MAX_PERIOD  2400        /* maximum allowable pitch period */

#if INT_MAX == 32767
#define MERGE_OFFSET    32768L      /* promote to long before offset */
#define abs32           labs        /* use long abs to avoid UB */
#else
#define MERGE_OFFSET    32768
#define abs32           abs
#endif

#define MAX_CORR    UINT32_MAX  /* maximum value for correlation ratios */

struct stretch_cnxt {
    int num_chans, inbuff_samples, shortest, longest, tail, head, fast_mode;
    int16_t *inbuff, *calcbuff;
    float outsamples_error;
    uint32_t *results;

    struct stretch_cnxt *next;
    int16_t *intermediate;
};

static void merge_blocks (int16_t *output, int16_t *input1, int16_t *input2, int samples);
static int find_period_fast (struct stretch_cnxt *cnxt, int16_t *samples);
static int find_period (struct stretch_cnxt *cnxt, int16_t *samples);

/*
 * Initialize a context of the time stretching code. The shortest and longest periods
 * are specified here. The longest period determines the lowest fundamental frequency
 * that can be handled correctly. Note that higher frequencies can be handled than the
 * shortest period would suggest because multiple periods can be combined, and the
 * worst-case performance will suffer if too short a period is selected.
 */

StretchHandle stretch_init (int shortest_period, int longest_period, int num_channels, int flags)
{
    struct stretch_cnxt *cnxt;

    if (flags & STRETCH_FAST_FLAG) {
        longest_period = (longest_period + 1) & ~1;
        shortest_period &= ~1;
    }

    if (longest_period <= shortest_period || shortest_period < MIN_PERIOD || longest_period > MAX_PERIOD) {
        fprintf (stderr, "stretch_init(): invalid periods!\n");
        return NULL;
    }

    cnxt = (struct stretch_cnxt *) calloc (1, sizeof (struct stretch_cnxt));

    if (cnxt) {
        cnxt->inbuff_samples = longest_period * num_channels * 6;
        cnxt->inbuff = calloc (cnxt->inbuff_samples, sizeof (*cnxt->inbuff));

        if (num_channels == 2 || (flags & STRETCH_FAST_FLAG))
            cnxt->calcbuff = calloc (longest_period * num_channels, sizeof (*cnxt->calcbuff));

        if ((flags & STRETCH_FAST_FLAG))
            cnxt->results = calloc (longest_period, sizeof (*cnxt->results));
    }

    if (!cnxt || !cnxt->inbuff || (num_channels == 2 && (flags & STRETCH_FAST_FLAG) && !cnxt->calcbuff) || ((flags & STRETCH_FAST_FLAG) && !cnxt->results)) {
        fprintf (stderr, "stretch_init(): out of memory!\n");
        return NULL;
    }

    cnxt->head = cnxt->tail = cnxt->longest = longest_period * num_channels;
    cnxt->fast_mode = (flags & STRETCH_FAST_FLAG) ? 1 : 0;
    cnxt->shortest = shortest_period * num_channels;
    cnxt->num_chans = num_channels;

    if (flags & STRETCH_DUAL_FLAG) {
        cnxt->next = stretch_init (shortest_period, longest_period, num_channels, flags & ~STRETCH_DUAL_FLAG);
        cnxt->intermediate = calloc (longest_period * num_channels * 4, sizeof (*cnxt->intermediate));
    }

    return (StretchHandle) cnxt;
}

/*
 * Re-Initialize a context of the time stretching code - as if freshly created
 * with stretch_init(). This drops all internal state.
 */

void stretch_reset (StretchHandle handle)
{
    struct stretch_cnxt *cnxt = (struct stretch_cnxt *) handle;
    cnxt->head = cnxt->tail = cnxt->longest;

    if (cnxt->next)
        cnxt->next->head = cnxt->next->tail = cnxt->next->longest;
}


/*
 * Process the specified samples with the given ratio (which is clipped to the
 * range 0.5 to 2.0). Note that the number of samples refers to total samples for
 * both channels in stereo and can be as large as desired (samples are buffered
 * here). The exact number of samples output is not possible to determine in
 * advance, but the maximum will be the number of input samples times the ratio
 * plus 3X the longest period (or 4X the longest period in "fast" mode).
 */

int stretch_samples (StretchHandle handle, const int16_t *samples, int num_samples, int16_t *output, float ratio)
{
    struct stretch_cnxt *cnxt = (struct stretch_cnxt *) handle;
    int out_samples = 0, next_samples = 0;
    int16_t *outbuf = output;

    if (cnxt->next) {
        outbuf = cnxt->intermediate;
        ratio = sqrt (ratio);
    }

    num_samples *= cnxt->num_chans;

    if (ratio < 0.5)
        ratio = 0.5;
    else if (ratio > 2.0)
        ratio = 2.0;

    /* while we have samples to do... */

    do {
        /* if there are more samples and room for them, copy in */

        if (num_samples && cnxt->head < cnxt->inbuff_samples) {
            int samples_to_copy = num_samples;

            if (samples_to_copy > cnxt->inbuff_samples - cnxt->head)
                samples_to_copy = cnxt->inbuff_samples - cnxt->head; 

            memcpy (cnxt->inbuff + cnxt->head, samples, samples_to_copy * sizeof (cnxt->inbuff [0]));
            num_samples -= samples_to_copy;
            samples += samples_to_copy;
            cnxt->head += samples_to_copy;
        }

        /* while there are enough samples to process, do so */

        while (cnxt->tail >= cnxt->longest && cnxt->head - cnxt->tail >= cnxt->longest * (cnxt->fast_mode ? 3 : 2)) {
            int period = cnxt->fast_mode ? find_period_fast (cnxt, cnxt->inbuff + cnxt->tail) :
                find_period (cnxt, cnxt->inbuff + cnxt->tail);
            float process_ratio;

            /*
             * Once we have calculated the best-match period, there are 4 possible transformations
             * available to convert the input samples to output samples. Obviously we can simply
             * copy the samples verbatim (1:1). Standard TDHS provides algorithms for 2:1 and
             * 1:2 scaling, and I have created an obvious extension for 2:3 scaling. To achieve
             * intermediate ratios we maintain a "error" term (in samples) and use that here to
             * calculate the actual transformation to apply.
             */

            if (cnxt->outsamples_error == 0.0)
                process_ratio = floor (ratio * 2.0 + 0.5) / 2.0;
            else if (cnxt->outsamples_error > 0.0)
                process_ratio = floor (ratio * 2.0) / 2.0;
            else
                process_ratio = ceil (ratio * 2.0) / 2.0;

            if (process_ratio == 0.5) {
                merge_blocks (outbuf + out_samples, cnxt->inbuff + cnxt->tail,
                    cnxt->inbuff + cnxt->tail + period, period);
                cnxt->outsamples_error += period - (period * 2.0 * ratio);
                out_samples += period;
                cnxt->tail += period * 2;
            }
            else if (process_ratio == 1.0) {
                memcpy (outbuf + out_samples, cnxt->inbuff + cnxt->tail, period * 2 * sizeof (cnxt->inbuff [0]));
                cnxt->outsamples_error += (period * 2.0) - (period * 2.0 * ratio);
                out_samples += period * 2;
                cnxt->tail += period * 2;
            }
            else if (process_ratio == 1.5) {
                memcpy (outbuf + out_samples, cnxt->inbuff + cnxt->tail, period * sizeof (cnxt->inbuff [0]));
                merge_blocks (outbuf + out_samples + period, cnxt->inbuff + cnxt->tail + period,
                    cnxt->inbuff + cnxt->tail, period);
                memcpy (outbuf + out_samples + period * 2, cnxt->inbuff + cnxt->tail + period, period * sizeof (cnxt->inbuff [0]));
                cnxt->outsamples_error += (period * 3.0) - (period * 2.0 * ratio);
                out_samples += period * 3;
                cnxt->tail += period * 2;
            }
            else if (process_ratio == 2.0) {
                merge_blocks (outbuf + out_samples, cnxt->inbuff + cnxt->tail,
                    cnxt->inbuff + cnxt->tail - period, period * 2);

                cnxt->outsamples_error += (period * 2.0) - (period * ratio);
                out_samples += period * 2;
                cnxt->tail += period;

                if (cnxt->fast_mode) {
                    merge_blocks (outbuf + out_samples, cnxt->inbuff + cnxt->tail,
                        cnxt->inbuff + cnxt->tail - period, period * 2);

                    cnxt->outsamples_error += (period * 2.0) - (period * ratio);
                    out_samples += period * 2;
                    cnxt->tail += period;
                }
            }
            else
                fprintf (stderr, "stretch_samples: fatal programming error: process_ratio == %g\n", process_ratio);

            if (cnxt->next) {
                next_samples += stretch_samples (cnxt->next, outbuf, out_samples / cnxt->num_chans, output + next_samples * cnxt->num_chans, ratio);
                out_samples = 0;
            }
        }

        /* if we're almost done with buffer, copy the rest back to beginning */

        if (cnxt->head == cnxt->inbuff_samples) {
            int samples_to_move = cnxt->inbuff_samples - cnxt->tail + cnxt->longest;

            memmove (cnxt->inbuff, cnxt->inbuff + cnxt->tail - cnxt->longest,
                samples_to_move * sizeof (cnxt->inbuff [0]));

            cnxt->head -= cnxt->tail - cnxt->longest;
            cnxt->tail = cnxt->longest;
        }

    } while (num_samples);

    return cnxt->next ? next_samples : out_samples / cnxt->num_chans;
}  

/* flush any leftover samples out at normal speed */

int stretch_flush (StretchHandle handle, int16_t *output)
{
    struct stretch_cnxt *cnxt = (struct stretch_cnxt *) handle;
    int samples_to_copy = (cnxt->head - cnxt->tail) / cnxt->num_chans;

    memcpy (output, cnxt->inbuff + cnxt->tail, samples_to_copy * cnxt->num_chans * sizeof (*output));
    cnxt->tail = cnxt->head;

    return samples_to_copy;
}

/* free handle */

void stretch_deinit (StretchHandle handle)
{
    struct stretch_cnxt *cnxt = (struct stretch_cnxt *) handle;

    free (cnxt->calcbuff);
    free (cnxt->results);
    free (cnxt->inbuff);

    if (cnxt->next)
        stretch_deinit (cnxt->next);

    free (cnxt);
}

/*
 * The pitch detection is done by finding the period that produces the
 * maximum value for the following correlation formula applied to two
 * consecutive blocks of the given period length:
 *
 *         sum of the absolute values of each sample in both blocks
 *   ---------------------------------------------------------------------
 *   sum of the absolute differences of each corresponding pair of samples
 *
 * This formula was chosen for two reasons.  First, it produces output values
 * that can directly compared regardless of the pitch period.  Second, the
 * numerator can be accumulated for successive periods, and only the
 * denominator need be completely recalculated.
 */

static int find_period (struct stretch_cnxt *cnxt, int16_t *samples)
{
    uint32_t sum, diff, factor, scaler, best_factor = 0;
    int16_t *calcbuff = samples;
    int period, best_period;
    int i, j;

    period = best_period = cnxt->shortest / cnxt->num_chans;

    // convert stereo to mono, and accumulate sum for longest period

    if (cnxt->num_chans == 2) {
        calcbuff = cnxt->calcbuff;

        for (sum = i = j = 0; i < cnxt->longest * 2; i += 2)
            sum += abs32 (calcbuff [j++] = ((int32_t) samples [i] + samples [i+1]) >> 1);
    }
    else
        for (sum = i = 0; i < cnxt->longest; ++i)
            sum += abs32 (calcbuff [i]) + abs32 (calcbuff [i+cnxt->longest]);

    // if silence return longest period, else calculate scaler based on largest sum

    if (sum)
        scaler = (MAX_CORR - 1) / sum;
    else
        return cnxt->longest;

    /* accumulate sum for shortest period size */

    for (sum = i = 0; i < period; ++i)
        sum += abs32 (calcbuff [i]) + abs32 (calcbuff [i+period]);

    /* this loop actually cycles through all period lengths */

    while (1) {
        int16_t *comp = calcbuff + period * 2;
        int16_t *ref = calcbuff + period;

        /* compute sum of absolute differences */

        diff = 0;

        while (ref != calcbuff)
            diff += abs32 ((int32_t) *--ref - *--comp);

        /*
         * Here we calculate and store the resulting correlation
         * factor.  Note that we must watch for a difference of
         * zero, meaning a perfect match.  Also, for increased
         * precision using integer math, we scale the sum.
         */

        factor = diff ? (sum * scaler) / diff : MAX_CORR;

        if (factor >= best_factor) {
            best_factor = factor;
            best_period = period;
        }

        /* see if we're done */

        if (period * cnxt->num_chans == cnxt->longest)
            break;

        /* update accumulating sum and current period */

        sum += abs32 (calcbuff [period * 2]) + abs32 (calcbuff [period * 2 + 1]);
        period++;
    }

    return best_period * cnxt->num_chans;
}

/*
 * This pitch detection function is similar to find_period() above, except that it
 * is optimized for speed. The audio data corresponding to two maximum periods is
 * averaged 2:1 into the calculation buffer, and then the calulations are done
 * for every other period length. Because the time is essentially proportional to
 * both the number of samples and the number of period lengths to try, this scheme
 * can reduce the time by a factor approaching 4x. The correlation results on either
 * side of the peak are compared to calculate a more accurate center of the period.
 */

static int find_period_fast (struct stretch_cnxt *cnxt, int16_t *samples)
{
    uint32_t sum, diff, scaler, best_factor = 0;
    int period, best_period;
    int i, j;

    best_period = period = cnxt->shortest / (cnxt->num_chans * 2);

    /* first step is compressing data 2:1 into calcbuff, and calculating maximum sum */

    if (cnxt->num_chans == 2)
        for (sum = i = j = 0; i < cnxt->longest * 2; i += 4)
            sum += abs32 (cnxt->calcbuff [j++] = ((int32_t) samples [i] + samples [i+1] + samples [i+2] + samples [i+3]) >> 2);
    else
        for (sum = i = j = 0; i < cnxt->longest * 2; i += 2)
            sum += abs32 (cnxt->calcbuff [j++] = ((int32_t) samples [i] + samples [i+1]) >> 1);

    // if silence return longest period, else calculate scaler based on largest sum

    if (sum)
        scaler = (MAX_CORR - 1) / sum;
    else
        return cnxt->longest;

    /* accumulate sum for shortest period */

    for (sum = i = 0; i < period; ++i)
        sum += abs32 (cnxt->calcbuff [i]) + abs32 (cnxt->calcbuff [i+period]);

    /* this loop actually cycles through all period lengths */

    while (1) {
        int16_t *comp = cnxt->calcbuff + period * 2;
        int16_t *ref = cnxt->calcbuff + period;

        /* compute sum of absolute differences */

        diff = 0;

        while (ref != cnxt->calcbuff)
            diff += abs32 ((int32_t) *--ref - *--comp);

        /*
         * Here we calculate and store the resulting correlation
         * factor.  Note that we must watch for a difference of
         * zero, meaning a perfect match.  Also, for increased
         * precision using integer math, we scale the sum.
         */

        cnxt->results [period] = diff ? (sum * scaler) / diff : MAX_CORR;

        if (cnxt->results [period] >= best_factor) {    /* check if best yet */
            best_factor = cnxt->results [period];
            best_period = period;
        }

        /* see if we're done */

        if (period * cnxt->num_chans * 2 == cnxt->longest)
            break;

        /* update accumulating sum and current period */

        sum += abs32 (cnxt->calcbuff [period * 2]) + abs32 (cnxt->calcbuff [period * 2 + 1]);
        period++;
    }

    if (best_period * cnxt->num_chans * 2 != cnxt->shortest && best_period * cnxt->num_chans * 2 != cnxt->longest) {
        uint32_t high_side_diff = cnxt->results [best_period] - cnxt->results [best_period+1];
        uint32_t low_side_diff = cnxt->results [best_period] - cnxt->results [best_period-1];

        if ((low_side_diff + 1) / 2 > high_side_diff)
            best_period = best_period * 2 + 1;
        else if ((high_side_diff + 1) / 2 > low_side_diff)
            best_period = best_period * 2 - 1;
        else
            best_period *= 2;
    }
    else
        best_period *= 2;           /* shortest or longest use as is */

    return best_period * cnxt->num_chans;
}

/*
 * To combine the two periods into one, each corresponding pair of samples
 * are averaged with a linearly sliding scale.  At the beginning of the period
 * the first sample dominates, and at the end the second sample dominates.  In
 * this way the resulting block blends with the previous and next blocks.
 *
 * The signed values are offset to unsigned for the calculation and then offset
 * back to signed.  This is done to avoid the compression around zero that occurs
 * with calculations of this type on C implementations that round division toward
 * zero.
 *
 * The maximum period handled here without overflow possibility is 65535 samples.
 * This corresponds to a maximum calculated period of 16383 samples (2x for stereo
 * and 2x for the "2.0" version of the stretch algorithm). Since the maximum
 * calculated period is currently set for 2400 samples, we have plenty of margin.
 */

static void merge_blocks (int16_t *output, int16_t *input1, int16_t *input2, int samples)
{
    int i;

    for (i = 0; i < samples; ++i)
        output [i] = (int32_t)(((uint32_t)(input1 [i] + MERGE_OFFSET) * (samples - i) +
            (uint32_t)(input2 [i] + MERGE_OFFSET) * i) / samples) - MERGE_OFFSET;
}
