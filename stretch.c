////////////////////////////////////////////////////////////////////////////
//                        **** AUDIO-STRETCH ****                         //
//                      Time Domain Harmonic Scaler                       //
//                    Copyright (c) 2019 David Bryant                     //
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

#include "stretch.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

typedef stretch_sample_t sample_t; /* input and output sample type - from header file */
typedef uint32_t         usum_t;   /* for sum of abs(samples) over period lengths */
typedef int32_t          isum_t;   /* for product in merge_blocks() */
typedef uint32_t         corr_t;   /* for correlation results */

#define N_SAMPLE_BITS    16
#define MERGE_SAMPLE_OFF 32768     /* = 1 << (N_SAMPLE_BITS -1) ; see merge_blocks() */
#define N_SUM_BITS       32
#define N_CORR_MAX_BITS  32
#define N_MIN_CORR_PREC  4         /* minimum precision bits for correlation with fixed-point integer division */
#define N_MAX_CORR_PREC  12        /* maximum precision bits to use */
#define N_DFLT_CORR_PREC 7         /* default/desired precision bits */

#define CORR_MAX         (UINT32_MAX)

#define MIN_PERIOD       24        /* minimum allowable pitch period */
#define MAX_PERIOD       1024      /* maximum allowable pitch period */

/* depending on above definitions, we can increase MAX_PERIOD */
/* #define MAX_PERIOD  ( 1 << MIN( N_SUM_BITS - N_SAMPLE_BITS, N_CORR_MAX_BITS - N_SAMPLE_BITS - N_MIN_CORR_PREC ) ) */


struct stretch_cnxt {
    int num_chans, inbuff_samples, shortest, longest, tail, head, fast_mode, prec_shift;
    int16_t *inbuff, *calcbuff;
    float outsamples_error;
};

static void merge_blocks (sample_t *output, const sample_t *input1, const sample_t *input2, int samples);
static int find_period_fast (struct stretch_cnxt *cnxt, const sample_t *samples);
static int find_period (struct stretch_cnxt *cnxt, const sample_t *samples);

/*
 * Initialize a context of the time stretching code. The shortest and longest periods
 * are specified here. The longest period determines the lowest fundamental frequency
 * that can be handled correctly. Note that higher frequencies can be handled than the
 * shortest period would suggest because multiple periods can be combined, and the
 * worst-case performance will suffer if too short a period is selected.
 */

StretchHandle stretch_init (int shortest_period, int longest_period, int num_channels, int fast_mode)
{
    struct stretch_cnxt *cnxt;

    if ( (8 * sizeof(sample_t)) != N_SAMPLE_BITS
         || (8 * sizeof(usum_t)) != N_SUM_BITS
         || (8 * sizeof(corr_t)) < N_CORR_MAX_BITS )
    {
        fprintf (stderr, "stretch_init(): invalid bit sizes of types!\n");
        return NULL;
    }

    if (num_channels < 1 || num_channels > 2) {
        fprintf (stderr, "stretch_init(): invalid number of channels! only 1 (=mono) or 2 (=stereo) supported.\n");
        return NULL;
    }

    if (fast_mode) {
        longest_period = (longest_period + 1) & ~1;
        shortest_period &= ~1;
    }

    if (longest_period < 2 * shortest_period || shortest_period < MIN_PERIOD || longest_period > MAX_PERIOD) {
        fprintf (stderr, "stretch_init(): invalid periods!\n");
        return NULL;
    }

    cnxt = (struct stretch_cnxt *) calloc (1, sizeof (struct stretch_cnxt));

    if (cnxt) {
        cnxt->inbuff_samples = longest_period * num_channels * 8;
        cnxt->inbuff = calloc (cnxt->inbuff_samples, sizeof (*cnxt->inbuff));
        cnxt->calcbuff = calloc (longest_period * num_channels, sizeof (*cnxt->calcbuff));
    }

    if (!cnxt || !cnxt->inbuff || !cnxt->calcbuff) {
        fprintf (stderr, "stretch_init(): out of memory!\n");
        return NULL;
    }

    cnxt->head = cnxt->tail = cnxt->longest = longest_period * num_channels;
    {
        int sum_bits, corr_bits, longest_lg2;
        longest_lg2 = 0;  /* longest_lg2 := next power of 2 for cnxt->longest */
        while ( (1 << longest_lg2) < cnxt->longest )
            ++longest_lg2;
        sum_bits = N_SAMPLE_BITS + longest_lg2;  /* bits in 'sum' of find_period() and find_period_fast() */
        assert( sum_bits <= N_SUM_BITS );
        corr_bits = (sum_bits <= N_CORR_MAX_BITS) ? (N_CORR_MAX_BITS - sum_bits) : 0;
        corr_bits = ( corr_bits > N_MAX_CORR_PREC ) ? N_MAX_CORR_PREC : corr_bits;
        assert( corr_bits >= N_MIN_CORR_PREC );
        cnxt->prec_shift = corr_bits;
    }
    cnxt->shortest = shortest_period * num_channels;
    cnxt->num_chans = num_channels;
    cnxt->fast_mode = fast_mode;
    cnxt->outsamples_error = 0.0F;

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
}


/*
 * Process the specified samples with the given ratio (which is clipped to the
 * range 0.5 to 2.0). Note that the number of samples refers to total samples for
 * both channels in stereo and can be as large as desired (samples are buffered
 * here). The exact number of samples output is not possible to determine in
 * advance, but it will be close to the number of input samples times the ratio
 * plus or minus 3X the longest period.
 */

int stretch_samples (StretchHandle handle, const sample_t *samples, int num_samples, sample_t *output, float ratio)
{
    struct stretch_cnxt *cnxt = (struct stretch_cnxt *) handle;
    int out_samples = 0;

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

        while (cnxt->tail >= cnxt->longest && cnxt->head - cnxt->tail >= cnxt->longest * 2) {
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
                merge_blocks (output + out_samples, cnxt->inbuff + cnxt->tail,
                    cnxt->inbuff + cnxt->tail + period, period);
                cnxt->outsamples_error += period - (period * 2.0 * ratio);
                out_samples += period;
                cnxt->tail += period * 2;
            }
            else if (process_ratio == 1.0) {
                memcpy (output + out_samples, cnxt->inbuff + cnxt->tail, period * 2 * sizeof (cnxt->inbuff [0]));
                cnxt->outsamples_error += (period * 2.0) - (period * 2.0 * ratio);
                out_samples += period * 2;
                cnxt->tail += period * 2;
            }
            else if (process_ratio == 1.5) {
                memcpy (output + out_samples, cnxt->inbuff + cnxt->tail, period * sizeof (cnxt->inbuff [0]));
                merge_blocks (output + out_samples + period, cnxt->inbuff + cnxt->tail + period,
                    cnxt->inbuff + cnxt->tail, period);
                memcpy (output + out_samples + period * 2, cnxt->inbuff + cnxt->tail + period, period * sizeof (cnxt->inbuff [0]));
                cnxt->outsamples_error += (period * 3.0) - (period * 2.0 * ratio);
                out_samples += period * 3;
                cnxt->tail += period * 2;
            }
            else if (process_ratio == 2.0) {
                merge_blocks (output + out_samples, cnxt->inbuff + cnxt->tail,
                    cnxt->inbuff + cnxt->tail - period, period * 2);

                cnxt->outsamples_error += (period * 2.0) - (period * ratio);
                out_samples += period * 2;
                cnxt->tail += period;
            }
            else
                fprintf (stderr, "stretch_samples: fatal programming error: process_ratio == %g\n", process_ratio);
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

    return out_samples / cnxt->num_chans;
}

/* flush any leftover samples out at normal speed */

int stretch_flush (StretchHandle handle, sample_t *output)
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
    free (cnxt->inbuff);
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

static int find_period (struct stretch_cnxt *cnxt, const sample_t *samples)
{
    const int prec_shift = cnxt->prec_shift;
    corr_t factor, best_factor = 0;
    usum_t sum, diff;
    const sample_t *calcbuff = samples;
    int period, best_period;
    int i, j;

    period = best_period = cnxt->shortest / cnxt->num_chans;

    if (cnxt->num_chans == 2) {
        sample_t *buff = cnxt->calcbuff;

        for (i = j = 0; i < cnxt->longest * 2; i += 2)
            buff [j++] = ((isum_t)samples [i] + samples [i+1]) >> 1;

        calcbuff = buff;
    }

    /* accumulate sum for shortest period size */

    for (sum = i = 0; i < period; ++i)
        sum += abs (calcbuff [i]) + abs (calcbuff [i+period]);

    /* this loop actually cycles through all period lengths */

    while (1) {
        const sample_t *comp = calcbuff + period * 2;
        const sample_t *ref = calcbuff + period;

        /* compute sum of absolute differences */

        diff = 0;

        while (ref != calcbuff)
            diff += abs (*--ref - *--comp);

        /*
         * Here we calculate and store the resulting correlation
         * factor.  Note that we must watch for a difference of
         * zero, meaning a perfect match.  Also, for increased
         * precision using integer math, we scale the sum.  Care
         * must be taken here to avoid possibility of overflow.
         * see prec_shift calculation to avoid overflow.
         */

        factor = diff ? (((corr_t)sum) << prec_shift) / diff : (sum ? CORR_MAX : 1);

        /* check with <= : prefer longer periods */
        if (factor >= best_factor) {
            best_factor = factor;
            best_period = period;
        }

        /* see if we're done */

        if (period * cnxt->num_chans == cnxt->longest)
            break;

        /* update accumulating sum and current period */

        sum += abs (calcbuff [period * 2]);
        sum += abs (calcbuff [period * 2 + 1]);
        ++period;
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

static int find_period_fast (struct stretch_cnxt *cnxt, const sample_t *samples)
{
    const int prec_shift = cnxt->prec_shift;
    corr_t best_factor = 0, best_factors_prev = 0, best_factors_next = 0;
    corr_t prev_result = 0, curr_result = 0;
    usum_t sum, diff;
    int period, best_period;
    int i, j;

    best_period = period = cnxt->shortest / (cnxt->num_chans * 2);

    /* first step is compressing data 2:1 into calcbuff */

    if (cnxt->num_chans == 2)
        for (i = j = 0; i < cnxt->longest * 2; i += 4)
            cnxt->calcbuff [j++] = ((isum_t)samples [i] + samples [i+1] + samples [i+2] + samples [i+3]) >> 2;
    else
        for (i = j = 0; i < cnxt->longest * 2; i += 2)
            cnxt->calcbuff [j++] = ((isum_t)samples [i] + samples [i+1]) >> 1;

    /* accumulate sum for shortest period */

    for (sum = i = 0; i < period; ++i)
        sum += abs (cnxt->calcbuff [i]) + abs (cnxt->calcbuff [i+period]);

    /* this loop actually cycles through all period lengths */

    while (1) {
        const sample_t *comp = cnxt->calcbuff + period * 2;
        const sample_t *ref = cnxt->calcbuff + period;

        /* compute sum of absolute differences */

        diff = 0;

        while (ref != cnxt->calcbuff)
            diff += abs (*--ref - *--comp);

        /*
         * Here we calculate and store the resulting correlation
         * factor.  Note that we must watch for a difference of
         * zero, meaning a perfect match.  Also, for increased
         * precision using integer math, we scale the sum.  Care
         * must be taken here to avoid possibility of overflow.
         * see prec_shift calculation to avoid overflow.
         */

        prev_result = curr_result;           /* shift results */

        curr_result = diff ? (((corr_t)sum) << prec_shift) / diff : (sum ? CORR_MAX : 1);  /* for current period */


        if ( period - 1 == best_period )
            best_factors_next = curr_result; /* complete necessary data for best_factor */

        /* check with <= : prefer longer periods */
        if (curr_result >= best_factor) {    /* check if best yet */
            best_factors_prev = prev_result;
            best_factors_next = 0;           /* in case we are at last iteration */
            best_factor = curr_result;
            best_period = period;
        }

        /* see if we're done */

        if (period * cnxt->num_chans * 2 == cnxt->longest)
            break;

        /* update accumulating sum and current period */

        sum += abs (cnxt->calcbuff [period * 2]);
        sum += abs (cnxt->calcbuff [period * 2 + 1]);
        ++period;
    }

    if (best_period * cnxt->num_chans * 2 != cnxt->shortest && best_period * cnxt->num_chans * 2 != cnxt->longest) {
        corr_t high_side_diff = best_factor - best_factors_next;
        corr_t low_side_diff = best_factor - best_factors_prev;

        if (low_side_diff / 2 > high_side_diff)
            best_period = best_period * 2 + 1;
        else if (high_side_diff / 2 > low_side_diff)
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
 */

static void merge_blocks (sample_t *output, const sample_t *input1, const sample_t *input2, int num_samples)
{
    /* check / predetermine possible temporary overflow!:
     * num_samples <= cnxt->longest  --> log2(MAX_PERIOD) <= 10
     * temporary product + sum requires  16 + 1 + 10 = 27 Bits => sufficient.
     */
    for (int i = 0; i < num_samples; ++i)
        output [i] = (
        (  (isum_t)(input1 [i] + MERGE_SAMPLE_OFF) * (num_samples - i)
         + (isum_t)(input2 [i] + MERGE_SAMPLE_OFF) * i
        ) / num_samples ) - MERGE_SAMPLE_OFF;
}

