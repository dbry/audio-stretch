////////////////////////////////////////////////////////////////////////////
//                        **** AUDIO-STRETCH ****                         //
//                      Time Domain Harmonic Scaler                       //
//                    Copyright (c) 2019 David Bryant                     //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// stretch.h

// Time Domain Harmonic Compression and Expansion
//
// This library performs time domain harmonic scaling with pitch detection
// to stretch the timing of a 16-bit PCM signal (either mono or stereo) from
// 1/2 to 2 times its original length. This is done without altering any of
// its tonal characteristics.
//
// Use stereo (num_chans = 2), when both channels are from same source
// and should contain approximately similar content.
// For independent channels, prefer using multiple StretchHandle-instances.
// see https://github.com/dbry/audio-stretch/issues/6
// Multiple instances, of course, wil consume more CPU load.
// In addition, different output amounts need to be handled.

#ifndef STRETCH_H
#define STRETCH_H

#ifdef __cplusplus
extern "C" {
#endif

typedef void *StretchHandle;

StretchHandle stretch_init (int shortest_period, int longest_period, int num_chans, int fast_mode);
int stretch_samples (StretchHandle handle, const short *samples, int num_samples, short *output, float ratio);
int stretch_flush (StretchHandle handle, short *output);
void stretch_reset (StretchHandle handle);
void stretch_deinit (StretchHandle handle);

#ifdef __cplusplus
}
#endif

#endif

