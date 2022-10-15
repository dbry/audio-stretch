#!/bin/bash

if [ ! -d output ]; then
  echo "creating directory output"
  mkdir output
fi

if [ ! -f samples/mono.wav ] || [ ! -f samples/stereo.wav ]; then
  WVUNPACK=$(which wvunpack)
  if [ -z "$WVUNPACK" ]; then
    echo "please build/install WavPack with wvunpack to convert .wv samples to .wav"
    exit 1
  fi
  $WVUNPACK -w samples/mono.wv
  $WVUNPACK -w samples/stereo.wv
fi

STARTER=""
if [ "$1" = "gdb" ]; then
  STARTER="gdb -q -ex run -ex quit --args"
  shift
fi

EXAMPLE="mono"
if [ "$1" = "mono" ]; then
  EXAMPLE="$1"
  shift
fi
if [ "$1" = "stereo" ]; then
  EXAMPLE="$1"
  shift
fi


if [ -z "$1" ] && [ -z "$2" ]; then
  echo "usage: $0 [mono|stereo] [f|n] [s|x]"
  echo "  'f': fast pitch detection"
  echo "  'n': normal pitch detection"
  echo "  's': simple range for ratio: 0.5 .. 2.0"
  echo "  'x': extended range for ratio: 0.25 .. 4.0"
  echo ""
fi

if [ -z "$1" ] || [ "$1" = "f" ]; then
  echo "testing with fast pitch detection"
  FO="-f"
  FN="f"
else
  echo "testing with normal pitch detection"
  FO="-n"
  FN="n"
fi


if [ -z "$2" ] || [ "$2" = "s" ]; then
  echo ""
  echo "testing normal range 0.5 .. 2.0"
  echo "x2.0"
  $STARTER ./audio-stretch -q -y $FO -r0.5   samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r050_x200.wav
  echo "x1.75"
  $STARTER ./audio-stretch -q -y $FO -r0.571 samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r057_x175.wav
  echo "x1.5"
  $STARTER ./audio-stretch -q -y $FO -r0.666 samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r066_x150.wav
  echo "x1.25"
  $STARTER ./audio-stretch -q -y $FO -r0.8   samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r080_x125.wav
  echo "x1.0"
  $STARTER ./audio-stretch -q -y $FO -r1.0   samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r100_x100.wav
  echo "x0.75"
  $STARTER ./audio-stretch -q -y $FO -r1.333 samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r133_x075.wav
  echo "x0.5"
  $STARTER ./audio-stretch -q -y $FO -r2.0   samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r200_x050.wav
fi
if [ -z "$2" ] || [ "$2" = "x" ]; then
  echo ""
  echo "testing extended range 0.25 .. 0.5 and 2.0 .. 4.0"
  echo "x4.0"
  $STARTER ./audio-stretch -q -y $FO -r0.25  samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r025_x400.wav
  echo "x3.5"
  $STARTER ./audio-stretch -q -y $FO -r0.285 samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r028_x350.wav
  echo "x3.0"
  $STARTER ./audio-stretch -q -y $FO -r0.333 samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r033_x300.wav
  echo "x2.5"
  $STARTER ./audio-stretch -q -y $FO -r0.4   samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r040_x250.wav
  echo "x0.4"
  $STARTER ./audio-stretch -q -y $FO -r2.5   samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r250_x040.wav
  echo "x0.333"
  $STARTER ./audio-stretch -q -y $FO -r3.0   samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r300_x033.wav
  echo "x0.285"
  $STARTER ./audio-stretch -q -y $FO -r3.5   samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r350_x028.wav
  echo "x0.25"
  $STARTER ./audio-stretch -q -y $FO -r4.0   samples/${EXAMPLE}.wav output/out_${EXAMPLE}_${FN}_r400_x025.wav
fi
