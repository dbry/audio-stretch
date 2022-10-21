#!/bin/bash

if [ -z "$1" ] || [ "$1" = "rel" ]; then
  echo "building release .."
  gcc -Ofast main.c stretch.c -lm -o audio-stretch
elif [ "$1" = "dbg" ]; then
  echo "building debug .."
  gcc -O0 -g main.c stretch.c -lm -o audio-stretch
elif [ "$1" = "ubsan" ]; then
  echo "building debug with undefined behaviour sanitizer .."
  gcc -O0 -g main.c stretch.c -fsanitize=undefined -lm -o audio-stretch
elif [ "$1" = "asan" ]; then
  echo "building debug with address sanitizer .."
  gcc -O0 -g main.c stretch.c -fsanitize=address -lm -o audio-stretch
else
  echo "error: unknown option '$1'"
fi

