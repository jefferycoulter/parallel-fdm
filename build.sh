#!/bin/bash

DIR="$( cd "$( dirname "$0" )" && pwd )"

cmake -S $DIR -B $DIR/build \
    -G "Unix Makefiles" \
    -DCMAKE_C_COMPILER="/usr/bin/clang" \
    -DCMAKE_BUILD_TYPE=Debug

if [[ -z "$MAKEFLAGS" ]]
then
      export MAKEFLAGS=-j$(sysctl -n hw.ncpu)
fi

cd $DIR/build; VERBOSE=1 make