#!/bin/bash

# directory where current shell script resides
PROJECTDIR=$(dirname "$BASH_SOURCE")

cd "$PROJECTDIR"

mode=$1 # debug/release mode
shift   # shift command-line arguments
        # the rest are cmake command-line arguments

mkdir -p build && cd build

# if debug directory does not exist, create it
mkdir -p debug
# if release directory does not exist, create it
mkdir -p release

# if debug mode
if [ "$mode" = "debug" ]; then
    cd debug
# if release mode
elif [ "$mode" = "release" ]; then
    cd release
else
    echo "Usage: ./clean-build.sh <debug/release> [cmake options]" 1>&2
    exit 1
fi

# directory with R include files
R_HOME=$(Rscript -e "cat(Sys.getenv(\"R_HOME\"))")
R_INCLUDE_DIR=$(Rscript -e "cat(Sys.getenv(\"R_INCLUDE_DIR\"))")
R_LIB_DIR="$R_HOME/lib"

rm -rf *
cmake ../../ \
    -DR_INCLUDE_DIR=$R_INCLUDE_DIR \
    -DR_LIB_DIR=$R_LIB_DIR \
    "$@"
cmake --build . -- -j12
