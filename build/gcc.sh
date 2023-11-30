#!/usr/bin/env bash

gfortran \
    -g -O3 \
    -fstack-arrays \
    -fopenmp \
    -foffload=nvptx-none \
    -foffload-options=nvptx-none="-O3 -lm -lgfortran -latomic -march=sm_80 -mptx=7.0" \
    -foffload-options=nvptx-none=-msoft-stack \
    -foffload-options=nvptx-none=-msoft-stack-reserve-local=1024 \
    -o microphys_driver.x \
    input.F90 output.F90 gfdl_cloud_microphys.F90 driver.F90

    # -D__GFORTRAN_TEST__ \

