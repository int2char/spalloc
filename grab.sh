#!/bin/bash
nvcc -O3 -pg -std=c++11 *.cpp *.cu --gpu-architecture=compute_35 --gpu-code=sm_35 -I ./include
./a.out
#gprof -b a.out gmon.out >report.txt
