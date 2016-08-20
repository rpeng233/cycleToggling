#!/bin/bash
rm exec/$1
if [ "$2" == '-g' ]; then
	g++ source/$1.cpp -o exec/$1 -Wall -Wextra -g -lmcheck -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC --std=c++11  -L/usr/local/lib   -lnsl -lz -lm -lquadmath
fi
if [ "$2" == '-O3' ]; then
	g++ source/$1.cpp -o exec/$1 -O3 --std=c++11 -L/home/kdeweese/build/gcc/4.8.2/lib64  -lnsl -lz -lm -lquadmath
fi


