# Builds kernel.so
#
# Reference: https://makefiletutorial.com/

CC = gcc-14
CFLAGS = -O3 -ffast-math -march=native -fopt-info-vec-missed

ckernel.so: kernel.c
	${CC} -shared ${CFLAGS} -o ckernel.so -fPIC kernel.c