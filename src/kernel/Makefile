# Builds kernel.so
#
# Reference: https://makefiletutorial.com/

CC = gcc
CFLAGS = -O3 -ffast-math -march=native

ckernel.so: kernel.c
	${CC} -shared ${CFLAGS} -o ckernel.so -fPIC kernel.c