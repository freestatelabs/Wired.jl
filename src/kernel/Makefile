# Builds kernel.so
#
# Reference: https://makefiletutorial.com/

CC = gcc
CFLAGS = -O3 -march=native -ffast-math 

ckernel.so: kernel.c
	${CC} -shared ${CFLAGS} -o ckernel.so -fPIC kernel.c