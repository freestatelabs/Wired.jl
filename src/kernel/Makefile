# Builds kernel.so
#
# Reference: https://makefiletutorial.com/

CC = gcc
CFLAGS = -O3 -march=native -ffast-math 

kernel.so:
	${CC} -shared ${CFLAGS} -o kernel.so -FPIC kernel.c