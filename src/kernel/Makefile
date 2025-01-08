# Builds kernel.so
#
# Reference: https://makefiletutorial.com/

CC = gcc-14 
CFLAGS = -O3 -march=native -ffast-math -shared -fPIC

kernel.so:
	${CC} ${CFLAGS} -o kernel.so kernel.c