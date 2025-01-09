CC = gcc
CFLAGS = -O3 -march=native -ffast-math 

test.so: test.c
	${CC} -shared ${CFLAGS} -o test.so -fPIC test.c
