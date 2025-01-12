# Builds kernel.so
#
# Reference: https://makefiletutorial.com/

all: wires_sp.so wires_dp.so rings_sp.so rings_dp.so
CC = gcc
CFLAGS = -O3 -ffast-math -march=native

wires_sp.so: wires_sp.c
	${CC} -shared ${CFLAGS} -o wires_sp.so -fPIC wires_sp.c

wires_dp.so: wires_dp.c
	${CC} -shared ${CFLAGS} -o wires_dp.so -fPIC wires_dp.c

rings_sp.so: rings_sp.c
	${CC} -shared ${CFLAGS} -o rings_sp.so -fPIC rings_sp.c

rings_dp.so: rings_dp.c
	${CC} -shared ${CFLAGS} -o rings_dp.so -fPIC rings_dp.c