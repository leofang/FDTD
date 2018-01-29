# Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
#
# This program is free software. It comes without any warranty,
# to the extent permitted by applicable law. You can redistribute
# it and/or modify it under the terms of the WTFPL, Version 2, as
# published by Sam Hocevar. See the accompanying LICENSE file or
# http://www.wtfpl.net/ for more details.

CC=gcc-5
OpenMP=
CFLAGS=-Wall -std=gnu99 -pedantic -O3 #-ggdb3 -Werror
SRCS=$(wildcard *.c)
OBJS=$(patsubst %.c, %.o, $(SRCS))
PROGRAM=FDTD
LDFLAGS=-lm

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(OpenMP) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.c 
	$(CC) -c $(CFLAGS) $(OpenMP) $<

clean:
	rm -f $(OBJS) $(PROGRAM) *~

depend:
	makedepend -Y -- $(CFLAGS) -- $(SRCS)


# DO NOT DELETE

dynamics.o: dynamics.h grid.h kv.h
grid.o: kv.h grid.h special_function.h dynamics.h NM_measure.h
kv.o: kv.h
main.o: grid.h kv.h dynamics.h NM_measure.h
NM_measure.o: NM_measure.h grid.h kv.h special_function.h dynamics.h
special_function.o: special_function.h
