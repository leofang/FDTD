# Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
#
# This program is free software. It comes without any warranty,
# to the extent permitted by applicable law. You can redistribute
# it and/or modify it under the terms of the WTFPL, Version 2, as
# published by Sam Hocevar. See the accompanying LICENSE file or
# http://www.wtfpl.net/ for more details.

CFLAGS=-Wall -std=gnu99 -pedantic -O3 #-ggdb3 -Werror
SRCS=$(wildcard *.c)
OBJS=$(patsubst %.c, %.o, $(SRCS))
PROGRAM=FDTD
LDFLAGS=-lm

$(PROGRAM): $(OBJS)
	gcc $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.c 
	gcc -c $(CFLAGS) $<

clean:
	rm -f $(OBJS) $(PROGRAM) *~

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)


# DO NOT DELETE

dynamics.o: dynamics.h grid.h kv.h
grid.o: kv.h grid.h
kv.o: kv.h
main.o: grid.h kv.h dynamics.h
