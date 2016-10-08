CFLAGS=-Wall -std=gnu99 -pedantic -Wno-comment -ggdb3 #-O3 -Werror
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


grid.o: grid.h kv.h
kv.o: kv.h
main.o: grid.h kv.h dynamics.h
