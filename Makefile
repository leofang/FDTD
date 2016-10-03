CFLAGS=-Wall -std=gnu99 -pedantic -ggdb3 #-Werror
SRCS=$(wildcard *.c)
OBJS=$(patsubst %.c, %.o, $(SRCS))
PROGRAM=FDTD

$(PROGRAM): $(OBJS)
	gcc $(CFLAGS) -o $@ $(OBJS)

%.o: %.c 
	gcc -c $(CFLAGS) $<

clean:
	rm -f $(OBJS) $(PROGRAM) *~

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)


grid.o: grid.h kv.h
kv.o: kv.h
main.o: grid.h kv.h
