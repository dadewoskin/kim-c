CC = gcc
LIBS = -lm
PROGRAM = kimFE
LFLAGS = -fopenmp
CFLAGS = -c -fopenmp
SOURCES = kimFE.c rhs.c
DEPS = rhs.h
OBJS = $(SOURCES:.c=.o)

$(PROGRAM): $(OBJS)
	$(CC) $(LFLAGS) $(LIBS) $(OBJS) -o $(PROGRAM)

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $<

clean:
	rm $(OBJS) $(PROGRAM)

tar:
	tar cfv $(PROGRAM).tar *.c *.h Makefile