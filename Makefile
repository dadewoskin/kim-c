CC = icc
#LIBS = -lm
PROGRAM = kimFE
LFLAGS = 
CFLAGS = -O3 -c
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
