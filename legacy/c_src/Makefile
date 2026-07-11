OPT := -O3
CFLAGS = $(OPT) -DHGCDTE_SUTR -DMOONLIGHT_
LDFLAGS = -lm

.c.o:
	$(CC) -c -o $*.o $(CFLAGS) $*.c

all: gsetc.x
	cp gsetc.x ../python/pfsspecsim/bin/.

gsetc.x: gsetc.o
	$(CC) gsetc.o $(LDFLAGS) -o gsetc.x

gsetc.o : modeldata.h

clean :
	$(RM) gsetc.x gsetc.o *~
