.c.o:
	$(CC) -c -o $*.o $(CFLAGS) $*.c

OPT := -O3
CFLAGS = $(OPT) -DHGCDTE_SUTR -DMOONLIGHT_
LDFLAGS = -lm

../python/pfsspecsim/bin/gsetc.x : gsetc.o
	$(CC) gsetc.o $(LDFLAGS) -o ../python/pfsspecsim/bin/gsetc.x
gsetc.o : modeldata.h

clean :
	$(RM) gsetc.o *~
