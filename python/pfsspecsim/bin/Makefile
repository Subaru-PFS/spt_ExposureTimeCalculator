.c.o:
	$(CC) -c -o $*.o $(CFLAGS) $*.c

OPT := -O3
CFLAGS = $(OPT) -DHGCDTE_SUTR -DMOONLIGHT_
LDFLAGS = -lm

gsetc.x :
	cd ../src; make ../bin/gsetc.x

clean :
	$(RM) ../bin/gsetc.x *~
