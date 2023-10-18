all :
	cd src && make -f Makefile all
	cd src && make -f Makefile.omp all

clean :
	cd src && make -f Makefile clean
	cd src && make -f Makefile.omp clean
	cd python/pfsspecsim/bin && make clean
