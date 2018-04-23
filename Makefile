all :
	cd src; make
#	./set_path.py
clean :
	cd src; make clean
	cd python/pfsspecsim/bin; make clean
