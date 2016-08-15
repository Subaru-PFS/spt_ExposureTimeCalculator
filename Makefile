all :
	cd src; make
	./set_path.py
clean :
	cd src; make clean
	cd bin; make clean
