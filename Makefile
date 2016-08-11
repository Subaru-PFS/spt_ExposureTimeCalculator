all :
	cd src; make
	./set_path.sh
clean :
	cd src; make clean
	cd bin; make clean
