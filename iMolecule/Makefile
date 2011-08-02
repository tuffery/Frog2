all: ancillary Config.py
	vi Config.py
	make ancillary
ancillary: 
	cd code_c; make all; cd ..
clean:
	cd code_c; make clean; cd ..
	\rm *~ *.pyc
