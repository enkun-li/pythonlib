#PYTHON3 = /home/ekli/.local/anaconda3/bin/python3.7
PYTHON3 = /home/pub/.local/bin/python3.7
PYTHON2 = python

.PHONY: default all clean

default: python3

all: python2 python3 cplib

python2:
	$(PYTHON2) setup.py build_ext --inplace

python3:
	$(PYTHON3) setup.py build_ext --inplace

cplib:
	cp *.so /home/ekli/myworks/pythonlib/

clean:
	-rm -rf build
	-rm -rf *.pyc *.c
	-rm -rf *.so
	-rm -rf lib/*
