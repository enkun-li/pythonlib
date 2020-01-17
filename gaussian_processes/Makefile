ifeq ($(shell hostname), ekcosm)
	PYTHON3 = python3.7
else ifeq ($(shell hostname), ek)
	PYTHON3 = pyana
endif
PYTHON2 = python

FC = gcc

GSLPATH ?= /home/ekli/lib/gsl-2.2
GSLLIB ?= $(GSLPATH)/lib
GSLINC ?= $(GSLPATH)/include

FLAGS = -Wall -lm -fpic 
FFLAGS = -shared

OBJECTS = subroutines.o

.PHONY: default all clean libsub

default: python3 cplib

all: python2 python3 cplib libsub

python2:
	$(PYTHON2) setup.py build_ext --inplace

python3:
	$(PYTHON3) setup.py build_ext --inplace

cplib:
	cp *.so /home/ekli/myworks/pythonlib/

libsub: $(OBJECTS)
	$(FC) $(FLAGS) $(FFLAGS) $(OBJECTS) -o $@.so

%.o: %.c
	$(FC) $(FLAGS) -c $*.c -o $*.o

clean:
	-rm -rf build
	-rm -rf *.pyc *.c
	-rm -rf *.so
	-rm -rf lib/*
