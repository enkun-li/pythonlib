FC = gcc

ifeq ($(shell hostname), ekcosm)
	GSLPATH ?= /home/ekli/lib/gsl-2.2
else ifeq ($(shell hostname), ek)
	GSLPATH ?= /home/ekli/.local/lib/gsl-2.5
endif
GSLLIB ?= $(GSLPATH)/lib
GSLINC ?= $(GSLPATH)/include

FLAGS = -Wall -fpic -g -fopenmp -lm -I$(GSLINC) 
FFLAGS = -L$(GSLLIB) -lgsl -lgslcblas

OBJECTS = subroutines.o

GP_OBJ = gaussian_process.o

OBJ_utils = utils.o kernel.o gapp.o gapp_rec.o gapp_rec_fs8.o cosmo_bk.o

.PHONY: clean libsub

default: gapp

all: libsub main test_utils gapp libgapp

## gcc subroutines.c -Wall -L /home/ekli/lib/gsl-2.2/lib -lgsl -lgslcblas -lm

libsub: $(OBJECTS)
	$(FC) $(FLAGS) $(FFLAGS) $(OBJECTS) -o $@

main: $(GP_OBJ)
	gcc -o $@ $(GP_OBJ) $(FLAGS) $(FFLAGS)

test_utils: $(OBJ_utils) test_utils.o
	$(FC) $(FLAGS) -I ./ $(OBJ_utils) $@.o -o $@

gapp: $(OBJ_utils) test_gapp.o
	$(FC) $(FLAGS) $(FFLAGS) $(OBJ_utils) test_gapp.o -o $@

pso: kernel.o gapp.o pso.o
	$(FC) $(FLAGS) $(FFLAGS) kernel.o gapp.o pso.o -o $@

getpars: kernel.o gapp.o optimization.o
	$(FC) $(FLAGS) $(FFLAGS) kernel.o gapp.o optimization.o -o $@

libgapp: $(OBJ_utils)
	$(FC) -shared -fPIC $(FFLAGS) $(OBJ_utils) -o $@.so
	cp $@.so ../../
	cp my_gapp.py ../../

test_libgapp:
	gcc test_gapp.c -L. -lgapp -L /home/ekli/lib/gsl-2.2/lib -lgsl -lgslcblas -lm -o $@

%.o: %.c
	$(FC) $(FLAGS) $(FFLAGS) -c $*.c -o $*.o

clean:
	-rm -rf *.so
	-rm -rf *.o
