FC = gcc

ifeq ($(shell hostname), ekcosm)
	GSLPATH ?= /home/ekli/lib/gsl-2.2
else ifeq ($(shell hostname), ek)
	GSLPATH ?= /home/ekli/.local/lib/gsl-2.5
endif
GSLLIB ?= $(GSLPATH)/lib
GSLINC ?= $(GSLPATH)/include

SRC = ./source
LIB = ./lib

FLAGS = -Wall -fpic -g -fopenmp -lm -I$(GSLINC) -I./include
FFLAGS = -L$(GSLLIB) -lgsl -lgslcblas -L$(LIB)

OBJ_utils = $(LIB)/utils.o $(LIB)/gsl_2d_integrate.o \
			$(LIB)/kernel.o $(LIB)/gapp.o \
			$(LIB)/gapp_rec.o $(LIB)/gapp_sn_lens.o \
			$(LIB)/cosmo_bk.o $(LIB)/cosmography.o \
		   	$(LIB)/pantheon.o $(LIB)/stronglens.o \
			$(LIB)/gapp_rec_SN.o $(LIB)/cosmo_gp.o
OBJ_opt = $(LIB)/utils.o $(LIB)/kernel.o $(LIB)/gapp.o

.PHONY: clean libsub

default: gapp

all: gapp libgapp

gapp: $(OBJ_utils) $(LIB)/test_gapp.o
	$(FC) $(FLAGS) $(FFLAGS) $(OBJ_utils) $(LIB)/test_gapp.o -o $@

libgapp: $(OBJ_utils)
	$(FC) -shared -fPIC $(FFLAGS) $(OBJ_utils) -o $(LIB)/$@.so
	cp $(LIB)/libgapp.so ../
	cp $(LIB)/my_gapp.py ../

$(LIB)/%.o: $(SRC)/%.c
	$(FC) $(FLAGS) $(FFLAGS) -c $(SRC)/$*.c -o $(LIB)/$*.o

getpars: $(OBJ_opt) $(LIB)/optimization.o
	$(FC) $(FLAGS) $(FFLAGS) $(OBJ_opt) $(LIB)/optimization.o -o $@

pso: $(OBJ_opt) $(LIB)/pso.o
	$(FC) $(FLAGS) $(FFLAGS) $(OBJ_opt) $(LIB)/pso.o -o $@

clean:
	-rm -rf $(LIB)/*.o 
	-rm -rf $(LIB)/*.so
	-rm -rf *.so
	-rm -rf *.o
