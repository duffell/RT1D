
INITIAL  = readtable
HYDRO    = euler
OUTPUT   = ascii

UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
H55 = /home/install/app/hdf5
endif
ifeq ($(UNAME),Darwin)
H55 = /opt/local
endif

CC = mpicc
FLAGS = -O3 -Wall -g

INC = -I$(H55)/include
LIB = -L$(H55)/lib -lm -lhdf5

OBJ = main.o mpisetup.o profiler.o readpar.o domain.o gridsetup.o geometry.o exchange.o misc.o timestep.o onestep.o riemann.o boundary.o plm.o $(INITIAL).o $(OUTPUT).o $(HYDRO).o #report.o

default: rt1d

%.o: %.c paul.h
	$(CC) $(FLAGS) $(INC) -c $<

$(TIMESTEP).o: Timestep/$(TIMESTEP).c paul.h
	$(CC) $(FLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c paul.h
	$(CC) $(FLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c paul.h
	$(CC) $(FLAGS) $(INC) -c Hydro/$(HYDRO).c

$(OUTPUT).o : Output/$(OUTPUT).c paul.h
	$(CC) $(FLAGS) $(INC) -c Output/$(OUTPUT).c

rt1d: $(OBJ) paul.h
	$(CC) $(FLAGS) $(LIB) -o rt1d $(OBJ)

clean:
	rm -f *.o rt1d
