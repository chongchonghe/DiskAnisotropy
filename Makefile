FC = gfortran
FFLAGS = -O2
#FFLAGS = -std=f2003
#FFLAGS = -std=f95 -Wall
#FFLAGS = -g -fcheck=all -std=f95 -Wall -Wextra

PROGRAMS = main

OBJECTS = vectors.o main.o

all: $(PROGRAMS)

main: $(OBJECTS)
	$(FC) -o main $(OBJECTS)

%.o : %.f95
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod

veryclean: clean
	rm -f *~ $(PROGRAMS)
