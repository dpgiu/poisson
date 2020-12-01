
FC = gfortran
TARGET = poisson 

SOURCES = $(wildcard *.f90)
OBJS = $(SOURCES:.f90=.o)

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o: 
	$(FC) -c  $<

all: $(OBJS)
	$(FC) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	rm *.o *.mod *~ $(TARGET)

distclean: clean
	rm $(TARGET)

poisson.o : precision.o

