############################################################################

#vzdeigen - find dominant eigenvalue of linear operator
#Original code by Vlad Zheligovsky
#Interface improvements by Manuel Baptista (2005)

############################################################################

F90 = ifort

F90FLAGS = -fast

AR = ar

ARFLAGS = rcs 

############################################################################

all: vzdeigen.o 

clean:
	@rm -rfv *.o *.a *.so

distclean: clean	
	@rm -rfv *~

vzdeigen.o: vzdeigen.f 
	$(F90) $(F90FLAGS) -c vzdeigen.f -o vzdeigen.o

