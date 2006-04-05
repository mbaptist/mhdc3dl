############################################################################

#vzdeigen - find dominant eigenvalue of linear operator
#Original code by Vlad Zheligovsky
#Interface improvements by Manuel Baptista (2005)

############################################################################

F90 = ifort

F90FLAGS = -xB -g 

AR = ar

ARFLAGS = rcs 

############################################################################

all: libvzdeigen.a libvzdeigen.so

clean:
	@rm -rfv *.o *~ *.a *.so

vzdeigen.o: vzdeigen.f 
	$(F90) $(F90FLAGS) -c vzdeigen.f

libvzdeigen.so: vzdeigen.o
	$(F90) $(F90FLAGS) -shared -o libvzdeigen.so vzdeigen.o

vzdeigen_static.o: vzdeigen.f
	$(F90) $(F90FLAGS) -o vzdeigen_static.o -c vzdeigen.f

libvzdeigen.a: vzdeigen_static.o
	$(AR) $(ARFLAGS) libvzdeigen.a vzdeigen_static.o

