############################################################################

#vzdeigen - find dominant eigenvalue of linear operator
#Original code by Vlad Zheligovsky
#Interface improvements by Manuel Baptista (2005)

############################################################################

F90 = ifort

F90FLAGS =

AR = ar

ARFLAGS = rcs 

############################################################################


all: libvzdeigen

clean:
	@rm -rfv *.o *~ *.a *.so

distclean: clean
	@rm -rfv *~

vzdeigen.o: vzdeigen.f
	$(F90) $(F90FLAGS) -static-libcxa -c vzdeigen.f

libvzdeigen: vzdeigen.o
	$(F90) $(F90FLAGS) -shared -o libvzdeigen.so vzdeigen.o
	$(AR) $(ARFLAGS) libvzdeigen.a vzdeigen.o 


#all: vzdeigen.o 
#
#clean:
#	@rm -rfv *.o *.a *.so
#
#distclean: clean	
#	@rm -rfv *~

#vzdeigen.o: vzdeigen.f 
#	$(F90) $(F90FLAGS) -c vzdeigen.f -o vzdeigen.o


