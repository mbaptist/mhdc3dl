################################################
#
# mhdcl Makefile
#
# MHD Convection in a layer
#
# Manuel Baptista (2004)
#
################################################

CXX = g++ -O3 -funroll-loops -fexpensive-optimizations   

DEBUG = 
#-g 

PYTHON_VERSION=$(shell python -c "import sys ; print sys.version[:3]")

INCLUDE = -I/home/mbaptist/work/codes/devel/cat -I/home/mbaptist/work/codes/devel/goops -I/home/mbaptist/work/codes/devel/lass -I/usr/include/python$(PYTHON_VERSION)
LIB = -L/home/mbaptist/work/codes/devel/cat -lcat -L/home/mbaptist/work/codes/devel/goops -lgoops -L/home/mbaptist/work/codes/devel/lass -llass -lfftw3 -lpython$(PYTHON_VERSION) -lm

INSTALL_ROOT = /usr/local

AR = ar

ARFLAGS = rcs


###############################################################################

IFLAGS = $(DEBUG) $(INCLUDE)
FLAGS = $(DEBUG) $(LIB) 

OBJECTS = globals.o input.o spectral.o gen_random.o basic.o linops.o lss.o sss.o

all:  libmhdc3dl mhdc3dl_python mhdc3dl_test mhdc3dl


clean:
	@rm -rfv mhdc3dl mhdc3dl_test *.o *.so && $(MAKE) -C vzdeigen clean

distclean: clean
	@rm -rfv *~ && $(MAKE) -C vzdeigen distclean


%.o: %.C *.h
	$(CXX) $(IFLAGS) -c $< 

mhdc3dl_python.o: mhdc3dl_python.C
	$(CXX) $(IFLAGS) -pthread -fPIC -c mhdc3dl_python.C


libmhdc3dl: $(OBJECTS) 
	$(CXX) $(FLAGS) -shared -o libmhdc3dl.so $(OBJECTS) 
	$(AR) $(ARFLAGS) libmhdc3dl.a $(OBJECTS)


mhdc3dl_python: mhdc3dl_python.o libmhdc3dl
	$(CXX)  $(FLAGS) -pthread -shared -o mhdc3dl_python.so mhdc3dl_python.o -L. -lmhdc3dl

mhdc3dl_test: mhdc3dl_test.o libmhdc3dl
	$(CXX) $(FLAGS) -o mhdc3dl_test mhdc3dl_test.o -L. -lmhdc3dl

mhdc3dl: mhdc3dl.o libmhdc3dl
	$(CXX) $(FLAGS) -o mhdc3dl mhdc3dl.o -lfftw3 -L. -lmhdc3dl


install:
	install -c -m 755 mhdc3dl $(INSTALL_ROOT)/bin

uninstall:
	@(cd $(INSTALL_ROOT)/bin && rm -rfv bin/mhdc3dl)



