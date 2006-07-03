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
LIB = -L/home/mbaptist/work/codes/devel/cat -lcat -L/home/mbaptist/work/codes/devel/goops -lgoops -lfftw3 -lpython$(PYTHON_VERSION) -lm

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

libmhdc3dl: $(OBJECTS) libvzdeigen
	$(CXX) $(FLAGS) -shared -o libmhdc3d.so $(OBJECTS) ./vzdeigen/libvzdeigen.a -L/opt/intel/fc/9.0/lib -lifcore -limf

libmhdc3dl_static: $(OBJECTS) libvzdeigen
	$(AR) $(ARFLAGS) libmhdc3d.a $(OBJECTS) ./vzdeigen/libvzdeigen.a -L/opt/intel/fc/9.0/lib -lifcore -limf

mhdc3dl_python: libmhdc3dl mhdc3dl_python.o
	$(CXX)  $(FLAGS) -pthread -shared -o mhdc3dl_python.so mhdc3dl_python.o -lmhdc3dl -lfftw3 -lm -lgoops -lcat

mhdc3dl_python.o: mhdc3dl_python.C
	$(CXX) $(IFLAGS) -pthread -fPIC -c mhdc3dl_python.C 

mhdc3dl_test: libmhdc3dl mhdc3dl_test.o
	$(CXX) $(FLAGS) -o mhdc3dl_test mhdc3dl_test.o -lfftw3 -lmhdc3dl

mhcd3dl: libmhdc3dl mhdc3dl.o
	$(CXX) $(FLAGS) -o mhdc3dl mhdc3dl.o -lfftw3 -lmhdc3dl

libvzdeigen:
	$(MAKE) -C vzdeigen

%.o: %.C *.h
	$(CXX) $(IFLAGS) -c $<

install:
	install -c -m 755 mhdc3dl $(INSTALL_ROOT)/bin

uninstall:
	@(cd $(INSTALL_ROOT)/bin && rm -rfv bin/mhdc3dl)



