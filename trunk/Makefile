################################################
#
# mhdcl Makefile
#
# MHD Convection in a layer
#
# Manuel Baptista (2004)
#
################################################

CXX = g++ -O3 -funroll-loops -fexpensive-optimizations -g  

DEBUG = 
#-g 

PYTHON_VERSION=$(shell python -c "import sys ; print sys.version[:3]")

INCLUDE = -I/home/mbaptist/work/codes/devel/cat -I/home/mbaptist/work/codes/devel/goops -I/home/mbaptist/work/codes/devel/lass -I/usr/include/python$(PYTHON_VERSION)
LIB = -L/home/mbaptist/work/codes/devel/cat -lcat -L/home/mbaptist/work/codes/devel/goops -lgoops -lfftw3 -lpython$(PYTHON_VERSION) -lm

INSTALL_ROOT = /usr/local

###############################################################################

IFLAGS = $(DEBUG) $(INCLUDE)
FLAGS = $(DEBUG) $(LIB) 

OBJECTS = globals.o input.o spectral.o gen_random.o basic.o linops.o lss.o 
#sss.o

all:  mhdc3dl mhdc3dl_python mhdc3dl_test

clean:
	@rm -rfv mhdc3dl mhdc3dl_test *.o *.so && $(MAKE) -C vzdeigen clean

distclean: clean
	@rm -rfv *~ && $(MAKE) -C vzdeigen distclean

mhdc3dl: $(OBJECTS) mhdc3dl.o sss.o libvzdeigen
	$(CXX) $(FLAGS) -o mhdc3dl $(OBJECTS) mhdc3dl.o sss.o ./vzdeigen/libvzdeigen.a -L/opt/intel/fc/9.0/lib -lifcore -limf

#mhdc3dl_sss: libvzdeigen $(OBJECTS)
#	$(CXX) $(FLAGS) -o mhdc3dl_sss $(OBJECTS) mhdc3dl_sss.o -lfftw3 ./vzdeigen/libvzdeigen.a -L/opt/intel_fc_80/lib -lifcore -lirc -lifport 

mhdc3dl_python: $(OBJECTS) mhdc3dl.o mhdc3dl_python.o
	$(CXX)  $(FLAGS) -pthread -shared -o mhdc3dl_python.so mhdc3dl_python.o $(OBJECTS) -lfftw3 -lm -lgoops -lcat

mhdc3dl_python.o: mhdc3dl_python.C
	$(CXX) $(IFLAGS) -pthread -fPIC -c mhdc3dl_python.C 

mhdc3dl_test: $(OBJECTS) mhdc3dl_test.o
	$(CXX) $(FLAGS) -o mhdc3dl_test $(OBJECTS) mhdc3dl_test.o -lfftw3

libvzdeigen:
	$(MAKE) -C vzdeigen

%.o: %.C *.h
	$(CXX) $(IFLAGS) -c $<

install:
	install -c -m 755 mhdc3dl $(INSTALL_ROOT)/bin

uninstall:
	@(cd $(INSTALL_ROOT)/bin && rm -rfv bin/mhdc3dl)



