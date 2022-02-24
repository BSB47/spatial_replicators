VER = ver4

CC = gcc
# CXX = /home/takeuchi/local/gstlfilt/gfilt /cand:L
CXX = g++

# Name of your program
PROJECT = dnaevo

# Your files
OBJMAIN = main.o para.o RandomLib/Random.o myca.o vesicles.o\
vesicle-info.o my_dist.o automaton.o cash-display.o vesicle-stat.o\
random_wrapper.o parallel-myca.o serial-myca.o upd_order.o\
atomic-vesicles.o atomic-vesicle-info.o

# You don't have to change here
OBJCASH = arithmetic.o basic.o color.o filter.o io.o logical.o \
margolus.o movie.o neighbors.o noise.o png.o ps.o random.o shift.o\
x11.o

# Options to compiler (both for CC and CXX)
#CCOPT = -O0 -ggdb -Wall -DTBB_USE_DEBUG# For debagging.
#CCOPT = -O3 -Wall -Winline  --param inline-unit-growth=150 # For weak debagging
COPT = -O3 -Wall -DNDEBUG # Most optimized.
CCOPT = -std=c++03 -O3 -Wall -DNDEBUG # Most optimized.
#CCOPT = -O3 -pg -Wall -Winline --param large-function-growth=450 --param inline-unit-growth=150 -DNDEBUG # For profiling

#CCOPT = -O3 -Wall -Winline -finline-limit=20000 --param inline-unit-growth=5000 --param large-function-growth=60000 -DNDEBUG # Most optimized.
#CCOPT = -O3 -pg -Wall -Winline -finline-limit=20000 --param inline-unit-growth=5000 --param large-function-growth=60000 -DNDEBUG # for profile

# Path to libraries and include files
LIBS = -lpng -lz -lX11 -lm -ltbb
LDFLAGS = -Wl,-rpath,/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_013oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_release/ -Wl,-rpath,/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_013oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_debug/
IDIR = -I. -I/home/takeuchi/local/include/ -I/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_013oss/include/
LDIR = -L/usr/X11R6/lib64/ -L/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_013oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_release/ -L/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_013oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_debug/

# Other version of TBB
#LDFLAGS = -Wl,-rpath,/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_20090809oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_release/ -Wl,-rpath,/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_20090809oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_debug/
#IDIR = -I. -I/home/takeuchi/local/include/ -I/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_20090809oss/include/
#LDIR = -L/usr/X11R6/lib64/ -L/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_20090809oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_release/ -L/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_20090809oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_debug/

# Flag used for RnadomLib
SFMTFLAGS = -DUSE_SFMT -DMTPREC64 -DHAVE_SSE2 -msse2\
-DHAVE_BOOST_SERIALIZATION=1 -funroll-loops -I.\
-Winline -finline-limit=5000 --param inline-unit-growth=2000\
--param large-function-growth=2000

# -fomit-frame-pointer -finline-functions 

# You don't have to change here
all: $(OBJCASH) $(OBJMAIN) source.tar.gz
	$(CXX) $(OBJCASH) $(OBJMAIN) $(CCOPT) -o $(PROJECT) $(LDFLAGS) $(LIBS) $(LDIR)

$(OBJCASH): Makefile cash.h

RandomLib/Random.o: Makefile RandomLib/Random.cpp RandomLib/Random.hpp
	$(CXX) -c $(CCOPT) $(SFMTFLAGS) $(IDIR) RandomLib/Random.cpp -o RandomLib/Random.o

# Dependency of files (add/modify if necessarly)

$(OBJMAIN): Makefile
main.o: para.hpp myca.hpp random_wrapper.hpp random_wrapper.hpp
para.o: RandomLib/Random.hpp para.hpp 
myca.o:          cellular-automata.hpp assert.hpp para.hpp automaton.hpp myca.hpp random_wrapper.hpp upd_order.hpp vesicles.hpp atomic-vesicles.hpp cash-display.hpp
parallel-myca.o: cellular-automata.hpp assert.hpp para.hpp automaton.hpp myca.hpp random_wrapper.hpp upd_order.hpp              atomic-vesicles.hpp
serial-myca.o:   cellular-automata.hpp assert.hpp para.hpp automaton.hpp myca.hpp random_wrapper.hpp upd_order.hpp vesicles.hpp
vesicles.o: vesicles.hpp assert.hpp para.hpp my_dist.hpp vesicle-info.hpp
vesicle-info.o: vesicle-info.hpp assert.hpp
my_dist.o: my_dist.hpp para.hpp
automaton.o: automaton.hpp assert.hpp para.hpp
cash-display.o: cash-display.hpp cash.h assert.hpp
vesicle-stat.o: para.hpp assert.hpp vesicle-stat.hpp
random_wrapper.o: RandomLib/Random.hpp para.hpp random_wrapper.hpp
upd_order.o: upd_order.hpp random_wrapper.hpp para.hpp
atomic-vesicles.o: atomic-vesicle-info.hpp atomic-vesicles.hpp  para.hpp
atomic-vesicle-info.o: atomic-vesicle-info.hpp assert.hpp

#CASH sources
C_SOURCE = $(addsuffix .c, $(basename $(OBJCASH)))
#CXX sources which also have a header
CXX_BASE = $(basename para.o myca.o vesicles.o vesicle-info.o vesicle-stat.o my_dist.o automaton.o cash-display.o random_wrapper.o atomic-vesicles.o atomic-vesicle-info.o upd_order.o)
#header only files, source only files, makefile and RandomLib
OTHERS = Makefile main.cpp cash.h assert.hpp \
cellular-automata.hpp RandomLib/ parallel-myca.cpp serial-myca.cpp


#everything needed to compile a program
EVERYTHING = $(C_SOURCE) $(addsuffix .cpp, $(CXX_BASE)) \
$(addsuffix .hpp, $(CXX_BASE)) $(OTHERS)

source.tar.gz: $(EVERYTHING)
	tar -zcf source.tar.gz $(EVERYTHING)

# Make rules

%.o : %.c
	$(CC) -c $(COPT) $(IDIR) $< -o $@

%.o : %.cpp
	$(CXX) -c $(CCOPT) $(IDIR) $< -o $@

clean:
	rm *.o
