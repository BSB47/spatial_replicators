VER = ver4

CC = gcc
# CXX = /home/takeuchi/local/gstlfilt/gfilt /cand:L
CXX = g++

# Name of your program
PROJECT = CenSur

# Your files
OBJMAIN = cash-display.o main.o newca.o

# You don't have to change here
OBJCASH = arithmetic.o basic.o color.o filter.o io.o logical.o \
margolus.o movie.o neighbors.o noise.o png.o ps.o random.o shift.o\
x11.o

# Options to compiler (both for CC and CXX)
#CCOPT = -O0 -ggdb -Wall -DTBB_USE_DEBUG# For debagging.
#CCOPT = -O3 -Wall -Winline  --param inline-unit-growth=150 # For weak debagging
COPT = -O3 -Wall -DNDEBUG # Most optimized.
CCOPT = -g -std=c++17 -O3 -Wall -DNDEBUG # Most optimized.
#CCOPT = -O3 -pg -Wall -Winline --param large-function-growth=450 --param inline-unit-growth=150 -DNDEBUG # For profiling

#CCOPT = -O3 -Wall -Winline -finline-limit=20000 --param inline-unit-growth=5000 --param large-function-growth=60000 -DNDEBUG # Most optimized.
#CCOPT = -O3 -pg -Wall -Winline -finline-limit=20000 --param inline-unit-growth=5000 --param large-function-growth=60000 -DNDEBUG # for profile

# Path to libraries and include files
LIBS = -lpng -lz -lX11 -lm

# Other version of TBB
#LDFLAGS = -Wl,-rpath,/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_20090809oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_release/ -Wl,-rpath,/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_20090809oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_debug/
#IDIR = -I. -I/home/takeuchi/local/include/ -I/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_20090809oss/include/
#LDIR = -L/usr/X11R6/lib64/ -L/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_20090809oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_release/ -L/home/takeuchi/study/dna-simple/program/$(VER)/tbb22_20090809oss/build/linux_intel64_gcc_cc4.1.2_libc2.3.5_kernel2.6.5_debug/

# -fomit-frame-pointer -finline-functions 

# You don't have to change here
all: $(OBJCASH) $(OBJMAIN) #source.tar.gz
	$(CXX) $(OBJCASH) $(OBJMAIN) $(CCOPT) -o $(PROJECT) $(LIBS) 

$(OBJCASH): Makefile cash.h

# Dependency of files (add/modify if necessarly)

$(OBJMAIN): Makefile
main.o: molecule.h newca.h para.h random.h
newca.o: newca.h cash-display.hpp cash.h molecule.h random.h
cash-display.o: cash-display.hpp cash.h assert.hpp
random.o: random.h

#CASH sources
C_SOURCE = $(addsuffix .c, $(basename $(OBJCASH)))
#CXX sources which also have a header
CXX_BASE = $(basename para.o newca.o cash-display.o)
#header only files, source only files, makefile and RandomLib
OTHERS = Makefile main.cpp cash.h assert.hpp \
cellular-automata.hpp

#everything needed to compile a program
EVERYTHING = $(C_SOURCE) $(addsuffix .cpp, $(CXX_BASE)) \
$(addsuffix .hpp, $(CXX_BASE)) $(OTHERS)

# source.tar.gz: $(EVERYTHING)
# 	tar -zcf source.tar.gz $(EVERYTHING)

# Make rules

%.o : %.c
	$(CC) -c $(COPT) $< -o $@

%.o : %.cpp
	$(CXX) -c $(CCOPT) $< -o $@

clean:
	rm *.o
