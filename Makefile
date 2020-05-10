#Makefile for Short Signature Scheme
#Ben Lynn
#
#Copyright (C) 2001 Benjamin Lynn (blynn@cs.stanford.edu)
#
#This file is part of the Stanford short signature system.
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

SSL_INCLUDE=/usr/local/ssl/include
SSL_LIB=/usr/local/ssl/lib
NTL_INCLUDE=/usr/local/include
LNTL=-lntl -lgmp
LSSL=-lssl -lcrypto
CXX=g++
CPPFLAGS= -O2 -mcpu=i686 -Wall -I$(NTL_INCLUDE) -I$(SSL_INCLUDE) -L$(SSL_LIB)

.PHONY: target dist clean

F3KLOS=f3kl.o f3k.o

WEILOS=weil.o $(F3KLOS)

BINARIES=benchmark f3kmaketable f3ktest f3kltest genminpoly geninitcode

target: $(BINARIES)

weil.o : weil.cc weil.h

benchmark.o : benchmark.cc

sss_lib.o : sss_lib.cc sss.h

test.o : test.cc

test : test.o sss_lib.o $(WEILOS)
	$(CXX) $(CPPFLAGS) -o test test.o sss_lib.o $(WEILOS) $(LSSL) $(LNTL)

geninitcode.o : geninitcode.cc

gen.o : gen.cc

genminpoly.o : genminpoly.cc

genminpoly: genminpoly.o $(WEILOS)
	$(CXX) $(CPPFLAGS) -o genminpoly genminpoly.o $(WEILOS) $(LNTL)

geninitcode: geninitcode.o $(WEILOS)
	$(CXX) $(CPPFLAGS) -o geninitcode geninitcode.o $(WEILOS) $(LNTL)

f3k.o : f3k.cc f3k.h

f3ktest.o : f3ktest.cc

f3ktest: f3ktest.o f3k.o
	$(CXX) $(CPPFLAGS) -o f3ktest f3ktest.o f3k.o $(LNTL)

f3kl.o : f3kl.cc f3kl.h f3k.h
f3kltest.o : f3kltest.cc

f3kltest: f3kltest.o $(F3KLOS)
	$(CXX) $(CPPFLAGS) -o f3kltest f3kltest.o $(F3KLOS) $(LNTL)

f3kmaketable.o : f3kmaketable.cc

f3kmaketable: f3kmaketable.o
	$(CXX) $(CPPFLAGS) -o f3kmaketable f3kmaketable.o $(LNTL)

benchmark: benchmark.o $(WEILOS)
	$(CXX) $(CPPFLAGS) -o benchmark benchmark.o $(WEILOS) $(LNTL)

ALLFILES=*.cc *.h Makefile HISTORY LICENSE

projname := $(shell awk '/SSS_VERSION/ { print $$3 }' version.h )

results: benchmark
	./benchmark 0 > results.0
	./benchmark 1 > results.1
	./benchmark 2 > results.2
	./benchmark 3 > results.3
	./benchmark 4 > results.4
	./benchmark 5 > results.5
	./benchmark 6 > results.6

dist: $(ALLFILES)
	-rm -rf $(projname)
	mkdir $(projname)
	cp -rl --parents $(ALLFILES) $(projname)
	tar chfz $(projname).tgz $(projname)
	-rm -rf $(projname)

clean:
	-rm *.o $(BINARIES) results.?
