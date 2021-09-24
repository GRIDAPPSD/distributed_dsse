

EIGENPATH=$(HOME)/git/eigen
INCLUDES=-I$(EIGENPATH)

# For SuiteSparse code only
SSPATH=$(HOME)/git/gridappsd-state-estimator/SuiteSparse

SSINCPATH=$(SSPATH)/include
# APR is part of the Linux distro by default, hopefully here
APRINCPATH=/usr/include/apr-1.0
SSINCLUDES=-I$(SSINCPATH) -I$(APRINCPATH)

# SparseSuite libraries
SSCONFLIBPATH=$(SSPATH)/SuiteSparse_config
CSLIBPATH=$(SSPATH)/CXSparse/Lib
AMDLIBPATH=$(SSPATH)/AMD/Lib
BTFLIBPATH=$(SSPATH)/BTF/Lib
COLAMDLIBPATH=$(SSPATH)/COLAMD/Lib
KLULIBPATH=$(SSPATH)/KLU/Lib
# use the first definition below for static libs, the second for shared libs
SSLIBPATHS=-L$(SSCONFLIBPATH) -L$(CSLIBPATH) -L$(AMDLIBPATH) -L$(BTFLIBPATH) -L$(COLAMDLIBPATH) -L$(KLULIBPATH)
#SSLIBPATHS=-L$(SSPATH)/lib
SSLIBS=-lklu -lamd -lbtf -lcolamd -lcxsparse -lsuitesparseconfig
# End SuiteSparse defs


# Changed from -Ofast to -O2 on 7/30/20 because under an Ubuntu 20 VM with
# g++ version 9.3.0 -Ofast was giving bad numerical results in terms of
# intermittent very high Supd condition numbers and incorrect state estimates
# Then on 9/11/20, when adding the Rmat code for uncertainty based on time
# since last measurement, it happened again where it would give KLU errors
# right away every few times SE was started.  Just some very innocuous code
# adding zlast to zary and not even using it caused the issue.  So now -O
# seems to be the one that produces reliable results.  Benchmarking on the 9500
# node model shows no discernable loss in performance fortunately.
CXXFLAGS=-O
#CXXFLAGS=-O2
#CXXFLAGS=-Ofast

all: cpp_sparse_benchmark cpp_suitesparse_benchmark cpp_dense_benchmark

cpp_suitesparse_benchmark: cpp_suitesparse_benchmark.o
	g++ $(CXXFLAGS) $? $(SSLIBPATHS) $(SSLIBS) -lstdc++ -o $@

cpp_suitesparse_benchmark.o: cpp_suitesparse_benchmark.cpp
	g++ $(CXXFLAGS) -Wno-deprecated-declarations -std=c++11 -c $(SSINCLUDES) $< -o $@

cpp_sparse_benchmark: cpp_sparse_benchmark.o
	g++ $(CXXFLAGS) $? -lstdc++ -o $@

cpp_sparse_benchmark.o: cpp_sparse_benchmark.cpp
	g++ $(CXXFLAGS) -Wno-deprecated-declarations -std=c++11 -c $(INCLUDES) $< -o $@

cpp_dense_benchmark: cpp_dense_benchmark.o
	g++ $(CXXFLAGS) $? -lstdc++ -o $@

cpp_dense_benchmark.o: cpp_dense_benchmark.cpp
	g++ $(CXXFLAGS) -Wno-deprecated-declarations -std=c++11 -c $(INCLUDES) $< -o $@

clean:
	rm -f cpp_sparse_benchmark cpp_sparse_benchmark.o cpp_suitesparse_benchmark cpp_suitesparse_benchmark.o cpp_dense_benchmark cpp_dense_benchmark.o
