# Compiler
CXX = g++

# Optimization and Debug flags
OPT = -O2 #-g -v

# Other compiler flags
FCFLAGS = -fPIC
CXXFLAGS = $(OPT) -fPIC -Wall -Wextra -g

# ROOT flags and libs
ROOT_FLAGS_INC = $(shell root-config --cflags --glibs)
ROOT_LIB = $(shell root-config --cflags --glibs)

# All Includes
INCS = $(ROOT_FLAGS_INC)

# All libs 
LIBS = -L. $(ROOT_LIB)


all: compile_storer run_storer 
#all: compile_storer run_storer compile_fitter run_fitter
#all: Storer
#
#Storer: Storer.o amp_func.o
#	$(CXX) -o $@ Storer.o amp_func.o $(LIBS)
#
#Storer.o: Storer.cpp
#	$(CXX) -c Storer.cpp $(CXXFLAGS) $(INCS)
#
#amp_func.o: amp_func.cpp
#	$(CXX) $(CXXFLAGS) $(INCS) -c amp_func.cpp -o amp_func.o


#compile_storer: Storer.cpp 
#	$(CXX) Storer.cpp `root-config --cflags --glibs` -o Storer
#amp_func.o: amp_func.cpp
#	$(CXX) $(CXXFLAGS) $(INCS) -c amp_func.cpp -o amp_func.o

compile_storer: Storer.o amp_func.o
	$(CXX) Storer.o amp_func.o `root-config --cflags --glibs` -o Storer
	#$(CXX) -o Storer Storer.o amp_func.o $(LIBS)
	#$(CXX) -o $@ Storer.o amp_func.o $(LIBS)
Storer.o: Storer.cpp
	$(CXX) -c -w Storer.cpp `root-config --cflags --glibs` -o Storer.o
	#$(CXX) -c -w Storer.cpp $(CXXFLAGS) $(INCS)
amp_func.o: amp_func.cpp
	$(CXX) -c -w amp_func.cpp `root-config --cflags --glibs` -o amp_func.o
	#$(CXX) $(CXXFLAGS) $(INCS) -c -w amp_func.cpp -o amp_func.o


run_storer: Storer
	./Storer

#compile_fitter: Fitter.cpp
#	$(CXX) Fitter.cpp `root-config --cflags --glibs` -o Fitter

#run_fitter: Fitter
#	./Fitter

clean:
	rm -rf Fitter Storer *.o *.so
