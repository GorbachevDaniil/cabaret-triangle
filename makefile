#dir:
ROOT_DIR = ./
PROJ_SRC = ${ROOT_DIR}src
INCLUDE_DIR = ${ROOT_DIR}include

#the c++ compiler to use:
CXX = g++

#compilation flags for c++ source files:
OPTS=-O3 -g0
CXXFLAGS = $(OPTS) -Wall -iquote $(INCLUDE_DIR)

#source and header files
vpath %.h $(INCLUDE_DIR)
vpath %.cpp $(PROJ_SRC)

objects = Main.o 

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@

all: cabaret_triangle

cabaret_triangle: $(objects)
	$(CXX) $(objects) -o cabaret_triangle 

clean:
	rm -rf *.o cabaret_triangle






