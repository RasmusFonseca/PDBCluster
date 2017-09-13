
# Primary flags for the compiler and linker
CPPFLAGS = -std=c++0x -Wall -g -O2
LDFLAGS = 
CXX = g++

# Source files
SOURCES = $(wildcard *.cpp)

# object file
OBJECTS = $(SOURCES:.cpp=.o)

PDBCluster: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $^ -o $@

# clean executable
clean:
	rm -f *.o PDBCluster

