# Start of makefile
# Defining variables
SOURCE = glm.cpp
OBJECTS = glm.o
CXX = g++

arghelper.o : ${OBJECTS}
	ld -r -o $@ $^

%.o : CXXFLAGS += -DDSFMT_MEXP=19937

clean:
	rm -f ${OBJECTS}
	@echo 'Clean done.'
