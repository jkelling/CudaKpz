OBJECTS = dSFMT.o

CXX = g++
CXXFLAGS = ${LINKFLAGS} -O3

all : libdSFMT.a

lib%.a : %.o
	ar cr $@ $<

dSFMT.o : dSFMT.cpp
	${CXX} ${CXXFLAGS} -c -o $@ $< -DDSFMT_MEXP=19937

clean :
	rm ${OBJECTS}
