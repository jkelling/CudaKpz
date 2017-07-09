CXXFLAGS = -std=c++0x -O3
CXXFLAGS += -DRLOG
#CXXFLAGS += -DSUMSIM
CXXFLAGS += -DDIMENSION=3
#CXXFLAGS += -DKMC_SWITCH_BENCHMARK
CXX = g++
EXTERNOBJECTS = ${KMCSVNDIR}/inc/dSFMT/libdSFMT.a

FILES = lattice.h lattice.cpp sumSim.h sumSim.cpp heightMap.cpp slopeSim.h slopeSim.cpp dim.h

scaling: scaling.cpp ${FILES}
	${CXX} -o $@ $< ${CXXFLAGS} ${EXTERNOBJECTS}

clean:
	rm scaling
