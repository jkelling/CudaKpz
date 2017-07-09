OPTFLAGS = -O2 -g
CXXFLAGS = -DMODE_CPU ${OPTFLAGS} -lpthread

EXTERN_OBJECTS = ${KMCSVNDIR}/inc/libKMC.a ${KMCSVNDIR}/inc/dSFMT/libdSFMT.a

schedulerRandom: schedulerRandom.cpp ../systemSize.o
	${CXX} -o $@ $^ ${CXXFLAGS} ${EXTERN_OBJECTS}
