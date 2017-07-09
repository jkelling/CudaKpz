all: acNoise sk2 sk

acNoise: acNoise.cpp
	${CXX} --std=c++0x -g -O1 -o $@ $< ../../inc/libKMC.a
sk: sk.cpp
	${CXX} --std=c++0x -g -O1 -o $@ $< ../../inc/libKMC.a
agingScaling: agingScaling.cpp
	${CXX} --std=c++1y -g -O1 -o $@ $< ../../inc/libKMC.a

sk2: sk2.c
	${CC} -g -O1 -lm -o $@ $<
