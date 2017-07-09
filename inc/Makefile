# Start of makefile
# Defining variables
SOURCE = arghelper.cpp filehelper.cpp zlibWrapper.cpp combineTables.cpp timer.cpp kmcThreadPool.cpp
LIBKMC_OBJECTS = arghelper.o filehelper.o combineTables.o timer.o kmcThreadPool.o
LIBKMC_CUDA_OBJECTS = kmcExceptCUDA.o
ARGHELPER_OBJECTS = arghelper/arghelper.o #arghelper/glm.o
LIBZWRAPPER_OBJECTS = zlibWrapper.o
CXX = g++
CXXFLAGS = --std=c++0x
LIBKMCSPLASH_OBJECTS = KMCsplash.o
LIBKMCH5_OBJECTS = KMCh5.o

all: libKMC.a libzWrapper.a libdSFMT.a libTinyMT.a
allall: all libKMCCUDA.a libKMCsplash.a libKMCall.a

.PHONY: arghelper/arghelper.o
arghelper/arghelper.o:
	${MAKE} -C arghelper
.PHONY: libdSFMT.a
libdSFMT.a:
	${MAKE} -C dSFMT
libTinyMT.a:
	${MAKE} -C tinyMT
.PHONY: libKpzCudaRng.a
libKpzCudaRng.a:
	${MAKE} -C kpzCudaRng

# KMCh5.o:KMCh5.cpp
# 	${CXX} -std=c++1y -I/home/kelling/checkout/hdf5-1.10.0/install/include/ -c $< -o $@

lib%.a : %.o

%.o : CXXFLAGS += -DDSFMT_MEXP=19937
%.o : %.cu
	nvcc -std=c++11 -c $< -o $@

libKMC.a : ${LIBKMC_OBJECTS} ${ARGHELPER_OBJECTS}
	ar cr $@ $^

libKMCCUDA.a : ${LIBKMC_CUDA_OBJECTS}
	ar cr $@ $^

libzWrapper.a : zlibWrapper.o
	ar cr $@ $^

libKMCsplash.a : ${LIBKMCSPLASH_OBJECTS}
	ar cr $@ $^

libKMCall.a : ${LIBKMC_OBJECTS} ${ARGHELPER_OBJECTS} ${LIBKMC_CUDA_OBJECTS} ${LIBZWRAPPER_OBJECTS} ${LIBKMCH5_OBJECTS} \
	dSFMT/dSFMT.o tinyMT/tinymt32.o kpzCudaRng/GPUTinyMT.o kpzCudaRng/GPUSLCG64.o ${LIBZWRAPPER_OBJECTS} 
	# kpzCudaRng/GPURandom123_Threefry.o
	ar cr $@ $^

clean:
	rm -f ${LIBKMC_OBJECTS} ${LIBZWRAPPER_OBJECTS} ${LIBKMCSPLASH_OBJECTS} ${LIBKMCH5_OBJECTS}
	@echo 'Clean done.'
	${MAKE} -C arghelper clean
	${MAKE} -C dSFMT clean

