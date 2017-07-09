CUDA_OBJECTS = GPUSLCG64.o GPUTinyMT.o GPURandom123_Threefry.o

ifndef CUDA_LIBRARY_PATH
CUDA_LIBRARY_PATH = ${CUDA_ROOT}/lib64
endif
LINKFLAGS = -lcudart -lcuda -L${CUDA_LIBRARY_PATH}

NXX = nvcc
NXXFLAGS = --ptxas-options=-v --std=c++11
# CXX = g++
# CXXFLAGS = ${LINKFLAGS} -O3

all : libKpzCudaRng.a

%.o : %.cu %.h
	nvcc ${NXXFLAGS} -c $< -o $@
# %.o : %.dat
# 	ld -r -b binary -o $@ $<
# 	objcopy --rename-section .data=.rodata,alloc,load,readonly,data,contents $@ $@

libKpzCudaRng.a : ${CUDA_OBJECTS}
	ar cr $@ $^

clean :
	rm ${CUDA_OBJECTS} libKpzCudaRng.a
