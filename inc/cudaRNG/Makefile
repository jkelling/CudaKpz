CUDA_OBJECTS = CudaRNG.o

ifndef CUDA_LIBRARY_PATH
CUDA_LIBRARY_PATH = /opt/pkg/devel/cuda/3.0/cuda/lib64
endif
LINKFLAGS = -lcudart -L ${CUDA_LIBRARY_PATH}

NXX = nvcc
NXXFLAGS = --ptxas-options=-v --host-compilation c++ # #-O3 #-deviceemu #-G0 #--compiler-bindir ~/bin/ 
CXX = g++
CXXFLAGS = ${LINKFLAGS} -O3

all : libCudaRNG.a

%.o : %.cu %.h
	nvcc ${NXXFLAGS} -c $< -o $@
%.o : %.dat
	ld -r -b binary -o $@ $<
	objcopy --rename-section .data=.rodata,alloc,load,readonly,data,contents $@ $@

libCudaRNG.a : CudaRNG.o MersenneTwister.o
	ar cr $@ $^
CudaRNG.o : MersenneTwister_NVGPUSDK.h MersenneTwister_kernel_NVGPUSDK.cu RC4.cu twofish/TwofishKernel.cu twofish/TwofishInternal.cu

test : test.cu libCudaRNG.a
	nvcc ${NXXFLAGS} -o $@ $< libCudaRNG.a ../fft.o ../dSFMT/libdSFMT.a

clean :
	rm ${CUDA_OBJECTS} MersenneTwister.o test
