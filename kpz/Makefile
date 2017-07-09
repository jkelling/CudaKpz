
BUILD_COUNTER = ${KMCSVNDIR}/tools/buildCounter
EXTERN_OBJECTS = ${KMCSVNDIR}/inc/libKMC.a ${KMCSVNDIR}/inc/dSFMT/libdSFMT.a ${KMCSVNDIR}/inc/libKMCsplash.a ${KMCSVNDIR}/inc/libKMCCUDA.a \
			${KMCSVNDIR}/inc/tinyMT/libTinyMT.a
CPP_SOURCE = schedulerCPUCell.cpp schedulerCPU.cpp systemSize.cpp threadLayout.cpp schedulerCL.cpp \
			 kpzScaling.cpp schedulerCPU1D.cpp roughness.cpp testKpz.cpp \
			schedulerCPUDB_instances.cpp schedulerCPUpThread.cpp CPURngIface.cpp schedulerCPUDT_instances.cpp \
			schedulerCPUDT_instances.cpp kpzFrontend.cpp schedulerService.cpp schedulerServiceStash.cpp schedulerCPUBase.cpp \
			correlator.cpp correlatorLines.cpp sampling.cpp \
			schedulerCPUBaseThreadedRng.cpp schedulerCPUDTr_instances.cpp schedulerCPUBaseThreaded.cpp \
			schedulerCPUCBBit_instances.cpp
CU_SOURCE = schedulerBase.cu scheduler.cu kpz.cu schedulerGPU.cu GPURandom.cu kpzSCA.cu schedulerSCA.cu GPUTinyMT.cu schedulerDT.cu kpzFrontendCuda.cu \
			kpzSCABit.cu schedulerSCABit.cu correlatorCuda.cu
CUDA_OBJECTS = schedulerBase.cu_o scheduler.cu_o kpz.cu_o schedulerGPU.cu_o GPURandom.cu_o kpzSCA.cu_o schedulerSCA.cu_o GPUTinyMT.cu_o \
			   schedulerDT.cu_o kpzFrontendCuda.cu_o schedulerSCABit.cu_o kpzSCABit.cu_o correlatorCuda.cu_o
OCL_OBJECTS = schedulerCL.o
CPU_SCHEDULER_OBJECTS = schedulerCPUCell.o schedulerCPU.o schedulerCPU1D.o schedulerCPUDB_instances.o schedulerCPUpThread.o CPURngIface.o \
			schedulerCPUDT_instances.o schedulerCPUCB_instances.o schedulerServiceStash.o schedulerCPUBase.o \
			schedulerCPUDTr_instances.o schedulerCPUBaseThreaded.o schedulerCPUCBBit_instances.o
OBJECTS = systemSize.o threadLayout.o roughness.o kpzFrontend.o schedulerService.o \
			correlator.o correlatorLines.o sampling.o
BASEDEPS = ${OBJECTS} ${CPU_SCHEDULER_OBJECTS} ${EXTERN_OBJECTS} 
DEPS = ${BASEDEPS}
OCLFLAGS = -ldl -lOpenCL
OPTFLAGS = -O3 -g -msse -msse2 -msse3 -msse4 -msse4.1 -msse4.2 #-mavx -mavx2

# ifndef CUDA_LIBRARY_PATH
# CUDA_LIBRARY_PATH = /opt/pkg/devel/cuda/5.0/lib64/
# endif
LINKFLAGS = -lhdf5-shared -lsplash

ifdef OPENCL
CPU_SCHEDULER_OBJECTS += schedulerCL.o
LINKFLAGS = -DMODE_OPENCL ${OCLFLAGS}
DEPS += ${OCL_OBJECTS} ${KMCSVNDIR}/inc/opencl/util.o
else
DEPS += ${CUDA_OBJECTS}
endif

NXX = nvcc
NXXFLAGS = -O3 -g --ptxas-options=-v -maxrregcount=32 -std=c++11 -DMODE_CUDA # #-O3 #-deviceemu #-G0 #--compiler-bindir ~/bin/
CXX = g++
CXXFLAGS = ${LINKFLAGS} ${OPTFLAGS} -lpthread -std=c++1y
MXX = mpic++

ifdef FERMI
NXXFLAGS += -arch=compute_20 -code=sm_20 -code=sm_35 -code=sm_37 -code=sm_50 -DKPZ_FERMI
CXXFLAGS += -DKPZ_FERMI
else ifdef KEPLER
NXXFLAGS += -arch=compute_30 -code=sm_35 -code=sm_37 -code=sm_50 -DKPZ_FERMI
CXXFLAGS += -DKPZ_FERMI
else
ifdef K80
	# do not use, using >48k smem does not work
NXXFLAGS += -arch=compute_37 -code=sm_37 -DKPZ_K80
CXXFLAGS += -DKPZ_K80
else
NXXFLAGS += -arch=compute_12 -code=sm_13 -code=sm_12
endif
endif

ifdef PRAND_PATH
CPPFLAGS = -DKPZ_SWITCH_PRAND -DKPZ_SWITCH_RNG=PrandMT
NXXFLAGS += "-I${PRAND_PATH}/include/" ${CPPFLAGS}
LINKFLAGS := -L${PRAND_PATH}/lib/ -lprand ${LINKFLAGS}
CXXFLAGS += ${CPPFLAGS}
CU_SOURCE += Prand.cu PrandDevice.cu
CUDA_OBJECTS += Prand.cu_o
endif

all: gpu cpu

kpzScaling: CXXFLAGS += -DMODE_CUDA
kpzScaling: LINKLAGS += -lcudart -lcuda -L ${CUDA_LIBRARY_PATH}
gpu: kpzScaling 

kpzScalingCpu: CXXFLAGS += -DMODE_CPU ${OPTFLAGS} -lpthread
cpu: kpzScalingCpu

cum: cum.cpp
	g++ -o cum cum.cpp ${KMCSVNDIR}/inc/libKMC.a -O2

testKpz: testKpz.cpp
	g++ -o testKpz testKpz.cpp ${KMCSVNDIR}/inc/dSFMT/libdSFMT.a -O2

# kpzBenchmark : ${DEPS}
kpzScaling: kpzScaling.cpp ${DEPS}
	${MXX} -o $@ $< ${DEPS} ${CXXFLAGS} ${EXTERN_OBJECTS} ${LINKLAGS}
kpzScalingCpu : kpzScaling.cpp ${BASEDEPS}
	${MXX} -o $@ $^ ${CXXFLAGS}
kpzConsistencyCheck: kpzConsistencyCheck.cpp schedulerService.o correlator.o systemSize.o roughness.o correlatorLines.o
	${MXX} -o $@ $^ ${CXXFLAGS} ${EXTERN_OBJECTS} ${LINKLAGS}
%.h : %.cpp ;

kpzConst.h : param.h
scheduler.h : schedulerService.h schedulerService.cpp
%.cu_o : %.cu kpzConst.h
	nvcc ${NXXFLAGS} -c $< -o $@
%.o : %.h kpzConst.h
kpz.cu : kpzCopy.cu kpzRandom.cu kpzSimulationPrimitivesTemplates.cu

doc:
	doxygen Doxyfile

counter:
	@${BUILD_COUNTER} update ${BUILDCOUNTERFILE}
clean:
	rm -f ${OBJECTS} ${CUDA_OBJECTS} ${CPU_SCHEDULER_OBJECTS} ${OCL_OBJECTS} .deps .cudeps
	@echo 'Clean done.'

.deps: ${CPP_SOURCE}
	g++ -E -M -std=c++0x ${CPP_SOURCE} > $@
ifneq ($(which nvcc),)
.cudeps:
	@echo -n '' > $@
	$(foreach f, ${CU_SOURCE}, nvcc -E -M $(f) >> $@; )
else
.cudeps:
	@echo -n '' > $@
endif

include .deps
include .cudeps

