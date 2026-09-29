# Builds the simulators into bin/.
#   make          CPU simulator (bin/fsm2d; needs a C++17 compiler with OpenMP)
#   make gpu      CUDA simulators (bin/fsm2d_gpu, bin/fsm3d_gpu; needs nvcc). Set CUDA_ARCH for your GPU, e.g.
#                 make gpu CUDA_ARCH=sm_89 (RTX 40xx, L4), sm_90 (H100), sm_120 (RTX 50xx / RTX PRO 6000, default)
#   make test     build everything that can be built and run the test suite (tests/)
CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -fopenmp -Wall -Wextra
NVCC ?= nvcc
CUDA_ARCH ?= sm_120
NVCCFLAGS ?= -O3 -arch=$(CUDA_ARCH) -std=c++17 --use_fast_math

COMMON := $(wildcard src/common/*.hpp)
FSM2D := $(wildcard src/fsm2d/*.hpp)

all: bin/fsm2d

gpu: bin/fsm2d_gpu bin/fsm3d_gpu

bin/fsm2d: src/fsm2d/cpu_main.cpp $(FSM2D) $(COMMON)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

bin/fsm2d_gpu: src/fsm2d/gpu_main.cu $(wildcard src/fsm2d/*.cuh) $(FSM2D) $(COMMON) src/common/cuda_util.cuh
	@mkdir -p bin
	$(NVCC) $(NVCCFLAGS) -o $@ $<

bin/fsm3d_gpu: src/fsm3d/main.cu $(wildcard src/fsm3d/*) $(COMMON) src/common/cuda_util.cuh
	@mkdir -p bin
	$(NVCC) $(NVCCFLAGS) -o $@ $<

test: all
	-$(MAKE) gpu
	python3 -m pytest -q tests

clean:
	rm -rf bin

.PHONY: all gpu test clean
