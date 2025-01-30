#
# Exact genetic sequence alignment
#
# Parallel computing (Degree in Computer Engineering)
# 2023/2024
#
# (c) 2024 Arturo Gonzalez-Escribano
# Grupo Trasgo, Universidad de Valladolid (Spain)
#

# Compilers
CC=gcc
OMPFLAG=-fopenmp
PFLAG=-pthread -D NUM_THREADS=$(num_threads)
MPICC=mpicc
CUDACC=nvcc

# Flags for optimization and external libs
LIBS=-lm
FLAGS=-O3 -Wall
CUDAFLAGS=-O3 -Xcompiler -Wall
GPROF=-pg

# Targets to build
OBJS=align_seq align_mpi align_pthread

# Rules. By default show help
help:
	@echo
	@echo "Exact genetic sequence alignment"
	@echo
	@echo "Group Trasgo, Universidad de Valladolid (Spain)"
	@echo
	@echo "make align_seq								Build only the sequential version"
	@echo "make align_omp num_threads = <num_threads>	Build only the PThreads version"
	@echo "make align_mpi								Build only the MPI version"
	@echo "make align_cuda								Build only the CUDA version"
	@echo
	@echo "make all	Build all versions (Sequential, OpenMP, MPI, CUDA)"
	@echo "make debug	Build all version with demo output for small sequences"
	@echo "make clean	Remove targets"
	@echo

all: $(OBJS)

align_seq: align.c rng.c
	$(CC) $(FLAGS) $(DEBUG) $< $(LIBS) $(GPROF) -o $@

align_mpi: align_mpi.c rng.c
	$(MPICC) $(FLAGS) $(DEBUG) $< $(LIBS) $(GPROF) -o $@

align_mpi_exp: align_mpi_exp.c rng.c
	$(MPICC) $(FLAGS) $(DEBUG) $< $(LIBS) $(GPROF) -o $@

align_mpi_bcast: align_mpi_bcast.c rng.c
	$(MPICC) $(FLAGS) $(DEBUG) $< $(LIBS) $(GPROF) -o $@

align_pthread: align_pthreads.c rng.c
	$(CC) $(FLAGS) $(DEBUG) $(PFLAG) $< $(LIBS) $(GPROF) -o $@

align_omp: align_omp.c rng.c
	$(CC) $(FLAGS) $(DEBUG) $(OMPFLAG) $< $(LIBS) $(GPROF) -o $@

align_mpipthread: align_mpipthreads.c rng.c
	$(MPICC) $(FLAGS) $(DEBUG) $(PFLAG) $< $(LIBS) $(GPROF) -o $@

align_mpipthread2: align_mpipthreads.c rng.c 
	$(MPICC) $(FLAGS) $(DEBUG) $(PFLAG) $< $(LIBS) $(GPROF) -o $@ 

align_mpipthread4: align_mpipthreads.c rng.c
	$(MPICC) $(FLAGS) $(DEBUG) $(PFLAG) $< $(LIBS) $(GPROF) -o $@

align_mpipthread8: align_mpipthreads.c rng.c
	$(MPICC) $(FLAGS) $(DEBUG) $(PFLAG) $< $(LIBS) $(GPROF) -o $@

align_cuda: align_cuda.cu rng.c
	$(CUDACC) $(CUDAFLAGS) $(DEBUG) $< $(LIBS) $(GPROF) -o $@

# Remove the target files
clean:
	rm -rf $(OBJS)

# Compile in debug mode
debug:
	make DEBUG="-DDEBUG -g" all
