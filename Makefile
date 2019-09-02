# defines
CC=gcc
MPICC=mpicc
CXX=g++
MPICXX=mpicxx
CFLAGS= -Wall -std=gnu99 -O3
CXXFLAGS= -Wall  -O3
LDFLAGS=  -lrt  -O3


TARGET = gen_valid_info validation gen_RMAT mst_reference mst_reference_mpi mst

all: $(TARGET)

# your own implementation, executable must called mst
mst: main.o mst.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

mst_reference_mpi: main_mpi.mpi.o graph_tools.mpi.o gen_RMAT_mpi.mpi.o mst_reference_mpi.mpi.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)


# reference implementation
mst_reference: main.o mst_reference.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gen_RMAT: gen_RMAT.o
	$(CXX) $^ -o $@ $(LDFLAGS)

validation: validation.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gen_valid_info: graph_tools.o mst_reference.o gen_valid_info.o
	$(CXX) $^ -o $@ $(LDFLAGS)

%.mpi.o: %.cpp
	$(MPICXX) -DUSE_MPI $(CXXFLAGS) -o $@ -c $<

%.mpi.o: %.c
	$(MPICC) -DUSE_MPI $(CFLAGS) -o $@ -c $<

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf *.o $(TARGET)
