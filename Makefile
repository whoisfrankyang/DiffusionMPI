CC = CC
CFLAGS = -lm -g
OPTFLAGS = -O3 -Ofast -funroll-loops -ftree-vectorize -march=native -mtune=native -flto -ffast-math -fassociative-math -fno-trapping-math -freciprocal-math -fprefetch-loop-arrays -fno-strict-aliasing -falign-functions=64 -falign-loops=16 -mavx2 -mfma -msse4.2 -fopt-info-vec -floop-parallelize-all -fno-exceptions -fno-rtti -funroll-all-loops -fipa-pta -fgcse-after-reload -fomit-frame-pointer -fipa-sra -fpredictive-commoning -freorder-blocks-and-partition -fstrict-aliasing -funsafe-math-optimizations -fno-math-errno -fipa-cp-clone -fwhole-program -foptimize-sibling-calls -floop-strip-mine -floop-interchange -floop-block -funswitch-loops
MPIFLAGS = -DMPI

all: diffusion

serial: main.cpp serial.cpp
	$(CC) $^ -o $@ $(CFLAGS) $(OPTFLAGS)

mpi: main.cpp mpi.cpp
	$(CC) $^ -o $@ $(MPIFLAGS) $(CFLAGS) $(OPTFLAGS)

clean:
	rm -f ./serial
	rm -f *.txt
	rm -f ./mpi
	rm -rf snapshots