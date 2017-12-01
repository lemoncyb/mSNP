export LC_CTYPE=C

CXX = /usr/local/mpi/bin/mpicxx
CXXFLAGS = -Wall -openmp -vec-report=0 -I/usr/local/mpi/include/

LFLAGS = -lz

all: msnp
.PHONY: all

objects: call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o main.o
$(objects): %.o: soap_snp.h makefile

msnp: call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o main.o makefile
	$(CXX) $(CXXFLAGS) call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o main.o -o msnp $(LFLAGS)

.PHONY: clean
clean:
	rm -f *.o msnp
