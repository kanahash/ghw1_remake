exec = ghw1
objs = ghw1.o 
#LFLAGS = -mcmodel=large -march=native -lm -O2 -fopenmp -lpthread -lfftw3_threads -lfftw3 #-pg
####LFLAGS = -mcmodel=large -march=native -lm -O2 -fopenmp -lfftw3_omp -lfftw3
LFLAGS = -mcmodel=large -march=native -lm -O2 -fopenmp -lfftw3_omp -lfftw3 
CFLAGS = -I.
LIBS = /usr/lib64 
all: $(exec)
$(exec): $(objs)
	g++  $(exec).cpp -o $@  $(LFLAGS) -L$(LIBS)
	rm *.o
ghw.o: ghw1.cpp 
	g++  -c $< -o $@
.PHONY : clean
clean :
	rm -f $(exec) $(objs) *~ 
