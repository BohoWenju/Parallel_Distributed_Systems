CC = gcc
CFLAGS = -o
CILKFLAG = -fcilkplus
OPENMPFLAG = -fopenmp
PTHREADFLAG = -pthread
RM = rm -f

all: V4_Serial V4_Cilk V4_PThread V4_OpenMP

V4_Cilk: 
	$(CC) V4_Cilk.c mmio.c $(CFLAGS) V4_Cilk $(CILKFLAG)

V4_OpenMP: 
	$(CC) V4_OpenMP.c mmio.c $(CFLAGS) V4_OpenMP $(OPENMPFLAG)

V4_PThread: 
	$(CC) V4_PThread.c mmio.c $(CFLAGS) V4_PThread $(PTHREADFLAG)

#clean temporary files 
clean: 
	$(RM) *.o V4_Serial V4_Cilk V4_PThread
