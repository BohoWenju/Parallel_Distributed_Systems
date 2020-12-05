CC = gcc
CFLAGS = -o
CILKFLAG = -fcilkplus
OPENMPFLAG = -fopenmp
PTHREADFLAG = -pthread
RM = rm -f

all: V3_Serial V3_Cilk V3_OpenMP V3_PThread 

V3_Serial:
	$(CC) V3_Serial.c mmio.c $(CFLAGS) V3_Serial

V3_Cilk:
	$(CC) V3_Cilk.c mmio.c $(CFLAGS) V3_Cilk $(CILKFLAG)

V3_OpenMP:
	$(CC) V3_OpenMP.c mmio.c $(CFLAGS) V3_OpenMP $(OPENMPFLAG)

V3_PThread:
	$(CC) V3_Pthread.c mmio.c $(CFLAGS) V3_PThread $(PTHREADFLAG)


#clean temporary files
clean:
	$(RM) *.o V3_Serial V3_Cilk V3_OpenMP V3_PThread 
