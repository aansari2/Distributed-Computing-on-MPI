#MAE 598 - HPC - Final Project - Adil Ansari
MAX_ITERATION = 4000
THRESHOLD = 1e-4
CC1 = mpicc
CC2 = ifort
CR = mpirun
CFLAGS1  = -lm
CFLAGS2 = -coarray=shared
CFRUN = -np
CORES = 16

TARGET0 = jacobi

TARGET1 = jacobimpi

TARGET2 = jacobicoarray

all: setcore clean $(TARGET0) $(TARGET1) $(TARGET2) run

setcore:
	export FOR_COARRAY_NUM_IMAGES=$(CORES); 

$(TARGET0): $(TARGET0).c
	icc -o $(TARGET0) $(TARGET0).c

$(TARGET1): $(TARGET1).c
	$(CC1) $(CFLAGS1) -o $(TARGET1) $(TARGET1).c
	
$(TARGET2): $(TARGET2).f90
	$(CC2) $(CFLAGS2) -o $(TARGET2) $(TARGET2).f90

clean:
	$(RM) $(TARGET1); $(RM) $(TARGET2);

run: $(TARGET0) $(TARGET1) $(TARGET2)
	./$(TARGET0) $(MAX_ITERATION) $(THRESHOLD); #Run Single Core Implement
	$(CR) $(CFRUN) $(CORES) ./$(TARGET1) $(MAX_ITERATION) $(THRESHOLD); #Run MPI Program
	./$(TARGET2) $(MAX_ITERATION) $(THRESHOLD); #Run Coarray Program
