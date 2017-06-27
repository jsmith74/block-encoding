CC = g++
CFLAGS = -O3 -funroll-loops -c
LFLAGS = -O3 -funroll-loops
OBJS = LinearOpticalTransform.o BFGS_Optimization.o AncillaAugment.o MeritFunction.o main.o
OMPFLAGS = -fopenmp

all: LinearOpticalSimulation Script DataProcess

DataProcess: dataProcess.cpp
	$(CC) dataProcess.cpp -o DataProcess

Script: script.cpp
	$(CC) $(OMPFLAGS) script.cpp -o Script

LinearOpticalSimulation: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o LinearOpticalSimulation

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

LinearOpticalTransform.o: LinearOpticalTransform.cpp
	$(CC) $(CFLAGS) LinearOpticalTransform.cpp

BFGS_Optimization.o: BFGS_Optimization.cpp
	$(CC) $(CFLAGS) BFGS_Optimization.cpp

AncillaAugment.o: AncillaAugment.cpp
	$(CC) $(CFLAGS) AncillaAugment.cpp

MeritFunction.o: MeritFunction.cpp
	$(CC) $(CFLAGS) MeritFunction.cpp

clean:
	rm *.o LinearOpticalSimulation *.dat Script DataProcess
