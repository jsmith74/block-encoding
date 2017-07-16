CC = g++
CFLAGS = -Ofast -funroll-loops -c
LFLAGS = -Ofast -funroll-loops
OBJS = LinearOpticalTransform.o BFGS_Optimization.o AncillaAugment.o MeritFunction.o main.o
OMPFLAGS = -fopenmp

all: LinearOpticalSimulation DataProcess

DataProcess: dataProcess.cpp
	$(CC) dataProcess.cpp -o DataProcess

LinearOpticalSimulation: $(OBJS)
	$(CC) $(OMPFLAGS) $(LFLAGS) $(OBJS) -o LinearOpticalSimulation

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

LinearOpticalTransform.o: LinearOpticalTransform.cpp
	$(CC) $(OMPFLAGS) $(CFLAGS) LinearOpticalTransform.cpp

BFGS_Optimization.o: BFGS_Optimization.cpp
	$(CC) $(CFLAGS) BFGS_Optimization.cpp

AncillaAugment.o: AncillaAugment.cpp
	$(CC) $(CFLAGS) AncillaAugment.cpp

MeritFunction.o: MeritFunction.cpp
	$(CC) $(CFLAGS) MeritFunction.cpp

clean:
	rm *.o LinearOpticalSimulation *.dat DataProcess
