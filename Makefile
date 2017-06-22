CC = g++
CFLAGS = -O3 -funroll-loops -c
LFLAGS = -O3 -funroll-loops
OBJS = LinearOpticalTransform.o BFGS_Optimization.o AncillaAugment.o MeritFunction.o main.o

all: LinearOpticalSimulation

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
	rm *.o LinearOpticalSimulation *.dat
