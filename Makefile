CC = g++ -Wall -g

poissonSolver: PoissonMatrix.o Algorithms.o PoissonVector.o MGMVector.o UpperMatrix.o LowerMatrix.o CGOperators.o PoissonSolver.o
	$(CC) -o $@ $+

PoissonMatrix.o: PoissonMatrix.cpp classes.h
	$(CC) -c -o $@ $<

Algorithms.o: Algorithms.cpp classes.h
	$(CC) -c -o $@ $<

PoissonVector.o: PoissonVector.cpp classes.h
	$(CC) -c -o $@ $<

MGMVector.o: MGMVector.cpp classes.h
	$(CC) -c -o $@ $<

UpperMatrix.o: UpperMatrix.cpp classes.h
	$(CC) -c -o $@ $<

LowerMatrix.o: LowerMatrix.cpp classes.h
	$(CC) -c -o $@ $<

CGOperators.o: CGOperators.cpp classes.h
	$(CC) -c -o $@ $<

PoissonSolver.o: PoissonSolver.cpp classes.h
	$(CC) -c -o $@ $<

clean:
	rm *.o poissonSolver