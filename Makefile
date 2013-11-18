CC = g++ -Wall -g

poissonSolver: PoissonMatrix.o Algorithms.o CGVectors.o UpperMatrix.o LowerMatrix.o Operators.o Preconditioner.o PoissonSolver.o
	$(CC) -o $@ $+

PoissonMatrix.o: PoissonMatrix.cpp classes.h
	$(CC) -c -o $@ $<

Algorithms.o: Algorithms.cpp classes.h
	$(CC) -c -o $@ $<

CGVectors.o: CGVectors.cpp classes.h
	$(CC) -c -o $@ $<

UpperMatrix.o: UpperMatrix.cpp classes.h
	$(CC) -c -o $@ $<

LowerMatrix.o: LowerMatrix.cpp classes.h
	$(CC) -c -o $@ $<

Operators.o: Operators.cpp classes.h
	$(CC) -c -o $@ $<

Preconditioner.o: Preconditioner.cpp classes.h
	$(CC) -c -o $@ $<

PoissonSolver.o: PoissonSolver.cpp classes.h
	$(CC) -c -o $@ $<

clean:
	rm *.o poissonSolver