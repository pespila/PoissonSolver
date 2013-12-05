CC = g++ -Wall -g

poissonSolver: PoissonMatrix.o Algorithms.o Vectors.o Operators.o PoissonSolver.o
	$(CC) -o $@ $+

PoissonMatrix.o: PoissonMatrix.cpp classes.h
	$(CC) -c -o $@ $<

Algorithms.o: Algorithms.cpp classes.h
	$(CC) -c -o $@ $<

Vectors.o: Vectors.cpp classes.h
	$(CC) -c -o $@ $<

Operators.o: Operators.cpp classes.h
	$(CC) -c -o $@ $<

PoissonSolver.o: PoissonSolver.cpp classes.h
	$(CC) -c -o $@ $<