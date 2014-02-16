CC = g++ -Wall -g

poissonSolver: PoissonMatrix.o Algorithms.o Vector.o LUMatrices.o PoissonSolver.o
	$(CC) -o $@ $+

PoissonMatrix.o: PoissonMatrix.cpp classes.h
	$(CC) -c -o $@ $<

Algorithms.o: Algorithms.cpp classes.h
	$(CC) -c -o $@ $<

Vector.o: Vector.cpp classes.h
	$(CC) -c -o $@ $<

LUMatrices.o: LUMatrices.cpp classes.h
	$(CC) -c -o $@ $<

PoissonSolver.o: PoissonSolver.cpp classes.h
	$(CC) -c -o $@ $<

clean:
	rm *.o