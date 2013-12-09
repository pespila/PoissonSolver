CC = g++ -Wall -g

poissonSolver: PoissonMatrix.o Algorithms.o PoissonVector.o UpperMatrix.o LowerMatrix.o PoissonOperators.o PoissonSolver.o
	$(CC) -o $@ $+

PoissonMatrix.o: PoissonMatrix.cpp classes.h
	$(CC) -c -o $@ $<

Algorithms.o: Algorithms.cpp classes.h
	$(CC) -c -o $@ $<

PoissonVector.o: PoissonVector.cpp classes.h
	$(CC) -c -o $@ $<

UpperMatrix.o: UpperMatrix.cpp classes.h
	$(CC) -c -o $@ $<

LowerMatrix.o: LowerMatrix.cpp classes.h
	$(CC) -c -o $@ $<

PoissonOperators.o: PoissonOperators.cpp classes.h
	$(CC) -c -o $@ $<

PoissonSolver.o: PoissonSolver.cpp classes.h
	$(CC) -c -o $@ $<

clean:
	rm *.o poissonSolver