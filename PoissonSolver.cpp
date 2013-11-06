#include "classes.h"

int main(int argc, char const *argv[]) {
	int arg, n;//, dim;
	if (argc > 1) {
        arg = atoi(&*argv[1]);
    } else {
        arg = 4;
    }
    n=arg-1;
	//dim = n*n;

    PoissonMatrix A(n);
    Operators O(n);
    CGVectors V(n,O);
    Algorithms Run(A);

    LowerPoissonMatrix L(n);
    UpperPoissonMatrix U(n);

    //Run.CG(A,O,V.x,V.b);
    Run.modifiedIncompleteLU(A,L,U);
    //Run.incompleteLU(A,L,U);
    printf("done\n");
    // A.printMat();
    //L.printMat();
    //U.printMat();

    Run.PCG(A,O,L,U,V.x,V.b);
    
    // for(int i=0;i<dim;i++)
    //     printf("%f ", V.x[i]);
    // printf("\n");

	return EXIT_SUCCESS;
}