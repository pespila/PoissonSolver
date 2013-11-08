#include "classes.h"

int main(int argc, char const *argv[]) {
    int arg, n, dim;
    if (argc > 1) {
        arg = atoi(&*argv[1]);
    } else {
        arg = 4;
    }
    n=arg-1;
    dim = n*n;

    PoissonMatrix A(n);
    Operators O(n);
    CGVectors V(n,O);
    Algorithms Run(A);
    LowerMatrix L(n);
    UpperMatrix U(n);

    // double time,start=0.0,end=0.0;
    // start = clock();

    //Run.modifiedIncompleteLU(A,L,U);
    Run.incompleteLU(A,L,U);

    // end = clock();
    // time=(end-start)/CLOCKS_PER_SEC;
    // printf("%f\n", time);

    //Run.CG(A,O,V);
    Run.PCG(A,O,L,U,V);

    V.PrintVector();
    V.WriteToFile(O);

    return EXIT_SUCCESS;
}