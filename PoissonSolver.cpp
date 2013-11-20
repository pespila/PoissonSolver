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
    Preconditioner W(n);
    Operators O(n);
    Vectors V(n,O);
    Algorithms Run(A);
    LowerMatrix L(n);
    UpperMatrix U(n);

    // double time,start=0.0,end=0.0;
    // start = clock();

    //Run.modifiedIncompleteLU(A,L,U);
    //Run.incompleteLU(A,L,U);
    //Run.incompleteCholesky(A,L,U);
    Run.modifiedIncompleteCholesky(W,L,U);
    A.PrintMatrix();
    L.PrintMatrix();
    U.PrintMatrix();

    // end = clock();
    // time=(end-start)/CLOCKS_PER_SEC;
    // printf("%f\n", time);

    //Run.CG(A,O,V);
    
    //V.PrintVector();

    //V.x.assign(dim,0);
    //Run.JacobiMethod(A,O,V);
    //Run.GaussSeidelMethod(A,O,V);
    //Run.SORMethod(A,O,V);

    Run.PCG(A,O,L,U,V);

    // A.PrintMatrix();
    // L.PrintMatrix();
    // U.PrintMatrix();

    V.PrintVector();
    V.WriteToFile(O);

    return EXIT_SUCCESS;
}