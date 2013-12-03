#include "classes.h"

int main(int argc, char const *argv[]) {
    int arg,n;//,m;
    if (argc > 1) {
        arg = atoi(&*argv[1]);
    } else {
        arg = 4;
    }
    n=arg;
    //m=n-1;

    PoissonMatrix A(n);
    Operators O;
    Vectors V(n,O);
    Algorithms Run(A);
    LowerMatrix L(n);
    UpperMatrix U(n);

    A.Preconditioning();

    printf("Started\n");
    double time,start=0.0,end=0.0;
    start = clock();

    //Run.LU(A,L,U);
    //Run.incompleteLU(A,L,U);
    //Run.incompleteCholesky(A,L,U);
    //Run.modifiedIncompleteLU(A,L,U);

    Run.CG(A,O,V);
    //Run.JacobiMethod(A,O,V.x,V.b,5000);
    //Run.GaussSeidelMethod(A,O,V.x,V.b,5000);
    //Run.SORMethod(A,O,V);
    //Run.SSORMethod(A,O,V);
    //Run.PCG(A,O,L,U,V);
    //Run.MultiGridMethod(V.x,V.b,O,n);

    end = clock();
    time=(end-start)/CLOCKS_PER_SEC;
    printf("%f\n", time);

    A.PrintMatrix();
    //L.PrintMatrix();
    //U.PrintMatrix();

    V.PrintVector();
    V.WriteToFile(O);

    return EXIT_SUCCESS;
}