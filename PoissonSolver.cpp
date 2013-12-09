#include "classes.h"

int main(int argc, char const *argv[]) {
    int arg,n;//,m;
    if (argc > 1) {
        arg = atoi(&*argv[1]);
    } else {
        arg = 8;
    }
    n=arg-1;
    //m=n-1;

    PoissonMatrix A(n);
    LowerMatrix L(n);
    UpperMatrix U(n);
    PoissonVector P(n);
    MGMVector M(n);
    PoissonOperators O;
    Algorithms Run;

    Run.InitHashMatrix(n);

    printf("Started\n");
    double time,start=0.0,end=0.0;
    start = clock();

    //Run.LU(A,L,U,O);
    //Run.incompleteLU(A,L,U,O);
    //Run.incompleteCholesky(A,L,U,O);
    //Run.modifiedIncompleteLU(A,L,U,O);

    //Run.JacobiMethod(A,O,V.x,V.b);
    //Run.GaussSeidelMethod(A,O,V.x,V.b,3);
    //Run.SORMethod(A,O,V.x,V.b);
    //Run.SSORMethod(A,O,V);
    //Run.CG(A,O,V);
    //Run.PCG(A,O,L,U,V);
    Run.MultiGridMethod(V.x,V.b,n,O);

    end = clock();
    time=(end-start)/CLOCKS_PER_SEC;
    printf("%f\n", time);

    A.PrintMatrix();
    //L.PrintMatrix();
    //U.PrintMatrix();

    V.PrintVector();
    //Vmgm.PrintVector();
    //V.WriteToFile(O);

    return EXIT_SUCCESS;
}