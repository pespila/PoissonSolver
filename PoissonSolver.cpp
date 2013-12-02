#include "classes.h"

int main(int argc, char const *argv[]) {
    int arg, n;
    if (argc > 1) {
        arg = atoi(&*argv[1]);
    } else {
        arg = 4;
    }
    n=arg-1;

    PoissonMatrix A(n);
    Operators O(n);
    Vectors V(n,O);
    Algorithms Run(A);
    LowerMatrix L(n);
    UpperMatrix U(n);
    //Preconditioner W(n);

    printf("Started\n");
    double time,start=0.0,end=0.0;
    start = clock();

    //Run.LU(A,L,U);
    //Run.incompleteLU(A,L,U);
    //Run.incompleteCholesky(A,L,U);
    Run.modifiedIncompleteLU(A,L,U);
    //Run.modifiedIncompleteCholesky(W,L,U);

    //Run.CG(A,O,V);
    //Run.JacobiMethod(A,O,V);
    //Run.GaussSeidelMethod(A,O,V);
    //Run.SORMethod(A,O,V);
    //Run.SSORMethod(A,O,V);
    Run.PCG(A,O,L,U,V);

    end = clock();
    time=(end-start)/CLOCKS_PER_SEC;
    printf("%f\n", time);

    A.PrintMatrix();
    L.PrintMatrix();
    U.PrintMatrix();

    V.PrintVector();
    //V.WriteToFile(O);

    printf("\n");
    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) {
            printf("o ");
        }
        printf("\n");
    }
    printf("\n");

    printf("\n");
    for(int i=0;i<n/2;i++) {
        printf("  ");
        for(int j=0;j<n/2;j++) {
            printf("o   ");
        }
        printf("\n");
        printf("\n");
    }
    printf("\n");

    printf("\n");
    for(int i=0;i<n/4;i++) {
        printf("    ");
        for(int j=0;j<n/4;j++) {
            printf("o    ");
        }
        printf("\n");
        printf("\n");
        printf("\n");
    }
    printf("\n");

    return EXIT_SUCCESS;
}