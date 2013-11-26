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

    double time,start=0.0,end=0.0;
    start = clock();

    //Run.modifiedIncompleteLU(A,L,U);
    Run.incompleteLU(A,L,U);
    //Run.incompleteCholesky(A,L,U);
    //Run.modifiedIncompleteCholesky(W,L,U);
    //Run.LU(A,L,U);

    //Run.CG(A,O,V);
    //Run.JacobiMethod(A,O,V);
    //Run.GaussSeidelMethod(A,O,V);
    //Run.SORMethod(A,O,V);
    Run.PCG(A,O,L,U,V);

    end = clock();
    time=(end-start)/CLOCKS_PER_SEC;
    printf("%f\n", time);

    // A.PrintMatrix();

    // for(int i=0;i<n*n;i++) {
    //     for(int j=0;j<5;j++) {
    //         printf("%d ", A.HashMatrix[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");


    // vector<vector<int> > Restriktion;

    // int i,j,k;
    // for(i=0;i<n*n;i++) {
    //     vector<int> push;
    //     push.assign(3,0);
    //     k=0;
    //     for(j=0;j<dim;j++) {
    //         if(k==0) push[k]=(j);
    //         if(k==1) push[k]=(j);
    //         if(k==2) push[k]=(j);
    //         k++;
    //     }
    //     Restriktion.push_back( push );
    // }




    // A.PrintMatrix();
    // L.PrintMatrix();
    // U.PrintMatrix();

    //V.PrintVector();
    V.WriteToFile(O);

    return EXIT_SUCCESS;
}