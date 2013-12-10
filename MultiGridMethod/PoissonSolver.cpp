#include "classes.h"

int main(int argc, char const *argv[]) {
    int n;
    if (argc>1) {
        n=atoi(&*argv[1])-1;
    } else {
        n=8-1;
    }

    PoissonMatrix A(n);
    Vectors V(n);
    Algorithms Run(n);

    A.PrintMatrix();

    //V.x=Run.MultiGridMethod(A,V,V.x,V.b);

    vector<double> r;
    r.assign(9,5);
    vector<double> b;
    b.assign(49,0);

    for(int i=0;i<49;i++) {
        printf("%.2f ", b[i]);
    }
    printf("\n");
    printf("\n");
    for(int i=0;i<9;i++) {
        printf("%.2f ", r[i]);
    }
    printf("\n");
    printf("\n");

    b=Run.Interpolation(r,7);

    for(int i=0;i<49;i++) {
        printf("%.2f ", b[i]);
    }
    printf("\n");
    printf("\n");
    for(int i=0;i<9;i++) {
        printf("%.2f ", r[i]);
    }
    printf("\n");
    printf("\n");

    V.PrintVector();
    V.WriteToFile();

    return EXIT_SUCCESS;
}