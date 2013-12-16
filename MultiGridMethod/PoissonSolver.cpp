#include "classes.h"

int main(int argc, char const *argv[]) {
    int n;
    if (argc>1) {
        n=atoi(&*argv[1])-1;
    } else {
        n=16-1;
    }

    PoissonMatrix A(n);
    Vectors V(n);
    Algorithms Run(n);

    A.PrintMatrix();
    A.InitHashMatrix();

    V.x=Run.MultiGridMethod(A,V,V.b,n);
    V.PrintVector();
    V.WriteToFile();

    return EXIT_SUCCESS;
}