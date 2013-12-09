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

    V.x=Run.MultiGridMethod(A,V,V.x,V.b);

    V.PrintVector();
    V.WriteToFile();

    return EXIT_SUCCESS;
}