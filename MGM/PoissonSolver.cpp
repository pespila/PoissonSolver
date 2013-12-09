#include "classes.h"

int main(int argc, char const *argv[]) {
    int arg,n;//,m;
    if (argc > 1) {
        arg = atoi(&*argv[1]);
    } else {
        arg = 4;
    }
    n=arg-1;
    //m=n-1;

    PoissonMatrix A(n);
    Operators O;
    Vectors V(n,O);
    Algorithms Run;

    Run.InitHashMatrix(n);

    //Run.JacobiMethod(A,O,V.x,V.b);
    //Run.GaussSeidelMethod(A,O,V.x,V.b,3);
    //Run.SORMethod(A,O,V.x,V.b);
    V.x=Run.MultiGridMethod(V.x,V.b,n,O);

    A.PrintMatrix();

    V.PrintVector();
    V.WriteToFile(O);

    return EXIT_SUCCESS;
}