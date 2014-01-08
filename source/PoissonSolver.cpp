#include "classes.h"

int main(int argc, char const *argv[]) {
    int arg,alg;
    if (argc==2) {
        arg=atoi(&*argv[1])-1;
    } else {
        arg=39;
    }

    PoissonMatrix A(arg);
    LowerMatrix L(arg);
    UpperMatrix U(arg);
    Operators O(arg);
    Vectors V(arg);
    Algorithms Run(arg);

    printf("Choose your Algorithm: \n[0] CG\n[1] JacobiMethod\n[2] GaussSeidelMethod\n[3] SORMethod\n[4] PCG\n");
    scanf("%d",&alg);

    printf("Started\n");
    double time,start=0.0,end=0.0;
    start = clock();

    if(alg==0) Run.CG(A,O,V);
    if(alg==1) Run.JacobiMethod(A,O,V);
    if(alg==2) Run.GaussSeidelMethod(A,O,V);
    if(alg==3) Run.SORMethod(A,O,V);
    if(alg==4) {
        Run.modifiedIncompleteLU(A,L,U);
        Run.PCG(A,O,L,U,V);
    }

    if(arg<=5) {
        A.PrintMatrix();
        L.PrintMatrix();
        U.PrintMatrix();
        V.PrintVector();
    }

    end = clock();
    time=(end-start)/CLOCKS_PER_SEC;
    printf("%f\n", time);

    V.WriteToFile();

    return EXIT_SUCCESS;
}

double f(double x,double y) {
    return -4.0;
}

double g(double x,double y) {
    return pow(x,2)+pow(y,2);
}