#include "classes.h"

int main(int argc, char const *argv[]) {
    int arg;
    if (argc==2) {
        arg=atoi(&*argv[1])-1;
    } else {
        arg=31;
    }

    PoissonMatrix A(arg);
    Vectors V(arg);
    Algorithms Run(arg);

    printf("Started\n");
    double time,start=0.0,end=0.0;
    start = clock();

    Run.MultiGridMethod(A,V);

    end = clock();
    time=(end-start)/CLOCKS_PER_SEC;
    printf("%f\n", time);

    if(arg<=5) {
        A.PrintMatrix();
        V.PrintVector();
    } else {
        printf("Couldn't print Matrix! Dimension is too high.\n");
    }
    V.WriteToFile();

    return EXIT_SUCCESS;
}

double f(double x,double y) {
    return -4.0;
}

double g(double x,double y) {
    return pow(x,2)+pow(y,2);
}