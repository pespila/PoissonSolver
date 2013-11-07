#include "classes.h"

void writeToFile(LowerPoissonMatrix&,UpperPoissonMatrix&);
void readFromFile(LowerPoissonMatrix&,UpperPoissonMatrix&);

int main(int argc, char const *argv[]) {
    int arg, n;//, dim;
    if (argc > 1) {
        arg = atoi(&*argv[1]);
    } else {
        arg = 4;
    }
    n=arg-1;
    //dim = n*n;

    PoissonMatrix A(n);
    Operators O(n);
    CGVectors V(n,O);
    Algorithms Run(A);

    LowerPoissonMatrix L(n);
    UpperPoissonMatrix U(n);

    //Run.CG(A,O,V.x,V.b);
    //Run.modifiedIncompleteLU(A,L,U);
    //Run.incompleteLU(A,L,U);
    //printf("done\n");
    //writeToFile(L,U);
    readFromFile(L,U);
    // A.printMat();
    //L.printMat();
    //U.printMat();

    Run.PCG(A,O,L,U,V.x,V.b);
    
    // for(int i=0;i<dim;i++)
    //     printf("%f ", V.x[i]);
    // printf("\n");

    return EXIT_SUCCESS;
}

void writeToFile(LowerPoissonMatrix& L, UpperPoissonMatrix& U) {
    FILE *file;
    file=fopen("LU.txt", "w");
    if(file==NULL)
        printf("ERROR!\n");
    for(int i=0;i<L.size();i++) {
        for(int j=0;j<L.size();j++) {
            fprintf(file, "%f\n", L.get(i,j));
            fprintf(file, "%f\n", U.get(i,j));
        }
    }
    fclose (file);
}

void readFromFile(LowerPoissonMatrix& L,UpperPoissonMatrix& U) {
    FILE *file;
    file = fopen("LU.txt", "r");
    double l,u;
    if(file==NULL)
        printf("ERROR!\n");
    for(int i=0;i<L.size();i++) {
        for(int j=0;j<L.size();j++) {
            fscanf(file, "%lf\n", &l);
            fscanf(file, "%lf\n", &u);
            L.set(i,j,l);
            U.set(i,j,u);
        }
    }
    fclose(file);
}