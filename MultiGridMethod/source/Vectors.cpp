#include "classes.h"

PoissonVector::PoissonVector(int n, double a) {
    this->n=n;
    this->h=1.0/(double)(n+1);
    x.assign(n*n,a);
    solved.assign(n*n,0);
    b.assign(n*n,0);

    for(int i=1,k=0;i<=n;i++) {
        for(int j=1;j<=n;j++,k++) {
            solved[k]=g(j*h,i*h);
        }
    }

    for(int i=1,k=0;i<=n;i++) {
        for(int j=1;j<=n;j++,k++) {
            b[k]=pow(h,2)*f(h,h);
            if(i==1) b[k]+=g(j*h,0);
            if(i==n) b[k]+=g(j*h,1);
            if(j==1) b[k]+=g(0,i*h);
            if(j==n) b[k]+=g(1,i*h);
        }
    }
}

PoissonVector::~PoissonVector() {
    vector<double>().swap(x);
    vector<double>().swap(b);
    vector<double>().swap(solved);
}

void PoissonVector::WriteToFile() {
    FILE *file;
    file=fopen("../Plot/plot.txt", "w");
    if(file==NULL) {
        printf("ERROR: Could not open file!\n");
    } else {
        for(int i=1,k=0;i<=n;i++) {
            for(int j=1;j<=n;j++,k++) {
                fprintf(file, "%f %f %f\n", i*h,j*h,x[k]);
            }
            fprintf(file, "\n");
        }
    }
    fclose (file);
}

void PoissonVector::PrintVector() {
    for(int i=0;i<n*n;i++) {
            printf("%f ", x[i]);
    }
    printf("\n");
}