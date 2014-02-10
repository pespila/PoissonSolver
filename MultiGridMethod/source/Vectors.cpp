#include "classes.h"

Vector::Vector(int n) {
    this->n=n;
    this->h=1.0/(double)(n+1);
    x.assign(n*n,0);
    b.resize(n*n);
    solved.resize(n*n);
    for(int i=1,k=0;i<=n;i++) {
        for(int j=1;j<=n;j++,k++) {
            solved[k]=g(j*h,i*h);
            b[k]=pow(h,2)*f(h,h);
            if(i==1) b[k]+=g(j*h,0);
            if(i==n) b[k]+=g(j*h,1);
            if(j==1) b[k]+=g(0,i*h);
            if(j==n) b[k]+=g(1,i*h);
        }
    }
}

Vector::~Vector() {
    vector<double>().swap(x);
    vector<double>().swap(b);
    vector<double>().swap(solved);
}

void Vector::WriteToFile(const vector<double>& x) {
    FILE *file;
    file=fopen("../Plot/plot.txt", "w");
    int n=sqrt(x.size());
    double h=1.0/(double)(n+1);
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

void Vector::PrintVector(const vector<double>& x) {
    for(int i=0;i<(int)x.size();i++) {
        printf("%f ", x[i]);
    }
    printf("\n");
}