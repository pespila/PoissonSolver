#include "classes.h"

PoissonVector::PoissonVector(int n) {
    this->n=n;
    this->dim=n*n;
    this->h=1.0/(double)(n+1);
}

PoissonVector::~PoissonVector() {
}

void PoissonVector::InitRightSide(vector<double>& b) {
    if((int)b.size()!=this->dim) printf("Error in InitRightSide : Size of vector f does not match!\n");
    else {
        for(int i=1,k=0;i<=this->n;i++) {
            for(int j=1;j<=this->n;j++,k++) {
                b[k]=pow(h,2)*f(h,h);
                if(i==1) {
                    b[k]+=g(j*h,0);
                }
                if(i==n) {
                    b[k]+=g(j*h,1);
                }
                if(j==1) {
                    b[k]+=g(0,i*h);
                }
                if(j==n) {
                    b[k]+=g(1,i*h);
                }
            }
        }
    }
}

void PoissonVector::InitSolution(vector<double>& solved) {
    if((int)solved.size()!=this->dim) printf("Error in InitSolution : Size of vector solved does not match!\n");
    else {
        for(int i=1,k=0;i<=n;i++) {
            for(int j=1;j<=n;j++,k++) {
                solved[k]=g(j*h,i*h);
            }
        }
    }
}

void PoissonVector::WriteToFile(const vector<double>& x) {
    if((int)x.size()!=this->dim) printf("Error in WriteToFile: Dimension doesn't match!\n");
    FILE *file;
    file=fopen("../Plot/plot.txt", "w");
    if(file==NULL) {
        printf("Error in WriteToFile: Could not open file!\n");
    } else {
        for(int i=1,k=0;i<=n;i++) {
            for(int j=1;j<=n;j++,k++) {
                fprintf(file, "%f %f %f\n", i*this->h,j*this->h,x[k]);
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