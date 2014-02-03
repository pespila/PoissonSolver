#include "classes.h"

PoissonVector1D::PoissonVector1D(int n) {
    this->n=n;
    this->dim=n;
    this->h=1.0/(double)(n+1);
}

PoissonVector1D::~PoissonVector1D() {
}

void PoissonVector1D::InitRightSide(vector<double>& b) {
    if((int)b.size()!=this->dim) printf("Size of Vector does not match!\n");
    else {
        for(int i=0;i<this->dim;i++) {
            b[i]=pow(h,2)*f(h);
            if(i==0) {
                b[i]+=g(0);
            }
            if(i==dim-1) {
                b[i]+=g(1);
            }
        }
    }
}

void PoissonVector1D::WriteToFile(const vector<double>& x) {
    FILE *file;
    file=fopen("../Plot/plot.txt", "w");
    if(file==NULL) {
        printf("ERROR: Could not open file!\n");
    } else if ((int)x.size()!=this->dim) {
        printf("Size of Vector does not match!\n");
    } else {
        for(int i=0;i<this->dim;i++) {
            if(i==0) {
                fprintf(file, "%d %f\n", i,g(i));
            }
            fprintf(file, "%f %f\n", (double)i*h,x[i]);
            if(i==dim-1) {
                fprintf(file, "%d %f\n", 1,g(1));
            }
        }
        fprintf(file, "\n");
    }
    fclose (file);
}

PoissonVector2D::PoissonVector2D(int n) {
    this->n=n;
    this->dim=n;
    this->h=1.0/(double)(n+1);
}

PoissonVector2D::~PoissonVector2D() {
}

void PoissonVector2D::InitRightSide(vector<double>& b) {
    if(sqrt(b.size())!=this->n) printf("Size of Vector b does not match!\n");
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

void PoissonVector2D::WriteToFile(const vector<double>& x) {
    FILE *file;
    file=fopen("../Plot/plot.txt", "w");
    if(file==NULL) {
        printf("ERROR: Could not open file!\n");
    } else if (sqrt(x.size())!=this->n) {
        printf("Size of Vector does not match!\n");
    }  else {
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