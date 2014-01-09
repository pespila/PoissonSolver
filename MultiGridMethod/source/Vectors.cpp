#include "classes.h"

Vectors::Vectors(int n) {
    this->n=n;
    this->dim=n*n;
    this->h=1.0/(double)(n+1);
    this->b.assign(dim,0);
    this->x.assign(dim,0);
    int i,j,k;

    for(i=1,k=0;i<=n;i++) {
        for(j=1;j<=n;j++,k++) {
            b[k]=f(i*h,j*h);
            if(i==1) {
                b[k]+=pow(1.0/h,2)*g(j*h,0);
            }
            if(i==n) {
                b[k]+=pow(1.0/h,2)*g(j*h,1);
            }
            if(j==1) {
                b[k]+=pow(1.0/h,2)*g(0,i*h);
            }
            if(j==n) {
                b[k]+=pow(1.0/h,2)*g(1,i*h);
            }
        }
    }
}

Vectors::~Vectors() {
    vector<double>().swap(x);
    vector<double>().swap(b);
}

void Vectors::PrintVector() {
    for(int i=0;i<dim;i++) {
        printf("%f ", x[i]);
    }
    printf("\n");
}

void Vectors::WriteToFile() {
    FILE *file;
    file=fopen("../Plot/plot.txt", "w");
    if(file==NULL) {
        printf("ERROR: Could not open file!\n");
    } else {
        for(int i=0,k=0;i<=n;i++) {
            for(int j=0;j<=n;j++) {
                if(i==0) {
                    fprintf(file, "%f %d %f\n", (double)j/(double)n,i,g(j*h,i));
                } else if(i!=0) {
                    if(j==0) {
                        fprintf(file, "%d %f %f\n", j,(double)i/(double)n,g(j,i*h));
                    } else if(j!=0){
                        fprintf(file, "%f %f %f\n", (double)j/(double)n,(double)i/(double)n,x[k]);
                        k++;
                    }
                }
            }
            fprintf(file, "\n");
        }
    }
    fclose (file);
}