#include "classes.h"

MGMVector::MGMVector(int n) {
    this->n=n;
    this->dim=n*n;
    this->h=1.0/(double)(n+1);
    b.resize(dim);
    x.resize(dim);

    vector<double> pushX;
    pushX.assign(n,0);
    x.assign(n,pushX);
    vector<double>().swap(pushX);

    vector<double> pushB;
    pushB.resize(n);
    int k,i,j;
    for(i=1;i<=n;i++) {
        for(j=1,k=0;j<=n;j++,k++) {
            if(i==1) {
                if (j==1) {
                    pushB[k]=O.f(h,h)+(n+1)*(n+1)*(O.g(0,h)+O.g(h,0));
                } else if (j==n) {
                    pushB[k]=O.f(1-h,1-h)+(n+1)*(n+1)*(O.g(1-h,0)+O.g(1,h));
                } else {
                    pushB[k]=O.f(j*h,h)+(n+1)*(n+1)*O.g(j*h,0);
                }
            } else if(i>1 && i<n) {
                if (j==1) {
                    pushB[k]=O.f(h,i*h)+(n+1)*(n+1)*O.g(0,i*h);            
                } else if (j==n) {
                    pushB[k]=O.f(1-h,i*h)+(n+1)*(n+1)*O.g(1,i*h);
                } else {
                    pushB[k]=O.f(j*h,i*h);
                }
            } else {
                if (j==1) {
                    pushB[k]=O.f(h,1-h)+(n+1)*(n+1)*(O.g(h,1)+O.g(0,1-h));
                } else if (j==n) {
                    pushB[k]=O.f(1-h,1-h)+(n+1)*(n+1)*(O.g(1-h,1)+O.g(1,1-h));
                } else {
                    pushB[k]=O.f(j*h,1-h)+(n+1)*(n+1)*O.g(j*h,1);
                }
            }
        }
        b.push_back(pushB);
    }
    vector<double>();swap(pushB);
}

MGMVector::~MGMVector() {
    vector<vector<double> >().swap(x);
}

double MGMVector::f(double x,double y) {
    return (-4.0);
}

double MGMVector::g(double x,double y) {
    return (x*x+y*y);
}

void MGMVector::PrintVector() {
    if(this->dim<=16) {
        for(int i=0;i<this->n;i++) {
            for(int j=0;j<this->n;j++) {
                printf("%f ",x[i][j]);
            }
        }
    }
    printf("\n");
}

void MGMVector::WriteToFile() {
    int i,j;
    FILE *file;
    file=fopen("./Plot/plot.txt", "w");
    if(file==NULL)
        printf("ERROR: Could not open file!\n");
    for(i=0;i<=n;i++) {
        for(j=0;j<=n;j++) {
            if(i==0) {
                fprintf(file, "%f %d %f\n", (double)j/(double)n,i,g(j*h,i));
            } else if(i!=0) {
                if(j==0) {
                    fprintf(file, "%d %f %f\n", j,(double)i/(double)n,g(j,i*h));
                } else if(j!=0){
                    fprintf(file, "%f %f %f\n", (double)j/(double)n,(double)i/(double)n,x[i-1][j-1]);
                }
            }
        }
        fprintf(file,"\n");
    }
    fclose (file);
}