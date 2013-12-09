#include "classes.h"

CGVector::CGVector(int m) {
    n=m;
    dim=n*n;
    h=1.0/(double)(n+1);
    b.resize(dim);
    x.resize(dim);

    int k,i,j;
    for(j=1,k=0;j<=n;j++,k++) {
        if(j==1) {
            b[k]=f(h,h)+pow((n+1),2)*(g(0,h)+g(h,0));
        } else if(j==n) {
            b[k]=f(1-h,1-h)+pow((n+1),2)*(g(1-h,0)+g(1,h));
        } else {
            b[k]=f(j*h,h)+pow((n+1),2)*g(j*h,0);
        }
    }
    for(i=2;i<=(n-1);i++,k++) {
        for (j=1;j<=n;j++) {
            if(j==1) {
                b[k]=f(h,i*h)+pow((n+1),2)*g(0,i*h);            
            } else if(j==n) {
                b[k]=f(1-h,i*h)+pow((n+1),2)*g(1,i*h);
            } else {
                b[k]=f(j*h,i*h);
            }
        }
    }
    for(j=1;j<=n;j++,k++) {
          if (j==1) {
            b[k]=f(h, 1-h)+pow((n+1),2)*(g(h,1)+g(0,1-h));
        } else if (j==n) {
            b[k]=f(1-h,1-h)+pow((n+1),2)*(g(1-h,1)+g(1,1-h));
        } else {
            b[k]=f(j*h,1-h)+pow((n+1),2)*g(j*h,1);
        }
    }
}

CGVector::~CGVector() {
    vector<double>().swap(x);
    vector<double>().swap(b);
}

double CGVector::f(double x,double y) {
    return (-4.0);
}

double CGVector::g(double x,double y) {
    return (x*x+y*y);
}

void CGVector::PrintVector() {
    if(this->dim<=20)
        for(int i=0;i<this->dim;i++) {
            printf("%f ",x[i]);
        }
    printf("\n");
}

void CGVector::WriteToFile() {
    int i,j,k;
    FILE *file;
    file=fopen("./Plot/plot.txt", "w");
    if(file==NULL)
        printf("ERROR: Could not open file!\n");
    for(i=0,k=0;i<=n;i++) {
        for(j=0;j<=n;j++) {
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
        fprintf(file,"\n");
    }
    fclose (file);
}