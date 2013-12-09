#include "classes.h"

Vectors::Vectors(int m,Operators& O) {
    n=m;
    dim=n*n;
    b.resize(dim);
    x.resize(dim);

    int k,i,j;
    double h=1.0/(double)(n+1);

    k=0;
    for (j=1;j<=n;j++) {
        if (j == 1) {
            b[k] = f(h, h) + (n+1)*(n+1) * (g(0, h) + g(h, 0));
        } else if (j == n) {
            b[k] = f(1-h, 1-h) + (n+1)*(n+1) * (g(1-h, 0) + g(1, h));
        } else {
            b[k] = f(j*h, h) + (n+1)*(n+1) * g(j*h, 0);
        }
        k++;
    }
    
    for (i=2; i<=(n-1);i++) {
        for (j = 1; j <= n; j++) {
            if (j == 1) {
                b[k] = f(h, i*h) + (n+1)*(n+1) * g(0, i*h);            
            } else if (j == n) {
                b[k] = f(1-h, i*h) + (n+1)*(n+1) * g(1, i*h);
            } else {
                b[k] = f(j*h, i*h);
            }
            k++;
        }
    }
    
    for (j=1;j<=n;j++) {
          if (j == 1) {
            b[k] = f(h, 1-h) + (n+1)*(n+1) * (g(h, 1)+g(0, 1-h));
        } else if (j == n) {
            b[k] = f(1-h, 1-h) + (n+1)*(n+1) * (g(1-h, 1)+g(1, 1-h));
        } else {
            b[k] = f(j*h, 1-h) + (n+1)*(n+1) * g(j*h, 1);
        }
        k++;
    }
}

Vectors::~Vectors() {
    vector<double>().swap(x);
    vector<double>().swap(b);
}

double Vectors::f(double x,double y) {
    return (-4.0);
}

double Vectors::g(double x,double y) {
    return (pow(x,2)+pow(y,2));
}

void Vectors::PrintVector() {
    if(dim<=16)
        for (int i=0;i<dim;i++) {
            printf("%f ",x[i]);
        }
    printf("\n");
}

void Vectors::WriteToFile(Operators& O) {
    int i,j,k;
    double h=1.0/(double)(n+1);
    FILE *file;
    file=fopen("./Plot/plot.txt", "w");
    k=0;
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
                    fprintf(file, "%f %f %f\n", (double)j/(double)n,(double)i/(double)n,x[k]);
                    k++;
                }
            }
        }
        fprintf(file, "\n");
    }
    fclose (file);
}