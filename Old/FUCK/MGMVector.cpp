#include "classes.h"

MGMVector::MGMVector(int m) {
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
    // this->n=n;
    // this->dim=n*n;
    // this->h=1.0/(double)(n+1);
    // b.resize(dim);
    // x.resize(dim);

    // vector<double> pushX;
    // pushX.assign(n,0);
    // x.assign(n,pushX);
    // vector<double>().swap(pushX);

    // vector<double> pushB;
    // pushB.resize(n);
    // int k,i,j;
    // for(i=1;i<=n;i++) {
    //     for(j=1,k=0;j<=n;j++,k++) {
    //         if(i==1) {
    //             if (j==1) {
    //                 pushB[k]=f(h,h)+(n+1)*(n+1)*(g(0,h)+g(h,0));
    //             } else if (j==n) {
    //                 pushB[k]=f(1-h,1-h)+(n+1)*(n+1)*(g(1-h,0)+g(1,h));
    //             } else {
    //                 pushB[k]=f(j*h,h)+(n+1)*(n+1)*g(j*h,0);
    //             }
    //         } else if(i>1 && i<n) {
    //             if (j==1) {
    //                 pushB[k]=f(h,i*h)+(n+1)*(n+1)*g(0,i*h);            
    //             } else if (j==n) {
    //                 pushB[k]=f(1-h,i*h)+(n+1)*(n+1)*g(1,i*h);
    //             } else {
    //                 pushB[k]=f(j*h,i*h);
    //             }
    //         } else {
    //             if (j==1) {
    //                 pushB[k]=f(h,1-h)+(n+1)*(n+1)*(g(h,1)+g(0,1-h));
    //             } else if (j==n) {
    //                 pushB[k]=f(1-h,1-h)+(n+1)*(n+1)*(g(1-h,1)+g(1,1-h));
    //             } else {
    //                 pushB[k]=f(j*h,1-h)+(n+1)*(n+1)*g(j*h,1);
    //             }
    //         }
    //     }
    //     b.push_back(pushB);
    // }
    // vector<double>().swap(pushB);
}

MGMVector::~MGMVector() {
    vector<double>().swap(x);
    vector<double>().swap(b);
    // vector<vector<double> >().swap(x);
}

double MGMVector::f(double x,double y) {
    return (-4.0);
}

double MGMVector::g(double x,double y) {
    return (x*x+y*y);
}

void MGMVector::PrintVector() {
    if(this->dim<=20)
        for(int i=0;i<this->dim;i++) {
            printf("%f ",x[i]);
        }
    printf("\n");
    // if(this->dim<=16) {
    //     for(int i=0;i<this->n;i++) {
    //         for(int j=0;j<this->n;j++) {
    //             printf("%f ",x[i][j]);
    //         }
    //     }
    // }
    // printf("\n");
}

void MGMVector::WriteToFile() {
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
    // int i,j;
    // FILE *file;
    // file=fopen("./Plot/plot.txt", "w");
    // if(file==NULL)
    //     printf("ERROR: Could not open file!\n");
    // for(i=0;i<=n;i++) {
    //     for(j=0;j<=n;j++) {
    //         if(i==0) {
    //             fprintf(file, "%f %d %f\n", (double)j/(double)n,i,g(j*h,i));
    //         } else if(i!=0) {
    //             if(j==0) {
    //                 fprintf(file, "%d %f %f\n", j,(double)i/(double)n,g(j,i*h));
    //             } else if(j!=0){
    //                 fprintf(file, "%f %f %f\n", (double)j/(double)n,(double)i/(double)n,x[i-1][j-1]);
    //             }
    //         }
    //     }
    //     fprintf(file,"\n");
    // }
    // fclose (file);
}