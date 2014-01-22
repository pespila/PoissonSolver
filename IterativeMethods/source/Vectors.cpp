#include "classes.h"

Boundary::Boundary(int n) {
    this->n=n;
    this->dim=n*n;
    this->h=1.0/(double)(n+1);
    this->b.resize(dim);

    for(int i=1,k=0;i<=n;i++) {
        for(int j=1;j<=n;j++,k++) {
            b[k]=pow(h,2)*f(i*h,j*h);
            if(i==1) {
                b[k]+=pow(h,2)*pow(1.0/h,2)*g(j*h,0);
            }
            if(i==n) {
                b[k]+=pow(h,2)*pow(1.0/h,2)*g(j*h,1);
            }
            if(j==1) {
                b[k]+=pow(h,2)*pow(1.0/h,2)*g(0,i*h);
            }
            if(j==n) {
                b[k]+=pow(h,2)*pow(1.0/h,2)*g(1,i*h);
            }
        }
    }
}

Boundary::~Boundary() {
    vector<double>().swap(b);
}

double Boundary::Get(int i) {
    return this->b[i];
}

int Boundary::Size() {
    return b.size();
}

Startvector::Startvector(int n,double value) {
    this->n=n;
    this->dim=n*n;
    this->value=value;
    this->x.assign(dim,value);
}

Startvector::~Startvector() {
    vector<double>().swap(x);
}

double Startvector::Get(int i) {
    return this->x[i];
}

int Startvector::Size() {
    return x.size();
}

void Vectors::PrintVector() {
    for (int i=0;i<Size();i++) {
        printf("%f ",Get(i));
    }
    printf("\n");
}

void Vectors::WriteToFile() {
    FILE *file;
    file=fopen("../Plot/plot.txt","w");
    if(file==NULL)
        printf("ERROR: Could not open file!\n");
    else {
        for(int i=0,k=0;i<=sqrt(Size());i++) {
            for(int j=0;j<=sqrt(Size());j++) {
                if(i==0) {
                    fprintf(file,"%f %d %f\n",(double)j/(double)(sqrt(Size())),i,g(j*1.0/(double)(sqrt(Size()+1)),i));
                } else if(i!=0) {
                    if(j==0) {
                        fprintf(file,"%d %f %f\n",j,(double)i/(double)(sqrt(Size())),g(j,i*1.0/(double)(sqrt(Size()+1))));
                    } else if(j!=0){
                        fprintf(file,"%f %f %f\n",(double)j/(double)(sqrt(Size())),(double)i/(double)(sqrt(Size())),Get(k));
                        k++;
                    }
                }
            }
            fprintf(file,"\n");
        }
    }
    fclose (file);
}