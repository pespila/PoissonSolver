#include "classes.h"

Operators::Operators(int n) {
    this->n=n;
    this->dim=n*n;
}

Operators::~Operators() {
}

double Operators::innerProduct(const vector<double>& x,const vector<double>& y) {
    double ip=0.0;
    for (int i=0;i<dim;i++)
        ip+=x[i]*y[i];
    return ip;
}

double Operators::vectorNorm(const vector<double>& x) {
    return sqrt(innerProduct(x,x));
}

void Operators::MatrixVectorMultiplyer(Matrix& M,const vector<double>& x,vector<double>& b) {
    for (int i=0;i<dim;i++) {
        b[i]+=x[i]*M.Get(i,i);
        if(i<(dim-n)){
            b[i]+=x[i+n]*M.Get(i,i+n);
        } 
        if(i>=n) {
            b[i]+=x[i-n]*M.Get(i,i-n);
        }
        if(i%n!=0) {
            b[i]+=x[i-1]*M.Get(i,i-1);
            b[i-1]+=x[i]*M.Get(i-1,i);
        }
    }
}

void Operators::LUsolverLower(Matrix& A,Matrix& L,vector<double>& z) {
    int m;
    for(int i=0;i<dim;i++) {
        for(int j=0;j<5;j++) {
            m=A.HashMatrix[i][j];
            if(m!=-1 && m<i) {
                z[i]-=L.Get(i,m)*z[m];
            }
        }
        z[i]/=L.Get(i,i);
    }
}

void Operators::LUsolverUpper(Matrix& A,Matrix& U,vector<double>& z) {
    int m;
    for(int i=dim-1;i>=0;i--) {
        for(int j=0;j<5;j++) {
            m=A.HashMatrix[i][j];
            if(m!=-1 && m>=i) {
                z[i]-=U.Get(i,m)*z[m];
            }
         }
        z[i]/=U.Get(i,i);
    }
}