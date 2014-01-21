#include "classes.h"

Operators::Operators(int n) {
    this->n=n;
    this->dim=n*n;
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