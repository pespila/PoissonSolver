#include "classes.h"

Operators::Operators() {
}

Operators::~Operators() {
}

double Operators::innerProduct(const vector<double>& x, const vector<double>& y) {
    double ip=0.0;
    if(x.size()!=y.size()) {
        printf("Returned 0.0, because dimension of x doesn't match dimension of y.\n");
        return 0.0;
    }
    int dim=x.size();
    for (int i=0;i<dim;i++) {
        ip+=x[i]*y[i];
    }
    return ip;
}

double Operators::vectorNorm(const vector<double>& x) {
    return sqrt(innerProduct(x,x));
}

void Operators::MatrixVectorMultiplyer(Matrix& M, const vector<double>& x, vector<double>& b) {
    int dim=x.size();
    int n=sqrt(dim);
    for (int i=0;i<dim;i++) {
        b[i]+=x[i]*M.Get(i,i);
        if (i<(dim-n)) {
            b[i]+=x[i+n]*M.Get(i,i+n);
        } 
        if (i>=n) {
            b[i]+=x[i-n]*M.Get(i,i-n);
        }
        if (i%n!=0) {
            b[i]+=x[i-1]*M.Get(i,i-1);
            b[i-1]+=x[i]*M.Get(i-1,i);
        }
    }
}

double Operators::f(double x, double y) {
    return (-4.0);
}

double Operators::g(double x, double y) {
    return (x*x + y*y);
}