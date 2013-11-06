#include "classes.h"

Operators::Operators(int m) {
    n=m;
    dim=n*n;
}

Operators::~Operators() {
}

double Operators::innerProduct(const vector<double>& x, const vector<double>& y) {
    double ip = 0.0;
    for (int i=0;i<dim;i++) {
        ip += x[i]*y[i];
    }
    return ip;
}

double Operators::vectorNorm(const vector<double>& x) {
    double ip, norm;
    ip=0.0;

    for (int i=0;i<dim;i++) {
        ip += x[i]*x[i];
    }
    norm = sqrt(ip);
    return norm;
}

void Operators::MatrixVectorMultiplyer(PoissonMatrix& Mat, const vector<double>& x, vector<double>& b) {
    for (int i = 0; i < dim; i++) {
        b[i] += x[i]*Mat.get(i,i);
        if (i < (dim-n)) {
            b[i] += x[i+n]*Mat.get(i,i+n);
        } 
        if (i >= n) {
            b[i] += x[i-n]*Mat.get(i,i-n);
        }
        if (i%n != 0) {
            b[i] += x[i-1]*Mat.get(i,i-1);
            b[i-1] += x[i]*Mat.get(i-1,i);
        }
    }
}

void Operators::LUsolverLower(LowerPoissonMatrix& L, vector<double>& z) {
    for(int i=0;i<dim;i++)
        for(int j=0;j<i;j++)
            z[i] -= L.get(i,j)*z[j];
}

void Operators::LUsolverUpper(UpperPoissonMatrix& U, vector<double>& z) {
    for(int i=dim-1;i>=0;i--) {
        for(int j=dim-1;j>=i+1;j--) {
            z[i] -= U.get(i,j)*z[j];
        }
        z[i] /= U.get(i,i);
    }
}

double Operators::f(double x, double y) {
    return (-4.0);
}

double Operators::g(double x, double y) {
    return (x*x + y*y);
}