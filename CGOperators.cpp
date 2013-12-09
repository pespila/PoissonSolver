#include "classes.h"

double CGOperators::innerProduct(const vector<double>& x, const vector<double>& y) {
    if(x.size()!=y.size()) {
        printf("Returned 0.0, because dimension of x doesn't match dimension of y.\n");
        return 0.0;
    }
    double ip=0.0;
    int dim=x.size();
    for(int i=0;i<dim;i++) {
        ip+=x[i]*y[i];
    }
    return ip;
}

double CGOperators::vectorNorm(const vector<double>& x) {
    return sqrt(innerProduct(x,x));
}

void CGOperators::MatrixVectorMultiplyer(PoissonMatrix& A, const vector<double>& x, vector<double>& b) {
    if(A.Size()!=x.size() || A.Size()!=b.size() || x.size!=b.size()) {
        printf("Filled b with 0.0, because the dimensions of A,x and b don't match!\n");
        b.assign(b.size(),0);
    } else {
        int dim=x.size();
        int n=sqrt(dim);
        for(int i=0;i<dim;i++) {
            b[i]+=x[i]*A.Get(i,i);
            if(i<(dim-n)) {
                b[i]+=x[i+n]*A.Get(i,i+n);
            } 
            if(i>=n) {
                b[i]+=x[i-n]*A.Get(i,i-n);
            }
            if(i%n!=0) {
                b[i]+=x[i-1]*A.Get(i,i-1);
                b[i-1]+=x[i]*A.Get(i-1,i);
            }
        }
    }
}

double CGOperators::f(double x, double y) {
    return (-4.0);
}

double CGOperators::g(double x, double y) {
    return (x*x + y*y);
}