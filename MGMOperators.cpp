#include "classes.h"

double MGMOperators::innerProduct(const vector<double>& x, const vector<double>& y) {
    if(x.size()!=y.size()) {
        printf("Returned 0.0, because dimension of x doesn't match dimension of y.\n");
        return 0.0;
    }
    double ip=0.0;
    int n=x.size();
    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) {
            ip+=x[i][j]*y[i][j];
        }
    }
    return ip;
}

double MGMOperators::vectorNorm(const vector<double>& x) {
    return sqrt(innerProduct(x,x));
}

void MGMOperators::MatrixVectorMultiplyer(PoissonMatrix& A, const vector<double>& x, vector<double>& b) {
    if(sqrt(A.Size())!=x.size() || sqrt(A.Size())!=b.size() || x.size!=b.size()) {
        printf("Filled b with 0.0, because the dimensions of A,x and b don't match!\n");
        vector<double> pushB;
        pushB.assign(b.size(),0);
        b.assign(b.size(),pushB);
        vector<double>().swap(pushB);
    } else  {
        int i,j,k;
        int n=sqrt(x.size());
        for(i=0,k=0;i<n;i++) {
            for(j=0;j<n;j++,k++) {
                b[i][j]+=x[i][j]*A.Get(k,k);
                if(j!=n) {
                    b[i][j]+=x[i][j+1]*A.Get(k,k+1);
                    b[i][j+1]+=x[i][j]*A.Get(k,k-1);
                }
                if(i>0) {
                    b[i][j]+=x[i][j-n]*A.Get(k,k-n);
                }
                if(i<n-1) {
                    b[i][j]+=x[i][j+n]*A.Get(k,k+n);
                }
            }
        }
    }
}

double MGMOperators::f(double x, double y) {
    return (-4.0);
}

double MGMOperators::g(double x, double y) {
    return (x*x + y*y);
}