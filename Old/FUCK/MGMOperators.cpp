#include "classes.h"

MGMOperators::MGMOperators(){
}

MGMOperators::~MGMOperators(){
}

double MGMOperators::innerProduct(const vector<double>& x, const vector<double>& y) {
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

double MGMOperators::vectorNorm(const vector<double>& x) {
    return sqrt(innerProduct(x,x));
}

void MGMOperators::MatrixVectorMultiplyer(Matrix& A, const vector<double>& x, vector<double>& b) {
    if(A.Size()!=(int)x.size() || A.Size()!=(int)b.size() || x.size()!=b.size()) {
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

// double MGMOperators::innerProduct(const vector<vector<double> >& x, const vector<vector<double> >& y) {
//     if(x.size()!=y.size()) {
//         printf("Returned 0.0, because dimension of x doesn't match dimension of y.\n");
//         return 0.0;
//     }
//     double ip=0.0;
//     int n=x.size();
//     for(int i=0;i<n;i++) {
//         for(int j=0;j<n;j++) {
//             ip+=x[i][j]*y[i][j];
//         }
//     }
//     return ip;
// }

// double MGMOperators::vectorNorm(const vector<vector<double> >& x) {
//     return sqrt(innerProduct(x,x));
// }

// void MGMOperators::MatrixVectorMultiplyer(Matrix& A, const vector<vector<double> >& x, vector<vector<double> >& b) {
//     if(sqrt(A.Size())!=(int)x.size() || sqrt(A.Size())!=(int)b.size() || x.size()!=b.size()) {
//         printf("Filled b with 0.0, because the dimensions of A,x and b don't match!\n");
//         vector<double> pushB;
//         pushB.assign(b.size(),0);
//         b.assign(b.size(),pushB);
//         vector<double>().swap(pushB);
//     } else  {
//         int i,j,k;
//         int n=sqrt(x.size());
//         for(i=0,k=0;i<n;i++) {
//             for(j=0;j<n;j++,k++) {
//                 b[i][j]+=x[i][j]*A.Get(k,k);
//                 if(j!=n) {
//                     b[i][j]+=x[i][j+1]*A.Get(k,k+1);
//                     b[i][j+1]+=x[i][j]*A.Get(k,k-1);
//                 }
//                 if(i>0) {
//                     b[i][j]+=x[i][j-n]*A.Get(k,k-n);
//                 }
//                 if(i<n-1) {
//                     b[i][j]+=x[i][j+n]*A.Get(k,k+n);
//                 }
//             }
//         }
//     }
// }

double MGMOperators::f(double x,double y) {
    return (-4.0);
}

double MGMOperators::g(double x,double y) {
    return (x*x+y*y);
}