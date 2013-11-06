#include "classes.h"

CGVectors::CGVectors(int m, Operators& O) {
    n=m;
    dim=n*n;
    b.resize(dim);
    x.resize(dim);

    int k,i,j;
    double h=1.0/(double)(n+1);

    k=0;
    for (j=1;j<=n;j++) {
        if (j == 1) {
            b[k] = O.f(h, h) + (n+1)*(n+1) * (O.g(0, h) + O.g(h, 0));
        } else if (j == n) {
            b[k] = O.f(1-h, 1-h) + (n+1)*(n+1) * (O.g(1-h, 0) + O.g(1, h));
        } else {
            b[k] = O.f(j*h, h) + (n+1)*(n+1) * O.g(j*h, 0);
        }
        k++;
    }
    
    for (i=2; i<=(n-1);i++) {
        for (j = 1; j <= n; j++) {
            if (j == 1) {
                b[k] = O.f(h, i*h) + (n+1)*(n+1) * O.g(0, i*h);            
            } else if (j == n) {
                b[k] = O.f(1-h, i*h) + (n+1)*(n+1) * O.g(1, i*h);
            } else {
                b[k] = O.f(j*h, i*h);
            }
            k++;
        }
    }
    
    for (j=1;j<=n;j++) {
          if (j == 1) {
            b[k] = O.f(h, 1-h) + (n+1)*(n+1) * (O.g(h, 1)+O.g(0, 1-h));
        } else if (j == n) {
            b[k] = O.f(1-h, 1-h) + (n+1)*(n+1) * (O.g(1-h, 1)+O.g(1, 1-h));
        } else {
            b[k] = O.f(j*h, 1-h) + (n+1)*(n+1) * O.g(j*h, 1);
        }
        k++;
    }
}

CGVectors::~CGVectors() {
}