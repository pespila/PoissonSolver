#include "classes.h"

PoissonMatrix::PoissonMatrix(int n,double a,double b,double c) {
	this->n=n;
	this->dim=n*n;
	this->diagonal=a;
	this->tridiagonal=b;
	this->identity=c;
}

PoissonMatrix::~PoissonMatrix() {
}

int PoissonMatrix::Size() {
	return dim;
}

double PoissonMatrix::Get(int i,int j) {
	if(i==j) {
		return diagonal;
	} else if((j==(i+1) && i%n!=(n-1)) || (j==(i-1) && i%n!=0)) {
		return tridiagonal;
	} else if((j==(i+n)) || (j==(i-n))) {
		return identity;
	} else {
		return 0.0;
	}
}

void Matrix::PrintMatrix() {
	for(int i=0;i<Size();i++) {
		for(int j=0;j<Size();j++) {
			printf(" %.2f ",Get(i,j));
		}
		printf("\n");
	}
	printf("\n");
}

vector<double> PoissonMatrix::operator*(const vector<double>& x) {
    vector<double> tmp;
    tmp.assign(x.size(),0);
    for(int i=0;i<(int)x.size();i++) {
        tmp[i]+=x[i]*diagonal;
        if(i<(dim-n)) {
            tmp[i]+=x[i+n]*identity;
            tmp[i+n]+=x[i]*identity;
        }
        if(i%n!=0) {
            tmp[i]+=x[i-1]*tridiagonal;
            tmp[i-1]+=x[i]*tridiagonal;
        }
    }
    return tmp;
}