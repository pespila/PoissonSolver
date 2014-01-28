#include "classes.h"

PoissonMatrix1D::PoissonMatrix1D(int n,double a,double b) {
	this->n=n;
	this->dim=n;
	this->diagonal=a;
	this->tridiagonal=b;
}

PoissonMatrix1D::~PoissonMatrix1D() {
}

int PoissonMatrix1D::Size() {
	return dim;
}

double PoissonMatrix1D::Get(int i,int j) {
	if(i==j) {
		return diagonal;
	} else if(j==(i+1)) {
		return tridiagonal;
	} else if(j==(i-1)) {
		return tridiagonal;
	} else {
		return 0.0;
	}
}

void PoissonMatrix1D::Resize(int n) {
    this->n=n;
    this->dim=n;
}

PoissonMatrix2D::PoissonMatrix2D(int n,double a,double b,double c) {
	this->n=n;
	this->dim=pow(n,2);
	this->diagonal=a;
	this->tridiagonal=b;
	this->identity=c;
}

PoissonMatrix2D::~PoissonMatrix2D() {
}

int PoissonMatrix2D::Size() {
	return dim;
}

double PoissonMatrix2D::Get(int i,int j) {
	if(i==j) {
		return diagonal;
	} else if(j==(i+1) && i%n!=(n-1)) {
		return tridiagonal;
	} else if(j==(i-1) && i%n!=0) {
		return tridiagonal;
	} else if(j==(i+n)) {
		return identity;
	} else if(j==(i-n)) {
		return identity;
	} else {
		return 0.0;
	}
}

void PoissonMatrix2D::Resize(int n) {
    this->n=n;
    this->dim=n*n;
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

vector<double> PoissonMatrix2D::operator*(const vector<double>& x) {
    vector<double> tmp;
    tmp.assign(x.size(),0);
    for(int i=0;i<dim;i++) {
        tmp[i]+=x[i]*4.0;
        if(i<(dim-n)) {
            tmp[i]+=x[i+n]*-1.0;
            tmp[i+n]+=x[i]*-1.0;
        }
        if(i%n!=0) {
            tmp[i]+=x[i-1]*-1.0;
            tmp[i-1]+=x[i]*-1.0;
        }
    }
    return tmp;
}

vector<double> PoissonMatrix1D::operator*(const vector<double>& x) {
    vector<double> tmp;
    tmp.assign(x.size(),0);
    for(int i=0;i<dim;i++) {
        tmp[i]+=x[i]*diagonal+x[i+1]*tridiagonal+x[i-1]*tridiagonal;
    }
    return tmp;
}