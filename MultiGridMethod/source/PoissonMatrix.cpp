#include "classes.h"

PoissonMatrix::PoissonMatrix(int n) {
	this->n=n;
	this->dim=n*n;
	this->diagonal=4.0*pow(n+1,2);
	this->tridiagonal=-1.0*pow(n+1,2);
	this->identity=tridiagonal;
}

PoissonMatrix::~PoissonMatrix() {
}

int PoissonMatrix::Size() {
	return dim;
}

double PoissonMatrix::Get(int i,int j) {
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

void PoissonMatrix::Resize(int n) {
    this->n=n;
    this->dim=n*n;
    this->diagonal=4.0*pow(n+1,2);
    this->tridiagonal=-1.0*pow(n+1,2);
    this->identity=tridiagonal;
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