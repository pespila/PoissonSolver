#include "classes.h"

PoissonMatrix::PoissonMatrix(int n) {
	this->n=n;
	this->dim=n*n;
	diagonal=4.0*(n+1)*(n+1);
	tridiagonal=-1.0*(n+1)*(n+1);
	identity=tridiagonal;
}

PoissonMatrix::~PoissonMatrix() {
}

int PoissonMatrix::Size() {
	return dim;
}

double PoissonMatrix::Get(int i,int j) {
	if (i==j) {
		return diagonal;
	} else if (j==(i+1) && i%n!=(n-1)) {
		return tridiagonal;
	} else if (j==(i-1) && i%n!=0) {
		return tridiagonal;
	} else if (j==(i+n)) {
		return identity;
	} else if (j==(i-n)) {
		return identity;
	} else {
		return 0.0;
	}
}

void PoissonMatrix::Resize(int n) {
	this->n=n;
	this->dim=n*n;
}

void Matrix::PrintMatrix(int n) {
	if(n*n>=25){
		printf("Couldn't print Matrix! Dimension is too high.\n");
	} else {
		for(int i=0;i<n*n;i++) {
			for(int j=0;j<n*n;j++) {
				printf(" %.2f ", Get(i,j));
			}
			printf("\n");
		}
		printf("\n");
	}
}