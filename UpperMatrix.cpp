#include "classes.h"

UpperPoissonMatrix::UpperPoissonMatrix(int m) : PoissonMatrix(m) {
	n=m;
	dim=n*n;
	Udiagonal.resize(dim);
	Utridiagonal.resize(dim-1);
	Uidentity.resize(dim-n);
	for(int i=0;i<dim;i++) {
		Udiagonal[i] = 0.0;
		if(i < (dim-1)) {
			Utridiagonal[i] = 0.0;
		}
		if(i < (dim-n)) {
			Uidentity[i] = 0.0;
		}
	}
}

UpperPoissonMatrix::~UpperPoissonMatrix() {
}

void UpperPoissonMatrix::set(int i,int j,double value) {
	if (i == j) {
		Udiagonal[i] = value;
	} else if (j == (i+1)) {
		Utridiagonal[i] = value;
	} else if (j == (i+n)) {
		Uidentity[i] = value;
	}
}

double UpperPoissonMatrix::get(int i,int j) {
	if (i == j) {
		return Udiagonal[i];
	} else if (j == (i+1)) {
		return Utridiagonal[i];
	} else if (j == (i+n)) {
		return Uidentity[i];
	} else {
		return 0.0;
	}
}