#include "classes.h"

LowerPoissonMatrix::LowerPoissonMatrix(int m) : PoissonMatrix(m) {
	n=m;
	dim=n*n;
	Ldiagonal.resize(dim);
	Ltridiagonal.resize(dim-1);
	Lidentity.resize(dim-n);
	for(int i=0;i<dim;i++) {
		Ldiagonal[i] = 1.0;
		if(i < (dim-1)) {
			Ltridiagonal[i] = 0.0;
		}
		if(i < (dim-n)) {
			Lidentity[i] = 0.0;
		}
	}
}

LowerPoissonMatrix::~LowerPoissonMatrix() {
}

void LowerPoissonMatrix::set(int i,int j,double value) {
	if (i == j) {
		Ldiagonal[i] = value;
	} else if (j == (i-1)) {
		Ltridiagonal[j] = value;
	} else if (j == (i-n)) {
		Lidentity[j] = value;
	}
}

double LowerPoissonMatrix::get(int i,int j) {
	if (i == j) {
		return Ldiagonal[i];
	} else if (j == (i-1)) {
		return Ltridiagonal[j];
	} else if (j == (i-n)) {
		return Lidentity[j];
	} else {
		return 0.0;
	}
}