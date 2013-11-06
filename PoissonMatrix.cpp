#include "classes.h"

PoissonMatrix::PoissonMatrix(int m) {
	n = m;
	dim = n*n;
	diagonal.resize(dim);
	tridiagonal.resize(dim-1);
	identity.resize(dim-n);

	for(int i=0;i<dim;i++) {
		diagonal[i] = 4.0*(n+1)*(n+1);
		if(i < (dim-1)) {
			if(i%n == (n-1))
				tridiagonal[i] = 0.0*(n+1)*(n+1);
			else
				tridiagonal[i] = -1.0*(n+1)*(n+1);
		}
		if(i < (dim-n)) {
			identity[i] = -1.0*(n+1)*(n+1);
		}
	}
}

PoissonMatrix::~PoissonMatrix() {
}

void PoissonMatrix::set(int i,int j,double value) {
	if (i == j) {
		diagonal[i] = value;
	} else if (j == (i+1)) {
		tridiagonal[i] = value;
	} else if (j == (i-1)) {
		tridiagonal[j] = value;
	} else if (j == (i+n)) {
		identity[i] = value;
	} else if (j == (i-n)) {
		identity[j] = value;
	}
}

double PoissonMatrix::get(int i, int j) {
	if (i == j) {
		return diagonal[i];
	} else if (j == (i+1)) {
		return tridiagonal[i];
	} else if (j == (i-1)) {
		return tridiagonal[j];
	} else if (j == (i+n)) {
		return identity[i];
	} else if (j == (i-n)) {
		return identity[j];
	} else {
		return 0.0;
	}
}

int PoissonMatrix::size() {
	return dim;
}

void PoissonMatrix::printMat() {
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			printf(" %.2f ", get(i,j));
		}
		printf("\n");
	}
	printf("\n");
}

PoissonMatrixLite::PoissonMatrixLite(int m) : PoissonMatrix(m) {
    n=m;
    dim=n*n;
    diag = 4*pow((n+1),2);
    tridiag = -1*pow((n+1),2);
    id = -1*pow((n+1),2);
}

PoissonMatrixLite::~PoissonMatrixLite(){
}
