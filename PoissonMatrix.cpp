#include "classes.h"

PoissonMatrix::PoissonMatrix(int m) {
	n = m;
	dim = n*n;
	diagonal = 4.0*(n+1)*(n+1);
	tridiagonal = -1.0*(n+1)*(n+1);
	identity = tridiagonal;

	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			if(Get(i,j)!=0){
				hash.i=i;
				hash.j=j;
				HashKeys.push_back(hash);
			}
		}
	}
}

PoissonMatrix::~PoissonMatrix() {
}

int PoissonMatrix::Size() {
	return dim;
}

double PoissonMatrix::Get(int i, int j) {
	if (i == j) {
		return diagonal;
	} else if (j == (i+1) && i%n != (n-1)) {
		return tridiagonal;
	} else if (j == (i-1) && i%n != 0) {
		return tridiagonal;
	} else if (j == (i+n)) {
		return identity;
	} else if (j == (i-n)) {
		return identity;
	} else {
		return 0.0;
	}
}

void Matrix::PrintMatrix() {
	if(Size()>25){
		printf("Couldn't print Matrix! Dimension is too high.\n");
	} else {
		for(int i = 0; i < Size(); i++) {
			for(int j = 0; j < Size(); j++) {
				printf(" %.2f ", Get(i,j));
			}
			printf("\n");
		}
		printf("\n");
	}
}