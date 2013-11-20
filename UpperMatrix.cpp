#include "classes.h"

UpperMatrix::UpperMatrix(int m) : LowerMatrix(m) {
	n=m;
	dim=n*n;
	diagonal.assign(dim,0);
	tridiagonal.resize(dim-1);
	for(int i=0;i<dim-1;i++) {
		if(i%n!=n-1) {
			tridiagonal[i]=-1.0*pow(n,2);
		}
	}
	identity.assign(dim-n,-1.0*pow(n,2));
}

UpperMatrix::~UpperMatrix() {
	vector<double>().swap(diagonal);
	vector<double>().swap(tridiagonal);
	vector<double>().swap(identity);
}

double UpperMatrix::Get(int i,int j) {
	return LowerMatrix::Get(j,i);
}

void UpperMatrix::Set(int i,int j,double value) {
	return LowerMatrix::Set(j,i,value);
}