#include "classes.h"

UpperMatrix::UpperMatrix(int m) : LowerMatrix(m) {
	n=m;
	dim=n*n;
	diagonal.assign(dim,0);
	tridiagonal.assign(dim-1,0);
	identity.assign(dim-n,0);
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