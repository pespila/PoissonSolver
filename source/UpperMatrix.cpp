#include "classes.h"

UpperMatrix::UpperMatrix(int n) : LowerMatrix(n) {
	this->n=n;
	this->dim=n*n;
	this->diagonal.assign(dim,0);
	this->tridiagonal.assign(dim-1,0);
	this->identity.assign(dim-n,0);
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