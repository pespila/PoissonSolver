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

vector<double> UpperMatrix::operator*(const vector<double>& x) {
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