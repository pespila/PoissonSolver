#include "classes.h"

Matrix::Matrix(int n,double a,double b,double c) {
	this->n=n;
	this->dim=pow(n,2);
	this->diagonal=a*pow(1.0/(n+1),2);
	this->tridiagonal=b*pow(1.0/(n+1),2);
	this->identity=c*pow(1.0/(n+1),2);
}

Matrix::~Matrix() {}

void Matrix::PrintMatrix() {
	for(int i=0;i<dim;i++) {
		for(int j=0;j<dim;j++) {
			if(i==j) printf(" %.2f ", diagonal);
			else if((j==(i+1) && i%n!=(n-1)) || (j==(i-1) && i%n!=0)) printf(" %.2f ", tridiagonal);
			else if((j==(i+n)) || (j==(i-n))) printf(" %.2f ", identity);
			else printf(" %.2f ", 0.0);
		}
		printf("\n");
	}
	printf("\n");
}

vector<double> Matrix::operator*(const vector<double>& x) {
	int dim=x.size(),n=sqrt(dim);
    vector<double> tmp(dim,0);
    for(int i=0;i<dim;i++) {
        tmp[i]+=x[i]*4.0*pow(n+1,2);
        if(i<(dim-n)) {
            tmp[i]+=x[i+n]*-1.0*pow(n+1,2);
            tmp[i+n]+=x[i]*-1.0*pow(n+1,2);
        }
        if(i%n!=0) {
            tmp[i]+=x[i-1]*-1.0*pow(n+1,2);
            tmp[i-1]+=x[i]*-1.0*pow(n+1,2);
        }
    }
    return tmp;
}