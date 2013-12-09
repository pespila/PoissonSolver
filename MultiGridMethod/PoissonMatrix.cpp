#include "classes.h"

PoissonMatrix::PoissonMatrix(int m) {
	n=m;
	dim=n*n;
	diagonal=4.0*(n+1)*(n+1);
	tridiagonal=-1.0*(n+1)*(n+1);
	identity=tridiagonal;

    int i;
    vector<int> push;
    push.assign(5,-1);
    for(i=0;i<dim;i++) {
        if(i<n) {
            if(i==0) {
                push[0]=i;
                push[1]=i+1;
                push[2]=i+n;
            } else if(i%n!=n-1 && i!=0) {
                push[0]=i-1;
                push[1]=i;
                push[2]=i+1;
                push[3]=i+n;
            } else {
                push[0]=i-1;
                push[1]=i;
                push[2]=i+n;
            }
        } else if(i>=n && i<=dim-n) {
            if(i%n==0) {
                push[0]=i-n;
                push[1]=i;
                push[2]=i+1;
                push[3]=i+n;
                push[4]=-1;
            } else if(i%n==n-1) {
                push[0]=i-n;
                push[1]=i-1;
                push[2]=i;
                push[3]=i+n;
                push[4]=-1;
            } else {
                push[0]=i-n;
                push[1]=i-1;
                push[2]=i;
                push[3]=i+1;
                push[4]=i+n;
            }
        } else if(i>dim-n) {
            if(i==dim-1) {
                push[0]=i-n;
                push[1]=i-1;
                push[2]=i;
                push[3]=-1;
                push[4]=-1;
            } else if(i%n!=0 && i!=dim-1) {
                push[0]=i-n;
                push[1]=i-1;
                push[2]=i;
                push[3]=i+1;
                push[4]=-1;
            } else {
                push[0]=i-n;
                push[1]=i;
                push[2]=i+1;
                push[3]=-1;
                push[4]=-1;
            }
        }
        HashMatrix.push_back(push);
    }
}

PoissonMatrix::~PoissonMatrix() {
}

int PoissonMatrix::Size() {
	return dim;
}

double PoissonMatrix::Get(int i,int j) {
	if(i==j) {
		return diagonal;
	} else if(j==(i+1) && i%n!=(n-1)) {
		return tridiagonal;
	} else if(j==(i-1) && i%n!=0) {
		return tridiagonal;
	} else if(j==(i+n)) {
		return identity;
	} else if(j==(i-n)) {
		return identity;
	} else {
		return 0.0;
	}
}

void PoissonMatrix::Resize(int m) {
    n=m;
    dim=n*n;
}

void Matrix::PrintMatrix() {
	if(Size()>25){
		printf("Couldn't print Matrix! Dimension is too high.\n");
	} else {
		for(int i=0;i<Size();i++) {
			for(int j=0;j<Size();j++) {
				printf(" %.2f ",Get(i,j));
			}
			printf("\n");
		}
		printf("\n");
	}
}