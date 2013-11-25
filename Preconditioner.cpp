#include "classes.h"

Preconditioner::Preconditioner(int m) {
	n=m;
	dim=n*n;
	diagonal.assign(dim,4);
	tridiagonal.assign(dim-n,-1);
	identity.assign(dim-n,-1);

	int i,j,k;
	for(i=0;i<dim;i++) {
		vector<int> push;
		push.assign(5,-1);
		k=0;
		for(j=0;j<dim;j++) {
			if(Get(i,j)!=0) {
				if(k==0) push[k]=(j);
				if(k==1) push[k]=(j);
				if(k==2) push[k]=(j);
				if(k==3) push[k]=(j);
				if(k==4) push[k]=(j);
				k++;
			}
		}
		HashMatrix.push_back( push );
	}
}

Preconditioner::~Preconditioner() {
	vector<double>().swap(diagonal);
	vector<double>().swap(tridiagonal);
	vector<double>().swap(identity);
}

int Preconditioner::Size() {
	return dim;
}

void Preconditioner::Set(int i,int j,double value) {
	if (i == j) {
		diagonal[i] = value;
	} else if (j == (i+1) && i%n != (n-1)) {
		tridiagonal[i] = value;
	} else if (j == (i-1) && i%n != 0) {
		tridiagonal[j] = value;
	} else if (j == (i+n)) {
		identity[i] = value;
	} else if (j == (i-n)) {
		identity[j] = value;
	}
}

double Preconditioner::Get(int i,int j) {
	if (i == j) {
		return diagonal[i];
	} else if (j == (i+1) && i%n != (n-1)) {
		return tridiagonal[i];
	} else if (j == (i-1) && i%n != 0) {
		return tridiagonal[j];
	} else if (j == (i+n)) {
		return identity[i];
	} else if (j == (i-n)) {
		return identity[j];
	} else {
		return 0.0;
	}
}