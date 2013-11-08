#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <vector>
using namespace std;

class PoissonMatrix
{
    private:
    	int dim;
    	int n;
    	vector<double> diagonal;
    	vector<double> tridiagonal;
    	vector<double> identity;
    public:
    	PoissonMatrix(int);
    	~PoissonMatrix();

    	virtual void set(int, int, double);
    	virtual double get(int, int);
    	virtual int size();
    	void printMat();
};

class PoissonMatrixLite : public PoissonMatrix
{
    private:
        int dim;
        int n;
        double diag;
        double tridiag;
        double id;
    public:
        PoissonMatrixLite(int);
        ~PoissonMatrixLite();
};

class LowerPoissonMatrix : public PoissonMatrix
{
    private:
        int dim;
        int n;
        vector<double> Ldiagonal;
        vector<double> Ltridiagonal;
        vector<double> Lidentity;
    public:
        LowerPoissonMatrix(int);
        ~LowerPoissonMatrix();
        void set(int, int, double);
        double get(int, int);
};

class UpperPoissonMatrix : public PoissonMatrix
{
	private:
		int dim;
		int n;
		vector<double> Udiagonal;
		vector<double> Utridiagonal;
		vector<double> Uidentity;
	public:
		UpperPoissonMatrix(int);
		~UpperPoissonMatrix();
		void set(int, int, double);
		double get(int, int);
};

class Operators
{
	private:
		int dim;
		int n;
	public:
		Operators(int);
		~Operators();
		double innerProduct(const vector<double>&,const vector<double>&);
		double vectorNorm(const vector<double>&);
		void MatrixVectorMultiplyer(PoissonMatrix&,const vector<double>&,vector<double>& y);
		void LUsolverLower(LowerPoissonMatrix&,vector<double>&);
		void LUsolverUpper(UpperPoissonMatrix&,vector<double>&);
		double f(double, double);
    	double g(double, double);
};

class CGVectors {
	private:
		int dim;
		int n;
    public:
    	vector<double> x;
    	vector<double> b;
    	CGVectors(int, Operators&);
    	~CGVectors();
};

class Algorithms {
	private:
		int dim;
		int n;
	public:
		Algorithms(PoissonMatrix&);
		~Algorithms();
		void modifiedIncompleteLU(PoissonMatrix&, LowerPoissonMatrix&, UpperPoissonMatrix&);
		void modifiedIncompleteCholesky(PoissonMatrix&, LowerPoissonMatrix&, UpperPoissonMatrix&);
		void incompleteLU(PoissonMatrix&, LowerPoissonMatrix&, UpperPoissonMatrix&);
		void incompleteCholesky(PoissonMatrix&, LowerPoissonMatrix&, UpperPoissonMatrix&);
		void CG(PoissonMatrix&, Operators&, vector<double>&, const vector<double>&);
		void PCG(PoissonMatrix&, Operators&, LowerPoissonMatrix&, UpperPoissonMatrix&, vector<double>&, const vector<double>&);
};

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

UpperPoissonMatrix::UpperPoissonMatrix(int m) : PoissonMatrix(m) {
	n=m;
	dim=n*n;
	Udiagonal.resize(dim);
	Utridiagonal.resize(dim-1);
	Uidentity.resize(dim-n);
	for(int i=0;i<dim;i++) {
		Udiagonal[i] = 0.0;
		if(i < (dim-1)) {
			Utridiagonal[i] = 0.0;
		}
		if(i < (dim-n)) {
			Uidentity[i] = 0.0;
		}
	}
}

UpperPoissonMatrix::~UpperPoissonMatrix() {
}

void UpperPoissonMatrix::set(int i,int j,double value) {
	if (i == j) {
		Udiagonal[i] = value;
	} else if (j == (i+1)) {
		Utridiagonal[i] = value;
	} else if (j == (i+n)) {
		Uidentity[i] = value;
	}
}

double UpperPoissonMatrix::get(int i,int j) {
	if (i == j) {
		return Udiagonal[i];
	} else if (j == (i+1)) {
		return Utridiagonal[i];
	} else if (j == (i+n)) {
		return Uidentity[i];
	} else {
		return 0.0;
	}
}

CGVectors::CGVectors(int m, Operators& O) {
	n=m;
	dim=n*n;
    b.resize(dim);
    x.resize(dim);

	int k,i,j;
    double h=1.0/(double)(n+1);

    k=0;
    for (j=1;j<=n;j++) {
        if (j == 1) {
            b[k] = O.f(h, h) + (n+1)*(n+1) * (O.g(0, h) + O.g(h, 0));
        } else if (j == n) {
            b[k] = O.f(1-h, 1-h) + (n+1)*(n+1) * (O.g(1-h, 0) + O.g(1, h));
        } else {
            b[k] = O.f(j*h, h) + (n+1)*(n+1) * O.g(j*h, 0);
        }
        k++;
    }
    
    for (i=2; i<=(n-1);i++) {
        for (j = 1; j <= n; j++) {
            if (j == 1) {
				b[k] = O.f(h, i*h) + (n+1)*(n+1) * O.g(0, i*h);            
            } else if (j == n) {
                b[k] = O.f(1-h, i*h) + (n+1)*(n+1) * O.g(1, i*h);
            } else {
                b[k] = O.f(j*h, i*h);
            }
            k++;
        }
    }
    
    for (j=1;j<=n;j++) {
    	  if (j == 1) {
           	b[k] = O.f(h, 1-h) + (n+1)*(n+1) * (O.g(h, 1)+O.g(0, 1-h));
        } else if (j == n) {
            b[k] = O.f(1-h, 1-h) + (n+1)*(n+1) * (O.g(1-h, 1)+O.g(1, 1-h));
        } else {
            b[k] = O.f(j*h, 1-h) + (n+1)*(n+1) * O.g(j*h, 1);
        }
        k++;
    }
}

CGVectors::~CGVectors() {
}

Operators::Operators(int m) {
	n=m;
	dim=n*n;
}

Operators::~Operators() {
}

double Operators::innerProduct(const vector<double>& x, const vector<double>& y) {
    double ip = 0.0;
    for (int i=0;i<dim;i++) {
        ip += x[i]*y[i];
    }
    return ip;
}

double Operators::vectorNorm(const vector<double>& x) {
	double ip, norm;
    ip=0.0;

    for (int i=0;i<dim;i++) {
        ip += x[i]*x[i];
    }
    norm = sqrt(ip);
    return norm;
}
//---------------------------------------------------------------------------------------------------------------------------
void Operators::MatrixVectorMultiplyer(PoissonMatrix& Mat, const vector<double>& x, vector<double>& b) {
    for (int i = 0; i < dim; i++) {
        b[i] += x[i]*Mat.get(i,i);
        if (i < (dim-n)) {
            b[i] += x[i+n]*Mat.get(i,i+n);
        } 
        if (i >= n) {
            b[i] += x[i-n]*Mat.get(i,i-n);
        }
        if (i%n != 0) {
            b[i] += x[i-1]*Mat.get(i,i-1);
            b[i-1] += x[i]*Mat.get(i-1,i);
        }
    }
}
//---------------------------------------------------------------------------------------------------------------------------
void Operators::LUsolverLower(LowerPoissonMatrix& L, vector<double>& z) {
    for(int i=0;i<dim;i++)
        for(int j=0;j<i;j++)
            z[i] -= L.get(i,j)*z[j];
}

void Operators::LUsolverUpper(UpperPoissonMatrix& U, vector<double>& z) {
    for(int i=dim-1;i>=0;i--) {
        for(int j=dim-1;j>=i+1;j--) {
            z[i] -= U.get(i,j)*z[j];
        }
        z[i] /= U.get(i,i);
    }
}

double Operators::f(double x, double y) {
	return (-4.0);
}

double Operators::g(double x, double y) {
	return (x*x + y*y);
}

Algorithms::Algorithms(PoissonMatrix& Mat) {
	dim = Mat.size();
	n = sqrt(dim);
}

Algorithms::~Algorithms() {
}

void Algorithms::modifiedIncompleteLU(PoissonMatrix& A, LowerPoissonMatrix& L, UpperPoissonMatrix& U) {
    int i,j,k;
    double s, tmp, drop;

    for(i=0;i<dim;i++) {
        drop = 0;
        for(k=0;k<i;k++) {
            s = 0;
            for(j=0;j<k;j++) {
                s += L.get(i,j) * U.get(j,k);
            }
            if(A.get(i,k) != 0) {
                tmp = (A.get(i,k) - s)/U.get(k,k);
                L.set(i,k,tmp);
            } else {
                drop += s;
            }
        }
        for(k=i;k<dim;k++) {
            s = 0;
            for(j=0;j<i;j++) {
                s += L.get(i,j) * U.get(j,k);
            }
            if(A.get(i,k) != 0) {
                tmp = A.get(i,k) - s;
                U.set(i,k,tmp);
            } else {
                drop += s;
            }
        }
        U.set(i,i,(U.get(i,i) - drop));
    }
}

void Algorithms::modifiedIncompleteCholesky(PoissonMatrix& A, LowerPoissonMatrix& L, UpperPoissonMatrix& Ltranspose) {
    int i,j,k;
    double sum, tmp;

    for(k=0;k<dim;k++) {
        sum=0;
        for(j=0;j<k;j++) {
            if (A.get(k,j) != 0) {
                sum+=L.get(k,j);
            }
        }
        tmp = sqrt(A.get(k,k)-sum);
        L.set(k,k,tmp);
        Ltranspose.set(k,k,tmp); //neu
        for(i=k+1;i<dim;i++) {
            if (A.get(i,k) != 0) {
                sum=0;
                for(j=0;j<k;j++) {
                    if ((A.get(i,j) != 0) && (A.get(k,j) != 0)) {
                        sum+=L.get(i,j)*L.get(k,j);
                    }
                }
                tmp = (A.get(i,k) - sum)/L.get(k,k);
                L.set(i,k,tmp);
                Ltranspose.set(k,i,tmp); //neu
            } else {
                sum=0;
                for(j=0;j<k;j++) {
                    if ((A.get(i,j) != 0) && (A.get(k,j) != 0)) {
                        sum+=L.get(i,j)*L.get(k,j);
                    }
                }
                tmp = A.get(i,i) - sum;
                A.set(i,i,tmp);
            }
        }
    }
}

void Algorithms::incompleteLU(PoissonMatrix& A, LowerPoissonMatrix& L, UpperPoissonMatrix& U) {
    int i,j,k;
    double s, tmp;

    for(i=0;i<dim;i++) {
        for(k=1;k<i;k++) {
            if (A.get(i,k) != 0) {
                s=0;
                for(j=0;j<k;j++) {
                    s+=L.get(i,j)*U.get(j,k);
                }
                tmp = (A.get(i,k)-s)/U.get(k,k);
                L.set(i,k,tmp);
            }
        }
        for(k=i;k<dim;k++) {
            if (A.get(i,k) != 0) {
                s=0;
                for(j=0;j<i;j++) {
                    s+=L.get(i,j)*U.get(j,k);
                }
                tmp = A.get(i,k)-s;
                U.set(i,k,tmp);
            }
        }
    }
}

void Algorithms::incompleteCholesky(PoissonMatrix& A, LowerPoissonMatrix& L, UpperPoissonMatrix& Ltranspose) {
    int i,j,k;
    double sum, tmp;

    for(k=0;k<dim;k++) {
        sum=0;
        for(j=0;j<k;j++) {
            if(A.get(k,j) != 0) {
                sum+=L.get(k,j)*L.get(k,j);
            }
        }
        tmp = sqrt(A.get(k,k)-sum);
        L.set(k,k,tmp);
        Ltranspose.set(k,k,tmp); //neu
        for(i=k+1;i<dim;i++) {
            if (A.get(i,k) != 0) {
                sum=0;
                for(j=0;j<k;j++) {
                    if(A.get(i,j) != 0 && A.get(k,j) != 0) {
                        sum+=L.get(i,j)*L.get(k,j);
                    }
                }
                tmp = (A.get(i,k) - sum)/L.get(k,k);
                L.set(i,k,tmp);
                Ltranspose.set(k,i,tmp);
            }
        }
    }
}

void Algorithms::CG(PoissonMatrix& A, Operators& O, vector<double>& x, const vector<double>& b) {
    double alpha,beta,eps,h,norm,skpOfRes,denom;
    int steps,i,j,k;

    steps=0;
    h=1.0/(double)(n+1);
    beta=0;

    vector<double> r;
    vector<double> solution;
    vector<double> stopNorm;
    vector<double> p;
    vector<double> Ap;
    vector<double> rTmp;
    r.resize(dim);
    solution.resize(dim);
    stopNorm.resize(dim);
    p.resize(dim);
    rTmp.resize(dim);

    vector<double> Ax;
    Ax.assign(dim,0);
    x.assign(dim,0);
    O.MatrixVectorMultiplyer(A,x,Ax);
    for(i=0;i<dim;i++)
        r[i]=b[i]-Ax[i];
    Ax.clear();

    k=0;
    for (i=1;i<(n+1);i++) {
        for (j=1;j<(n+1);j++) {
            solution[k] = O.g(j*h,i*h);
            stopNorm[k] = x[k]-solution[k];
            k++;
        }
    }

    eps=pow(10,-3)*O.vectorNorm(stopNorm);
    norm=1.0;

    p=r;
    rTmp=r;

    while (eps <= norm) {
        steps++;

        Ap.assign(dim,0);
        O.MatrixVectorMultiplyer(A,p,Ap);

        skpOfRes=O.innerProduct(rTmp,rTmp);
        denom = O.innerProduct(p,Ap);
        alpha=skpOfRes/denom;

        for(i=0;i<dim;i++) {
            x[i] += alpha*p[i];
            r[i] -= alpha*Ap[i];
        }

        denom = O.innerProduct(r,r);
        beta = denom/skpOfRes;

        for(i=0;i<dim;i++) {
            p[i]=r[i]+beta*p[i];
            stopNorm[i]=x[i]-solution[i];
        }
        rTmp = r;
        Ap.clear();
        norm=O.vectorNorm(stopNorm);
    }
    printf("CGSteps: %d\n", steps);
}


//Ax=b -> Wz=r
//Ly=b -> Ly=r
//Rx=y -> Rz=y
void Algorithms::PCG(PoissonMatrix& A, Operators& O, LowerPoissonMatrix& L, UpperPoissonMatrix& U, vector<double>& x, const vector<double>& b) {
    double alpha,beta,eps,h,norm,skpOfRes,denom;
    int steps,i,j,k;

    steps=0;
    h=1.0/(double)(n+1);

    vector<double> r;
    vector<double> z;
    vector<double> solution;
    vector<double> stopNorm;
    vector<double> p;
    vector<double> Ap;
    vector<double> rTmp;
    vector<double> zTmp;
    r.resize(dim);
    z.resize(dim);
    solution.resize(dim);
    stopNorm.resize(dim);
    p.resize(dim);
    rTmp.resize(dim);
    zTmp.resize(dim);

    vector<double> Ax;
    Ax.assign(dim,0);
    x.assign(dim,0);
    O.MatrixVectorMultiplyer(A,x,Ax);
    for(i=0;i<dim;i++)
        r[i]=b[i]-Ax[i];

    z=r;

    O.LUsolverLower(L,z);
    O.LUsolverUpper(U,z);

    p=z;
    rTmp=r;
    zTmp=z;

    k=0;
    for (i=1;i<(n+1);i++) {
        for (j=1;j<(n+1);j++) {
            solution[k] = O.g(j*h,i*h);
            stopNorm[k] = x[k]-solution[k];
            k++;
        }
    }

    eps=pow(10,-3)*O.vectorNorm(stopNorm);
    norm=1.0;

    while (eps <= norm) {
        steps++;

        Ap.assign(dim,0);
        O.MatrixVectorMultiplyer(A,p,Ap);

        skpOfRes=O.innerProduct(zTmp,rTmp);
        denom = O.innerProduct(p,Ap);
        alpha=skpOfRes/denom;

        for(i=0;i<dim;i++) {
            x[i] += alpha*p[i];
            r[i] -= alpha*Ap[i];
        }

        z=r;
        O.LUsolverLower(L,z);
        O.LUsolverUpper(U,z);
        zTmp=z;

        denom = O.innerProduct(z,r);
        beta = denom/skpOfRes;

        for(i=0;i<dim;i++) {
            p[i]=z[i]+beta*p[i];
            stopNorm[i]=x[i]-solution[i];
        }
        rTmp = r;
        Ap.clear();
        norm=O.vectorNorm(stopNorm);
        if (steps > 100)
            break;
    }
    printf("PCGSteps: %d\n", steps);
}

int main(int argc, char const *argv[]) {
	int arg, n, dim;
	if (argc > 1) {
        arg = atoi(&*argv[1]);
    } else {
        arg = 4;
    }
    n=arg-1;
	dim = n*n;

    PoissonMatrix A(n);
    Operators O(n);
    CGVectors V(n,O);
    Algorithms Run(A);

    LowerPoissonMatrix L(n);
    UpperPoissonMatrix U(n);

    //Run.CG(A,O,V.x,V.b);

    //Run.modifiedIncompleteLU(A,L,U);
    Run.incompleteLU(A,L,U);
    printf("done\n");
    // A.printMat();
    //L.printMat();
    //U.printMat();

    // double matrix[dim][dim];

    // double tmp;
    // for(int i=0;i<dim;i++){
    //     for(int j=0;j<dim;j++){
    //         tmp=0;
    //         for(int k=0;k<dim;k++){
    //             tmp += L.get(i,k)*U.get(k,j);
    //         }
    //     matrix[i][j] = tmp;
    //     }
    // }

    // for(int i=0;i<dim;i++){
    //     for(int j=0;j<dim;j++){
    //         printf("%f ", matrix[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    //Run.CG(L,O,V.x,V.b);

    Run.PCG(A,O,L,U,V.x,V.b);
    
    // for(int i=0;i<dim;i++)
    //     printf("%f ", V.x[i]);
    // printf("\n");

	return EXIT_SUCCESS;
}