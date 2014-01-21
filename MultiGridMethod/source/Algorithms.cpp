#include "classes.h"

Algorithms::Algorithms(int n) {
    this->n=n;
    this->dim=n*n;
    this->h=1.0/(double)(n+1);
}

Algorithms::~Algorithms() {
}

double Algorithms::vectorNorm(const vector<double>& x) {
    double norm=0.0;
    for(int i=0;i<dim;i++) {
        norm+=pow(x[i],2);
    }
    return sqrt(norm);
}

void Algorithms::JacobiMethod(Matrix& A, vector<double>& x, const vector<double>& b, int steps) {
    int dim=A.Size();
    vector<double> r(dim);
    for(int j=0;j<steps;j++) {
        r=A*x;
        for(int i=0;i<dim;i++) {
            r[i]=1.0/4.0*(b[i]-r[i]);
            x[i]+=r[i];
        }
    }
}

void Algorithms::JacobiRelaxationMethod(Matrix& A, vector<double>& x, const vector<double>& b, int steps) {
    int dim=A.Size();
    vector<double> r(dim);
    for(int j=0;j<steps;j++) {
        r=A*x;
        for(int i=0;i<dim;i++) {
            r[i]=1.0/5.0*(b[i]-r[i]);
            x[i]+=r[i];
        }
    }
}

void Algorithms::JacobiRelaxationSolver(Matrix& A,vector<double>& x,const vector<double>& b) {
    int dim=A.Size(),steps=0;
    vector<double> r(dim,1),solved(dim);

    for(int i=1,k=0;i<(n+1);i++) {
        for(int j=1;j<(n+1);j++,k++) {
            solved[k]=g(j*h,i*h);
        }
    }
    double TOL=pow(10,-3)*vectorNorm(solved);
    while(TOL<vectorNorm(r)) {
        r=A*x;
        for(int i=0;i<dim;i++) {
            r[i]=1.0/5.0*(b[i]-r[i]);
            x[i]+=r[i];
            r[i]=x[i]-solved[i];
        }
        // x+=1.0/4.0*(b-A*x);
        steps++;
    }
    printf("JacobianSteps: %d\n", steps);
}

void Algorithms::GaussSeidelMethod(PoissonMatrix& A, vector<double>& x, const vector<double>& b, int steps) {
    double sum;
    int dim=A.Size();

    for(int j=0;j<steps;j++) {
        for(int i=0;i<dim;i++) {
            sum=0.0;
            if(i==0) {
                sum=x[i+1]*A.Get(i,i+1)+x[i+n]*A.Get(i,i+n);
            } else if(i!=0 && i<n) {
                sum=x[i-1]*A.Get(i,i-1)+x[i+1]*A.Get(i,i+1)+x[i+n]*A.Get(i,i+n);
            } else if(i>=n && i<dim-n) {
                sum=x[i-n]*A.Get(i,i-n)+x[i-1]*A.Get(i,i-1)+x[i+1]*A.Get(i,i+1)+x[i+n]*A.Get(i,i+n);
            } else if(i!=dim-1 && i>=dim-n) {
                sum=x[i-n]*A.Get(i,i-n)+x[i-1]*A.Get(i,i-1)+x[i+1]*A.Get(i,i+1);
            } else if(i==dim-1) {
                sum=x[i-n]*A.Get(i,i-n)+x[i-1]*A.Get(i,i-1);
            }
            x[i]=1.0/A.Get(i,i)*(b[i]-sum);
        }
    }
}

void Algorithms::SORMethod(PoissonMatrix& A, vector<double>& x, const vector<double>& b, int steps) {
    double sum,Pi=3.141592654,omega=2/(1+sqrt(1-pow(cos(Pi*h),2)));
    int dim=A.Size();

    for(int j=0;j<steps;j++) {
        for(int i=0;i<dim;i++) {
            sum=0.0;
            if(i==0) {
                sum=x[i+1]*A.Get(i,i+1)+x[i+n]*A.Get(i,i+n);
            } else if(i!=0 && i<n) {
                sum=x[i-1]*A.Get(i,i-1)+x[i+1]*A.Get(i,i+1)+x[i+n]*A.Get(i,i+n);
            } else if(i>=n && i<dim-n) {
                sum=x[i-n]*A.Get(i,i-n)+x[i-1]*A.Get(i,i-1)+x[i+1]*A.Get(i,i+1)+x[i+n]*A.Get(i,i+n);
            } else if(i!=dim-1 && i>=dim-n) {
                sum=x[i-n]*A.Get(i,i-n)+x[i-1]*A.Get(i,i-1)+x[i+1]*A.Get(i,i+1);
            } else if(i==dim-1) {
                sum=x[i-n]*A.Get(i,i-n)+x[i-1]*A.Get(i,i-1);
            }
            x[i]=omega*1.0/A.Get(i,i)*(b[i]-sum)+(1-omega)*x[i];
        }
    }
}

void Algorithms::Restriction(const vector<double>& r,vector<double>& r2h,int n) {
    for(int i=1,l=0,k=0;i<=n;i++) {
        for(int j=1;j<=n;j++,k++) {
            if(i%2==0 && j%2==0) {
                r2h[l]=1.0/16.0*(4.0*r[k]+2.0*(r[k-1]+r[k+1]+r[k-n]+r[k+n])+r[k+n-1]+r[k+n+1]+r[k-n-1]+r[k-n+1]);
                l++;
            }
        }
    }
}

void Algorithms::Interpolation(const vector<double>& E2h,vector<double>& E,Vectors& V,int n) {
    for(int i=1,k=0,l=0;i<=n;i++) {
        for(int j=1;j<=n;j++,k++) {
            if(i%2==0 && j%2==0) {
                E[k]=E2h[l];
                l++;
            }
        }
    }

    for(int i=1,k=0;i<=n;i++) {
        for(int j=1;j<=n;j++,k++) {
            if(i%2==0 && j%2!=0) {
                if(j!=1 && j!=n) {
                    E[k]=(double)1/(double)2*(E[k-1]+E[k+1]);
                }
                if(j==1) {
                    E[k]=(double)1/(double)2*(E[k+1]);
                }
                if(j==n) {
                    E[k]=(double)1/(double)2*(E[k-1]);
                }
            }
            if(i%2!=0 && j%2==0) {
                if(i!=1 && i!=n) {
                    E[k]=(double)1/(double)2*(E[k-n]+E[k+n]);
                }
                if(i==1) {
                    E[k]=(double)1/(double)2*(E[k+n]);
                }
                if(i==n) {
                    E[k]=(double)1/(double)2*(E[k-n]);
                }
            }
            if(i%2!=0 && j%2!=0) {
                if(i!=1 && j!=1 && i!=n && j!=n) {
                    E[k]=(double)1/(double)4*(E[k+n-1]+E[k+n+1]+E[k-n+1]+E[k-n-1]);
                }
                if(i==1 && j==1) {
                    E[k]=(double)1/(double)4*(E[k+n+1]);
                }
                if(i==n && j==n) {
                    E[k]=(double)1/(double)4*(E[k-n-1]);
                }
                if(i==1 && j==n) {
                    E[k]=(double)1/(double)4*(E[k-n+1]);
                }
                if(i==n && j==1) {
                    E[k]=(double)1/(double)4*(E[k+n-1]);
                }
                if(i==1 && j!=1 && j!=n) {
                    E[k]=(double)1/(double)4*(E[k+n+1]+E[k+n-1]);
                }
                if(i==n && j!=1 && j!=n) {
                    E[k]=(double)1/(double)4*(E[k-n-1]+E[k-n+1]);
                }
                if(j==1 && i!=1 && i!=n) {
                    E[k]=(double)1/(double)4*(E[k+n+1]+E[k-n+1]);
                }
                if(j==n && i!=1 && i!=n) {
                    E[k]=(double)1/(double)4*(E[k+n-1]+E[k-n-1]);
                }
            }

        }
    }
}

vector<double> Algorithms::MultiGridAlgorithm(PoissonMatrix& A,Vectors& V,const vector<double>& b,int n) {
    int dim=pow(n,2),N2h=(n+1)/2-1,dim2h=pow(N2h,2);
    vector<double> x(dim,0),r(dim,0),E(dim,0),r2h(dim2h,0),E2h(dim2h,0);
    if(n==this->n) {
        x=V.x;
    }
    if(n==1) {
        double a;
        a=4.0;
        x[0]=b[0]/a;
        return x;
    } else {
        JacobiRelaxationMethod(A,x,b,3);
        r=b-A*x;
        Restriction(r,r2h,n);
        A.Resize(N2h);
        E2h=MultiGridAlgorithm(A,V,r2h,N2h);
        E2h=MultiGridAlgorithm(A,V,r2h,N2h);
        Interpolation(E2h,E,V,n);
        x+=E;
        A.Resize(n);
        JacobiRelaxationMethod(A,x,b,3);
        return x;
    }
}

void Algorithms::MultiGridMethod(PoissonMatrix& A,Vectors& V,vector<double>& x,vector<double>& b) {
    int steps=0;
    vector<double> r(dim,1),solved(dim);

    for(int i=1,k=0;i<(n+1);i++) {
        for(int j=1;j<(n+1);j++,k++) {
            solved[k]=g(j*h,i*h);
        }
    }
    r=b-A*x;
    double TOL=pow(10,-3)*vectorNorm(r);
    while(TOL<=vectorNorm(r)) {
        V.x=MultiGridAlgorithm(A,V,b,n);
        r=b-A*x;
        // r=V.x-solved;
        steps++;
        if(steps==200) {
            r.assign(dim,0);
        }
    }
    printf("MultiGridSteps: %d\n", steps);
}