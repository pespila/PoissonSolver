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

void Algorithms::MatrixVectorMultiplyer(PoissonMatrix& A, const vector<double>& x, vector<double>& b) {
    int dim=A.Size();
    for(int i=0;i<dim;i++) {
        b[i]+=x[i]*A.Get(i,i);
        if(i<(dim-n)) {
            b[i]+=x[i+n]*A.Get(i,i+n);
        }
        if(i>=n) {
            b[i]+=x[i-n]*A.Get(i,i-n);
        }
        if(i%n!=0) {
            b[i]+=x[i-1]*A.Get(i,i-1);
            b[i-1]+=x[i]*A.Get(i-1,i);
        }
    }
}

void Algorithms::JacobiMethod(PoissonMatrix& A, vector<double>& x, const vector<double>& b, int steps) {
    double sum;
    int dim=A.Size();
    vector<double> tmp;
    tmp.resize(dim);

    for(int j=0;j<steps;j++) {
        tmp=x;
        for(int i=0;i<dim;i++) {
            sum=0.0;
            if(i==0) {
                sum=tmp[i+1]*A.Get(i,i+1)+tmp[i+n]*A.Get(i,i+n);
            } else if(i!=0 && i<n) {
                sum=tmp[i-1]*A.Get(i,i-1)+tmp[i+1]*A.Get(i,i+1)+tmp[i+n]*A.Get(i,i+n);
            } else if(i>=n && i<dim-n) {
                sum=tmp[i-n]*A.Get(i,i-n)+tmp[i-1]*A.Get(i,i-1)+tmp[i+1]*A.Get(i,i+1)+tmp[i+n]*A.Get(i,i+n);
            } else if(i!=dim-1 && i>=dim-n) {
                sum=tmp[i-n]*A.Get(i,i-n)+tmp[i-1]*A.Get(i,i-1)+tmp[i+1]*A.Get(i,i+1);
            } else if(i==dim-1) {
                sum=tmp[i-n]*A.Get(i,i-n)+tmp[i-1]*A.Get(i,i-1);
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
    int dim=n*n;
    vector<double> x;
    x.resize(dim);
    if(n==this->n) {
        x=V.x;
    }
    if(n==1) {
        double a;
        a=4.0*pow(n+1,2);
        x[0]=b[0]/a;
        return x;
    } else {
        int smallerN=(n+1)/2-1;
        SORMethod(A,x,b,4);
        vector<double> Ax,r;
        Ax.assign(dim,0);
        r.assign(dim,0);
        MatrixVectorMultiplyer(A,x,Ax);
        for(int i=0;i<dim;i++) {
            r[i]=b[i]-Ax[i];
        }
        vector<double> r2h;
        r2h.assign(smallerN*smallerN,0);
        Restriction(r,r2h,n);
        A.Resize(smallerN);
        vector<double> E2h;
        E2h.assign(smallerN*smallerN,0);
        E2h=MultiGridAlgorithm(A,V,r2h,smallerN);
        vector<double> E;
        E.assign(n*n,0);
        Interpolation(E2h,E,V,n);
        for(int i=0;i<dim;i++) {
            x[i]+=E[i];
        }
        A.Resize(n);
        SORMethod(A,x,b,4);
        return x;
    }
}

void Algorithms::MultiGridMethod(PoissonMatrix& A,Vectors& V) {
    int steps=0;
    vector<double> Ax,r;
    Ax.assign(dim,0);
    r.assign(dim,0);
    MatrixVectorMultiplyer(A,V.x,Ax);
    for(int i=0;i<dim;i++) {
        r[i]=V.b[i]-Ax[i];
    }
    double TOL=pow(10,-6);
    double eps=vectorNorm(r)*TOL;
    double norm=100.0;
    while(eps<=norm) {
        steps++;
        V.x=MultiGridAlgorithm(A,V,V.b,n);
        Ax.assign(dim,0);
        MatrixVectorMultiplyer(A,V.x,Ax);
        for(int i=0;i<dim;i++) {
            r[i]=V.b[i]-Ax[i];
        }
        norm=vectorNorm(r);
        printf("%f\n", norm);
        if(steps==200) {
            norm=0.0;
        }
    }
    printf("MultiGridSteps: %d\n", steps);
}