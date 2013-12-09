#include "classes.h"

Algorithms::Algorithms(int m) {
    n=m;
    dim=n*n;
    h=1.0/(double)(n+1);
}

Algorithms::~Algorithms() {
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
    int i,j,k,m;
    for(j=0;j<steps;j++) {
        for(i=0;i<dim;i++) {
            sum=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m!=i) {
                    sum+=A.Get(i,m)*x[m];
                }
            }
            x[i]=1/A.Get(i,i)*(b[i]-sum);
        }
    }
}

void Algorithms::GaussSeidelMethod(PoissonMatrix& A, vector<double>& x, const vector<double>& b, int steps) {
    double sum1,sum2;
    int i,j,k,m;
    for(j=0;j<steps;j++) {
        for(i=0;i<dim;i++) {
            sum1=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m<i) {
                    sum1+=A.Get(i,m)*x[m];
                }
            }
            sum2=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m>i) {
                    sum2+=A.Get(i,m)*x[m];
                }
            }
            x[i] = 1/A.Get(i,i)*(b[i]-sum1-sum2);
        }
    }
}

void Algorithms::SORMethod(PoissonMatrix& A, vector<double>& x, const vector<double>& b, int steps) {
    double sum1,sum2,Pi=3.141592654,omega=2/(1+sqrt(1-pow(cos(Pi*h),2)));
    int i,j,k,m;
    for(j=0;j<steps;j++) {
        for(i=0;i<dim;i++) {
            sum1=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m<i) {
                    sum1+=A.Get(i,m)*x[m];
                }
            }
            sum2=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m>=i) {
                    sum2+=A.Get(i,m)*x[m];
                }
            }
            x[i] = x[i]-omega*1/A.Get(i,i)*(sum1+sum2-b[i]);
        }
    }
}

vector<double> Algorithms::Restriction(vector<double>& r, int n) {
    vector<double> r2h;
    r2h.resize(pow((((n+1)/2)-1),2));
    int k=0;
    int l=0;
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            if(i%2==0 && j%2==0) {
                r2h[l]=r[k];
                l++;
            }
            k++;
        }
    }
    return r2h;
}

vector<double> Algorithms::Interpolation(vector<double>& E2h, int n) {
    int i,j,k,l;
    vector<double> E;
    E.resize(n*n);
    // for(i=0,k=0;i<n,k<E2h.size();i++,k++) {
    //     for(j=0,l=0;j<n,l<E2h.size();j++,l++) {

    //     }
    // }
    for(i=1,k=0;i<=n;i++) {
        for(j=1;j<=n;j++,k++) {
            if(i%2==0 && j%2==0) {
                E[k]=E2h[l];
                l++;
            }
        }
    }
    k=0;
    for(i=0;i<=n;i++) {
        for(j=1;j<=n;j++) {
            if(i%2!=0 && j%2==0) {
                if(i==0) {
                    E[k]=1/2*E[k+1];
                } else if(i==n) {
                    E[k]=1/2*E[k-1];
                } else {
                    E[k]=1/2*(E[k-1]+E[k+1]);
                }
            }
            if(i>=1 && j>=1) k++;
        }
    }
    k=0;
    for(i=1;i<=n;i++) {
        for(j=0;j<=n;j++) {
            if(i%2==0 && j%2!=0) {
                if(j==0) {
                    E[k]=1/2*E[k+n];
                } else if(j==n) {
                    E[k]=1/2*E[k-n];
                } else {
                    E[k]=1/2*(E[k-n]+E[k+n]);
                }
            }
            if(i>=1 && j>=1) k++;
        }
    }
    k=0;
    for(i=1;i<=n;i++) {
        for(j=1;j<=n;j++) {
            if(i%2!=0 && j%2!=0) {
                E[k]=1/4*(E[k-n-1]+E[k-n+1]+E[k+n-1]+E[k+n+1]);
            } else if(i==1 && j==1) {
                E[k]=1/4*E[k+n+1];
            } else if(i==n && j==1) {
                E[k]=1/4*E[k+n-1];
            } else if(i==1 && j==n) {
                E[k]=1/4*E[k-n+1];
            } else if(i==n && j==n) {
                E[k]=1/4*E[k-n-1];
            } else if(i==1 && j%2!=0 && j!=1 && j!=n) {
                E[k]=1/4*(E[k+n-1]+E[k+n-1]);
            } else if(i==n && j%2!=0 && j!=1 && j!=n) {
                E[k]=1/4*(E[k-n-1]+E[k-n+1]);
            } else if(j==1 && i%2!=0 && i!=1 && i!=n) {
                E[k]=1/4*(E[k+n+1]+E[k-n+1]);
            } else if(j==n && i%2!=0 && i!=1 && i!=n) {
                E[k]=1/4*(E[k+n-1]+E[k-n-1]);
            }
            if(i>=0 && j>=0) {
                k++;
            }
        }
    }
    return E;
}

vector<double> Algorithms::MultiGridMethod(PoissonMatrix& A,Vectors& V,vector<double>& x,vector<double>& b) {
    int i,dim=A.Size(),n=sqrt(dim);
    if(n==1) {
        x[0]=1;
        return x;
    } else {
        vector<double> Ax,r,r2h,E2h,E;
        Ax.assign(dim,0);
        r.resize(dim);
        r2h.resize(pow(((n+1)/2)-1,2));
        E2h.resize(pow(((n+1)/2)-1,2));
        GaussSeidelMethod(A,x,b,3);
        MatrixVectorMultiplyer(A,V.x,Ax);
        for(i=0;i<dim;i++) {
            r[i]=b[i]-Ax[i];
        }
        r2h=Restriction(r,n);
        A.Resize(((n+1)/2)-1);
        E2h=MultiGridMethod(A,V,E2h,r2h);
        E.resize(dim);
        E=Interpolation(E2h,n);
        for(i=0;i<dim;i++) {
            x[i]=x[i]+E[i];
        }
        GaussSeidelMethod(A,x,b,3);
        return x;
    }
}