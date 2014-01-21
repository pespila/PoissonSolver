#include "classes.h"

Algorithms::Algorithms(int n) {
    this->n=n;
    this->dim=n*n;
    this->h=1.0/(double)(n+1);
}

Algorithms::~Algorithms() {
}

double Algorithms::innerProduct(const vector<double>& x,const vector<double>& y) {
    double ip=0.0;
    for (int i=0;i<dim;i++)
        ip+=x[i]*y[i];
    return ip;
}

double Algorithms::vectorNorm(const vector<double>& x) {
    double norm=0.0;
    for(int i=0;i<dim;i++) {
        norm+=pow(x[i],2);
    }
    return sqrt(norm);
}

void Algorithms::modifiedIncompleteLU(Matrix& A,WriteableMatrix& L,WriteableMatrix& U) {
    int i,j,k,m,u;
    double sum, drop;

    for(i=0;i<dim;i++) {
        drop=0;
        for(k=0;k<5;k++) {
            m=A.HashMatrix[i][k];
            if(m!=-1 && m<i) {
                sum=0;
                for(j=0;j<5;j++) {
                    u=A.HashMatrix[i][j];
                    if(u!=-1 && u<k) {
                        sum+=L.Get(i,u)*U.Get(u,m);
                    }
                }
                L.Set(i,m,(A.Get(i,m)-sum)/U.Get(m,m));
                drop+=sum;
            } else if(m!=-1 && m>=i) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m>=i) {
                    sum=0;
                    for(j=0;j<5;j++) {
                        u=A.HashMatrix[i][j];
                        if(u!=-1 && u<i) {
                            sum+=L.Get(i,u)*U.Get(u,m);
                        }
                    }
                    U.Set(i,m,(A.Get(i,m)-sum));
                    drop+=sum;
                }
            }
        }
        U.Set(i,i,(U.Get(i,i)-drop));
    }
}

void Algorithms::incompleteLU(Matrix& A, WriteableMatrix& L, WriteableMatrix& U) {
    int m,u;
    double sum;

    for(int i=0;i<dim;i++) {
        for(int k=0;k<5;k++) {
            m=A.HashMatrix[i][k];
            if(m!=-1 && m<i) {
                sum=0;
                for(int j=0;j<5;j++) {
                    u=A.HashMatrix[i][j];
                    if(u!=-1 && u<k) {
                        sum+=L.Get(i,u)*U.Get(u,m);
                    }
                }
                L.Set(i,m,((A.Get(i,m)-sum)/U.Get(m,m)));
            } else if(m!=-1 && m>=i) {
                sum=0;
                for(int j=0;j<5;j++) {
                    u=A.HashMatrix[i][j];
                    if(u!=-1 && u<i) {
                        sum+=L.Get(i,u)*U.Get(u,m);
                    }
                }
                U.Set(i,m,(A.Get(i,m)-sum));
            }
        }
    }
}

void Algorithms::JacobiMethod(Matrix& A,vector<double>& x,const vector<double>& b) {
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
            r[i]=1.0/4.0*(b[i]-r[i]);
            x[i]+=r[i];
            r[i]=x[i]-solved[i];
        }
        // x+=1.0/4.0*(b-A*x);
        steps++;
    }
    printf("JacobianSteps: %d\n", steps);
}

void Algorithms::JacobiRelaxationMethod(Matrix& A,vector<double>& x,const vector<double>& b) {
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

void Algorithms::GaussSeidelMethod(Matrix& A,vector<double>& x,const vector<double>& b) {
    double sum;
    int steps=0;

    vector<double> solved(dim),r(dim,1);
    for(int i=1,k=0;i<(n+1);i++) {
        for(int j=1;j<(n+1);j++,k++) {
            solved[k]=g(j*h,i*h);
        }
    }

    double TOL=pow(10,-3)*vectorNorm(solved);

    while(TOL<vectorNorm(r)) {
        for(int i=0;i<dim;i++) {
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
            r[i]=x[i]-solved[i];
        }
        steps++;
    }
    printf("GaussSeidelSteps: %d\n", steps);
}

void Algorithms::SORMethod(Matrix& A,vector<double>& x,const vector<double>& b) {
    double sum,Pi=3.141592654,omega=2/(1+sqrt(1-pow(cos(Pi*h),2)));
    int steps=0;

    vector<double> solved(dim),r(dim,1);
    for(int i=1,k=0;i<(n+1);i++) {
        for(int j=1;j<(n+1);j++,k++) {
            solved[k]=g(j*h,i*h);
        }
    }

    double TOL=pow(10,-3)*vectorNorm(solved);

    while(TOL<vectorNorm(r)) {
        for(int i=0;i<dim;i++) {
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
            r[i]=x[i]-solved[i];
        }
        steps++;
    }
    printf("GaussSeidelSteps: %d\n", steps);
}

void Algorithms::CG(Matrix& A, vector<double>& x,const vector<double>& b) {
    double alpha,beta=0.0,num1,num2,denom;
    int steps=0;

    vector<double> solved(dim),r(dim),Ap(dim),tmp(dim,1);
    for(int i=1,k=0;i<(n+1);i++) {
        for(int j=1;j<(n+1);j++,k++) {
            solved[k]=g(j*h,i*h);
        }
    }
    r=b-A*x;
    vector<double> p(r),rTmp(r);

    num2=innerProduct(rTmp,rTmp);
    num1=num2;
    double TOL=pow(10,-3)*vectorNorm(solved);
    while(TOL<vectorNorm(tmp)) {
        Ap=A*p;
        denom=innerProduct(p,Ap);
        alpha=num2/denom;

        for(int i=0;i<dim;i++) {
            x[i]+=alpha*p[i];
            r[i]-=alpha*Ap[i];
        }

        num2=innerProduct(r,r);
        beta=num2/num1;

        for(int i=0;i<dim;i++) {
            p[i]=r[i]+beta*p[i];
            tmp[i]=x[i]-solved[i];
        }
        rTmp=r;
        num1=num2;
        steps++;
    }
    printf("CGSteps: %d\n", steps);
}

void Algorithms::PCG(Matrix& A,WriteableMatrix& L,WriteableMatrix& U,vector<double>& x,const vector<double>& b) {
    double alpha,beta=0.0,num1,num2,denom;
    int steps=0;

    vector<double> solved(dim),r(dim),Ap(dim),tmp(dim,1);
    for(int i=1,k=0;i<(n+1);i++) {
        for(int j=1;j<(n+1);j++,k++) {
            solved[k]=g(j*h,i*h);
        }
    }
    r=b-A*x;
    vector<double> z(r),rTmp(r);
    LUsolverLower(A,L,z);
    LUsolverUpper(A,U,z);
    vector<double> p(z),zTmp(z);

    num2=innerProduct(zTmp,rTmp);
    num1=num2;
    double TOL=pow(10,-3)*vectorNorm(solved);
    while(TOL<vectorNorm(tmp)) {
        Ap=A*p;
        denom=innerProduct(p,Ap);
        alpha=num2/denom;
        
        for(int i=0;i<dim;i++) {
            x[i]+=alpha*p[i];
            r[i]-=alpha*Ap[i];
        }

        z=r;
        LUsolverLower(A,L,z);
        LUsolverUpper(A,U,z);
        zTmp=z;

        num2=innerProduct(z,r);
        beta=num2/num1;

        for(int i=0;i<dim;i++) {
            p[i]=z[i]+beta*p[i];
            tmp[i]=x[i]-solved[i];
        }
        rTmp=r;
        num1=num2;
        steps++;
    }
    printf("PCGSteps: %d\n", steps);
}

void Algorithms::LUsolverLower(Matrix& A,Matrix& L,vector<double>& z) {
    int m;
    for(int i=0;i<dim;i++) {
        for(int j=0;j<5;j++) {
            m=A.HashMatrix[i][j];
            if(m!=-1 && m<i) {
                z[i]-=L.Get(i,m)*z[m];
            }
        }
        z[i]/=L.Get(i,i);
    }
}

void Algorithms::LUsolverUpper(Matrix& A,Matrix& U,vector<double>& z) {
    int m;
    for(int i=dim-1;i>=0;i--) {
        for(int j=0;j<5;j++) {
            m=A.HashMatrix[i][j];
            if(m!=-1 && m>=i) {
                z[i]-=U.Get(i,m)*z[m];
            }
         }
        z[i]/=U.Get(i,i);
    }
}