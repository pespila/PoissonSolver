#include "classes.h"

Algorithms::Algorithms(int n) {
    this->n=n;
    this->dim=n*n;
    this->h=1.0/(double)(n+1);
    this->Vcounter=0;
}

Algorithms::~Algorithms() {}

void Algorithms::PCGdirect(Matrix& A,WriteableMatrix& L,WriteableMatrix& U,vector<double>& x,const vector<double>& b) {
    int dim=x.size();
    double alpha,beta=0.0,num1,num2,denom;

    vector<double> r(dim),Ap(dim);
    r=b-A*x;
    vector<double> z(r),rTmp(r);
    z=L*z;
    LUsolverUpper(A,U,z);
    vector<double> p(z),zTmp(z);

    num2=zTmp*rTmp;
    num1=num2;
    double TOL=pow(10,-2)*(r|r);
    while(TOL<(r|r)) {
        Ap=A*p;
        denom=p*Ap;
        alpha=num2/denom;
        
        for(int i=0;i<dim;i++) {
            x[i]+=alpha*p[i];
            r[i]-=alpha*Ap[i];
        }

        z=r;
        z=L*z;
        LUsolverUpper(A,U,z);
        zTmp=z;

        num2=z*r;
        beta=num2/num1;

        for(int i=0;i<dim;i++) {
            p[i]=z[i]+beta*p[i];
        }
        rTmp=r;
        num1=num2;
    }
}

void Algorithms::CGdirect(Matrix& A,vector<double>& x,const vector<double>& b) {
    double alpha,beta=0.0,num1,num2,denom;
    int dim=x.size();

    vector<double> r(dim),Ap(dim);
    r=b-A*x;
    vector<double> p(r),rTmp(r);

    num2=rTmp*rTmp;
    num1=num2;
    double TOL=pow(10,-3)*(r|r);
    while(TOL<(r|r)) {
        Ap=A*p;
        denom=p*Ap;
        alpha=num2/denom;

        for(int i=0;i<dim;i++) {
            x[i]+=alpha*p[i];
            r[i]-=alpha*Ap[i];
        }

        num2=r*r;
        beta=num2/num1;

        for(int i=0;i<dim;i++) {
            p[i]=r[i]+beta*p[i];
        }
        rTmp=r;
        num1=num2;
    }
}

void Algorithms::JacobiRelaxation(Matrix& A, vector<double>& x, const vector<double>& b, int steps) {
    vector<double> tmp;
    tmp=x;
    double sum,omega=4.0/5.0;
    int n=sqrt(x.size()),dim=n*n;
    for(int k=0;k<steps;k++) {
        for(int i=0;i<n*n;i++) {
            sum=0.0;
            if(i>=n && i<dim-n) sum+=tmp[i-n]+tmp[i+n];
            if(i<n) sum+=tmp[i+n];
            if(i>=dim-n) sum+=tmp[i-n];
            if(i%n!=0) sum+=tmp[i-1];
            if(i%n!=n-1) sum+=tmp[i+1];
            x[i]=omega*1.0/(4.0*pow(1.0/h,2))*(b[i]+pow(1.0/h,2)*sum)+tmp[i]*(1-omega);
            tmp[i]=x[i];
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

void Algorithms::Interpolation(const vector<double>& E2h,vector<double>& E,int n) {
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
                if(j!=1 && j!=n) E[k]=1.0/2.0*(E[k-1]+E[k+1]);
                if(j==1) E[k]=1.0/1.0*(E[k+1]);
                if(j==n) E[k]=1.0/1.0*(E[k-1]);
            }
            if(i%2!=0 && j%2==0) {
                if(i!=1 && i!=n) E[k]=1.0/2.0*(E[k-n]+E[k+n]);
                if(i==1) E[k]=1.0/1.0*(E[k+n]);
                if(i==n) E[k]=1.0/1.0*(E[k-n]);
            }
            if(i%2!=0 && j%2!=0) {
                if(i!=1 && j!=1 && i!=n && j!=n) E[k]=1.0/4.0*(E[k+n-1]+E[k+n+1]+E[k-n+1]+E[k-n-1]);
                if(i==1 && j==1) E[k]=1.0/1.0*(E[k+n+1]);
                if(i==n && j==n) E[k]=1.0/1.0*(E[k-n-1]);
                if(i==1 && j==n) E[k]=1.0/1.0*(E[k+n-1]);
                if(i==n && j==1) E[k]=1.0/1.0*(E[k-n+1]);
                if(i==1 && j!=1 && j!=n) E[k]=1.0/2.0*(E[k+n+1]+E[k+n-1]);
                if(i==n && j!=1 && j!=n) E[k]=1.0/2.0*(E[k-n-1]+E[k-n+1]);
                if(j==1 && i!=1 && i!=n) E[k]=1.0/2.0*(E[k+n+1]+E[k-n+1]);
                if(j==n && i!=1 && i!=n) E[k]=1.0/2.0*(E[k+n-1]+E[k-n-1]);
            }

        }
    }
}

vector<double> Algorithms::Cycle(Matrix& A,vector<double>& x,const vector<double>& b,int lambda,int theta,Matrix& B,WriteableMatrix& L,WriteableMatrix& U) {
    int dim=x.size(),n=sqrt(dim),N2h=(n+1)/2-1,dim2h=pow(N2h,2);
    vector<double> r(dim,0),E(dim,0),r2h(dim2h,0),E2h(dim2h,0);
    if(this->Vcounter==lambda) {
        PCGdirect(B,L,U,x,b);
        // CGdirect(A,x,b);
        return x;
    } else {
        this->Vcounter++;
        JacobiRelaxation(A,x,b,3);
        r=b-A*x;
        Restriction(r,r2h,n);
        E2h=Cycle(A,E2h,r2h,lambda,theta,B,L,U);
        if(theta==1) E2h=Cycle(A,E2h,r2h,lambda,theta,B,L,U);
        Interpolation(E2h,E,n);
        x+=E;
        JacobiRelaxation(A,x,b,3);
        this->Vcounter--;
        return x;
    }
}

int Algorithms::MultiGridMethod(Matrix& A,vector<double>& x,const vector<double>& b,const vector<double>& solved,int alg) {
    int steps=0,dim=x.size();
    vector<double> r(dim);
    r=x-solved;
    int numberOfGrids=2,VW=0,n=0;
    if(alg==4) {
        numberOfGrids=1;
        VW=0;
        n=((int)sqrt(dim)+1)/2-1;
    }
    if(alg==5) {
        numberOfGrids=2;
        VW=0;
        n=((int)sqrt(dim)+1)/4-1;
    }
    if(alg==6) {
        numberOfGrids=2;
        VW=1;
        n=((int)sqrt(dim)+1)/4-1;
    }
    PoissonMatrix B(n);
    LowerMatrix L(n);
    UpperMatrix U(n);
    B.InitHashMatrix();
    modifiedIncompleteLU(B,L,U);
    double TOL=pow(10,-3)*(r|r);
    while(TOL<=(r|r)) {
        x=Cycle(A,x,b,numberOfGrids,VW,B,L,U);
        r=x-solved;
        steps++;
    }
    return steps;
}

int Algorithms::JacobiMethod(Matrix& A,vector<double>& x,const vector<double>& b,const vector<double>& solved) {
    vector<double> tmp;
    double sum;
    int n=sqrt(x.size()),dim=n*n,steps=0;
    double h=1.0/(double)(n+1);
    vector<double> r(x.size());
    r=x-solved;
    double TOL=(r|r)*pow(10,-3);
    while(TOL<=(r|r)) {
        tmp=x;
        for(int i=0;i<n*n;i++) {
            sum=0.0;
            if(i>=n && i<dim-n) sum+=tmp[i-n]+tmp[i+n];
            if(i<n) sum+=tmp[i+n];
            if(i>=dim-n) sum+=tmp[i-n];
            if(i%n!=0) sum+=tmp[i-1];
            if(i%n!=n-1) sum+=tmp[i+1];
            x[i]=1.0/(4.0*pow(1.0/h,2))*(b[i]+pow(1.0/h,2)*sum);
            r[i]=x[i]-solved[i];
        }
        steps++;
    }
    return steps;
}

int Algorithms::JacobiRelaxationMethod(Matrix& A,vector<double>& x,const vector<double>& b,const vector<double>& solved) {
    vector<double> tmp;
    double sum,omega=4.0/5.0;
    int n=sqrt(x.size()),dim=n*n,steps=0;
    double h=1.0/(double)(n+1);
    vector<double> r(x.size());
    r=x-solved;
    double TOL=(r|r)*pow(10,-3);
    while(TOL<=(r|r)) {
        tmp=x;
        for(int i=0;i<n*n;i++) {
            sum=0.0;
            if(i>=n && i<dim-n) sum+=tmp[i-n]+tmp[i+n];
            if(i<n) sum+=tmp[i+n];
            if(i>=dim-n) sum+=tmp[i-n];
            if(i%n!=0) sum+=tmp[i-1];
            if(i%n!=n-1) sum+=tmp[i+1];
            x[i]=omega*1.0/(4.0*pow(1.0/h,2))*(b[i]+pow(1.0/h,2)*sum)+tmp[i]*(1-omega);
            r[i]=x[i]-solved[i];
        }
        steps++;
    }
    return steps;
}

int Algorithms::GaussSeidelMethod(Matrix& A,vector<double>& x,const vector<double>& b,const vector<double>& solved) {
    double sum;
    int n=sqrt(x.size()),dim=n*n,steps=0;
    double h=1.0/(double)(n+1);
    vector<double> r(x.size());
    r=x-solved;
    double TOL=(r|r)*pow(10,-3);
    while(TOL<=(r|r)) {
        for(int i=0;i<n*n;i++) {
            sum=0.0;
            if(i>=n && i<dim-n) sum+=x[i-n]+x[i+n];
            if(i<n) sum+=x[i+n];
            if(i>=dim-n) sum+=x[i-n];
            if(i%n!=0) sum+=x[i-1];
            if(i%n!=n-1) sum+=x[i+1];
            x[i]=1.0/(4.0*pow(1.0/h,2))*(b[i]+pow(1.0/h,2)*sum);
            r[i]=x[i]-solved[i];
        }
        steps++;
    }
    return steps;
}

int Algorithms::SORMethod(Matrix& A,vector<double>& x,const vector<double>& b,const vector<double>& solved) {
    double sum,Pi=3.141592654,omega=2/(1+sqrt(1-pow(cos(Pi*h),2)));
    int n=sqrt(x.size()),dim=n*n,steps=0;
    double h=1.0/(double)(n+1);
    vector<double> r(x.size());
    r=x-solved;
    double TOL=(r|r)*pow(10,-3);
    while(TOL<=(r|r)) {
        for(int i=0;i<n*n;i++) {
            sum=0.0;
            if(i>=n && i<dim-n) sum+=x[i-n]+x[i+n];
            if(i<n) sum+=x[i+n];
            if(i>=dim-n) sum+=x[i-n];
            if(i%n!=0) sum+=x[i-1];
            if(i%n!=n-1) sum+=x[i+1];
            x[i]=omega*1.0/(4.0*pow(1.0/h,2))*(b[i]+pow(1.0/h,2)*sum)+x[i]*(1-omega);
            r[i]=x[i]-solved[i];
        }
        steps++;
    }
    return steps;
}

int Algorithms::CG(Matrix& A, vector<double>& x,const vector<double>& b,const vector<double>& solved) {
    double alpha,beta=0.0,num1,num2,denom;
    int steps=0;

    vector<double> r(dim),Ap(dim),tmp(dim,1);
    r=b-A*x;
    vector<double> p(r),rTmp(r);

    num2=rTmp*rTmp;
    num1=num2;
    double TOL=pow(10,-3)*((x-solved)|(x-solved));
    while(TOL<(tmp|tmp)) {
        Ap=A*p;
        denom=p*Ap;
        alpha=num2/denom;

        for(int i=0;i<dim;i++) {
            x[i]+=alpha*p[i];
            r[i]-=alpha*Ap[i];
        }

        num2=r*r;
        beta=num2/num1;

        for(int i=0;i<dim;i++) {
            p[i]=r[i]+beta*p[i];
            tmp[i]=x[i]-solved[i];
        }
        rTmp=r;
        num1=num2;
        steps++;
    }
    return steps;
}

int Algorithms::PCG(Matrix& A,WriteableMatrix& L,WriteableMatrix& U,vector<double>& x,const vector<double>& b,const vector<double>& solved) {
    double alpha,beta=0.0,num1,num2,denom;
    int steps=0;

    vector<double> r(dim),Ap(dim),tmp(dim,1);
    r=b-A*x;
    vector<double> z(r),rTmp(r);
    z=L*z;
    LUsolverUpper(A,U,z);
    vector<double> p(z),zTmp(z);

    num2=zTmp*rTmp;
    num1=num2;
    double TOL=pow(10,-3)*((x-solved)|(x-solved));
    while(TOL<(tmp|tmp)) {
        Ap=A*p;
        denom=p*Ap;
        alpha=num2/denom;
        
        for(int i=0;i<dim;i++) {
            x[i]+=alpha*p[i];
            r[i]-=alpha*Ap[i];
        }

        z=r;
        z=L*z;
        LUsolverUpper(A,U,z);
        zTmp=z;

        num2=z*r;
        beta=num2/num1;

        for(int i=0;i<dim;i++) {
            p[i]=z[i]+beta*p[i];
            tmp[i]=x[i]-solved[i];
        }
        rTmp=r;
        num1=num2;
        steps++;
    }
    return steps;
}

void Algorithms::modifiedIncompleteLU(Matrix& A,WriteableMatrix& L,WriteableMatrix& U) {
    int i,j,k,m,u,dim=A.Size();
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

void Algorithms::LUsolverUpper(Matrix& A,Matrix& U,vector<double>& z) {
    int m;
    int dim=z.size();
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