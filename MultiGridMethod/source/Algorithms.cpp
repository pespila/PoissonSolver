#include "classes.h"

Algorithms::Algorithms(int n) {
    this->n=n;
    this->dim=n*n;
    this->h=1.0/(double)(n+1);
}

Algorithms::~Algorithms() {}

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
            x[i]=omega*1.0/4.0*(b[i]+sum)+tmp[i]*(1-omega);
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
                if(j==1) E[k]=1.0/2.0*(E[k+1]);
                if(j==n) E[k]=1.0/2.0*(E[k-1]);
            }
            if(i%2!=0 && j%2==0) {
                if(i!=1 && i!=n) E[k]=1.0/2.0*(E[k-n]+E[k+n]);
                if(i==1) E[k]=1.0/2.0*(E[k+n]);
                if(i==n) E[k]=1.0/2.0*(E[k-n]);
            }
            if(i%2!=0 && j%2!=0) {
                if(i!=1 && j!=1 && i!=n && j!=n) E[k]=1.0/4.0*(E[k+n-1]+E[k+n+1]+E[k-n+1]+E[k-n-1]);
                if(i==1 && j==1) E[k]=1.0/4.0*(E[k+n+1]);
                if(i==n && j==n) E[k]=1.0/4.0*(E[k-n-1]);
                if(i==1 && j==n) E[k]=1.0/4.0*(E[k-n+1]);
                if(i==n && j==1) E[k]=1.0/4.0*(E[k+n-1]);
                if(i==1 && j!=1 && j!=n) E[k]=1.0/4.0*(E[k+n+1]+E[k+n-1]);
                if(i==n && j!=1 && j!=n) E[k]=1.0/4.0*(E[k-n-1]+E[k-n+1]);
                if(j==1 && i!=1 && i!=n) E[k]=1.0/4.0*(E[k+n+1]+E[k-n+1]);
                if(j==n && i!=1 && i!=n) E[k]=1.0/4.0*(E[k+n-1]+E[k-n-1]);
            }

        }
    }
}

vector<double> Algorithms::MultiGridAlgorithm(Matrix& A,vector<double>& X,const vector<double>& b,int n) {
    int dim=pow(n,2),N2h=(n+1)/2-1,dim2h=pow(N2h,2);
    vector<double> x(dim,0),r(dim,0),E(dim,0),r2h(dim2h,0),E2h(dim2h,0);
    if(n==this->n) x=X;
    if(n==1) {
        x[0]=b[0]/4.0;
        return x;
    } else {
        JacobiRelaxation(A,x,b,6);        
        r=b-A*x;
        Restriction(r,r2h,n);
        E2h=MultiGridAlgorithm(A,X,r2h,N2h);
        Interpolation(E2h,E,n);
        x+=E;
        JacobiRelaxation(A,x,b,6);        
        return x;
    }
}

void Algorithms::MultiGridMethod(Matrix& A,vector<double>& x,const vector<double>& b,const vector<double>& solved) {
    int steps=0;
    vector<double> r(dim);
    r=x-solved;
    double TOL=pow(10,-3)*(r|r);
    while(TOL<=(r|r)) {
        x=MultiGridAlgorithm(A,x,b,n);
        r=x-solved;
        steps++;
        if(steps==500) r.assign(dim,0);
    }
    printf("MultiGridSteps: %d\n", steps);
}