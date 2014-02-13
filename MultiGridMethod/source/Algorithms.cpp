#include "classes.h"

Algorithms::Algorithms(int n) {
    this->n=n;
    this->dim=n*n;
    this->h=1.0/(double)(n+1);
    this->Vcounter=0;
    this->Wcounter1=0;
    this->Wcounter2=0;
}

Algorithms::~Algorithms() {}

void Algorithms::CG(Matrix& A,vector<double>& x,const vector<double>& b) {
    double alpha,beta=0.0,num1,num2,denom;
    int dim=x.size();
    // int steps=0;

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
        // steps++;
    }
    // printf("CGSteps: %d\n", steps);
}

void Algorithms::JacobiSolver(Matrix& A, vector<double>& x, const vector<double>& b) {
    vector<double> tmp;
    tmp=x;
    double sum;//,omega=4.0/5.0;
    int n=sqrt(x.size()),dim=n*n;
    double h=1.0/(double)(n+1);
    vector<double> r(x.size(),0);
    r=b-A*x;
    double TOL=(r|r)*pow(10,-2);
    while(TOL<=(r|r)) {
        for(int i=0;i<n*n;i++) {
            sum=0.0;
            if(i>=n && i<dim-n) sum+=tmp[i-n]+tmp[i+n];
            if(i<n) sum+=tmp[i+n];
            if(i>=dim-n) sum+=tmp[i-n];
            if(i%n!=0) sum+=tmp[i-1];
            if(i%n!=n-1) sum+=tmp[i+1];
            x[i]=1.0/(4.0*pow(1.0/h,2))*(b[i]+pow(1.0/h,2)*sum);
            tmp[i]=x[i];
        }
        r=b-A*x;
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
                // r2h[l]=1.0/8.0*(4.0*r[k]+r[k-1]+r[k+1]+r[k+n]+r[k-n]);
                r2h[l]=1.0/16.0*(4.0*r[k]+2.0*(r[k-1]+r[k+1]+r[k-n]+r[k+n])+r[k+n-1]+r[k+n+1]+r[k-n-1]+r[k-n+1]);
                // r2h[l]=r[k];
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
                if(j==1) E[k]=1.0/2.0*(E[k+1]);//+g(0,i*h));
                if(j==n) E[k]=1.0/2.0*(E[k-1]);//+g(1,i*h));
            }
            if(i%2!=0 && j%2==0) {
                if(i!=1 && i!=n) E[k]=1.0/2.0*(E[k-n]+E[k+n]);
                if(i==1) E[k]=1.0/2.0*(E[k+n]);//+g(j*h,0));
                if(i==n) E[k]=1.0/2.0*(E[k-n]);//+g(j*h,1));
            }
            if(i%2!=0 && j%2!=0) {
                if(i!=1 && j!=1 && i!=n && j!=n) E[k]=1.0/4.0*(E[k+n-1]+E[k+n+1]+E[k-n+1]+E[k-n-1]);
                if(i==1 && j==1) E[k]=1.0/1.0*(E[k+n+1]);//+g(0,0)+g(0,i*h)+g(j*h,0));
                if(i==n && j==n) E[k]=1.0/1.0*(E[k-n-1]);//+g(1,1)+g(0,i*h)+g(j*h,0));
                if(i==1 && j==n) E[k]=1.0/1.0*(E[k+n-1]);//+g(1,0)+g(0,i*h)+g(j*h,0));
                if(i==n && j==1) E[k]=1.0/1.0*(E[k-n+1]);//+g(0,1)+g(0,i*h)+g(j*h,0));
                if(i==1 && j!=1 && j!=n) E[k]=1.0/2.0*(E[k+n+1]+E[k+n-1]);//+g((j-1)*h,0)+g((j+1)*h,0));
                if(i==n && j!=1 && j!=n) E[k]=1.0/2.0*(E[k-n-1]+E[k-n+1]);//+g((j-1)*h,1)+g((j+1)*h,1));
                if(j==1 && i!=1 && i!=n) E[k]=1.0/2.0*(E[k+n+1]+E[k-n+1]);//+g(0,(i+1)*h)+g(0,(i-1)*h));
                if(j==n && i!=1 && i!=n) E[k]=1.0/2.0*(E[k+n-1]+E[k-n-1]);//+g(1,(i+1)*h)+g(1,(i-1)*h));
            }

        }
    }
}

void Algorithms::TwoGrid(Matrix& A,vector<double>& x,const vector<double>& b) {
    int dim=pow(n,2),N2h=(n+1)/2-1,dim2h=pow(N2h,2);
    vector<double> r(dim,0),E(dim,0),r2h(dim2h,0),E2h(dim2h,0);
    JacobiRelaxation(A,x,b,4);
    r=b-A*x;
    Restriction(r,r2h,n);
    CG(A,E2h,r2h);
    Interpolation(E2h,E,n);
    x+=E;
    JacobiRelaxation(A,x,b,4);
}

// vector<double> Algorithms::V_Cycle(Matrix& A,vector<double>& x,const vector<double>& b,int n) {
//     int dim=pow(n,2),N2h=(n+1)/2-1,dim2h=pow(N2h,2);
//     vector<double> r(dim,0),E(dim,0),r2h(dim2h,0),E2h(dim2h,0);
//     if(n==15) {
//         CG(A,x,b);
//         // x[0]=b[0]/4.0;
//         return x;
//     } else {
//         JacobiRelaxation(A,x,b,3);        
//         r=b-A*x;
//         Restriction(r,r2h,n);
//         E2h=V_Cycle(A,E2h,r2h,N2h);
//         E2h=V_Cycle(A,E2h,r2h,N2h);
//         Interpolation(E2h,E,n);
//         x+=E;
//         JacobiRelaxation(A,x,b,3);        
//         return x;
//     }
// }

vector<double> Algorithms::V_Cycle(Matrix& A,vector<double>& x,const vector<double>& b,int lambda,int theta) {
    int dim=x.size(),n=sqrt(dim),N2h=(n+1)/2-1,dim2h=pow(N2h,2);
    vector<double> r(dim,0),E(dim,0),r2h(dim2h,0),E2h(dim2h,0);
    if(this->Vcounter==lambda) {
        CG(A,x,b);
        return x;
    } else {
        this->Vcounter++;
        JacobiRelaxation(A,x,b,3);        
        r=b-A*x;
        Restriction(r,r2h,n);
        E2h=V_Cycle(A,E2h,r2h,lambda,theta);
        if(theta==1) E2h=V_Cycle(A,E2h,r2h,lambda,theta);
        Interpolation(E2h,E,n);
        x+=E;
        JacobiRelaxation(A,x,b,3);
        this->Vcounter--;
        return x;
    }
}

// vector<double> Algorithms::W_Cycle(Matrix& A,vector<double>& x,const vector<double>& b,int lambda,int theta,int kappa) {
//     this->Vcounter++;
//     this->Wcounter1++;
//     this->Wcounter2++;
//     int dim=x.size(),n=sqrt(dim),N2h=(n+1)/2-1,dim2h=pow(N2h,2);
//     vector<double> r(dim,0),E(dim,0),r2h(dim2h,0),E2h(dim2h,0);
//     if(this->Vcounter==lambda) {
//         printf("Solve direct\n");
//         CG(A,x,b);
//         return x;
//     } else if(this->Wcounter1<theta && this->Vcounter<theta) {
//         printf("First node\n");
//         JacobiRelaxation(A,x,b,3);
//         r=b-A*x;
//         Restriction(r,r2h,n);
//         E2h=W_Cycle(A,E2h,r2h,lambda,theta,kappa);
//         Interpolation(E2h,E,n);
//         x+=E;
//         JacobiRelaxation(A,x,b,3);
//         return x;
//     } else if(this->Wcounter2<kappa && this->Vcounter<kappa) {
//         printf("Second node\n");
//         JacobiRelaxation(A,x,b,3);
//         r=b-A*x;
//         Restriction(r,r2h,n);
//         E2h=W_Cycle(A,E2h,r2h,lambda,theta,kappa);
//         Interpolation(E2h,E,n);
//         x+=E;
//         JacobiRelaxation(A,x,b,3);
//         return x;
//     } else {
//         printf("ois easy...\n");
//         JacobiRelaxation(A,x,b,3);        
//         r=b-A*x;
//         Restriction(r,r2h,n);
//         E2h=W_Cycle(A,E2h,r2h,lambda,theta,kappa);
//         Interpolation(E2h,E,n);
//         x+=E;
//         JacobiRelaxation(A,x,b,3);      
//         return x;
//     }
// }

void Algorithms::MultiGridMethod(Matrix& A,vector<double>& x,const vector<double>& b,const vector<double>& solved) {
    int steps=0,dim=x.size();
    vector<double> r(dim),tmp(dim);
    r=x-solved;
    double TOL=pow(10,-3)*(r|r);
    while(TOL<=(r|r)) {
        x=V_Cycle(A,x,b,2,0);
        // x=W_Cycle(A,x,b,4,2,3);
        r=x-solved;
        steps++;
        if(steps==500) r.assign(dim,0);
    }
    printf("MultiGridSteps: %d\n", steps);
}