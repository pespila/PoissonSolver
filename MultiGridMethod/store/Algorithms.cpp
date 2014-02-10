#include "classes.h"

Algorithms::Algorithms() {
}

Algorithms::~Algorithms() {
}

void Algorithms::JacobiMethod(Matrix& A, vector<double>& x, const vector<double>& b, int steps) {
    if((int)x.size()!=(int)b.size()) {
        printf("Error in JacobiMethod: Dimension doesn't match!\n");
    } else {
        for(int j=0;j<steps;j++) {
            x+=1.0/4.0*(b-A*x);
        }
    }
}

void Algorithms::JacobiRelaxationMethod(Matrix& A, vector<double>& x, const vector<double>& b, int steps) {
    if((int)x.size()!=(int)b.size()) {
        printf("Error in JacobiRelaxationMethod: Dimension doesn't match!\n");
    } else {
        double omega=4.0/5.0*1.0/4.0;
        for(int j=0;j<steps;j++) {
            x+=omega*(b-A*x);
        }
    }
}

void Algorithms::JacobiRelaxationSolver(Matrix& A,vector<double>& x,const vector<double>& b) {
    if((int)x.size()!=(int)b.size()) {
        printf("Error in JacobiRelaxationSolver: Dimension doesn't match!\n");
    } else {
        vector<double> r((int)x.size(),0);
        r=b-A*x;
        double TOL=(r|r);
        while(TOL<=(r|r)) {
            x+=1.0/4.0*r;
            r=b-A*x;
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
                if(j!=1 && j!=n) {
                    E[k]=1.0/2.0*(E[k-1]+E[k+1]);
                }
                if(j==1) {
                    E[k]=1.0/2.0*(E[k+1]);
                }
                if(j==n) {
                    E[k]=1.0/2.0*(E[k-1]);
                }
            }
            if(i%2!=0 && j%2==0) {
                if(i!=1 && i!=n) {
                    E[k]=1.0/2.0*(E[k-n]+E[k+n]);
                }
                if(i==1) {
                    E[k]=1.0/2.0*(E[k+n]);
                }
                if(i==n) {
                    E[k]=1.0/2.0*(E[k-n]);
                }
            }
            if(i%2!=0 && j%2!=0) {
                if(i!=1 && j!=1 && i!=n && j!=n) {
                    E[k]=1.0/4.0*(E[k+n-1]+E[k+n+1]+E[k-n+1]+E[k-n-1]);
                }
                if(i==1 && j==1) {
                    E[k]=1.0/4.0*(E[k+n+1]);
                }
                if(i==n && j==n) {
                    E[k]=1.0/4.0*(E[k-n-1]);
                }
                if(i==1 && j==n) {
                    E[k]=1.0/4.0*(E[k-n+1]);
                }
                if(i==n && j==1) {
                    E[k]=1.0/4.0*(E[k+n-1]);
                }
                if(i==1 && j!=1 && j!=n) {
                    E[k]=1.0/4.0*(E[k+n+1]+E[k+n-1]);
                }
                if(i==n && j!=1 && j!=n) {
                    E[k]=1.0/4.0*(E[k-n-1]+E[k-n+1]);
                }
                if(j==1 && i!=1 && i!=n) {
                    E[k]=1.0/4.0*(E[k+n+1]+E[k-n+1]);
                }
                if(j==n && i!=1 && i!=n) {
                    E[k]=1.0/4.0*(E[k+n-1]+E[k-n-1]);
                }
            }

        }
    }
}

vector<double> Algorithms::MultiGridAlgorithm(Matrix& A,vector<double>& u,const vector<double>& f) {
    if((int)u.size()!=(int)f.size()) {
        printf("Error in MultiGridAlgorithm: Dimensions doesn't match!\n");
        return u;
    }
    else {
        int dim=(int)u.size(),n,n2h,dim2h;
        if(A.Get(0,0)==4) {
            n=sqrt(dim);n2h=(n+1)/2-1;dim2h=pow(n2h,2);   
        }
        if(A.Get(0,0)==2) {
            n=dim;n2h=(n+1)/2-1;dim2h=n2h;
        }
        vector<double> E(dim,0),r2h(dim2h,0),E2h(dim2h,0);
        if(n==1) {
            JacobiRelaxationSolver(A,u,f);
            return u;
        } else {
            JacobiRelaxationMethod(A,u,f,6);
            Restriction(f-A*u,r2h,n);
            E2h=MultiGridAlgorithm(A,E2h,r2h);
            Interpolation(E2h,E,n);
            u+=E;
            JacobiRelaxationMethod(A,u,f,6);        
            return u;
        }
    }
}

void Algorithms::MultiGridMethod(Matrix& A,vector<double>& u,const vector<double>& f,const vector<double>& solved) {
    if((int)u.size()!=(int)f.size()) {
        printf("Error in MultiGridMethod: Dimensions doesn't match!\n");
    }
    else {
        int steps=0,dim=(int)u.size();
        vector<double> r(dim);
        r=u-solved;
        double TOL=pow(10,-3)*(r|r);
        while(TOL<=(r|r)) {
            u=MultiGridAlgorithm(A,u,f);
            r=u-solved;
            steps++;
            if(steps==200) {
                r.assign(dim,0);
            }
        }
        printf("MultiGridSteps: %d\n", steps);
    }
}