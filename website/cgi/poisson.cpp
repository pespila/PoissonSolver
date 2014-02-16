#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <vector>
#include <assert.h>
#include <time.h>
#include <string>
using namespace std;

void print_header(void) {
    printf("Content-Type:text/html;charset=utf-8\n\n");
}

void print_table(int arg,int steps,double timer,const char* method) {
    printf("<td>%s</td><td>%d</td><td>%d</td><td>%f</td>",method,arg,steps,timer);
}

vector<double> operator-(const vector<double>& lhs, const vector<double>& rhs) {
    vector<double> tmp(lhs);
    for(int i=0;i<(int)lhs.size();i++) {
        tmp[i]=lhs[i]-rhs[i];
    }
    return tmp;
}

void operator+=(vector<double>& lhs, const vector<double>& rhs) {
    for(int i=0;i<(int)rhs.size();i++) {
        lhs[i]+=rhs[i];
    }
}

vector<double> operator*(double x, vector<double> rhs) {
    vector<double> tmp(rhs);
    for(int i=0;i<(int)rhs.size();i++) {
        tmp[i]=x*rhs[i];
    }
    return tmp;
}

double operator|(const std::vector<double>& x,const std::vector<double>& y) {
    double norm=0.0;
    for(int i=0;i<(int)x.size();i++) {
        norm+=x[i]*y[i];
    }
    return sqrt(norm);
}

double operator*(const std::vector<double>& x,const std::vector<double>& y) {
    double ip=0.0;
    for (int i=0;i<(int)x.size();i++)
        ip+=x[i]*y[i];
    return ip;
}

class Matrix {
    public:
        vector<vector<int> > HashMatrix;
        virtual int Size()=0;
        virtual double Get(int,int)=0;
        virtual vector<double> operator*(const vector<double>&)=0;
};

class WriteableMatrix : public Matrix {
    public:
        virtual void Set(int,int,double)=0;
};

class PoissonMatrix : public Matrix {
    private:
        int dim;
        int n;
        double diagonal;
        double tridiagonal;
        double identity;
    public:
        PoissonMatrix(int);
        ~PoissonMatrix();
        int Size();
        double Get(int, int);
        void InitHashMatrix();
        vector<double> operator*(const vector<double>&);
};

PoissonMatrix::PoissonMatrix(int n) {
   this->n=n;
   this->dim=n*n;
   this->diagonal=4.0*pow((n+1),2);
   this->tridiagonal=-1.0*pow((n+1),2);
   this->identity=-1.0*pow((n+1),2);
}

void PoissonMatrix::InitHashMatrix() {
    vector<int> tmp(5,-1);
    for(int i=0;i<this->dim;i++) {
        if(i<this->n) {
            if(i==0) {
                tmp[0]=i;
                tmp[1]=i+1;
                tmp[2]=i+this->n;
                tmp[3]=-1;
                tmp[4]=-1;
            } else if(i%this->n!=this->n-1 && i!=0) {
                tmp[0]=i-1;
                tmp[1]=i;
                tmp[2]=i+1;
                tmp[3]=i+this->n;
                tmp[4]=-1;
            } else {
                tmp[0]=i-1;
                tmp[1]=i;
                tmp[2]=i+this->n;
                tmp[3]=-1;
                tmp[4]=-1;
            }
        } else if(i>=this->n && i<this->dim-this->n) {
            if(i%this->n==0) {
                tmp[0]=i-this->n;
                tmp[1]=i;
                tmp[2]=i+1;
                tmp[3]=i+this->n;
                tmp[4]=-1;
            } else if(i%this->n==this->n-1) {
                tmp[0]=i-this->n;
                tmp[1]=i-1;
                tmp[2]=i;
                tmp[3]=i+this->n;
                tmp[4]=-1;
            } else {
                tmp[0]=i-this->n;
                tmp[1]=i-1;
                tmp[2]=i;
                tmp[3]=i+1;
                tmp[4]=i+this->n;
            }
        } else if(i>=this->dim-this->n) {
            if(i==this->dim-1) {
                tmp[0]=i-this->n;
                tmp[1]=i-1;
                tmp[2]=i;
                tmp[3]=-1;
                tmp[4]=-1;
            } else if(i%this->n!=0 && i!=this->dim-1) {
                tmp[0]=i-this->n;
                tmp[1]=i-1;
                tmp[2]=i;
                tmp[3]=i+1;
                tmp[4]=-1;
            } else {
                tmp[0]=i-this->n;
                tmp[1]=i;
                tmp[2]=i+1;
                tmp[3]=-1;
                tmp[4]=-1;
            }
        }
        HashMatrix.push_back(tmp);
    }
}

PoissonMatrix::~PoissonMatrix() {
    vector<vector<int> >().swap(HashMatrix);
}

int PoissonMatrix::Size() {
   return dim;
}

double PoissonMatrix::Get(int i,int j) {
   if(i==j) {
      return diagonal;
   } else if((j==(i+1) && i%n!=(n-1)) || (j==(i-1) && i%n!=0)) {
      return tridiagonal;
   } else if((j==(i+n)) || (j==(i-n))) {
      return identity;
   } else {
      return 0.0;
   }
}

vector<double> PoissonMatrix::operator*(const vector<double>& x) {
    int dim=x.size(),n=sqrt(dim);
    vector<double> tmp(dim,0);
    for(int i=0;i<dim;i++) {
        tmp[i]+=x[i]*4.0*pow(n+1,2);
        if(i<(dim-n)) {
            tmp[i]+=x[i+n]*-1.0*pow(n+1,2);
            tmp[i+n]+=x[i]*-1.0*pow(n+1,2);
        }
        if(i%n!=0) {
            tmp[i]+=x[i-1]*-1.0*pow(n+1,2);
            tmp[i-1]+=x[i]*-1.0*pow(n+1,2);
        }
    }
    return tmp;
}

class LowerMatrix : public WriteableMatrix {
    private:
        int dim;
        int n;
    public:
        LowerMatrix(int);
        ~LowerMatrix();
        vector<double> diagonal;
        vector<double> tridiagonal;
        vector<double> identity;
        int Size();
        double Get(int,int);
        void Set(int,int,double);
        vector<double> operator*(const vector<double>&);
};

LowerMatrix::LowerMatrix(int n) {
    this->n=n;
    this->dim=n*n;
    this->diagonal.assign(dim,1);
    this->tridiagonal.assign(dim-1,0);
    this->identity.assign(dim-n,0);
}

LowerMatrix::~LowerMatrix() {
    vector<double>().swap(diagonal);
    vector<double>().swap(tridiagonal);
    vector<double>().swap(identity);
}

int LowerMatrix::Size() {
    return dim;
}

void LowerMatrix::Set(int i,int j,double value) {
    if(i==j) {
        diagonal[i]=value;
    } else if(j==(i-1)) {
        tridiagonal[j]=value;
    } else if(j==(i-n)) {
        identity[j]=value;
    }
}

double LowerMatrix::Get(int i,int j) {
    if(i==j) {
        return diagonal[i];
    } else if(j==(i-1)) {
        return tridiagonal[j];
    } else if(j==(i-n)) {
        return identity[j];
    } else {
        return 0.0;
    }
}

std::vector<double> LowerMatrix::operator*(const std::vector<double>& z) {
    dim=z.size();
    std::vector<double> tmp(z);
    for(int i=0;i<dim;i++) {
        if(i>=n) {
            tmp[i]-=tmp[i-n]*Get(i,i-n);
        }
        if(i%n!=0) {
            tmp[i]-=tmp[i-1]*Get(i,i-1);
        }
        tmp[i]/=Get(i,i);
    }
    return tmp;
}

class UpperMatrix : public LowerMatrix {
   private:
      int dim;
      int n;
   public:
      UpperMatrix(int);
      ~UpperMatrix();
        double Get(int,int);
      void Set(int,int,double);
        vector<double> operator*(const vector<double>&);
};

UpperMatrix::UpperMatrix(int n) : LowerMatrix(n) {
   this->n=n;
   this->dim=n*n;
   this->diagonal.assign(dim,0);
   this->tridiagonal.assign(dim-1,0);
   this->identity.assign(dim-n,0);
}

UpperMatrix::~UpperMatrix() {
   vector<double>().swap(diagonal);
   vector<double>().swap(tridiagonal);
   vector<double>().swap(identity);
}

double UpperMatrix::Get(int i,int j) {
   return LowerMatrix::Get(j,i);
}

void UpperMatrix::Set(int i,int j,double value) {
   return LowerMatrix::Set(j,i,value);
}

std::vector<double> UpperMatrix::operator*(const std::vector<double>& z) {
    std::vector<double> tmp(z);
    for(int i=dim-1;i>=0;i--) {
        tmp[i]-=tmp[i]*Get(i,i);
        if(i<=(dim-n)) {
            tmp[i]-=tmp[i+n]*Get(i,i+n);
        }
        if(i%n!=0) {
            tmp[i]-=tmp[i+1]*Get(i,i+1);
        }
        tmp[i]/=Get(i,i);
    }
    return tmp;
}

class Vectors {
    public:
        virtual double Get(int)=0;
        virtual int Size()=0;
        void WriteToFile();
        double f(double,double,int);
        double g(double,double,int);
};

double Vectors::f(double x,double y,int k) {
    double val;
    if(k==1) val=-4.0;
    if(k==2) val=0.0;
    return val;
}

double Vectors::g(double x,double y,int k) {
    double val;
    if(k==1) val=pow(x,2)+pow(y,2);
    if(k==2) val=1.0;
    return val;
}

void Vectors::WriteToFile() {
    FILE *file;
    file=fopen("../Plot/plot.dat","w");
    if(file==NULL)
        printf("<td colspan=2>ERROR: Could not open file!</td>");
    else {
        for(int i=0,k=0;i<=sqrt(Size());i++) {
            for(int j=0;j<=sqrt(Size());j++) {
                if(i==0) {
                    continue;
                } else if(i!=0) {
                    if(j==0) {
                        continue;
                    } else if(j!=0){
                        fprintf(file,"%f %f %f\n",(double)j/(double)(sqrt(Size())),(double)i/(double)(sqrt(Size())),Get(k));
                        k++;
                    }
                }
            }
            fprintf(file,"\n");
        }
    }
    fclose (file);
}

class Boundary : public Vectors {
    private:
        int dim;
        int n;
        double h;
        int k;
    public:
        Boundary(int,int);
        ~Boundary();
        vector<double> b;
        vector<double> solved;
        double Get(int);
        int Size();
};

Boundary::Boundary(int n,int k) {
    this->n=n;
    this->dim=n*n;
    this->h=1.0/(double)(n+1);
    this->k=k;
    this->b.resize(dim);
    this->solved.resize(dim);

    for(int i=1,k=0;i<=n;i++) {
        for(int j=1;j<=n;j++,k++) {
            solved[k]=g(j*h,i*h,this->k);
            b[k]=f(i*h,j*h,this->k);
            if(i==1) b[k]+=pow(1.0/h,2)*g(j*h,0,this->k);
            if(i==n) b[k]+=pow(1.0/h,2)*g(j*h,1,this->k);
            if(j==1) b[k]+=pow(1.0/h,2)*g(0,i*h,this->k);
            if(j==n) b[k]+=pow(1.0/h,2)*g(1,i*h,this->k);
        }
    }
}

Boundary::~Boundary() {
    vector<double>().swap(b);
}

double Boundary::Get(int i) {
    return this->b[i];
}

int Boundary::Size() {
    return b.size();
}

class Startvector : public Vectors {
    private:
        int dim;
        int n;
        double value;
    public:
        Startvector(int,double);
        ~Startvector();
        vector<double> x;
        double Get(int);
        int Size();
};

Startvector::Startvector(int n,double value) {
    this->n=n;
    this->dim=n*n;
    this->value=value;
    this->x.assign(dim,value);
}

Startvector::~Startvector() {
    vector<double>().swap(x);
}

double Startvector::Get(int i) {
    return this->x[i];
}

int Startvector::Size() {
    return x.size();
}

class Algorithms {
   private:
        int n;
        int dim;
        double h;
        int Vcounter;
   public:
        Algorithms(int);
        ~Algorithms();
        int JacobiMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        int JacobiRelaxationMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        int GaussSeidelMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        int SORMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        int CG(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        int PCG(Matrix&,WriteableMatrix&,WriteableMatrix&,vector<double>&,const vector<double>&,const vector<double>&);
        void modifiedIncompleteLU(Matrix&, WriteableMatrix&, WriteableMatrix&);
        void incompleteLU(Matrix&,WriteableMatrix&,WriteableMatrix&);
        void LUsolverUpper(Matrix&,Matrix&,vector<double>&);

        void JacobiRelaxation(Matrix&,vector<double>&,const vector<double>&,int);
        void Restriction(const vector<double>&,vector<double>&,int);
        void Interpolation(const vector<double>&,vector<double>&,int);
        vector<double> Cycle(Matrix&,vector<double>&,const vector<double>&,int,int,Matrix&,WriteableMatrix&,WriteableMatrix&);
        int MultiGridMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&,int);
        void PCGdirect(Matrix&,WriteableMatrix&,WriteableMatrix&,vector<double>&,const vector<double>&);
        void CGdirect(Matrix&,vector<double>&,const vector<double>&);
};

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
                if(j==1) E[k]=1.0/1.0*(E[k+1]);//+g(0,i*h));
                if(j==n) E[k]=1.0/1.0*(E[k-1]);//+g(1,i*h));
            }
            if(i%2!=0 && j%2==0) {
                if(i!=1 && i!=n) E[k]=1.0/2.0*(E[k-n]+E[k+n]);
                if(i==1) E[k]=1.0/1.0*(E[k+n]);//+g(j*h,0));
                if(i==n) E[k]=1.0/1.0*(E[k-n]);//+g(j*h,1));
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

int main(void) {
    char *data;
    const char *method;
    int arg,alg,func,steps=0;
    print_header();
    data=getenv("QUERY_STRING");
    if(data==NULL)
        printf("<td colspan=2>Error! Error in passing data from form to script.</td>");
    else if(sscanf(data,"arg=%d&alg=%d&func=%d",&arg,&alg,&func)!=3)
        printf("<td colspan=2>Error! Invalid data. Data must be numeric.</td>");
    else {
        PoissonMatrix A(arg);
        LowerMatrix L(arg);
        UpperMatrix U(arg);
        Boundary B(arg,func);
        Startvector X(arg,0.0);
        Algorithms Run(arg);

        double timer,start=0.0,end=0.0;
        start=clock();
        
        if(alg==1) {
            method="Conjugate Gradient";
            if(arg<512) steps=Run.CG(A,X.x,B.b,B.solved);
            else steps=-1;
        }
        if(alg==2) {
            method="P. Conjugate Gradient (ICG)";
            if(arg<1024) {
                A.InitHashMatrix();
                Run.incompleteLU(A,L,U);
                steps=Run.PCG(A,L,U,X.x,B.b,B.solved);
            } else steps=-1;
        }
        if(alg==3) {
            method="P. Conjugate Gradient (MICG)";
            if(arg<1024) {
                A.InitHashMatrix();
                Run.modifiedIncompleteLU(A,L,U);
                steps=Run.PCG(A,L,U,X.x,B.b,B.solved);
            } else steps=-1;
        }
        if(alg==4) {
            method="Two Grid";
            if(arg<2048) steps=Run.MultiGridMethod(A,X.x,B.b,B.solved,alg);
            else steps=-1;
        }
        if(alg==5) {
            method="V-Cycle";
            if(arg<2048) steps=Run.MultiGridMethod(A,X.x,B.b,B.solved,alg);
            else steps=-1;
        }
        if(alg==6) {
            method="W-Cycle";
            if(arg<2048) steps=Run.MultiGridMethod(A,X.x,B.b,B.solved,alg);
            else steps=-1;
        }
        if(alg==7) {
            method="Jacobi Method";
            if(arg<128) steps=Run.JacobiMethod(A,X.x,B.b,B.solved);
            else steps=-1;
        }
        if(alg==8) {
            method="Jacobi Relaxation Method";
            if(arg<128) steps=Run.JacobiRelaxationMethod(A,X.x,B.b,B.solved);
            else steps=-1;
        }
        if(alg==9) {
            method="Gauss-Seidel-Method";
            if(arg<128) steps=Run.GaussSeidelMethod(A,X.x,B.b,B.solved);
            else steps=-1;
        }
        if(alg==10) {
            method="SOR Method";
            if(arg<256) steps=Run.SORMethod(A,X.x,B.b,B.solved);
            else steps=-1;
        }

        end=clock();
        timer=(end-start)/CLOCKS_PER_SEC;
        print_table(arg,steps,timer,method);
        X.WriteToFile();
    }
    return EXIT_SUCCESS;
}