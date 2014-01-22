#include "classes.h"

int main(int argc, char const *argv[]) {
    int arg,alg;
    if (argc==2) {
        arg=atoi(&*argv[1])-1;
    } else {
        arg=39;
    }

    PoissonMatrix A(arg);
    LowerMatrix L(arg);
    UpperMatrix U(arg);
    Boundary B(arg);
    Startvector X(arg,0.0);
    Algorithms Run(arg);

    printf("Choose your Algorithm: \n[1] CG\n[2] ICCG\n[3] MICCG\n[4] JacobiMethod\n[5] JacobiRelaxationMethod\n[6] GaussSeidelMethod\n[7] SORMethod\n");
    scanf("%d",&alg);

    printf("Started\n");
    double time,start=0.0,end=0.0;
    start=clock();

    if(alg==1) Run.CG(A,X.x,B.b);
    if(alg==2) {
        Run.incompleteLU(A,L,U);
        Run.PCG(A,L,U,X.x,B.b);
    }
    if(alg==3) {
        Run.modifiedIncompleteLU(A,L,U);
        Run.PCG(A,L,U,X.x,B.b);
    }
    if(alg==4) Run.JacobiMethod(A,X.x,B.b);
    if(alg==5) Run.JacobiRelaxationMethod(A,X.x,B.b);
    if(alg==6) Run.GaussSeidelMethod(A,X.x,B.b);
    if(alg==7) Run.SORMethod(A,X.x,B.b);

    end=clock();
    time=(end-start)/CLOCKS_PER_SEC;
    printf("%f\n", time);

    if(arg<=5) {
        A.PrintMatrix();
        L.PrintMatrix();
        U.PrintMatrix();
        X.PrintVector();
    }
    X.WriteToFile();

    return EXIT_SUCCESS;
}

double f(double x,double y) {
    return -4.0;
}

double g(double x,double y) {
    return pow(x,2)+pow(y,2);
}

std::vector<double> operator-(const std::vector<double>& lhs, const std::vector<double>& rhs) {
    std::vector<double> tmp(lhs);
    for(int i=0;i<(int)lhs.size();i++) {
        tmp[i]=lhs[i]-rhs[i];
    }
    return tmp;
}

std::vector<double> operator*(double x, std::vector<double> rhs) {
    std::vector<double> tmp(rhs);
    for(int i=0;i<(int)rhs.size();i++) {
        tmp[i]=x*rhs[i];
    }
    return tmp;
}

void operator+=(std::vector<double>& lhs, const std::vector<double>& rhs) {
    for(int i=0;i<(int)rhs.size();i++) {
        lhs[i]+=rhs[i];
    }
}

double operator*(const std::vector<double>& x,const std::vector<double>& y) {
    double ip=0.0;
    for (int i=0;i<(int)x.size();i++)
        ip+=x[i]*y[i];
    return ip;
}

double operator|(const std::vector<double>& x,const std::vector<double>& y) {
    double norm=0.0;
    for(int i=0;i<(int)x.size();i++) {
        norm+=x[i]*y[i];
    }
    return sqrt(norm);
}