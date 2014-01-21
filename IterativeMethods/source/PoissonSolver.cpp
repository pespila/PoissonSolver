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
    Vectors V(arg);
    Algorithms Run(arg);

    printf("Choose your Algorithm: \n[1] CG\n[2] ICCG\n[3] MICCG\n[4] JacobiMethod\n[5] JacobiRelaxationMethod\n[6] GaussSeidelMethod\n[7] SORMethod\n");
    scanf("%d",&alg);

    printf("Started\n");
    double time,start=0.0,end=0.0;
    start = clock();

    if(alg==1) Run.CG(A,V.x,V.b);
    if(alg==2) {
        Run.incompleteLU(A,L,U);
        Run.PCG(A,L,U,V.x,V.b);
    }
    if(alg==3) {
        Run.modifiedIncompleteLU(A,L,U);
        Run.PCG(A,L,U,V.x,V.b);
    }
    if(alg==4) Run.JacobiMethod(A,V.x,V.b);
    if(alg==5) Run.JacobiRelaxationMethod(A,V.x,V.b);
    if(alg==6) Run.GaussSeidelMethod(A,V.x,V.b);
    if(alg==7) Run.SORMethod(A,V.x,V.b);

    end = clock();
    time=(end-start)/CLOCKS_PER_SEC;
    printf("%f\n", time);

    if(arg<=5) {
        A.PrintMatrix();
        L.PrintMatrix();
        U.PrintMatrix();
        V.PrintVector();
    }
    V.WriteToFile();

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
    //tmp.insert(tmp.end(),rhs.begin(),rhs.end());
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