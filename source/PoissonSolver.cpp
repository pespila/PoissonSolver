#include "classes.h"

int main(int argc, char const *argv[]) {
    const char *method;
    int arg,alg,func,steps=0;
    if(argc<4) {
        printf("\n");
        printf("Choose your prefered algorithm:\n[1] Conjugate Gradient\n[2] P. Conjugate Gradient (ICG)\n[3] P.Conjugate Gradient (MICG)\n[4] Two Grid (square numbers of 2 recommended!)\n[5] V-Cycle (square numbers of 2 recommended!)\n[6] W-Cycle (square numbers of 2 recommended!)\n[7] Jacobi Method\n[8] Jacobi Relaxation Method \n[9] Gauss-Seidel-Method\n[10] SOR Method\n");
        scanf("%d",&alg);
        printf("Choose your m, where h=1/m and N=m-1 number of grid points in x and y direction:\n");
        scanf("%d",&arg);
        printf("Choose your prefered equation:\n[1] -nabla u(x,y) = -4, and u(x,y) = x^2 + y^2.\n[2] -nabla u(x,y) = 0, and u(x,y) = 1.\n");
        scanf("%d",&func);
        arg=arg-1;
    } else {
	    arg=atoi(&*argv[1])-1;alg=atoi(&*argv[2]);func=atoi(&*argv[3]);
	}

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
        if(arg<4096) steps=Run.MultiGridMethod(A,X.x,B.b,B.solved,alg);
        else steps=-1;
    }
    if(alg==5) {
        method="V-Cycle";
        if(arg<4096) steps=Run.MultiGridMethod(A,X.x,B.b,B.solved,alg);
        else steps=-1;
    }
    if(alg==6) {
        method="W-Cycle";
        if(arg<4096) steps=Run.MultiGridMethod(A,X.x,B.b,B.solved,alg);
        else steps=-1;
    }
    if(alg==7) {
        method="Jacobi Method";
        if(arg<256) steps=Run.JacobiMethod(A,X.x,B.b,B.solved);
        else steps=-1;
    }
    if(alg==8) {
        method="Jacobi Relaxation Method";
        if(arg<256) steps=Run.JacobiRelaxationMethod(A,X.x,B.b,B.solved);
        else steps=-1;
    }
    if(alg==9) {
        method="Gauss-Seidel-Method";
        if(arg<256) steps=Run.GaussSeidelMethod(A,X.x,B.b,B.solved);
        else steps=-1;
    }
    if(alg==10) {
        method="SOR Method";
        if(arg<512) steps=Run.SORMethod(A,X.x,B.b,B.solved);
        else steps=-1;
    }

    end=clock();
    timer=(end-start)/CLOCKS_PER_SEC;
    printf("%s with %d steps in %f seconds.\n", method,steps,timer);
    X.WriteToFile();

    return EXIT_SUCCESS;
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