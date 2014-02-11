#include "classes.h"

int main(int argc, char const *argv[]) {
    int arg;
    if (argc==2) arg=atoi(&*argv[1])-1;
    else arg=7;

    Matrix A(arg,4.0,-1.0,-1.0);
    Vector V(arg);
    Algorithms Run(arg);

    printf("Started\n");
    double time,start=0.0,end=0.0;
    start=clock();

    Run.MultiGridMethod(A,V.x,V.b,V.solved);
    // int steps=0;
    // vector<double> r(arg*arg);
    // r=V.x-V.solved;
    // double TOL=(r|r)*pow(10,-3);
    // while(TOL<(r|r)) {
    //     Run.TwoGrid(A,V.x,V.b);
    //     r=V.x-V.solved;
    //     steps++;
    // }
    // printf("%d\n", steps);

    // vector<double> E(49,0),E2h(9,0);
    // E2h=V.b;
    // Run.Interpolation(E2h,E,7);

    // for(int i=0,k=0;i<3;i++) {
    //     for(int j=0;j<3;j++,k++) {
    //         if(E2h[k]>0) printf("%f  ", E2h[k]);
    //         else printf("%f ", E2h[k]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // printf("\n");

    // for(int i=0,k=0;i<7;i++) {
    //     for(int j=0;j<7;j++,k++) {
    //         if(E[k]>0) printf("%f  ", E[k]);
    //         else printf("%f ", E[k]);
    //     }
    //     printf("\n");
    // }



    // vector<double> r(49),r2h(9,0);
    // r=V.b;
    // Run.Restriction(r,r2h,7);
    // // V.PrintVector(r2h);

    // for(int i=0,k=0;i<3;i++) {
    //     for(int j=0;j<3;j++,k++) {
    //         if(r2h[k]>0) printf("%f  ", r2h[k]);
    //         else printf("%f ", r2h[k]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // printf("\n");

    // for(int i=0,k=0;i<7;i++) {
    //     for(int j=0;j<7;j++,k++) {
    //         if(r[k]>0) printf("%f  ", r[k]);
    //         else printf("%f ", r[k]);
    //     }
    //     printf("\n");
    // }

    end=clock();
    time=(end-start)/CLOCKS_PER_SEC;
    printf("%f\n", time);

    if(arg<=5) {
        A.PrintMatrix();
        V.PrintVector(V.x);
    } else {
        printf("Couldn't print Matrix! Dimension is too high.\n");
    }
    V.WriteToFile(V.x);

    return EXIT_SUCCESS;
}

double f(double x,double y) {
    return -4.0;
}

double g(double x,double y) {
    return pow(x,2)+pow(y,2);
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