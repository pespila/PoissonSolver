#include "classes.h"

int main(int argc, char const *argv[]) {
    int arg;
    if (argc==2) {
        arg=atoi(&*argv[1])-1;
    } else {
        arg=31;
    }

    PoissonMatrix2D A(arg,4.0,-1.0,-1.0);
    PoissonVector2D V(arg);
    Algorithms Run(arg);
    vector<double> x(arg*arg,0),b(arg*arg,0);
    V.InitRightSide(b);

    printf("Started\n");
    double time,start=0.0,end=0.0;
    start=clock();

    Run.MultiGridMethod(A,x,b);
    
    end=clock();
    time=(end-start)/CLOCKS_PER_SEC;
    printf("%f\n", time);

    if(arg<=5) {
        A.PrintMatrix();
        V.PrintVector(x);
    } else {
        printf("Couldn't print Matrix! Dimension is too high.\n");
    }
    V.WriteToFile(x);

    return EXIT_SUCCESS;
}

double f(double x,double y) {
    return -4.0;
}

double g(double x,double y) {
    return pow(x,2)+pow(y,2);
}

double f(double x) {
    return -2.0;
}

double g(double x) {
    return pow(x,2);
}

vector<double> operator-(const vector<double>& lhs, const vector<double>& rhs) {
    vector<double> tmp(lhs);
    for(int i=0;i<(int)lhs.size();i++) {
        tmp[i]=lhs[i]-rhs[i];
    }
    //tmp.insert(tmp.end(),rhs.begin(),rhs.end());
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