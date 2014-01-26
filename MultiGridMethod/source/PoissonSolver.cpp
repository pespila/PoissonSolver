#include "classes.h"

int main(int argc, char const *argv[]) {
    int arg;
    if (argc==2) {
        arg=atoi(&*argv[1])-1;
    } else {
        arg=31;
    }

    PoissonMatrix A(arg);
    Vectors V(arg);
    Algorithms Run(arg);

    printf("Started\n");
    double time,start=0.0,end=0.0;
    start=clock();

    // std::vector<double> E(arg*arg,0);
    // double h=1.0/(double)(arg+1);
    // for(int i=1,k=0;i<arg;i++) {
    //     for(int j=1;j<arg;j++,k++) {
    //         E[k]=V.x[k]-g(i*h,j*h);
    //     }
    // }
    // FILE *file;
    // file=fopen("../Plot/plot.txt", "w");
    // if(file==NULL) {
    //     printf("ERROR: Could not open file!\n");
    // } else {
    //     for(int i=1,k=0;i<arg;i++) {
    //         for(int j=1;j<arg;j++,k++) {
    //             fprintf(file, "%f %f %f\n", (double)j/(double)arg,(double)i/(double)arg,E[k]);
    //         }
    //         fprintf(file, "\n");
    //     }
    // }
    // fclose (file);
    Run.MultiGridMethod(A,V,V.x,V.b);
    // Run.JacobiMethod(A,V.x,V.b,100);
    // for(int i=1,k=0;i<arg;i++) {
    //     for(int j=1;j<arg;j++,k++) {
    //         E[k]=V.x[k]-g(i*h,j*h);
    //     }
    // }
    // FILE *file;
    // file=fopen("../Plot/plot.txt", "w");
    // if(file==NULL) {
    //     printf("ERROR: Could not open file!\n");
    // } else {
    //     for(int i=0;i<=arg;i++) {
    //             fprintf(file, "%d %f\n",i,sin(3.14*(double)i/(double)arg));
    //     }
    // }
    // fclose (file);

    end=clock();
    time=(end-start)/CLOCKS_PER_SEC;
    printf("%f\n", time);

    if(arg<=5) {
        A.PrintMatrix();
        V.PrintVector();
    } else {
        printf("Couldn't print Matrix! Dimension is too high.\n");
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

void operator+=(std::vector<double>& lhs, const std::vector<double>& rhs) {
    for(int i=0;i<(int)rhs.size();i++) {
        lhs[i]+=rhs[i];
    }
}