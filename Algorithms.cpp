#include "classes.h"

Algorithms::Algorithms(int n) {
    this->n=n;
    this->dim=n*n;
    this->h=1.0/(double)(n+1);
}

Algorithms::~Algorithms() {
}

void Algorithms::modifiedIncompleteLU(Matrix& A,WriteableMatrix& L,WriteableMatrix& U) {
    int i,j,k,m,u;
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
            }
        }
        for(k=0;k<5;k++) {
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
        U.Set(i,i,(U.Get(i,i)-drop));
    }
}

//Mit Residuum
// void Algorithms::JacobiMethod(Matrix& A,Operators& O,Vectors& V) {
//     double eps,norm,sum;
//     int steps=0;

//     vector<double> Ax,r,tmp;
//     r.resize(dim);
//     Ax.assign(dim,0);
//     tmp.resize(dim);
//     O.MatrixVectorMultiplyer(A,V.x,Ax);

//     for(int i=0;i<dim;i++) {
//         r[i]=V.b[i]-Ax[i];
//     }

//     eps=pow(10,-6)*O.vectorNorm(r);
//     norm=1.0;

//     while(eps<=norm) {
//         tmp=V.x;
//         Ax.assign(dim,0);
//         O.MatrixVectorMultiplyer(A,V.x,Ax);
//         for(int i=0;i<dim;i++) {
//             if(i==0) {
//                 sum=tmp[i+1]*A.Get(i,i+1)+tmp[i+n]*A.Get(i,i+n);
//             } else if(i!=0 && i<n) {
//                 sum=tmp[i-1]*A.Get(i,i-1)+tmp[i+1]*A.Get(i,i+1)+tmp[i+n]*A.Get(i,i+n);
//             } else if(i>=n && i<dim-n) {
//                 sum=tmp[i-n]*A.Get(i,i-n)+tmp[i-1]*A.Get(i,i-1)+tmp[i+1]*A.Get(i,i+1)+tmp[i+n]*A.Get(i,i+n);
//             } else if(i!=dim-1 && i>=dim-n) {
//                 sum=tmp[i-n]*A.Get(i,i-n)+tmp[i-1]*A.Get(i,i-1)+tmp[i+1]*A.Get(i,i+1);
//             } else if(i==dim-1) {
//                 sum=tmp[i-n]*A.Get(i,i-n)+tmp[i-1]*A.Get(i,i-1);
//             }
//             V.x[i]=1.0/A.Get(i,i)*(V.b[i]-sum);
//             r[i]=V.b[i]-Ax[i];
//         }
//         norm=O.vectorNorm(r);
//         steps++;
//     }
//     printf("JacobianSteps: %d\n", steps);
// }

void Algorithms::JacobiMethod(Matrix& A,Operators& O,Vectors& V) {
    double eps,norm,sum;
    int steps=0;

    vector<double> solution;
    vector<double> stopNorm;
    solution.resize(dim);
    stopNorm.resize(dim);
    vector<double> tmp;
    tmp.resize(dim);

    int k=0;
    for (int i=1;i<(n+1);i++) {
        for (int j=1;j<(n+1);j++) {
            solution[k] = g(j*h,i*h);
            stopNorm[k] = V.x[k]-solution[k];
            k++;
        }
    }

    eps=pow(10,-3)*O.vectorNorm(stopNorm);
    norm=1.0;

    while(eps<=norm) {
        tmp=V.x;
        for(int i=0;i<dim;i++) {
            if(i==0) {
                sum=tmp[i+1]*A.Get(i,i+1)+tmp[i+n]*A.Get(i,i+n);
            } else if(i!=0 && i<n) {
                sum=tmp[i-1]*A.Get(i,i-1)+tmp[i+1]*A.Get(i,i+1)+tmp[i+n]*A.Get(i,i+n);
            } else if(i>=n && i<dim-n) {
                sum=tmp[i-n]*A.Get(i,i-n)+tmp[i-1]*A.Get(i,i-1)+tmp[i+1]*A.Get(i,i+1)+tmp[i+n]*A.Get(i,i+n);
            } else if(i!=dim-1 && i>=dim-n) {
                sum=tmp[i-n]*A.Get(i,i-n)+tmp[i-1]*A.Get(i,i-1)+tmp[i+1]*A.Get(i,i+1);
            } else if(i==dim-1) {
                sum=tmp[i-n]*A.Get(i,i-n)+tmp[i-1]*A.Get(i,i-1);
            }
            V.x[i]=1.0/A.Get(i,i)*(V.b[i]-sum);
            stopNorm[i] = V.x[i]-solution[i];
        }
        norm=O.vectorNorm(stopNorm);
        steps++;
    }
    printf("JacobianSteps: %d\n", steps);
}

void Algorithms::GaussSeidelMethod(Matrix& A,Operators& O,Vectors& V) {
    double eps,norm,sum;
    int steps=0;

    vector<double> solution;
    vector<double> stopNorm;
    solution.resize(dim);
    stopNorm.resize(dim);
    vector<double> tmp;
    tmp.resize(dim);

    int k=0;
    for (int i=1;i<(n+1);i++) {
        for (int j=1;j<(n+1);j++) {
            solution[k] = g(j*h,i*h);
            stopNorm[k] = V.x[k]-solution[k];
            k++;
        }
    }

    eps=pow(10,-3)*O.vectorNorm(stopNorm);
    norm=1.0;

    while(eps<=norm) {
        tmp=V.x;
        for(int i=0;i<dim;i++) {
            if(i==0) {
                sum=V.x[i+1]*A.Get(i,i+1)+V.x[i+n]*A.Get(i,i+n);
            } else if(i!=0 && i<n) {
                sum=V.x[i-1]*A.Get(i,i-1)+V.x[i+1]*A.Get(i,i+1)+V.x[i+n]*A.Get(i,i+n);
            } else if(i>=n && i<dim-n) {
                sum=V.x[i-n]*A.Get(i,i-n)+V.x[i-1]*A.Get(i,i-1)+V.x[i+1]*A.Get(i,i+1)+V.x[i+n]*A.Get(i,i+n);
            } else if(i!=dim-1 && i>=dim-n) {
                sum=V.x[i-n]*A.Get(i,i-n)+V.x[i-1]*A.Get(i,i-1)+V.x[i+1]*A.Get(i,i+1);
            } else if(i==dim-1) {
                sum=V.x[i-n]*A.Get(i,i-n)+V.x[i-1]*A.Get(i,i-1);
            }
            V.x[i]=1.0/A.Get(i,i)*(V.b[i]-sum);
            stopNorm[i] = V.x[i]-solution[i];
        }
        norm=O.vectorNorm(stopNorm);
        steps++;
    }
    printf("GaussSeidelSteps: %d\n", steps);
}

void Algorithms::SORMethod(Matrix& A,Operators& O,Vectors& V) {
    double eps,norm,sum,Pi=3.141592654,omega=2/(1+sqrt(1-pow(cos(Pi*h),2)));
    int steps=0;

    vector<double> solution;
    vector<double> stopNorm;
    solution.resize(dim);
    stopNorm.resize(dim);
    vector<double> tmp;
    tmp.resize(dim);

    int k=0;
    for (int i=1;i<(n+1);i++) {
        for (int j=1;j<(n+1);j++) {
            solution[k] = g(j*h,i*h);
            stopNorm[k] = V.x[k]-solution[k];
            k++;
        }
    }

    eps=pow(10,-3)*O.vectorNorm(stopNorm);
    norm=1.0;

    while(eps<=norm) {
        tmp=V.x;
        for(int i=0;i<dim;i++) {
            if(i==0) {
                sum=V.x[i+1]*A.Get(i,i+1)+V.x[i+n]*A.Get(i,i+n);
            } else if(i!=0 && i<n) {
                sum=V.x[i-1]*A.Get(i,i-1)+V.x[i+1]*A.Get(i,i+1)+V.x[i+n]*A.Get(i,i+n);
            } else if(i>=n && i<dim-n) {
                sum=V.x[i-n]*A.Get(i,i-n)+V.x[i-1]*A.Get(i,i-1)+V.x[i+1]*A.Get(i,i+1)+V.x[i+n]*A.Get(i,i+n);
            } else if(i!=dim-1 && i>=dim-n) {
                sum=V.x[i-n]*A.Get(i,i-n)+V.x[i-1]*A.Get(i,i-1)+V.x[i+1]*A.Get(i,i+1);
            } else if(i==dim-1) {
                sum=V.x[i-n]*A.Get(i,i-n)+V.x[i-1]*A.Get(i,i-1);
            }
            V.x[i]=omega*1.0/A.Get(i,i)*(V.b[i]-sum)+(1-omega)*V.x[i];
            stopNorm[i] = V.x[i]-solution[i];
        }
        norm=O.vectorNorm(stopNorm);
        steps++;
    }
    printf("GaussSeidelSteps: %d\n", steps);
}

void Algorithms::CG(Matrix& A, Operators& O, Vectors& V) {
    double alpha,beta,eps,norm,skpOfRes,denom;
    int steps,i,j,k;

    steps=0;
    beta=0;

    vector<double> r;
    vector<double> solution;
    vector<double> stopNorm;
    vector<double> p;
    vector<double> Ap;
    vector<double> rTmp;
    r.resize(dim);
    solution.resize(dim);
    stopNorm.resize(dim);
    p.resize(dim);
    rTmp.resize(dim);

    vector<double> Ax;
    Ax.assign(dim,0);
    V.x.assign(dim,0);
    O.MatrixVectorMultiplyer(A,V.x,Ax);
    for(i=0;i<dim;i++)
        r[i]=V.b[i]-Ax[i];
    Ax.clear();

    k=0;
    for (i=1;i<(n+1);i++) {
        for (j=1;j<(n+1);j++) {
            solution[k] = g(j*h,i*h);
            stopNorm[k] = V.x[k]-solution[k];
            k++;
        }
    }

    eps=pow(10,-3)*O.vectorNorm(stopNorm);
    norm=1.0;

    p=r;
    rTmp=r;

    while (eps <= norm) {
        steps++;

        Ap.assign(dim,0);
        O.MatrixVectorMultiplyer(A,p,Ap);

        skpOfRes=O.innerProduct(rTmp,rTmp);
        denom = O.innerProduct(p,Ap);
        alpha=skpOfRes/denom;

        for(i=0;i<dim;i++) {
            V.x[i] += alpha*p[i];
            r[i] -= alpha*Ap[i];
        }

        denom = O.innerProduct(r,r);
        beta = denom/skpOfRes;

        for(i=0;i<dim;i++) {
            p[i]=r[i]+beta*p[i];
            stopNorm[i]=V.x[i]-solution[i];
        }
        rTmp = r;
        Ap.clear();
        norm=O.vectorNorm(stopNorm);
    }
    printf("CGSteps: %d\n", steps);
}

void Algorithms::PCG(Matrix& A, Operators& O, WriteableMatrix& L, WriteableMatrix& U, Vectors& V) {
    double alpha,beta,eps,norm,skpOfRes,denom;
    int steps,i,j,k;

    steps=0;

    vector<double> r;
    vector<double> z;
    vector<double> solution;
    vector<double> stopNorm;
    vector<double> p;
    vector<double> Ap;
    vector<double> rTmp;
    vector<double> zTmp;
    r.resize(dim);
    z.resize(dim);
    solution.resize(dim);
    stopNorm.resize(dim);
    p.resize(dim);
    rTmp.resize(dim);
    zTmp.resize(dim);

    vector<double> Ax;
    Ax.assign(dim,0);
    V.x.assign(dim,0);
    O.MatrixVectorMultiplyer(A,V.x,Ax);
    for(i=0;i<dim;i++)
        r[i]=V.b[i]-Ax[i];

    z=r;

    O.LUsolverLower(A,L,z);
    O.LUsolverUpper(A,U,z);

    p=z;
    rTmp=r;
    zTmp=z;

    k=0;
    for (i=1;i<(n+1);i++) {
        for (j=1;j<(n+1);j++) {
            solution[k] = g(j*h,i*h);
            stopNorm[k] = V.x[k]-solution[k];
            k++;
        }
    }

    eps=pow(10,-3)*O.vectorNorm(stopNorm);
    norm=1.0;

    while (eps <= norm) {
        steps++;

        Ap.assign(dim,0);
        O.MatrixVectorMultiplyer(A,p,Ap);

        skpOfRes=O.innerProduct(zTmp,rTmp);
        denom = O.innerProduct(p,Ap);
        alpha=skpOfRes/denom;

        for(i=0;i<dim;i++) {
            V.x[i] += alpha*p[i];
            r[i] -= alpha*Ap[i];
        }

        z=r;
        O.LUsolverLower(A,L,z);
        O.LUsolverUpper(A,U,z);
        zTmp=z;

        denom = O.innerProduct(z,r);
        beta = denom/skpOfRes;

        for(i=0;i<dim;i++) {
            p[i]=z[i]+beta*p[i];
            stopNorm[i]=V.x[i]-solution[i];
        }
        rTmp = r;
        Ap.clear();
        norm=O.vectorNorm(stopNorm);
    }
    printf("PCGSteps: %d\n", steps);
}