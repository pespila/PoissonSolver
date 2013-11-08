#include "classes.h"

Algorithms::Algorithms(Matrix& M) {
    dim = M.Size();
    n = sqrt(dim);
}

Algorithms::~Algorithms() {
}

int Algorithms::HashTable(int i, int j) {
    if (i == j) {
        return 1;
    } else if (j == (i+1) && i%n != (n-1)) {
        return 1;
    } else if (j == (i-1) && i%n != 0) {
        return 1;
    } else if (j == (i+n)) {
        return 1;
    } else if (j == (i-n)) {
        return 1;
    } else {
        return 0;
    }
}

void Algorithms::modifiedIncompleteLU(Matrix& A, WriteableMatrix& L, WriteableMatrix& U) {
    int i,j,k;
    double s, drop;

    for(i=0;i<dim;i++) {
        drop = 0;
        for(k=0;k<i;k++) {
            s = 0;
            for(j=0;j<k;j++) {
                if(HashTable(i,j) && HashTable(j,k))
                    s += L.Get(i,j) * U.Get(j,k);
            }
            if(HashTable(i,k)) {
                L.Set(i,k,((A.Get(i,k) - s)/U.Get(k,k)));
            } else {
                drop += s;
            }
        }
        for(k=i;k<dim;k++) {
            s = 0;
            for(j=0;j<i;j++) {
                if(HashTable(i,j) && HashTable(j,k))
                    s += L.Get(i,j) * U.Get(j,k);
            }
            if(HashTable(i,k)) {
                U.Set(i,k,(A.Get(i,k) - s));
            } else {
                drop += s;
            }
        }
        U.Set(i,i,(U.Get(i,i) - drop));
    }
}

void Algorithms::incompleteLU(Matrix& A, WriteableMatrix& L, WriteableMatrix& U) {
    int i,j,k;
    double s;

    for(i=0;i<dim;i++) {
        for(k=1;k<i;k++) {
            if (HashTable(i,k)) {
                s=0;
                for(j=0;j<k;j++) {
                    if(HashTable(i,j) && HashTable(j,k))
                        s+=L.Get(i,j)*U.Get(j,k);
                }
                L.Set(i,k,((A.Get(i,k)-s)/U.Get(k,k)));
            }
        }
        for(k=i;k<dim;k++) {
            if (HashTable(i,k) != 0) {
                s=0;
                for(j=0;j<i;j++) {
                    if(HashTable(i,j) && HashTable(j,k))
                        s+=L.Get(i,j)*U.Get(j,k);
                }
                U.Set(i,k,(A.Get(i,k)-s));
            }
        }
    }
}

// void Algorithms::modifiedIncompleteCholesky(PoissonMatrix& A, LowerPoissonMatrix& L, UpperPoissonMatrix& Ltranspose) {
//     int i,j,k;
//     double sum, tmp;

//     for(k=0;k<dim;k++) {
//         sum=0;
//         for(j=0;j<k;j++) {
//             if (A.Get(k,j) != 0) {
//                 sum+=L.Get(k,j);
//             }
//         }
//         tmp = sqrt(A.Get(k,k)-sum);
//         L.Set(k,k,tmp);
//         Ltranspose.Set(k,k,tmp); //neu
//         for(i=k+1;i<dim;i++) {
//             if (A.Get(i,k) != 0) {
//                 sum=0;
//                 for(j=0;j<k;j++) {
//                     if ((A.Get(i,j) != 0) && (A.Get(k,j) != 0)) {
//                         sum+=L.Get(i,j)*L.Get(k,j);
//                     }
//                 }
//                 tmp = (A.Get(i,k) - sum)/L.Get(k,k);
//                 L.Set(i,k,tmp);
//                 Ltranspose.Set(k,i,tmp); //neu
//             } else {
//                 sum=0;
//                 for(j=0;j<k;j++) {
//                     if ((A.Get(i,j) != 0) && (A.Get(k,j) != 0)) {
//                         sum+=L.Get(i,j)*L.Get(k,j);
//                     }
//                 }
//                 tmp = A.Get(i,i) - sum;
//                 A.Set(i,i,tmp);
//             }
//         }
//     }
// }

// void Algorithms::incompleteCholesky(PoissonMatrix& A, LowerPoissonMatrix& L, UpperPoissonMatrix& Ltranspose) {
//     int i,j,k;
//     double sum, tmp;

//     for(k=0;k<dim;k++) {
//         sum=0;
//         for(j=0;j<k;j++) {
//             if(A.Get(k,j) != 0) {
//                 sum+=L.Get(k,j)*L.Get(k,j);
//             }
//         }
//         tmp = sqrt(A.Get(k,k)-sum);
//         L.Set(k,k,tmp);
//         Ltranspose.Set(k,k,tmp); //neu
//         for(i=k+1;i<dim;i++) {
//             if (A.Get(i,k) != 0) {
//                 sum=0;
//                 for(j=0;j<k;j++) {
//                     if(A.Get(i,j) != 0 && A.Get(k,j) != 0) {
//                         sum+=L.Get(i,j)*L.Get(k,j);
//                     }
//                 }
//                 tmp = (A.Get(i,k) - sum)/L.Get(k,k);
//                 L.Set(i,k,tmp);
//                 Ltranspose.Set(k,i,tmp);
//             }
//         }
//     }
// }

void Algorithms::CG(Matrix& A, Operators& O, CGVectors& V) {
    double alpha,beta,eps,h,norm,skpOfRes,denom;
    int steps,i,j,k;

    steps=0;
    h=1.0/(double)(n+1);
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
            solution[k] = O.g(j*h,i*h);
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

void Algorithms::PCG(Matrix& A, Operators& O, WriteableMatrix& L, WriteableMatrix& U, CGVectors& V) {
    double alpha,beta,eps,h,norm,skpOfRes,denom;
    int steps,i,j,k;

    steps=0;
    h=1.0/(double)(n+1);

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

    O.LUsolverLower(L,z);
    O.LUsolverUpper(U,z);

    p=z;
    rTmp=r;
    zTmp=z;

    k=0;
    for (i=1;i<(n+1);i++) {
        for (j=1;j<(n+1);j++) {
            solution[k] = O.g(j*h,i*h);
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
        O.LUsolverLower(L,z);
        O.LUsolverUpper(U,z);
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
        if (steps > 100)
            break;
    }
    printf("PCGSteps: %d\n", steps);
}