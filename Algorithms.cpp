#include "classes.h"

Algorithms::Algorithms(PoissonMatrix& Mat) {
    dim = Mat.size();
    n = sqrt(dim);
}

Algorithms::~Algorithms() {
}

void Algorithms::modifiedIncompleteLU(PoissonMatrix& A, LowerPoissonMatrix& L, UpperPoissonMatrix& U) {
    int i,j,k;
    double s, tmp, drop;

    for(i=0;i<dim;i++) {
        drop = 0;
        for(k=0;k<i;k++) {
            s = 0;
            for(j=0;j<k;j++) {
                s += L.get(i,j) * U.get(j,k);
            }
            if(A.get(i,k) != 0) {
                tmp = (A.get(i,k) - s)/U.get(k,k);
                L.set(i,k,tmp);
            } else {
                drop += s;
            }
        }
        for(k=i;k<dim;k++) {
            s = 0;
            for(j=0;j<i;j++) {
                s += L.get(i,j) * U.get(j,k);
            }
            if(A.get(i,k) != 0) {
                tmp = A.get(i,k) - s;
                U.set(i,k,tmp);
            } else {
                drop += s;
            }
        }
        U.set(i,i,(U.get(i,i) - drop));
    }
}

void Algorithms::modifiedIncompleteCholesky(PoissonMatrix& A, LowerPoissonMatrix& L, UpperPoissonMatrix& Ltranspose) {
    int i,j,k;
    double sum, tmp;

    for(k=0;k<dim;k++) {
        sum=0;
        for(j=0;j<k;j++) {
            if (A.get(k,j) != 0) {
                sum+=L.get(k,j);
            }
        }
        tmp = sqrt(A.get(k,k)-sum);
        L.set(k,k,tmp);
        Ltranspose.set(k,k,tmp); //neu
        for(i=k+1;i<dim;i++) {
            if (A.get(i,k) != 0) {
                sum=0;
                for(j=0;j<k;j++) {
                    if ((A.get(i,j) != 0) && (A.get(k,j) != 0)) {
                        sum+=L.get(i,j)*L.get(k,j);
                    }
                }
                tmp = (A.get(i,k) - sum)/L.get(k,k);
                L.set(i,k,tmp);
                Ltranspose.set(k,i,tmp); //neu
            } else {
                sum=0;
                for(j=0;j<k;j++) {
                    if ((A.get(i,j) != 0) && (A.get(k,j) != 0)) {
                        sum+=L.get(i,j)*L.get(k,j);
                    }
                }
                tmp = A.get(i,i) - sum;
                A.set(i,i,tmp);
            }
        }
    }
}

void Algorithms::incompleteLU(PoissonMatrix& A, LowerPoissonMatrix& L, UpperPoissonMatrix& U) {
    int i,j,k;
    double s, tmp;

    for(i=0;i<dim;i++) {
        for(k=1;k<i;k++) {
            if (A.get(i,k) != 0) {
                s=0;
                for(j=0;j<k;j++) {
                    s+=L.get(i,j)*U.get(j,k);
                }
                tmp = (A.get(i,k)-s)/U.get(k,k);
                L.set(i,k,tmp);
            }
        }
        for(k=i;k<dim;k++) {
            if (A.get(i,k) != 0) {
                s=0;
                for(j=0;j<i;j++) {
                    s+=L.get(i,j)*U.get(j,k);
                }
                tmp = A.get(i,k)-s;
                U.set(i,k,tmp);
            }
        }
    }
}

void Algorithms::incompleteCholesky(PoissonMatrix& A, LowerPoissonMatrix& L, UpperPoissonMatrix& Ltranspose) {
    int i,j,k;
    double sum, tmp;

    for(k=0;k<dim;k++) {
        sum=0;
        for(j=0;j<k;j++) {
            if(A.get(k,j) != 0) {
                sum+=L.get(k,j)*L.get(k,j);
            }
        }
        tmp = sqrt(A.get(k,k)-sum);
        L.set(k,k,tmp);
        Ltranspose.set(k,k,tmp); //neu
        for(i=k+1;i<dim;i++) {
            if (A.get(i,k) != 0) {
                sum=0;
                for(j=0;j<k;j++) {
                    if(A.get(i,j) != 0 && A.get(k,j) != 0) {
                        sum+=L.get(i,j)*L.get(k,j);
                    }
                }
                tmp = (A.get(i,k) - sum)/L.get(k,k);
                L.set(i,k,tmp);
                Ltranspose.set(k,i,tmp);
            }
        }
    }
}

void Algorithms::CG(PoissonMatrix& A, Operators& O, vector<double>& x, const vector<double>& b) {
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
    x.assign(dim,0);
    O.MatrixVectorMultiplyer(A,x,Ax);
    for(i=0;i<dim;i++)
        r[i]=b[i]-Ax[i];
    Ax.clear();

    k=0;
    for (i=1;i<(n+1);i++) {
        for (j=1;j<(n+1);j++) {
            solution[k] = O.g(j*h,i*h);
            stopNorm[k] = x[k]-solution[k];
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
            x[i] += alpha*p[i];
            r[i] -= alpha*Ap[i];
        }

        denom = O.innerProduct(r,r);
        beta = denom/skpOfRes;

        for(i=0;i<dim;i++) {
            p[i]=r[i]+beta*p[i];
            stopNorm[i]=x[i]-solution[i];
        }
        rTmp = r;
        Ap.clear();
        norm=O.vectorNorm(stopNorm);
    }
    printf("CGSteps: %d\n", steps);
}

void Algorithms::PCG(PoissonMatrix& A, Operators& O, LowerPoissonMatrix& L, UpperPoissonMatrix& U, vector<double>& x, const vector<double>& b) {
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
    x.assign(dim,0);
    O.MatrixVectorMultiplyer(A,x,Ax);
    for(i=0;i<dim;i++)
        r[i]=b[i]-Ax[i];

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
            stopNorm[k] = x[k]-solution[k];
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
            x[i] += alpha*p[i];
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
            stopNorm[i]=x[i]-solution[i];
        }
        rTmp = r;
        Ap.clear();
        norm=O.vectorNorm(stopNorm);
        if (steps > 100)
            break;
    }
    printf("PCGSteps: %d\n", steps);
}