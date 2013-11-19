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
                if(HashTable(i,j) && HashTable(j,k)) {
                    s += L.Get(i,j) * U.Get(j,k);
                    printf("Matrix L mit %.2f: %d,%d && %d,%d\n", s, i, j, j, k);
                }
            }
            //printf("in L   s_ik: %.2f_%d%d\n", s,i,k);
            if(HashTable(i,k)) {
                L.Set(i,k,((A.Get(i,k) - s)/U.Get(k,k)));
            } else {
                drop += s;
            }
        }
        //printf("in L   drop_ik: %.2f_%d%d\n", drop,i,k);
        for(k=i;k<dim;k++) {
            s = 0;
            for(j=0;j<i;j++) {
                if(HashTable(i,j) && HashTable(j,k)) {
                    s += L.Get(i,j) * U.Get(j,k);
                    printf("Matrix R mit %.2f: %d,%d && %d,%d\n", s, i, j, j, k);
                }
            }
            //printf("in R   s_ik: %.2f_%d%d\n", s,i,k);
            if(HashTable(i,k)) {
                U.Set(i,k,(A.Get(i,k) - s));
            } else {
                drop += s;
            }
        }
        printf("Drop: %.2f\n", drop);
        //printf("in R   drop_ik: %.2f_%d%d\n", 2*drop,i,k);
        U.Set(i,i,(U.Get(i,i) - drop));
        printf("\n");
    }
}

void Algorithms::LUNEW(Matrix& A, WriteableMatrix& L, WriteableMatrix& U) {
    int i;
    double sum, drop;
    int stop = A.HashKeys.size();

    for(i=0;i<stop;i++) {
        if(i==0){
            continue;
        } else if(i>=n&&HashTable(i,i+1)&&i!=(dim-1)) {
            printf("Calculate s(L_i,i-n+1)\n");
        } else if(HashTable(i,i)) {
            printf("Calculate s(U_i,i)\n");
        } else if(i<(dim-n)&&HashTable(i,i-1)&&i!=0) {
            printf("Calculate s(U_i,i+n-1\n");
        }

        if(i>n&&!HashTable(i,i+1)) {
            drop=0; //s(U_i,i) + s(U_i,i+n-1)
        } else if(i>(dim-n)&&i!=(dim-1)) {
            drop=0; //s(L_i,i-n+1)
        } else {
            drop=0; //s(U_i,i)
        }
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

void Algorithms::JacobiMethod(Matrix& A, Operators& O, Vectors& V) {
    double eps,h,norm,sum;
    int steps,i,j,k;

    steps=0;
    h=1.0/(double)(n+1);

    vector<double> solution;
    vector<double> stopNorm;
    solution.resize(dim);
    stopNorm.resize(dim);

    V.x.assign(dim,0);

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

    while(eps<norm) {
        for(i=0;i<dim;i++) {
            sum=0;
            for(j=0;j<dim;j++) {
                if(i==j) continue;
                if(HashTable(i,j)) {
                    sum+=A.Get(i,j)*V.x[j];
                } 
            }
                V.x[i]=1/A.Get(i,i)*(V.b[i]-sum);
        }

        for(i=0;i<dim;i++)
            stopNorm[i]=V.x[i]-solution[i];
        norm=O.vectorNorm(stopNorm);
        steps++;
    }
    printf("JacobianSteps: %d\n", steps);
}

void Algorithms::GaussSeidelMethod(Matrix& A, Operators& O, Vectors& V) {
    double eps,h,norm,sum1,sum2;
    int steps,i,j,k;

    steps=0;
    h=1.0/(double)(n+1);

    vector<double> solution;
    vector<double> stopNorm;
    solution.resize(dim);
    stopNorm.resize(dim);

    V.x.assign(dim,0);

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
    
    while(eps<norm) {
        for(i=0;i<dim;i++) {
            sum1=0;
            for(j=0;j<i;j++) {
                if(HashTable(i,j)) {
                    sum1+=A.Get(i,j)*V.x[j];
                }
            }
            sum2=0;
            for(j=i+1;j<dim;j++) {
                if(HashTable(i,j)) {
                    sum2+=A.Get(i,j)*V.x[j];
                }
            }
            V.x[i] = 1/A.Get(i,i)*(V.b[i]-sum1-sum2);
        }

        for(i=0;i<dim;i++)
            stopNorm[i]=V.x[i]-solution[i];
        norm=O.vectorNorm(stopNorm);
        steps++;
    }

    printf("GaussSeidelSteps: %d\n", steps);
}

void Algorithms::SORMethod(Matrix& A, Operators& O, Vectors& V) {
    double eps,h,norm,sum1,sum2,omega;
    int steps,i,j,k;

    steps=0;
    h=1.0/(double)(n+1);

    vector<double> solution;
    vector<double> stopNorm;
    solution.resize(dim);
    stopNorm.resize(dim);

    V.x.assign(dim,0);

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
    double Pi=3.141592654;
    omega=2/(1+sqrt(1-pow(cos(Pi*h),2)));

    while(eps<norm) {
        for(i=0;i<dim;i++) {
            sum1=0;
            for(j=0;j<i;j++) {
                if(HashTable(i,j)) {
                    sum1+=A.Get(i,j)*V.x[j];
                }
            }
            sum2=0;
            for(j=i;j<dim;j++) {
                if(HashTable(i,j)) {
                    sum2+=A.Get(i,j)*V.x[j];
                }
            }
            V.x[i] = V.x[i]-omega*1/A.Get(i,i)*(sum1+sum2-V.b[i]);
        }

        for(i=0;i<dim;i++)
            stopNorm[i]=V.x[i]-solution[i];
        norm=O.vectorNorm(stopNorm);
        steps++;
    }

    printf("SORSteps: %d\n", steps);
}

void Algorithms::SSORMethod(Matrix& A, Operators& O, Vectors& V) {
    double eps,h,norm,sum1,sum2,omega;
    int steps,i,j,k;

    steps=0;
    h=1.0/(double)(n+1);

    vector<double> solution;
    vector<double> stopNorm;
    solution.resize(dim);
    stopNorm.resize(dim);

    V.x.assign(dim,0);

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
    double Pi=3.141592654;
    omega=2/(1+sin(Pi*h));

    //while(eps<norm) {
        for(i=0;i<dim;i++) {
            sum1=0;
            for(j=0;j<i;j++) {
                if(HashTable(i,j)) {
                    sum1+=A.Get(i,j)*V.x[j];
                }
            }
            sum2=0;
            for(j=i;j<dim;j++) {
                if(HashTable(i,j)) {
                    sum2+=A.Get(i,j)*V.x[j];
                }
            }
            V.x[i] = V.x[i]-omega*1/A.Get(i,i)*(sum1+sum2-V.b[i]);
        }

        // for(i=0;i<dim;i++)
        //     stopNorm[i]=V.x[i]-solution[i];
        // norm=O.vectorNorm(stopNorm);
        // steps++;
    //}

    //printf("SORSteps: %d\n", steps);
}

void Algorithms::modifiedIncompleteCholesky(WriteableMatrix& A, WriteableMatrix& L, WriteableMatrix& Ltranspose) {
    int i,j,k;
    double sum;

    // for(j=0;j<dim;j++) {
    //     sum=0;
    //     for(k=0;k<j;k++) {
    //         if(HashTable(j,k)) {
    //             sum+=pow(L.Get(j,k),2);
    //         }
    //     }
    //     L.Set(j,j,(A.Get(j,j)-sum));
    //     Ltranspose.Set(j,j,(A.Get(j,j)-sum));
    //     for(i=j+1;i<dim;i++) {
    //         if(HashTable(i,j)) {
    //             sum=0;
    //             for(k=0;k<j;k++) {
    //                 if(HashTable(j,k) && HashTable(i,k)) {
    //                     sum+=L.Get(j,k)*L.Get(i,k);
    //                 }
    //             }
    //             L.Set(i,j,((A.Get(i,j)-sum)/L.Get(j,j)));
    //             Ltranspose.Set(i,j,((A.Get(i,j)-sum)/L.Get(j,j)));
    //         }
    //     }
    // }

    for(k=0;k<dim;k++) {
        sum=0;
        for(j=0;j<k;j++) {
            if (HashTable(k,j)) {
                sum+=L.Get(k,j);
            }
        }
        L.Set(k,k,(sqrt(A.Get(k,k)-sum)));
        Ltranspose.Set(k,k,(sqrt(A.Get(k,k)-sum)));
        for(i=k+1;i<dim;i++) {
            if (HashTable(i,k)) {
                sum=0;
                for(j=0;j<k;j++) {
                    if (HashTable(i,j) && HashTable(k,j)) {
                        sum+=L.Get(i,j)*L.Get(k,j);
                    }
                }
                L.Set(i,k,((A.Get(i,k) - sum)/L.Get(k,k)));
                Ltranspose.Set(k,i,((A.Get(i,k) - sum)/L.Get(k,k)));
            } else {
                sum=0;
                for(j=0;j<k;j++) {
                    if (HashTable(i,j) && HashTable(k,j)) {
                        sum+=L.Get(i,j)*L.Get(k,j);
                    }
                }
                A.Set(i,i,(A.Get(i,i) - sum));
            }
        }
    }
}

void Algorithms::incompleteCholesky(PoissonMatrix& A, LowerMatrix& L, UpperMatrix& Ltranspose) {
    int i,j,k;
    double sum;

    for(k=0;k<dim;k++) {
        sum=0;
        for(j=0;j<k;j++) {
            if(HashTable(k,j)) {
                sum+=L.Get(k,j)*L.Get(k,j);
            }
        }
        L.Set(k,k,(sqrt(A.Get(k,k)-sum)));
        Ltranspose.Set(k,k,(sqrt(A.Get(k,k)-sum)));
        for(i=k+1;i<dim;i++) {
            if (HashTable(i,k)) {
                sum=0;
                for(j=0;j<k;j++) {
                    if(HashTable(i,j) && HashTable(k,j)) {
                        sum+=L.Get(i,j)*L.Get(k,j);
                    }
                }
                L.Set(i,k,((A.Get(i,k) - sum)/L.Get(k,k)));
                Ltranspose.Set(k,i,((A.Get(i,k) - sum)/L.Get(k,k)));
            }
        }
    }
}

void Algorithms::CG(Matrix& A, Operators& O, Vectors& V) {
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

void Algorithms::PCG(Matrix& A, Operators& O, WriteableMatrix& L, WriteableMatrix& U, Vectors& V) {
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