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
    int i,j,k,m,u,v,tmp=0;
    double sum, drop;

    for(i=0;i<dim;i++) {
        drop=0;
        for(k=0;k<5;k++) {
            m=A.HashMatrix[i][k];
            if(m!=-1 && m<i) {
                sum=0;
                for(j=0;j<5;j++) {
                    u=A.HashMatrix[i][j];
                    if(u!=-1 && u<k) {// && v!=-1 && v<k) {
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
                    if(u!=-1 && u<i) {// && v!=-1 && v<i) {
                        sum+=L.Get(i,u)*U.Get(u,m);
                    }
                }
                U.Set(i,m,(A.Get(i,m)-sum));
                drop+=sum;
            }
        }
        //printf("%f\n", drop);
        U.Set(i,i,(U.Get(i,i)-drop));
    }

    // for(i=0;i<dim;i++) {
    //     drop = 0;
    //     for(k=0;k<i;k++) {
    //         sum = 0;
    //         for(j=0;j<k;j++) {
    //             if(HashTable(i,j) && HashTable(j,k)) {
    //                 sum += L.Get(i,j) * U.Get(j,k);
    //             }
    //         }
    //         if(HashTable(i,k)) {
    //             L.Set(i,k,((A.Get(i,k) - sum)/U.Get(k,k)));
    //         } else {
    //             drop += sum;
    //         }
    //     }
    //     for(k=i;k<dim;k++) {
    //         sum = 0;
    //         for(j=0;j<i;j++) {
    //             if(HashTable(i,j) && HashTable(j,k)) {
    //                 sum += L.Get(i,j) * U.Get(j,k);
    //             }
    //         }
    //         if(HashTable(i,k)) {
    //             U.Set(i,k,(A.Get(i,k) - sum));
    //         } else {
    //             drop += sum;
    //         }
    //     }
    //     printf("%f\n", drop);
    //     U.Set(i,i,(U.Get(i,i) - drop));
    // }
}

// void Algorithms::LU(Matrix& A, WriteableMatrix& L, WriteableMatrix& U) {
//     int i,j,k,m,u;
//     double sum;
//     for(i=0;i<dim;i++) {
//         for(k=0;k<5;k++) {
//             m = A.HashMatrix[i][k];
//             if(m!=-1 && m<=i-1) {
//                 sum=0;
//                 for(j=0;j<m;j++) {
//                     u = A.HashMatrix[i][j];
//                     if(u!=-1) {
//                         sum+=L.Get(i,j)*U.Get(j,m);
//                     }
//                 }
//                 L.Set(i,m,(A.Get(i,m)-sum)/U.Get(m,m));
//             }
//         }
//         for(k=0;k<5;k++) {
//             m = A.HashMatrix[i][k];
//             if(m!=-1 && m>=i-1) {
//                 sum=0;
//                 for(j=0;j<i;j++) {
//                     u = A.HashMatrix[i][j];
//                     if(u!=-1) {
//                         sum+=L.Get(i,j)*U.Get(j,m);
//                     }
//                 }
//                 U.Set(i,m,(A.Get(i,m)-sum));
//             }
//         }
//     }
// }

void Algorithms::incompleteLU(Matrix& A, WriteableMatrix& L, WriteableMatrix& U) {
    int i,j,k;
    double sum;

    for(i=0;i<dim;i++) {
        for(k=1;k<i;k++) {
            if (HashTable(i,k)) {
                sum=0;
                for(j=0;j<k;j++) {
                    if(HashTable(i,j) && HashTable(j,k))
                        sum+=L.Get(i,j)*U.Get(j,k);
                }
                L.Set(i,k,((A.Get(i,k)-sum)/U.Get(k,k)));
            }
        }
        for(k=i;k<dim;k++) {
            if (HashTable(i,k) != 0) {
                sum=0;
                for(j=0;j<i;j++) {
                    if(HashTable(i,j) && HashTable(j,k))
                        sum+=L.Get(i,j)*U.Get(j,k);
                }
                U.Set(i,k,(A.Get(i,k)-sum));
            }
        }
    }
}

void Algorithms::JacobiMethod(Matrix& A, Operators& O, Vectors& V) {
    double eps,h,norm,sum;
    int steps,i,j,k,m;

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

    while(eps<=norm) {
        for(i=0;i<dim;i++) {
            sum=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m!=i) {
                    sum+=A.Get(i,m)*V.x[m];
                }
            }
            V.x[i]=1/A.Get(i,i)*(V.b[i]-sum);
            stopNorm[i]=V.x[i]-solution[i];
        }
        norm=O.vectorNorm(stopNorm);
        steps++;
    }
    printf("JacobianSteps: %d\n", steps);
}

void Algorithms::GaussSeidelMethod(Matrix& A, Operators& O, Vectors& V) {
    double eps,h,norm,sum1,sum2;
    int steps,i,j,k,m;

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
    
    while(eps<=norm) {
        for(i=0;i<dim;i++) {
            sum1=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m<i) {
                    sum1+=A.Get(i,m)*V.x[m];
                }
            }
            sum2=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m>i) {
                    sum2+=A.Get(i,m)*V.x[m];
                }
            }
            V.x[i] = 1/A.Get(i,i)*(V.b[i]-sum1-sum2);
            stopNorm[i]=V.x[i]-solution[i];
        }
        norm=O.vectorNorm(stopNorm);
        steps++;
    }

    printf("GaussSeidelSteps: %d\n", steps);
}

void Algorithms::SORMethod(Matrix& A, Operators& O, Vectors& V) {
    double eps,h,norm,sum1,sum2,omega;
    int steps,i,j,k,m;

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

    while(eps<=norm) {
        for(i=0;i<dim;i++) {
            sum1=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m<i) {
                    sum1+=A.Get(i,m)*V.x[m];
                }
            }
            sum2=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m>=i) {
                    sum2+=A.Get(i,m)*V.x[m];
                }
            }
            V.x[i] = V.x[i]-omega*1/A.Get(i,i)*(sum1+sum2-V.b[i]);
            stopNorm[i]=V.x[i]-solution[i];
        }
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
    steps++;
    eps++;
    norm++;

        // for(i=0;i<dim;i++)
        //     stopNorm[i]=V.x[i]-solution[i];
        // norm=O.vectorNorm(stopNorm);
        // steps++;
    //}

    //printf("SORSteps: %d\n", steps);
}

void Algorithms::modifiedIncompleteCholesky(WriteableMatrix& A, WriteableMatrix& L, WriteableMatrix& Ltranspose) {
    int i,j,k,m,u,v;
    double sum;

    // for(k=0;k<dim;k++) {
    //     sum=0;
    //     for(j=0;j<k;j++) {
    //         sum+=L.Get(k,j);
    //     }
    //     L.Set(k,k,sqrt(A.Get(k,k)-sum));
    //     Ltranspose.Set(k,k,sqrt(A.Get(k,k)-sum));
    //     for(i=k+1;i<dim;i++) {
    //         if(A.Get(i,k)!=0) {
    //             sum=0;
    //             for(j=0;j<k;j++) {
    //                 sum+=L.Get(i,j)*L.Get(k,j);
    //             }
    //             L.Set(i,k,(A.Get(i,k)-sum)/L.Get(k,k));
    //             Ltranspose.Set(k,i,(A.Get(i,k)-sum)/L.Get(k,k));
    //         } else {
    //             sum=0;
    //             for(j=0;j<k;j++) {
    //                 sum+=L.Get(i,j)*L.Get(k,j);
    //             }
    //             A.Set(i,i,A.Get(i,i)-sum);
    //         }
    //     }
    // }

    for(k=0;k<dim;k++) {
        sum=0;
        for(j=0;j<5;j++) {
            m=A.HashMatrix[k][j];
            if(m!=-1 && m<k) {
                sum+=L.Get(k,m);
            }
        }
        L.Set(k,k,sqrt(A.Get(k,k)-sum));
        Ltranspose.Set(k,k,sqrt(A.Get(k,k)-sum));
        for(j=0;j<5;j++) {
            m=A.HashMatrix[k][j];
            if(m!=-1 && m>=k+1) {
                sum=0;
                for(i=0;i<5;i++) {
                    u=A.HashMatrix[m][i];
                    v=A.HashMatrix[k][i];
                    if(u!=-1 && u<k && v!=-1 && v<k) {
                        sum+=L.Get(m,u)*L.Get(m,v);
                    }
                }
                L.Set(m,k,(A.Get(m,k)-sum)/L.Get(k,k));
                Ltranspose.Set(k,m,(A.Get(m,k)-sum)/L.Get(k,k));
            } else if(m!=-1 && m<k+1) {
                sum=0;
                for(i=0;i<5;i++) {
                    u=A.HashMatrix[m][i];
                    v=A.HashMatrix[k][i];
                    if(u!=-1 && u<k && v!=-1 && v<k) {
                        sum+=L.Get(m,u)*L.Get(m,v);
                    }
                }
                A.Set(m,m,A.Get(m,m)-sum);
            }
        }
    }
}

void Algorithms::incompleteCholesky(PoissonMatrix& A, LowerMatrix& L, UpperMatrix& Ltranspose) {
    int i,j,k,m,u,v;
    double sum;

    for(k=0;k<dim;k++) {
        sum=0;
        for(j=0;j<5;j++) {
            m=A.HashMatrix[k][j];
            if(m!=-1 && m<k) {
                sum+=pow(L.Get(k,m),2);
            }
        }
        L.Set(k,k,sqrt(A.Get(k,k)-sum));
        Ltranspose.Set(k,k,sqrt(A.Get(k,k)-sum));
        for(j=0;j<5;j++) {
            m=A.HashMatrix[k][j];
            if(m!=-1 && m>=k+1) {
                sum=0;
                for(i=0;i<5;i++) {
                    u=A.HashMatrix[m][i];
                    v=A.HashMatrix[k][i];
                    if(u!=-1 && u<k && v!=-1 && v<k) {
                        sum+=L.Get(m,u)*L.Get(m,v);
                    }
                }
                L.Set(m,k,(A.Get(m,k)-sum)/L.Get(k,k));
                Ltranspose.Set(k,m,(A.Get(m,k)-sum)/L.Get(k,k));
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

    O.LUsolverLower(A,L,z);
    O.LUsolverUpper(A,U,z);

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
        if (steps > 100)
            break;
    }
    printf("PCGSteps: %d\n", steps);
}