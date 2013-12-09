#include "classes.h"

CGAlgorithm::CGAlgorithm() {
}

CGAlgorithm::~CGAlgorithm() {
}

void CGAlgorithm::InitHashMatrix(int n) {
    //vector<vector<int> >().swap(HashMatrix);
    int i,dim;
    dim=n*n;
    vector<int> push;
    push.assign(5,-1);
    for(i=0;i<dim;i++) {
        if(i<n) {
            if(i==0) {
                push[0]=i;
                push[1]=i+1;
                push[2]=i+n;
            } else if(i%n!=n-1 && i!=0) {
                push[0]=i-1;
                push[1]=i;
                push[2]=i+1;
                push[3]=i+n;
            } else {
                push[0]=i-1;
                push[1]=i;
                push[2]=i+n;
            }
        } else if(i>=n && i<=dim-n) {
            if(i%n==0) {
                push[0]=i-n;
                push[1]=i;
                push[2]=i+1;
                push[3]=i+n;
                push[4]=-1;
            } else if(i%n==n-1) {
                push[0]=i-n;
                push[1]=i-1;
                push[2]=i;
                push[3]=i+n;
                push[4]=-1;
            } else {
                push[0]=i-n;
                push[1]=i-1;
                push[2]=i;
                push[3]=i+1;
                push[4]=i+n;
            }
        } else if(i>dim-n) {
            if(i==dim-1) {
                push[0]=i-n;
                push[1]=i-1;
                push[2]=i;
                push[3]=-1;
                push[4]=-1;
            } else if(i%n!=0 && i!=dim-1) {
                push[0]=i-n;
                push[1]=i-1;
                push[2]=i;
                push[3]=i+1;
                push[4]=-1;
            } else {
                push[0]=i-n;
                push[1]=i;
                push[2]=i+1;
                push[3]=-1;
                push[4]=-1;
            }
        }
        HashMatrix.push_back(push);
    }
}

void CGAlgorithm::modifiedIncompleteLU(Matrix& A, WriteableMatrix& L, WriteableMatrix& U, PoissonOperators& O) {
    int i,j,k,m,u;
    double sum, drop;
    int dim=A.Size()*A.Size();

    for(i=0;i<dim;i++) {
        drop=0;
        for(k=0;k<5;k++) {
            m=HashMatrix[i][k];
            if(m!=-1 && m<i) {
                sum=0;
                for(j=0;j<5;j++) {
                    u=HashMatrix[i][j];
                    if(u!=-1 && u<k) {
                        sum+=L.Get(i,u)*U.Get(u,m);
                    }
                }
                L.Set(i,m,(A.Get(i,m)-sum)/U.Get(m,m));
                drop+=sum;
            }
        }
        for(k=0;k<5;k++) {
            m=HashMatrix[i][k];
            if(m!=-1 && m>=i) {
                sum=0;
                for(j=0;j<5;j++) {
                    u=HashMatrix[i][j];
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

void CGAlgorithm::modifiedIncompleteCholesky(WriteableMatrix& A, WriteableMatrix& L, WriteableMatrix& Ltranspose, PoissonOperators& O) {
    int i,j,k,m,u,v;
    double sum;
    int dim=A.Size()*A.Size();

    for(k=0;k<dim;k++) {
        sum=0;
        for(j=0;j<5;j++) {
            m=HashMatrix[k][j];
            if(m!=-1 && m<k) {
                sum+=L.Get(k,m);
            }
        }
        L.Set(k,k,sqrt(A.Get(k,k)-sum));
        Ltranspose.Set(k,k,sqrt(A.Get(k,k)-sum));
        for(j=0;j<5;j++) {
            m=HashMatrix[k][j];
            if(m!=-1 && m>=k+1) {
                sum=0;
                for(i=0;i<5;i++) {
                    u=HashMatrix[m][i];
                    v=HashMatrix[k][i];
                    if(u!=-1 && u<k && v!=-1 && v<k) {
                        sum+=L.Get(m,u)*L.Get(m,v);
                    }
                }
                L.Set(m,k,(A.Get(m,k)-sum)/L.Get(k,k));
                Ltranspose.Set(k,m,(A.Get(m,k)-sum)/L.Get(k,k));
            } else if(m!=-1 && m<k+1) {
                sum=0;
                for(i=0;i<5;i++) {
                    u=HashMatrix[m][i];
                    v=HashMatrix[k][i];
                    if(u!=-1 && u<k && v!=-1 && v<k) {
                        sum+=L.Get(m,u)*L.Get(m,v);
                    }
                }
                A.Set(m,m,A.Get(m,m)-sum);
            }
        }
    }
}

void CGAlgorithm::incompleteCholesky(PoissonMatrix& A, LowerMatrix& L, UpperMatrix& Ltranspose, PoissonOperators& O) {
    int i,j,k,m,u,v;
    double sum;
    int dim=A.Size()*A.Size();

    for(k=0;k<dim;k++) {
        sum=0;
        for(j=0;j<5;j++) {
            m=HashMatrix[k][j];
            if(m!=-1 && m<k) {
                sum+=pow(L.Get(k,m),2);
            }
        }
        L.Set(k,k,sqrt(A.Get(k,k)-sum));
        Ltranspose.Set(k,k,sqrt(A.Get(k,k)-sum));
        for(j=0;j<5;j++) {
            m=HashMatrix[k][j];
            if(m!=-1 && m>=k+1) {
                sum=0;
                for(i=0;i<5;i++) {
                    u=HashMatrix[m][i];
                    v=HashMatrix[k][i];
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

// void CGAlgorithm::incompleteLU(Matrix& A, WriteableMatrix& L, WriteableMatrix& U, PoissonOperators& O) {
//     int i,j,k;
//     double sum;

//     for(i=0;i<dim;i++) {
//         for(k=1;k<i;k++) {
//             if (HashTable(i,k)) {
//                 sum=0;
//                 for(j=0;j<k;j++) {
//                     if(HashTable(i,j) && HashTable(j,k))
//                         sum+=L.Get(i,j)*U.Get(j,k);
//                 }
//                 L.Set(i,k,((A.Get(i,k)-sum)/U.Get(k,k)));
//             }
//         }
//         for(k=i;k<dim;k++) {
//             if (HashTable(i,k) != 0) {
//                 sum=0;
//                 for(j=0;j<i;j++) {
//                     if(HashTable(i,j) && HashTable(j,k))
//                         sum+=L.Get(i,j)*U.Get(j,k);
//                 }
//                 U.Set(i,k,(A.Get(i,k)-sum));
//             }
//         }
//     }
// }

void CGAlgorithm::LUsolverLower(Matrix& A, Matrix& L, vector<double>& z, PoissonOperators& O) {
    int m;
    int dim=z.size();
    for(int i=0;i<dim;i++) {
        for(int j=0;j<5;j++) {
            m=HashMatrix[i][j];
            if(m!=-1 && m < i) {
                z[i]-=L.Get(i,m)*z[m];
            }
        }
        z[i]/=L.Get(i,i);
    }
}

void CGAlgorithm::LUsolverUpper(Matrix& A, Matrix& U, vector<double>& z, PoissonOperators& O) {
    int m;
    int dim=z.size();
    for(int i=dim-1;i>=0;i--) {
        for(int j=0;j<5;j++) {
            m=HashMatrix[i][j];
            if(m!=-1 && m >= i) {
                z[i]-=U.Get(i,m)*z[m];
            }
         }
        z[i]/=U.Get(i,i);
    }
}

// void CGAlgorithm::SSORMethod(Matrix& A, PoissonOperators& O, PoissonVector& V) {
//     double eps,h,norm,sum1,sum2,omega;
//     int steps,i,j,k;
//     int dim=A.Size();
//     int n=sqrt(dim);

//     steps=0;
//     h=1.0/(double)(n+1);

//     vector<double> solution;
//     vector<double> stopNorm;
//     solution.resize(dim);
//     stopNorm.resize(dim);

//     V.x.assign(dim,0);

//     k=0;
//     for (i=1;i<(n+1);i++) {
//         for (j=1;j<(n+1);j++) {
//             solution[k] = O.g(j*h,i*h);
//             stopNorm[k] = V.x[k]-solution[k];
//             k++;
//         }
//     }

//     eps=pow(10,-3)*O.vectorNorm(stopNorm);
//     norm=1.0;
//     double Pi=3.141592654;
//     omega=2/(1+sin(Pi*h));

//     //while(eps<norm) {
//         for(i=0;i<dim;i++) {
//             sum1=0;
//             for(j=0;j<i;j++) {
//                 //if(HashTable(i,j)) {
//                     sum1+=A.Get(i,j)*V.x[j];
//                 //}
//             }
//             sum2=0;
//             for(j=i;j<dim;j++) {
//                 //if(HashTable(i,j)) {
//                     sum2+=A.Get(i,j)*V.x[j];
//                 //}
//             }
//             V.x[i] = V.x[i]-omega*1/A.Get(i,i)*(sum1+sum2-V.b[i]);
//         }
//     steps++;
//     eps++;
//     norm++;

//         // for(i=0;i<dim;i++)
//         //     stopNorm[i]=V.x[i]-solution[i];
//         // norm=O.vectorNorm(stopNorm);
//         // steps++;
//     //}

//     //printf("SORSteps: %d\n", steps);
// }

void CGAlgorithm::CG(Matrix& A, PoissonOperators& O, PoissonVector& V) {
    double alpha,beta,eps,h,norm,skpOfRes,denom;
    int steps,i,j,k;
    int dim=A.Size();
    int n=sqrt(dim);

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

void CGAlgorithm::PCG(Matrix& A, PoissonOperators& O, WriteableMatrix& L, WriteableMatrix& U, PoissonVector& V) {
    double alpha,beta,eps,h,norm,skpOfRes,denom;
    int steps,i,j,k;
    int dim=A.Size();
    int n=sqrt(dim);

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

    LUsolverLower(A,L,z,O);
    LUsolverUpper(A,U,z,O);

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
        LUsolverLower(A,L,z,O);
        LUsolverUpper(A,U,z,O);
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