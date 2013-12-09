#include "classes.h"

Algorithms::Algorithms(Matrix& M) {
    dim = M.Size();
    n = sqrt(dim);
}

Algorithms::~Algorithms() {
}

void Algorithms::modifiedIncompleteLU(Matrix& A, WriteableMatrix& L, WriteableMatrix& U) {
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
}

// void Algorithms::incompleteLU(Matrix& A, WriteableMatrix& L, WriteableMatrix& U) {
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

void Algorithms::JacobiMethod(Matrix& A, Operators& O, vector<double>& x, const vector<double>& b, int maxIterations = 5000) {
    double eps,h,norm,sum;
    int steps,i,j,k,m;

    steps=0;
    h=1.0/(double)(n+1);

    vector<double> solution;
    vector<double> stopNorm;
    solution.resize(dim);
    stopNorm.resize(dim);

    x.assign(dim,0);

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

    while(eps<=norm) {
        for(i=0;i<dim;i++) {
            sum=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m!=i) {
                    sum+=A.Get(i,m)*x[m];
                }
            }
            x[i]=1/A.Get(i,i)*(b[i]-sum);
            stopNorm[i]=x[i]-solution[i];
        }
        norm=O.vectorNorm(stopNorm);
        steps++;
        if(steps==maxIterations) {
            norm=0.0;
        }
    }
    printf("JacobianSteps: %d\n", steps);
}

void Algorithms::GaussSeidelMethod(Matrix& A, Operators& O, vector<double>& x, const vector<double>& b, int maxIterations = 5000) {
    double eps,h,norm,sum1,sum2;
    int steps,i,j,k,m;

    steps=0;
    h=1.0/(double)(n+1);

    vector<double> solution;
    vector<double> stopNorm;
    solution.resize(dim);
    stopNorm.resize(dim);

    x.assign(dim,0);

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
    
    while(eps<=norm) {
        for(i=0;i<dim;i++) {
            sum1=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m<i) {
                    sum1+=A.Get(i,m)*x[m];
                }
            }
            sum2=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m>i) {
                    sum2+=A.Get(i,m)*x[m];
                }
            }
            x[i] = 1/A.Get(i,i)*(b[i]-sum1-sum2);
            stopNorm[i]=x[i]-solution[i];
        }
        norm=O.vectorNorm(stopNorm);
        steps++;
        if(steps==maxIterations) {
            norm=0.0;
        }
    }

    printf("GaussSeidelSteps: %d\n", steps);
}

void Algorithms::SORMethod(Matrix& A, Operators& O, vector<double>& x, const vector<double>& b, int maxIterations = 1000) {
    double eps,h,norm,sum1,sum2,omega;
    int steps,i,j,k,m;

    steps=0;
    h=1.0/(double)(n+1);

    vector<double> solution;
    vector<double> stopNorm;
    solution.resize(dim);
    stopNorm.resize(dim);

    x.assign(dim,0);

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
    double Pi=3.141592654;
    omega=2/(1+sqrt(1-pow(cos(Pi*h),2)));

    while(eps<=norm) {
        for(i=0;i<dim;i++) {
            sum1=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m<i) {
                    sum1+=A.Get(i,m)*x[m];
                }
            }
            sum2=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m>=i) {
                    sum2+=A.Get(i,m)*x[m];
                }
            }
            x[i] = x[i]-omega*1/A.Get(i,i)*(sum1+sum2-b[i]);
            stopNorm[i]=x[i]-solution[i];
        }
        norm=O.vectorNorm(stopNorm);
        steps++;
        if(steps==maxIterations) {
            norm=0.0;
        }
    }
    printf("SORSteps: %d\n", steps);
}

// void Algorithms::SSORMethod(Matrix& A, Operators& O, Vectors& V) {
//     double eps,h,norm,sum1,sum2,omega;
//     int steps,i,j,k;

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
//                 if(HashTable(i,j)) {
//                     sum1+=A.Get(i,j)*V.x[j];
//                 }
//             }
//             sum2=0;
//             for(j=i;j<dim;j++) {
//                 if(HashTable(i,j)) {
//                     sum2+=A.Get(i,j)*V.x[j];
//                 }
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

void Algorithms::modifiedIncompleteCholesky(WriteableMatrix& A, WriteableMatrix& L, WriteableMatrix& Ltranspose) {
    int i,j,k,m,u,v;
    double sum;

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
    }
    printf("PCGSteps: %d\n", steps);
}

void Algorithms::Restriction(const vector<double>& r, vector<double>& r2h, int n) {
    int k=0;
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            if(i%2==0 && j%2==0) {
                r2h[k]=r[k];
            }
            k++;
        }
    }
}

void Algorithms::Interpolation(const vector<double>& r2h, vector<double>& E, int n) {
    //int k=0;
    for(int i=1;i<=n;i++) {
        for(int j=1;j<=n;j++) {

        }
    }
}

void Algorithms::MultiGridMethod(vector<double>& x, const vector<double>& b, Operators& O, int m) {
    int i;
    if(m==1) {
        double a,b,h;
        
        h=1.0/(double)(m+1);
        a=4*pow(n+1,2);
        b=O.f(h,h)+pow(n+1,2)*(O.g(0,h)+O.g(h,0)+O.g(1-h,1)+O.g(1,1-h));
        x.resize(m,b/a);
    } else {
        PoissonMatrix A(m);
        GaussSeidelMethod(A,O,x,b,3);
        vector<double> Ax,r,E,r2h,xTmp;
        Ax.resize(m*m);
        r.assign(m*m,0);
        E.resize(m*m);
        r2h.resize((m-1)*(m-1));
        xTmp.assign((m-1)*(m-1),0);
        O.MatrixVectorMultiplyer(A,x,Ax);
        for(i=0;i<m*m;i++) {
            r[i]=b[i]-Ax[i];
        }
        Restriction(r,r2h,m);
        MultiGridMethod(xTmp,r2h,O,m-1);
        Interpolation(r2h,E,m);
        //Compute Interpolation!!!
        for(i=0;i<m*m;i++) {
            x[i]=x[i]+E[i];
        }
        GaussSeidelMethod(A,O,x,b,3);
        vector<double>().swap(Ax);
        vector<double>().swap(r);
        vector<double>().swap(E);
        vector<double>().swap(r2h);
        vector<double>().swap(xTmp);
    }
}

// void Algorithms::MultiGridMethod() {
//     GaussSeidelMethod(A,O,V,3);
//     vector<double> Au;
//     Au.resize(dim);
//     O.MatrixVectorMultiplyer(A,V.xh,Au);
//     int k=0;
//     for(int i=1;i<=n;i++) {
//         for(int j=1;j<=n;j++) {
//             if(i%2==0 && j%2==0) {
//                 V.r2h[k]=V.rh[k];
//             }
//             k++;
//         }
//     }
//     // for(int i=0;i<dim;i++) {
//     //     V.rh[i] = V.bh[i] - Au[i];
//     //     if(i%2!=0) {
//     //         V.r2h[i] = V.rh[i];
//     //     }
//     // }
//     vector<double> E2h;
//     E2h.assign(V.r2h.Size(),0);
//     int N = n+1;
//     int n2h = N/2;
//     PoissonMatrix A2h(n2h);
//     GaussSeidelMethod(A2h,O,V,3);
//     vector<double> Eh;
//     Eh.resize(dim);
//     int k=0;
//     //int l=0;
//     for(int i=0;i<=n;i++) {
//         for(int j=0;j<=n;j++) {
//             if(i%2==0 && j%2==0 && i>0 && j>0) {
//                 V.rh[k]=V.r2h[k];
//                 //l++;
//             }
//             if(i%2!=0 && j%2==0 && j>0) {
//                 V.rh[k]=1/2*(V.r2h[]+V.r2h[]);
//             }
//             if(i%2==0 && j%2!=0 && i>0) {
//                 V.rh[k]=1/2*
//             }
//             if(i%2!=0 && j%2!=0) {
//                 V.rh[k]=1/4*
//             }
//             k++;
//         }
//     }
//     GaussSeidelMethod(A,O,V,3);
// }