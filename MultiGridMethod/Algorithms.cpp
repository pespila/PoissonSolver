#include "classes.h"

Algorithms::Algorithms(int m) {
    n=m;
    dim=n*n;
    h=1.0/(double)(n+1);
}

Algorithms::~Algorithms() {
}

void Algorithms::MatrixVectorMultiplyer(PoissonMatrix& A, const vector<double>& x, vector<double>& b) {
    int dim=A.Size();
    for(int i=0;i<dim;i++) {
        b[i]+=x[i]*A.Get(i,i);
        if(i<(dim-n)) {
            b[i]+=x[i+n]*A.Get(i,i+n);
        }
        if(i>=n) {
            b[i]+=x[i-n]*A.Get(i,i-n);
        }
        if(i%n!=0) {
            b[i]+=x[i-1]*A.Get(i,i-1);
            b[i-1]+=x[i]*A.Get(i-1,i);
        }
    }
}

void Algorithms::JacobiMethod(PoissonMatrix& A, vector<double>& x, const vector<double>& b, int steps) {
    double sum;
    int i,j,k,m;
    int dim=A.Size();
    for(j=0;j<steps;j++) {
        for(i=0;i<dim;i++) {
            sum=0;
            for(k=0;k<5;k++) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m!=i) {
                    sum+=A.Get(i,m)*x[m];
                }
            }
            x[i]=1/A.Get(i,i)*(b[i]-sum);
        }
    }
}

void Algorithms::GaussSeidelMethod(PoissonMatrix& A, vector<double>& x, const vector<double>& b, int steps) {
    double sum1,sum2;
    int i,j,k,m;
    int dim=A.Size();
    for(j=0;j<steps;j++) {
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
            x[i]=1/A.Get(i,i)*(b[i]-sum1-sum2);
        }
    }
}

void Algorithms::SORMethod(PoissonMatrix& A, vector<double>& x, const vector<double>& b, int steps) {
    double sum1,sum2,Pi=3.141592654,omega=2/(1+sqrt(1-pow(cos(Pi*h),2)));
    int i,j,k,m;
    int dim=A.Size();
    for(j=0;j<steps;j++) {
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
            x[i]=x[i]-omega*1/A.Get(i,i)*(sum1+sum2-b[i]);
        }
    }
}

void Algorithms::Restriction(const vector<double>& r,vector<double>& r2h,int n) {
    //vector<double> r2h;
    // int tmp=sqrt(r.size());
    // int n=((tmp+1)/2)-1;
    //int tmp=n*n;
    //r2h.resize(n*n);
    int i,j,k,l;
    k=0;
    l=0;
    // for(i=1;i<=n;i++) {
    //     for(j=1;j<=n;j++) {
    //         if(i%2==0 && j%2==0) {
    //             r2h[l]=r[k];
    //             l++;
    //         }
    //         k++;
    //     }
    // }
    for(i=1;i<=n;i++) {
        for(j=1;j<=n;j++) {
            if(i%2==0 && j%2==0) {
                r2h[l]=1.0/16.0*(4.0*r[k]+2.0*(r[k-1]+r[k+1]+r[k-n]+r[k+n])+r[k+n-1]+r[k+n+1]+r[k-n-1]+r[k-n+1]);
                r[k];
                l++;
            }
            k++;
        }
    }
    //return r2h;
}

void Algorithms::Interpolation(const vector<double>& E2h,vector<double>& E,Vectors& V,int n) {
    // int tmp=sqrt(E2h.size());
    // int n=((tmp+1)*2)-1;
    int k,l,i,j;
    // vector<double> E;
    // E.resize(n*n);
    k=0;
    l=0;
    for(i=1;i<=n;i++) {
        for(j=1;j<=n;j++) {
            if(i%2==0 && j%2==0) {
                E[k]=E2h[l];
                l++;
            }
            k++;
        }
    }
    // k=0;
    // for(int i=1;i<=n;i++) {
    //     for(int j=1;j<=n;j++) {
    //         if(i%2==0 && j%2!=0) {
    //             if(j!=1 && j!=n) {
    //                 E[k]=(double)1/(double)2*(E[k-1]+E[k+1]);
    //             }
    //             if(j==1) {
    //                 E[k]=(double)1/(double)2*(E[k+1]+V.g(0,(double)i/(double)(n+1)));
    //             }
    //             if(j==n) {
    //                 E[k]=(double)1/(double)2*(E[k-1]+V.g(1,(double)i/(double)(n+1)));
    //             }
    //         }
    //         if(i%2!=0 && j%2==0) {
    //             if(i!=1 && i!=n) {
    //                 E[k]=(double)1/(double)2*(E[k-n]+E[k+n]);
    //             }
    //             if(i==1) {
    //                 E[k]=(double)1/(double)2*(E[k+n]+V.g((double)j/(double)(n+1),0));
    //             }
    //             if(i==n) {
    //                 E[k]=(double)1/(double)2*(E[k-n]+V.g((double)j/(double)(n+1),1));
    //             }
    //         }
    //         if(i%2!=0 && j%2!=0) {
    //             if(i!=1 && j!=1 && i!=n && j!=n) {
    //                 E[k]=(double)1/(double)4*(E[k+n-1]+E[k+n+1]+E[k-n+1]+E[k-n-1]);
    //             }
    //             if(i==1 && j==1) {
    //                 E[k]=(double)1/(double)4*(E[k+n+1]+V.g(0,0)+V.g(0,(double)(i+1)/(double)(n+1))+V.g((double)(i+1)/(double)(n+1),0));
    //             }
    //             if(i==n && j==n) {
    //                 E[k]=(double)1/(double)4*(E[k-n-1]+V.g(1,1)+V.g(1,(double)(i-1)/(double)(n+1))+V.g((double)(i-1)/(double)(n+1),1));
    //             }
    //             if(i==1 && j==n) {
    //                 E[k]=(double)1/(double)4*(E[k-n+1]+V.g(0,1)+V.g(0,(double)(i-1)/(double)(n+1))+V.g((double)(i+1)/(double)(n+1),1));
    //             }
    //             if(i==n && j==1) {
    //                 E[k]=(double)1/(double)4*(E[k+n-1]+V.g(1,0)+V.g(1,(double)(i+1)/(double)(n+1))+V.g((double)(i-1)/(double)(n+1),0));
    //             }
    //             if(i==1 && j!=1 && j!=n) {
    //                 E[k]=(double)1/(double)4*(E[k+n+1]+E[k+n-1]+V.g((double)(i-1)/(double)(n+1),0)+V.g((double)(i+1)/(double)(n+1),0));
    //             }
    //             if(i==n && j!=1 && j!=n) {
    //                 E[k]=(double)1/(double)4*(E[k-n-1]+E[k-n+1]+V.g((double)(i-1)/(double)(n+1),1)+V.g((double)(i+1)/(double)(n+1),1));
    //             }
    //             if(j==1 && i!=1 && i!=n) {
    //                 E[k]=(double)1/(double)4*(E[k+n+1]+E[k-n+1]+V.g(0,(double)(j+1)/(double)(n+1))+V.g(0,(double)(j-1)/(double)(n+1)));
    //             }
    //             if(j==n && i!=1 && i!=n) {
    //                 E[k]=(double)1/(double)4*(E[k+n-1]+E[k-n-1]+V.g(1,(double)(j+1)/(double)(n+1))+V.g(1,(double)(j-1)/(double)(n+1)));
    //             }
    //         }
    //         k++;
    //     }
    // }
    k=0;
    for(int i=1;i<=n;i++) {
        for(int j=1;j<=n;j++) {
            if(i%2==0 && j%2!=0) {
                if(j!=1 && j!=n) {
                    E[k]=(double)1/(double)2*(E[k-1]+E[k+1]);
                }
                if(j==1) {
                    E[k]=(double)1/(double)2*(E[k+1]);//+V.g(0,(double)i/(double)(n+1)));
                }
                if(j==n) {
                    E[k]=(double)1/(double)2*(E[k-1]);//+V.g(1,(double)i/(double)(n+1)));
                }
            }
            if(i%2!=0 && j%2==0) {
                if(i!=1 && i!=n) {
                    E[k]=(double)1/(double)2*(E[k-n]+E[k+n]);
                }
                if(i==1) {
                    E[k]=(double)1/(double)2*(E[k+n]);//+V.g((double)j/(double)(n+1),0));
                }
                if(i==n) {
                    E[k]=(double)1/(double)2*(E[k-n]);//+V.g((double)j/(double)(n+1),1));
                }
            }
            if(i%2!=0 && j%2!=0) {
                if(i!=1 && j!=1 && i!=n && j!=n) {
                    E[k]=(double)1/(double)4*(E[k+n-1]+E[k+n+1]+E[k-n+1]+E[k-n-1]);
                }
                if(i==1 && j==1) {
                    E[k]=(double)1/(double)4*(E[k+n+1]);//+V.g(0,0)+V.g(0,(double)(i+1)/(double)(n+1))+V.g((double)(i+1)/(double)(n+1),0));
                }
                if(i==n && j==n) {
                    E[k]=(double)1/(double)4*(E[k-n-1]);//+V.g(1,1)+V.g(1,(double)(i-1)/(double)(n+1))+V.g((double)(i-1)/(double)(n+1),1));
                }
                if(i==1 && j==n) {
                    E[k]=(double)1/(double)4*(E[k-n+1]);//+V.g(0,1)+V.g(0,(double)(i-1)/(double)(n+1))+V.g((double)(i+1)/(double)(n+1),1));
                }
                if(i==n && j==1) {
                    E[k]=(double)1/(double)4*(E[k+n-1]);//+V.g(1,0)+V.g(1,(double)(i+1)/(double)(n+1))+V.g((double)(i-1)/(double)(n+1),0));
                }
                if(i==1 && j!=1 && j!=n) {
                    E[k]=(double)1/(double)4*(E[k+n+1]+E[k+n-1]);//+V.g((double)(i-1)/(double)(n+1),0)+V.g((double)(i+1)/(double)(n+1),0));
                }
                if(i==n && j!=1 && j!=n) {
                    E[k]=(double)1/(double)4*(E[k-n-1]+E[k-n+1]);//+V.g((double)(i-1)/(double)(n+1),1)+V.g((double)(i+1)/(double)(n+1),1));
                }
                if(j==1 && i!=1 && i!=n) {
                    E[k]=(double)1/(double)4*(E[k+n+1]+E[k-n+1]);//+V.g(0,(double)(j+1)/(double)(n+1))+V.g(0,(double)(j-1)/(double)(n+1)));
                }
                if(j==n && i!=1 && i!=n) {
                    E[k]=(double)1/(double)4*(E[k+n-1]+E[k-n-1]);//+V.g(1,(double)(j+1)/(double)(n+1))+V.g(1,(double)(j-1)/(double)(n+1)));
                }
            }
            k++;
        }
    }
    // printf("E: ");
    // for(i=0;i<E.size();i++) {
    //     printf("%.3f ", E[i]);
    // }
    // printf("\n");
    //return E;
}

vector<double> Algorithms::MultiGridMethod(PoissonMatrix& A,Vectors& V,const vector<double>& b,int n) {
    int i,dim=n*n;
    vector<double> x,solution;
    x.assign(dim,0);
    // solution.assign(dim,0);
    // int k=0;
    // for (i=1;i<=n;i++) {
    //     for (int j=1;j<=n;j++) {
    //         solution[k] = V.g(j*h,i*h);
    //         k++;
    //     }
    // }
    if(n==1) {
        double a;
        //h=(double)1/(double)(n+1);
        a=4.0*(double)pow(n+1,2),
        //d=V.f(h,h)+pow(n+1,2)*(V.g(0,h)+V.g(h,0)+V.g(1-h,1)+V.g(1,1-h));
        x[0]=b[0]/a;
        return x;
    } else {
        int smallerN=(n+1)/2-1;
        GaussSeidelMethod(A,x,b,4);
        vector<double> Ax,r;
        Ax.assign(dim,0);
        r.assign(dim,0);
        MatrixVectorMultiplyer(A,x,Ax);
        for(i=0;i<dim;i++) {
            r[i]=b[i]-Ax[i];
        }
        vector<double> r2h;
        r2h.assign(smallerN*smallerN,0);
        Restriction(r,r2h,n);
        A.Resize(smallerN);
        A.InitHashMatrix();
        vector<double> E2h;
        E2h.assign(smallerN*smallerN,0);
        E2h=MultiGridMethod(A,V,r2h,smallerN);
        vector<double> E;
        E.assign(n*n,0);
        Interpolation(E2h,E,V,n);
        for(i=0;i<dim;i++) {
            x[i]+=E[i];
        }
        A.Resize(n);
        A.InitHashMatrix();
        GaussSeidelMethod(A,x,b,4);
        return x;
    }
}

// vector<double> Algorithms::MultiGridMethod(PoissonMatrix& A,Vectors& V,const vector<double>& b) {
//     int i,dim=b.size(),n=sqrt(dim);
//     vector<double> x;
//     x.assign(dim,0);
//     if(n==1) {
//         double a,b,h;
//         h=(double)1/(double)(n+1);
//         a=4.0*(double)pow(n+1,2),
//         b=V.f(h,h)+pow(n+1,2)*(V.g(0,h)+V.g(h,0)+V.g(1-h,1)+V.g(1,1-h));
//         x[0]=b/a;
//         return x;
//     } else {
//         int tmp=((n+1)/2)-1;
//         int biggerN=tmp*tmp;
//         vector<double> Ax,r,r2h,E2h,E;
//         Ax.assign(dim,0);
//         r.assign(dim,0);
//         r2h.resize(smallerN);
//         E2h.resize(newN);
//         E.assign(dim,0);
//         GaussSeidelMethod(A,x,b,3);
//         MatrixVectorMultiplyer(A,V.x,Ax);
//         for(i=0;i<dim;i++) {
//             r[i]=b[i]-Ax[i];
//         }
//         vector<double>().swap(Ax);
//         r2h=Restriction(r);
//         vector<double>().swap(r);
//         A.Resize(tmp);
//         E2h=MultiGridMethod(A,V,r2h);
//         vector<double>().swap(r2h);
//         // vector<double> E;
//         // E.resize(dim);
//         E=Interpolation(E2h,V);
//         vector<double>().swap(E2h);
//         for(i=0;i<dim;i++) {
//             x[i]=x[i]+E[i];
//         }
//         vector<double>().swap(E);
//         A.Resize(n);
//         GaussSeidelMethod(A,x,b,3);
//         return x;
//     }
// }

// vector<double> Algorithms::MultiGridMethod(PoissonMatrix& A,Vectors& V,vector<double>& x,vector<double>& b) {
//     int i,dim=A.Size(),n=sqrt(dim);
//     if(n==1) {
//         x[0]=1;
//         return x;
//     } else {
//         vector<double> Ax,r,r2h,E2h,E;
//         Ax.assign(dim,0);
//         r.resize(dim);
//         r2h.resize(pow(((n+1)/2)-1,2));
//         E2h.resize(pow(((n+1)/2)-1,2));
//         GaussSeidelMethod(A,x,b,3);
//         MatrixVectorMultiplyer(A,V.x,Ax);
//         for(i=0;i<dim;i++) {
//             r[i]=b[i]-Ax[i];
//         }
//         r2h=Restriction(r);
//         A.Resize(((n+1)/2)-1);
//         E2h=MultiGridMethod(A,V,E2h,r2h);
//         E.resize(dim);
//         E=Interpolation(E2h,V);
//         for(i=0;i<dim;i++) {
//             x[i]=x[i]+E[i];
//         }
//         GaussSeidelMethod(A,x,b,3);
//         return x;
//     }
// }

// vector<double> Algorithms::Restriction(vector<double>& r, int n) {
//     vector<double> r2h;
//     r2h.resize(pow((((n+1)/2)-1),2));
//     int k=0;
//     int l=0;
//     for(int i=1;i<=n;i++){
//         for(int j=1;j<=n;j++){
//             if(i%2==0 && j%2==0) {
//                 r2h[l]=r[k];
//                 l++;
//             }
//             k++;
//         }
//     }
//     return r2h;
// }

// vector<double> Algorithms::Interpolation(vector<double>& E2h, int n) {
//     int i,j,k,l;
//     vector<double> E;
//     //int dim=pow(((n+1)/2)-1,2);
//     E.resize(n*n);
//     // for(i=0,k=0;i<n,k<E2h.size();i++,k++) {
//     //     for(j=0,l=0;j<n,l<E2h.size();j++,l++) {

//     //     }
//     // }
//     k=0;
//     l=0;
//     for(i=1;i<=n;i++) {
//         for(j=1;j<=n;j++) {
//             if(i%2==0 && j%2==0) {
//                 E[k]=E2h[l];
//                 l++;
//             }
//             k++;
//         }
//     }
//     k=0;
//     for(i=1;i<=n;i++) {
//         for(j=0;j<=n;j++) {
//             if(i%2!=0 && j%2==0) {
//                 if(i==0) {
//                     E[k]=(double)1/(double)2*E[k+n];
//                 } else if(i==n) {
//                     E[k]=(double)1/(double)2*E[k-n];
//                 } else {
//                     E[k]=(double)1/(double)2*(E[k-1]+E[k+1]);
//                 }
//             }
//             if(i>=1 && j>=1) k++;
//         }
//     }
//     k=0;
//     for(i=0;i<=n;i++) {
//         for(j=1;j<=n;j++) {
//             if(i%2==0 && j%2!=0) {
//                 if(j==0) {
//                     E[k]=(double)1/(double)2*E[k+1];
//                 } else if(j==n) {
//                     E[k]=(double)1/(double)2*E[k-1];
//                 } else {
//                     E[k]=(double)1/(double)2*(E[k-n]+E[k+n]);
//                 }
//             }
//             if(i>=1 && j>=1) k++;
//         }
//     }
//     k=0;
//     for(i=1;i<=n;i++) {
//         for(j=1;j<=n;j++) {
//             if(i%2!=0 && j%2!=0) {
//                 E[k]=(double)1/(double)4*(E[k-n-1]+E[k-n+1]+E[k+n-1]+E[k+n+1]);
//             } else if(i==1 && j==1) {
//                 E[k]=(double)1/(double)4*E[k+n+1];
//             } else if(i==n && j==1) {
//                 E[k]=(double)1/(double)4*E[k+n-1];
//             } else if(i==1 && j==n) {
//                 E[k]=(double)1/(double)4*E[k-n+1];
//             } else if(i==n && j==n) {
//                 E[k]=(double)1/(double)4*E[k-n-1];
//             } else if(i==1 && j%2!=0 && j!=1 && j!=n) {
//                 E[k]=(double)1/(double)4*(E[k+n-1]+E[k+n-1]);
//             } else if(i==n && j%2!=0 && j!=1 && j!=n) {
//                 E[k]=(double)1/(double)4*(E[k-n-1]+E[k-n+1]);
//             } else if(j==1 && i%2!=0 && i!=1 && i!=n) {
//                 E[k]=(double)1/(double)4*(E[k+n+1]+E[k-n+1]);
//             } else if(j==n && i%2!=0 && i!=1 && i!=n) {
//                 E[k]=(double)1/(double)4*(E[k+n-1]+E[k-n-1]);
//             }
//             if(i>=0 && j>=0) {
//                 k++;
//             }
//         }
//     }
//     return E;
// }