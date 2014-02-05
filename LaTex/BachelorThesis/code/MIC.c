void Algorithms::modifiedIncompleteLU(Matrix& A,WriteableMatrix& L,WriteableMatrix& U) {
    int m,u;
    double sum, drop;

    for(int i=0;i<dim;i++) {
        drop=0;
        for(int k=0;k<5;k++) {
            m=A.HashMatrix[i][k];
            if(m!=-1 && m<i) {
                sum=0;
                for(int j=0;j<5;j++) {
                    u=A.HashMatrix[i][j];
                    if(u!=-1 && u<k) {
                        sum+=L.Get(i,u)*U.Get(u,m);
                    }
                }
                L.Set(i,m,(A.Get(i,m)-sum)/U.Get(m,m));
                drop+=sum;
            } else if(m!=-1 && m>=i) {
                m=A.HashMatrix[i][k];
                if(m!=-1 && m>=i) {
                    sum=0;
                    for(int j=0;j<5;j++) {
                        u=A.HashMatrix[i][j];
                        if(u!=-1 && u<i) {
                            sum+=L.Get(i,u)*U.Get(u,m);
                        }
                    }
                    U.Set(i,m,(A.Get(i,m)-sum));
                    drop+=sum;
                }
            }
        }
        U.Set(i,i,(U.Get(i,i)-drop));
    }
}