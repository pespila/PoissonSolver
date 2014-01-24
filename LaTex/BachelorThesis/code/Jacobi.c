void Algorithms::JacobiMethod(Matrix& A, vector<double>& x, const vector<double>& b, int steps) {
    int dim=A.Size();
    vector<double> r(dim);
    for(int j=0;j<steps;j++) {
        r=A*x;
        for(int i=0;i<dim;i++) {
            r[i]=1.0/4.0*(b[i]-r[i]);
            x[i]+=r[i];
        }
    }
}
