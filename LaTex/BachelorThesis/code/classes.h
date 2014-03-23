#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <vector>
#include <assert.h>
#include <time.h>
using namespace std;

vector<double> operator-(const vector<double>&,const vector<double>&);
void operator+=(vector<double>&,const vector<double>&);
vector<double> operator*(double,vector<double>);
double operator|(const vector<double>&,const vector<double>&);
double operator*(const vector<double>&,const vector<double>&);

class Matrix {
    public:
        vector<vector<int> > HashMatrix;
        virtual int Size()=0;
        virtual double Get(int,int)=0;
        virtual vector<double> operator*(const vector<double>&)=0;
};

class WriteableMatrix : public Matrix {
    public:
        virtual void Set(int,int,double)=0;
};

class PoissonMatrix : public Matrix {
    private:
        int dim;
        int n;
        double diagonal;
        double tridiagonal;
        double identity;
    public:
        PoissonMatrix(int);
        ~PoissonMatrix();
        int Size();
        double Get(int, int);
        void InitHashMatrix();
        vector<double> operator*(const vector<double>&);
};

class LowerMatrix : public WriteableMatrix {
    private:
        int dim;
        int n;
    public:
        LowerMatrix(int);
        ~LowerMatrix();
        vector<double> diagonal;
        vector<double> tridiagonal;
        vector<double> identity;
        int Size();
        double Get(int,int);
        void Set(int,int,double);
        vector<double> operator*(const vector<double>&);
};

class UpperMatrix : public LowerMatrix {
   private:
      int dim;
      int n;
   public:
      UpperMatrix(int);
      ~UpperMatrix();
        double Get(int,int);
      void Set(int,int,double);
        vector<double> operator*(const vector<double>&);
};

class Vectors {
    public:
        virtual double Get(int)=0;
        virtual int Size()=0;
        void WriteToFile();
        double f(double,double,int);
        double g(double,double,int);
};

class Boundary : public Vectors {
    private:
        int dim;
        int n;
        double h;
        int k;
    public:
        Boundary(int,int);
        ~Boundary();
        vector<double> b;
        vector<double> solved;
        double Get(int);
        int Size();
};

class Startvector : public Vectors {
    private:
        int dim;
        int n;
        double value;
    public:
        Startvector(int,double);
        ~Startvector();
        vector<double> x;
        double Get(int);
        int Size();
};

class Algorithms {
   private:
        int n;
        int dim;
        double h;
        int Vcounter;
   public:
        Algorithms(int);
        ~Algorithms();
        int JacobiMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        int JacobiRelaxationMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        int GaussSeidelMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        int SORMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        int CG(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        int PCG(Matrix&,WriteableMatrix&,WriteableMatrix&,vector<double>&,const vector<double>&,const vector<double>&);
        void modifiedIncompleteCholesky(Matrix&, WriteableMatrix&, WriteableMatrix&);
        void incompleteLU(Matrix&,WriteableMatrix&,WriteableMatrix&);
        void LUsolverUpper(Matrix&,Matrix&,vector<double>&);

        void JacobiRelaxation(Matrix&,vector<double>&,const vector<double>&,int);
        void Restriction(const vector<double>&,vector<double>&,int);
        void Interpolation(const vector<double>&,vector<double>&,int);
        vector<double> Cycle(Matrix&,vector<double>&,const vector<double>&,int,int,Matrix&,WriteableMatrix&,WriteableMatrix&);
        int MultiGridMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&,int);
        void PCGdirect(Matrix&,WriteableMatrix&,WriteableMatrix&,vector<double>&,const vector<double>&);
        void CGdirect(Matrix&,vector<double>&,const vector<double>&);
};

#endif