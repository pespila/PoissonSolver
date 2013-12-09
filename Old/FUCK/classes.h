#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <vector>
#include <time.h>
using namespace std;

class Matrix
{
    public:
        virtual int Size()=0;
        virtual double Get(int,int)=0;
        void PrintMatrix(int);
};

class WriteableMatrix : public Matrix
{
    public:
        virtual void Set(int,int,double)=0;
};

class PoissonMatrix : public Matrix
{
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
    	double Get(int,int);
        void Resize(int);
};

class LowerMatrix : public WriteableMatrix
{
    private:
        int dim;
        int n;
        vector<double> diagonal;
        vector<double> tridiagonal;
        vector<double> identity;
    public:
        LowerMatrix(int);
        ~LowerMatrix();
        int Size();
        double Get(int,int);
        void Set(int,int,double);
};

class UpperMatrix : public LowerMatrix
{
	private:
		int dim;
		int n;
        vector<double> diagonal;
        vector<double> tridiagonal;
        vector<double> identity;
	public:
		UpperMatrix(int);
		~UpperMatrix();
        double Get(int,int);
		void Set(int,int,double);
        void PrintMatrix(int);
};

class PoissonVector {
    private:
        int dim;
        int n;
        double h;
    public:
        PoissonVector(int);
        ~PoissonVector();
        void PrintVector();
        void WriteToFile();
        double f(double,double);
        double g(double,double);
        vector<double> x;
        vector<double> b;
};

class PoissonOperators {
    public:
        double innerProduct(const vector<double>&,const vector<double>&);
        double vectorNorm(const vector<double>&);
        void MatrixVectorMultiplyer(Matrix&,const vector<double>&,vector<double>&);
        double f(double,double);
        double g(double,double);
};

class Algorithms {
    public:
        Algorithms();
        ~Algorithms();
        vector<vector<int> > HashMatrix;
        void InitHashMatrix(int);
        void incompleteLU(Matrix&,WriteableMatrix&,WriteableMatrix&,PoissonOperators&);
        void modifiedIncompleteLU(Matrix&,WriteableMatrix&,WriteableMatrix&,PoissonOperators&);
        void LU(Matrix&,WriteableMatrix&,WriteableMatrix&,PoissonOperators&);
        void modifiedIncompleteCholesky(WriteableMatrix&,WriteableMatrix&,WriteableMatrix&,PoissonOperators&);
        void incompleteCholesky(PoissonMatrix&,LowerMatrix&,UpperMatrix&,PoissonOperators&);
        void LUsolverLower(Matrix&,Matrix&,vector<double>&,PoissonOperators&);
        void LUsolverUpper(Matrix&,Matrix&,vector<double>&,PoissonOperators&);
        void JacobiMethod(Matrix&,PoissonOperators&,vector<double>&,const vector<double>&,int maxIterations=5000);
        void GaussSeidelMethod(Matrix&,PoissonOperators&,vector<double>&,const vector<double>&,int maxIterations=5000);
        void SORMethod(Matrix&,PoissonOperators&,vector<double>&,const vector<double>&,int maxIterations=5000);
        void SSORMethod(Matrix&,PoissonOperators&,PoissonVector&);
        void CG(Matrix&,PoissonOperators&,PoissonVector&);
        void PCG(Matrix&,PoissonOperators&,WriteableMatrix&,WriteableMatrix&,PoissonVector&);
        void MultiGridMethod(vector<double>&,const vector<double>&,int,PoissonOperators&);
        vector<double> Restriction(vector<double>&,int);
        vector<double> Interpolation(vector<double>&,int);
};

#endif