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
        virtual void PrintMatrix();
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
        void PrintMatrix();
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
        void PrintMatrix();
};

class UpperMatrix : public LowerMatrix
{
	private:
		int dim;
		int n;
	public:
		UpperMatrix(int);
		~UpperMatrix();
        double Get(int,int);
		void Set(int,int,double);
        void PrintMatrix();
};

class Vectors {
    public:
        virtual int Size()=0;
        virtual void PrintVector()=0;
        virtual void WriteToFile()=0;
};

class PoissonVector : public Vectors {
    private:
        int dim;
        int n;
        double h;
    public:
        PoissonVector(int m);
        ~PoissonVector();
        int Size();
        void PrintVector();
        void WriteToFile();
        double f(double,double);
        double g(double,double);
        vector<double> x;
        vector<double> b;
};

class MGMVector : public Vectors {
    private:
        int dim;
        int n;
        double h;
    public:
        MGMVector(int m);
        ~MGMVector();
        int Size();
        void PrintVector();
        void WriteToFile();
        double f(double,double);
        double g(double,double);
        vector<vector<double> > b;
        vector<vector<double> > x;
};

class Operators {
    public:
        virtual double innerProduct(const vector<double>&,const vector<double>&)=0;
        virtual double vectorNorm(const vector<double>&)=0;
        virtual void MatrixVectorMultiplyer(Matrix&,const vector<double>&,vector<double>&)=0;
};

class CGOperators : public Operators
{
    public:
        PoissonOperators();
        ~PoissonOperators();
        double innerProduct(const vector<double>&,const vector<double>&);
        double vectorNorm(const vector<double>&);
        void MatrixVectorMultiplyer(Matrix&,const vector<double>&,vector<double>& y);
        double f(double, double);
        double g(double, double);
};

class MGMOperators : public Operators
{
    public:
        MGMOperators();
        ~MGMOperators();
        double innerProduct(const vector<double>&,const vector<double>&);
        double vectorNorm(const vector<double>&);
        void MatrixVectorMultiplyer(Matrix&,const vector<double>&,vector<double>& y);
        double f(double, double);
        double g(double, double);
};

class Algorithms {
	public:
		Algorithms();
		~Algorithms();
        vector<vector<int> > HashMatrix;
        void InitHashMatrix(int);
        void incompleteLU(Matrix&,WriteableMatrix&,WriteableMatrix&,Operators&);
        void modifiedIncompleteLU(Matrix&,WriteableMatrix&,WriteableMatrix&,Operators&);
        void LU(Matrix&,WriteableMatrix&,WriteableMatrix&,Operators&);
		void modifiedIncompleteCholesky(WriteableMatrix&,WriteableMatrix&,WriteableMatrix&,Operators&);
		void incompleteCholesky(PoissonMatrix&,LowerMatrix&,UpperMatrix&,Operators&);
        void LUsolverLower(Matrix&,Matrix&,vector<double>&,Operators&);
        void LUsolverUpper(Matrix&,Matrix&,vector<double>&,Operators&);
        void JacobiMethod(Matrix&,Operators&,vector<double>&,const vector<double>&,int maxIterations=5000);
        void GaussSeidelMethod(Matrix&,Operators&,vector<double>&,const vector<double>&,int maxIterations=5000);
        void SORMethod(Matrix&,Operators&,vector<double>&,const vector<double>&,int maxIterations=5000);
        void SSORMethod(Matrix&,Operators&,Vectors&);
		void CG(Matrix&,Operators&,Vectors&);
		void PCG(Matrix&,Operators&,WriteableMatrix&,WriteableMatrix&,Vectors&);
        void MultiGridMethod(vector<double>&,const vector<double>&,int,Operators&);
        vector<double> Restriction(vector<double>&,int);
        vector<double> Interpolation(vector<double>&,int);
};

#endif