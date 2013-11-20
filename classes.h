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

struct Tuple
{
    int i;
    int j;
};

class Matrix
{
    public:
        Tuple hash;
        vector<Tuple> HashKeys;
        virtual int Size()=0;
        virtual double Get(int,int)=0;
        void PrintMatrix();
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
    	double Get(int, int);
};

class Preconditioner : public WriteableMatrix
{
    private:
        int dim;
        int n;
        vector<double> diagonal;
        vector<double> tridiagonal;
        vector<double> identity;
    public:
        Preconditioner(int);
        ~Preconditioner();

        int Size();
        double Get(int,int);
        void Set(int,int,double);
};

class LowerMatrix : public WriteableMatrix
{
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
        double Get(int, int);
        void Set(int, int, double);
};

class UpperMatrix : public LowerMatrix
{
	private:
		int dim;
		int n;
	public:
		UpperMatrix(int);
		~UpperMatrix();
        double Get(int, int);
		void Set(int, int, double);
};

class Operators
{
	private:
		int dim;
		int n;
	public:
		Operators(int);
		~Operators();
		double innerProduct(const vector<double>&,const vector<double>&);
		double vectorNorm(const vector<double>&);
		void MatrixVectorMultiplyer(Matrix&,const vector<double>&,vector<double>& y);
		void LUsolverLower(Matrix&,vector<double>&);
		void LUsolverUpper(Matrix&,vector<double>&);
		double f(double, double);
    	double g(double, double);
};

class Vectors {
    private:
        int dim;
        int n;
    public:
        vector<double> x;
        vector<double> b;
        Vectors(int, Operators&);
        ~Vectors();
        void PrintVector();
        void WriteToFile(Operators&);
};

class Algorithms {
	private:
		int dim;
		int n;
	public:
		Algorithms(Matrix&);
		~Algorithms();
        int HashTable(int,int);
        void incompleteLU(Matrix&, WriteableMatrix&, WriteableMatrix&);
        void modifiedIncompleteLU(Matrix&, WriteableMatrix&, WriteableMatrix&);
		void modifiedIncompleteCholesky(WriteableMatrix&, WriteableMatrix&, WriteableMatrix&);
		void incompleteCholesky(PoissonMatrix&, LowerMatrix&, UpperMatrix&);
        void JacobiMethod(Matrix&, Operators&, Vectors&);
        void GaussSeidelMethod(Matrix&, Operators&, Vectors&);
        void SORMethod(Matrix&, Operators&, Vectors&);
        void SSORMethod(Matrix&, Operators&, Vectors&);
		void CG(Matrix&, Operators&, Vectors&);
		void PCG(Matrix&, Operators&, WriteableMatrix&, WriteableMatrix&, Vectors&);
};

#endif