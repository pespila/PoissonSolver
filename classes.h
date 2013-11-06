#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <vector>
using namespace std;

class PoissonMatrix
{
    private:
    	int dim;
    	int n;
    	vector<double> diagonal;
    	vector<double> tridiagonal;
    	vector<double> identity;
    public:
    	PoissonMatrix(int);
    	~PoissonMatrix();

    	virtual void set(int, int, double);
    	virtual double get(int, int);
    	virtual int size();
    	void printMat();
};

class PoissonMatrixLite : public PoissonMatrix
{
    private:
        int dim;
        int n;
        double diag;
        double tridiag;
        double id;
    public:
        PoissonMatrixLite(int);
        ~PoissonMatrixLite();
};

class LowerPoissonMatrix : public PoissonMatrix
{
    private:
        int dim;
        int n;
        vector<double> Ldiagonal;
        vector<double> Ltridiagonal;
        vector<double> Lidentity;
    public:
        LowerPoissonMatrix(int);
        ~LowerPoissonMatrix();
        void set(int, int, double);
        double get(int, int);
};

class UpperPoissonMatrix : public PoissonMatrix
{
	private:
		int dim;
		int n;
		vector<double> Udiagonal;
		vector<double> Utridiagonal;
		vector<double> Uidentity;
	public:
		UpperPoissonMatrix(int);
		~UpperPoissonMatrix();
		void set(int, int, double);
		double get(int, int);
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
		void MatrixVectorMultiplyer(PoissonMatrix&,const vector<double>&,vector<double>& y);
		void LUsolverLower(LowerPoissonMatrix&,vector<double>&);
		void LUsolverUpper(UpperPoissonMatrix&,vector<double>&);
		double f(double, double);
    	double g(double, double);
};

class CGVectors {
	private:
		int dim;
		int n;
    public:
    	vector<double> x;
    	vector<double> b;
    	CGVectors(int, Operators&);
    	~CGVectors();
};

class Algorithms {
	private:
		int dim;
		int n;
	public:
		Algorithms(PoissonMatrix&);
		~Algorithms();
		void modifiedIncompleteLU(PoissonMatrix&, LowerPoissonMatrix&, UpperPoissonMatrix&);
		void modifiedIncompleteCholesky(PoissonMatrix&, LowerPoissonMatrix&, UpperPoissonMatrix&);
		void incompleteLU(PoissonMatrix&, LowerPoissonMatrix&, UpperPoissonMatrix&);
		void incompleteCholesky(PoissonMatrix&, LowerPoissonMatrix&, UpperPoissonMatrix&);
		void CG(PoissonMatrix&, Operators&, vector<double>&, const vector<double>&);
		void PCG(PoissonMatrix&, Operators&, LowerPoissonMatrix&, UpperPoissonMatrix&, vector<double>&, const vector<double>&);
};

#endif