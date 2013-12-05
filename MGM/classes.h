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
        vector<vector<int> > HashMatrix;
        virtual int Size()=0;
        virtual double Get(int,int)=0;
        void PrintMatrix();
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
        void Resize(int);
};

class Operators
{
	public:
		Operators();
		~Operators();
		double innerProduct(const vector<double>&,const vector<double>&);
		double vectorNorm(const vector<double>&);
		void MatrixVectorMultiplyer(Matrix&,const vector<double>&,vector<double>& y);
		double f(double, double);
    	double g(double, double);
};

class Vectors {
    private:
        int dim;
        int n;
    public:
        vector<vector<double> > x;
        vector<vector<double> > b;
        Vectors(int, Operators&);
        ~Vectors();
        void PrintVector();
        void WriteToFile(Operators&);
};

class Algorithms {
	public:
		Algorithms();
		~Algorithms();
        vector<vector<int> > HashMatrix;
        void InitHashMatrix(int);
        void JacobiMethod(Matrix&,Operators&,vector<double>&,const vector<double>&,int maxIterations=5000);
        void GaussSeidelMethod(Matrix&,Operators&,vector<double>&,const vector<double>&,int maxIterations=5000);
        void SORMethod(Matrix&,Operators&,vector<double>&,const vector<double>&,int maxIterations=5000);
        vector<double> MultiGridMethod(vector<double>&,const vector<double>&,int,Operators&);
        vector<double> Restriction(vector<double>&,int);
        vector<double> Interpolation(vector<double>&,int);
};

#endif