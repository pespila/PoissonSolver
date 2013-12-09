#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <vector>
using namespace std;

class Matrix
{
    public:
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
        vector<vector<int> > HashMatrix;
        int Size();
    	double Get(int, int);
        void Resize(int);
};

class Vectors {
    private:
        int dim;
        int n;
        double h;
    public:
        vector<double> x;
        vector<double> b;
        Vectors(int);
        ~Vectors();
        void PrintVector();
        void WriteToFile();
        void Resize(int);
        double f(double,double);
        double g(double,double);
};

class Algorithms {
	private:
		int dim;
		int n;
        double h;
	public:
		Algorithms(int);
		~Algorithms();
        void MatrixVectorMultiplyer(PoissonMatrix&,const vector<double>&,vector<double>& y);
        void JacobiMethod(PoissonMatrix&,vector<double>&,const vector<double>&,int);
        void GaussSeidelMethod(PoissonMatrix&,vector<double>&,const vector<double>&,int);
        void SORMethod(PoissonMatrix&,vector<double>&,const vector<double>&,int);
        vector<double> MultiGridMethod(PoissonMatrix&,Vectors&,vector<double>&,vector<double>&);
        vector<double> Restriction(vector<double>&,int);
        vector<double> Interpolation(vector<double>&,int);
};

#endif