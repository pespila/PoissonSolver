#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <vector>
using namespace std;

double f(double,double);
double g(double,double);
std::vector<double> operator-(const std::vector<double>&,const std::vector<double>&);
void operator+=(std::vector<double>&,const std::vector<double>&);

class Matrix
{
    public:
        virtual int Size()=0;
        virtual double Get(int,int)=0;
        void PrintMatrix();
        virtual vector<double> operator*(const vector<double>&)=0;
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
        vector<double> operator*(const vector<double>&);
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
};

class Algorithms {
	private:
		int dim;
		int n;
        double h;
	public:
		Algorithms(int);
		~Algorithms();
        double vectorNorm(const vector<double>&);
        void JacobiMethod(Matrix&,vector<double>&,const vector<double>&,int);
        void JacobiRelaxationMethod(Matrix&,vector<double>&,const vector<double>&,int);
        void JacobiRelaxationSolver(Matrix&,vector<double>&,const vector<double>&);
        void GaussSeidelMethod(PoissonMatrix&,vector<double>&,const vector<double>&,int);
        void SORMethod(PoissonMatrix&,vector<double>&,const vector<double>&,int);
        vector<double> MultiGridAlgorithm(PoissonMatrix&,Vectors&,const vector<double>&,int);
        void MultiGridMethod(PoissonMatrix&,Vectors&,vector<double>&,vector<double>&);
        void Restriction(const vector<double>&,vector<double>&,int);
        void Interpolation(const vector<double>&,vector<double>&,Vectors&,int);
};

#endif