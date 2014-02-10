#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <vector>
#include <omp.h>
using namespace std;

double f(double,double);
double g(double,double);
vector<double> operator-(const vector<double>&,const vector<double>&);
void operator+=(vector<double>&,const vector<double>&);
vector<double> operator*(double, vector<double>);
double operator|(const std::vector<double>&,const std::vector<double>&);

class Matrix
{
    public:
        virtual int Size()=0;
        virtual double Get(int,int)=0;
        virtual vector<double> operator*(const vector<double>&)=0;
        void PrintMatrix();
        vector<vector<int> > HashMatrix;
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
        PoissonMatrix(int,double,double,double);
        ~PoissonMatrix();
        int Size();
        double Get(int, int);
        vector<double> operator*(const vector<double>&);
        // vector<vector<int> > HashMatrix;
};

class PoissonVector {
    private:
        int n;
        double h;
    public:
        vector<double> x;
        vector<double> b;
        vector<double> solved;
        PoissonVector(int,double);
        ~PoissonVector();
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
        void JacobiRelaxationMethod(Matrix&,vector<double>&,const vector<double>&,int);
        void JacobiSolver(Matrix&,vector<double>&,const vector<double>&);
        vector<double> MultiGridAlgorithm(Matrix&,vector<double>&,const vector<double>&,int);
        void MultiGridMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        void Restriction(const vector<double>&,vector<double>&,int);
        void Interpolation(const vector<double>&,vector<double>&,int);
};

#endif