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
};

class Vector {
    public:
        void PrintVector(const vector<double>&);
        virtual void WriteToFile(const vector<double>&)=0;
};

class PoissonVector : public Vector  {
    private:
        int dim;
        int n;
        double h;
    public:
        PoissonVector(int);
        ~PoissonVector();
        void InitRightSide(vector<double>&);
        void WriteToFile(const vector<double>&);
        void InitSolution(vector<double>&);
};

class Algorithms {
	private:
		int dim;
		int n;
        double h;
	public:
		Algorithms();
		~Algorithms();
        void JacobiMethod(Matrix&,vector<double>&,const vector<double>&,int);
        void JacobiRelaxationMethod(Matrix&,vector<double>&,const vector<double>&,int);
        void JacobiRelaxationSolver(Matrix&,vector<double>&,const vector<double>&);
        vector<double> MultiGridAlgorithm(Matrix&,vector<double>&,const vector<double>&);
        void MultiGridMethod(Matrix&,vector<double>&,const vector<double>&,const vector<double>&);
        void Restriction(const vector<double>&,vector<double>&,int);
        void Interpolation(const vector<double>&,vector<double>&,int);
};

#endif