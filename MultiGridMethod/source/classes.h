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
double f(double);
double g(double);
vector<double> operator-(const vector<double>&,const vector<double>&);
void operator+=(vector<double>&,const vector<double>&);
vector<double> operator*(double, vector<double>);
double operator|(const std::vector<double>&,const std::vector<double>&);

class Matrix
{
    public:
        virtual int Size()=0;
        virtual double Get(int,int)=0;
        virtual void Resize(int)=0;
        virtual vector<double> operator*(const vector<double>&)=0;
        void PrintMatrix();
};

class PoissonMatrix1D : public Matrix
{
    private:
        int dim;
        int n;
        double diagonal;
        double tridiagonal;
        double identity;
    public:
        PoissonMatrix1D(int,double,double);
        ~PoissonMatrix1D();
        int Size();
        double Get(int, int);
        void Resize(int);
        vector<double> operator*(const vector<double>&);
};

class PoissonMatrix2D : public Matrix
{
    private:
        int dim;
        int n;
        double diagonal;
        double tridiagonal;
        double identity;
    public:
        PoissonMatrix2D(int,double,double,double);
        ~PoissonMatrix2D();
        int Size();
        double Get(int, int);
        void Resize(int);
        vector<double> operator*(const vector<double>&);
};

class Vector {
    public:
        void PrintVector(const vector<double>&);
        virtual void WriteToFile(const vector<double>&)=0;
};

class PoissonVector1D : public Vector {
    private:
        int dim;
        int n;
        double h;
    public:
        PoissonVector1D(int);
        ~PoissonVector1D();
        void InitRightSide(vector<double>&);
        void WriteToFile(const vector<double>&);
};

class PoissonVector2D : public Vector  {
    private:
        int dim;
        int n;
        double h;
    public:
        PoissonVector2D(int);
        ~PoissonVector2D();
        void InitRightSide(vector<double>&);
        void WriteToFile(const vector<double>&);
};

class Algorithms {
	private:
		int dim;
		int n;
        double h;
	public:
		Algorithms(int);
		~Algorithms();
        void JacobiMethod(Matrix&,vector<double>&,const vector<double>&,int);
        void JacobiRelaxationMethod(Matrix&,vector<double>&,const vector<double>&,int);
        void JacobiRelaxationSolver(Matrix&,vector<double>&,const vector<double>&);
        vector<double> MultiGridAlgorithm(Matrix&,vector<double>&,const vector<double>&,int);
        void MultiGridMethod(Matrix&,vector<double>&,vector<double>&);
        void Restriction(const vector<double>&,vector<double>&,int);
        void Interpolation(const vector<double>&,vector<double>&,int);
};

#endif