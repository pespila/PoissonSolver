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

double f(double,double);
double g(double,double);
std::vector<double> operator-(const std::vector<double>&, const std::vector<double>&);
std::vector<double> operator*(double, std::vector<double>);
void operator+=(std::vector<double>&,const std::vector<double>&);
double operator*(const std::vector<double>&,const std::vector<double>&);
double operator|(const std::vector<double>&,const std::vector<double>&);

class Matrix {
    public:
        vector<vector<int> > HashMatrix;
        virtual int Size()=0;
        virtual double Get(int,int)=0;
        void PrintMatrix();
        virtual vector<double> operator*(const vector<double>&)=0;
};

class WriteableMatrix : public Matrix {
    public:
        virtual void Set(int,int,double)=0;
};

class PoissonMatrix : public Matrix {
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
        vector<double> operator*(const vector<double>&);
};

class LowerMatrix : public WriteableMatrix {
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
        double Get(int,int);
        void Set(int,int,double);
        vector<double> operator*(const vector<double>&);
};

class UpperMatrix : public LowerMatrix {
	private:
		int dim;
		int n;
	public:
		UpperMatrix(int);
		~UpperMatrix();
        double Get(int,int);
		void Set(int,int,double);
        vector<double> operator*(const vector<double>&);
};

class Vectors {
    public:
        virtual double Get(int)=0;
        virtual int Size()=0;
        void PrintVector();
        void WriteToFile();
};

class Boundary : public Vectors {
    private:
        int dim;
        int n;
        double h;
    public:
        Boundary(int);
        ~Boundary();
        vector<double> b;
        double Get(int);
        int Size();
};

class Startvector : public Vectors {
    private:
        int dim;
        int n;
        double value;
    public:
        Startvector(int,double);
        ~Startvector();
        vector<double> x;
        double Get(int);
        int Size();
};

class Algorithms {
	private:
		int dim;
		int n;
        double h;
	public:
		Algorithms(int);
		~Algorithms();
        void JacobiMethod(Matrix&,vector<double>&,const vector<double>&);
        void JacobiRelaxationMethod(Matrix&,vector<double>&,const vector<double>&);
        void GaussSeidelMethod(Matrix&,vector<double>&,const vector<double>&);
        void SORMethod(Matrix&,vector<double>&,const vector<double>&);
		void CG(Matrix&,vector<double>&,const vector<double>&);
		void PCG(Matrix&,WriteableMatrix&,WriteableMatrix&,vector<double>&,const vector<double>&);
        void modifiedIncompleteLU(Matrix&, WriteableMatrix&, WriteableMatrix&);
        void incompleteLU(Matrix&,WriteableMatrix&,WriteableMatrix&);
        void LUsolverUpper(Matrix&,Matrix&,vector<double>&);
};

#endif