#pragma once
#ifndef WAVESPECTRAL_H
#define WAVESPECTRAL_H
#include<string>

class data
{
protected:
	double sigma_a;
	double sigma_b;
	double beta_j;
	double pi;
	double gama;
	double C2F;
	double temp;
	

public:
	data();
	void print();
};


class wavespectral: public data
{
private:
	double Hs;
	double Tp;
	int	n;
	double f_p;
	double* f;
	double* spec;
	


public:
	wavespectral(double Hs, double Tp, int n, double *f);
	double * Spectral();
	void print2();


};

class readmyfile
{
public:
	double** data;
	int line;
	double* freq;
	double* Angle;
	int F_n;
	int Angle_n;
	readmyfile();
	void print();
};

class force
{
private:
	double* Angle;
	int Angle_n;

public:
	force(double* Angle, int Angle_n);
	double theta;//»Î…‰Ω«
	double outforce;
	void find_data();

};

#endif // !WAVESPECTRAL_H



