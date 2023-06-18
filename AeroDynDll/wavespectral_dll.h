#pragma once
#ifndef WAVESPECTRAL_H
#define WAVESPECTRAL_H
#include<iostream>
#include<cmath>
#include<string>
#include<stdlib.h>
#include<Windows.h>
#include <fstream>
#include<algorithm>
#include<random>
#include<time.h>
#include"pch.h"
#include"wavespectral_dll.h"
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
	double F2W;
	double TwoPi;
	

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
	double* ec1;
	double* ec2;

};

class wavespectral2 :public data
{
private:
	double Hs;
	double Tp;
	double DT;
	double sigma;
	double S;

public:
	wavespectral2(double Hs, double Tp, double DT);
	double Spectral(double omega);

};

class readmyfile
{
public:
	double** data;
	int line;
	double* freq;
	double* period;
	double* Angle;
	int F_n;
	int Angle_n;
	readmyfile();
	void print();
};

class readrandom
{
public:
	int line;
	double** data;
	readrandom();
};

class Find
{
private:
	double* Angle;
	int Angle_n;
	double theta;//»Î…‰Ω«
public:
	Find (double* Angle, int Angle_n, double theta);
	void find_data();
	int p;
	int q;
	double percent;
};

class waveNumber
{
private:
	double g = 9.81;
	double h;
	double omega;
	double ec=0.000001;//æ´∂»

public:
	waveNumber(double h);
	double newton(double w);
	double F(double k);
	double F1(double k, double dk);
};

class writeOut
{
public:
	writeOut(double time, double motion[6], double RoterSpeed, double GenTorq, double Fairlead_Tension[9], double BaseLoads[6], double TowerTopLoads[6]);

};

class writeOut_Blade_force
{
private:
	double time;
	int numblade;
	double blade1_x[19];
	double blade1_y[19];
	double blade1_z[19];

public:
	writeOut_Blade_force(double time, int num_blade, double blade1_x[19], double blade1_y[19], double blade1_z[19]);
	int write();
	
};

class writeOut_AeroDyn_force
{
private:
	double time;
	double blade_a[3][19];
	double blade_b[3][19];
	double blade_c[3][19];

public:
	writeOut_AeroDyn_force(double time, double blade_a[3][19], double blade_b[3][19], double blade_c[3][19]);
	int write();
};




class morison_velocity
{
private:
	double pi = 3.1415926;
	double delta_omega;
	double omega;
	double S;
	double AW;
	double f1;
	double f2;
	double f3;
	double f4;
	double f5;
	double zeta = 0;
	double zeta1 = 0;
	double velocity[3] = { 0 };
	double distance;
	double k_w;
	double temp1;
	double temp2;
	int N;
	double incide_angle;
	double motion[3] = { 0 };

public:
	morison_velocity(double time, double motion[3], double WaveTmax, double WaveDT, double waterDepth, double Hs, double Tp, double U1[15000], double U2[15000], double incide_angle);
	~morison_velocity();
	double* Velocity();
	int checkValue(double value);
	void morison_velocity2(double time, double motion[3], double WaveTmax, double WaveDT, double waterDepth, double Hs, double Tp, double U1[15000], double U2[15000], double incide_angle);
};


#endif // !WAVESPECTRAL_H1



