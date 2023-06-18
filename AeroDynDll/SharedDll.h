#pragma once
#include<iostream>
#include<cmath>
#include<string>
#include<stdlib.h>
#include<Windows.h>
#include"pch.h"
#include"SharedDll.h"
#include"wavespectral_dll.h"
using namespace std;

extern "C" _declspec(dllexport) int WF_GT1(double input_theta, double t, double* towerRootMotion, double* towerRootAngle, double* towerRootVel, double* towerRootAngleVel, double* fa, double* fb, double* fc, double* ft, double* m);
//extern "C" _declspec(dllexport) int WF_GT1(double a, double t, double *f1, double* f2, double* f3, double* m);
extern "C" _declspec(dllexport) int GETTHETA(double* theta2, double* towerRootMotion, double* towerRootAngle, double* towerRootVel, double* towerRootAngleVel);
//extern "C" _declspec(dllexport) int GETTHETA(double* theta2);
extern "C" _declspec(dllexport) int WRFO1(double* time, double Fa[19][3], double Fb[19][3], double Fc[19][3], double FT[10][3], double* moment);
extern "C" _declspec(dllexport) int Start();
extern "C" _declspec(dllexport) int retard_cal(double acce, double t, double* val_retard);
extern "C" _declspec(dllexport) int retard_int_cal(double time, double tp, double acce[6], double* val);
extern "C" _declspec(dllexport) int storeI();
extern "C" _declspec(dllexport) int excitforce(double forwards_angle, double time, double motion[6], double* F);
extern "C" _declspec(dllexport) int fluitVelocity(double time, double motion[3], double* Velocity);
extern "C" _declspec(dllexport) int PrintOut(double time, double motion[6], double RoterSpeed, double GenTorq, double Fairlead_Tension[9], double BaseLoads[6], double TowerTopLoads[6], double *xx);
extern "C" _declspec(dllexport) int PrintBladeForce(double time, int blade_num, double blade_x[19], double blade_y[19], double blade_z[19], int *xx);

//extern "C" _declspec(dllimport) int RETARDA(double* val);