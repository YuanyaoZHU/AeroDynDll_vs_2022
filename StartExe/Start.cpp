#include <iostream>
#include <Windows.h>
#include "Start.h"
#include"wavespectral.h"
using namespace std;

//#pragma comment (lib,"../x64/Release/AeroDynDll.lib")

int main()
{
	HANDLE event1 = CreateEvent(NULL, FALSE, FALSE, L"写入角度");
	HANDLE event2 = CreateEvent(NULL, FALSE, FALSE, L"写入外力");
	//HANDLE event3 = CreateEvent(NULL, FALSE, FALSE, L"结束程序");
	cout << "启动程序..." << endl;

	double a = 1.2;
	double b = 1.3;
	double c ;

	HINSTANCE h = LoadLibrary(L"AeroDynDll.dll");
	typedef int (*START)();
	typedef int (*RETARD_CAL)(double a, double b, double* c);
	START start = (START)GetProcAddress(h, "Start");
	RETARD_CAL retard_cal = (RETARD_CAL)GetProcAddress(h, "retard_cal");


	retard_cal(a, b, &c);
/*
    /*--------------------------调用生成波浪谱--------------------------*
    /*-----------------调用动态链接库（显示链接）---------------------*
    typedef void(_stdcall* WAVESPECTRALGENERATOR)(double* H_s, double* T_p, double S[201],double S2[200], double Omega[201], double* dm);
    HINSTANCE hLibrary = LoadLibrary(L"timeDomainSolution_by_fortran.dll");
    if (hLibrary == NULL)
    {
        cout << "can't find the dll file" << endl;
        return -1;
    }
    else
    {
        cout << "the dll file has been found!" << endl;
    }
    WAVESPECTRALGENERATOR waveSpectralGenerator = (WAVESPECTRALGENERATOR)GetProcAddress(hLibrary, "WAVESPECTRALGENERATOR");
    if (waveSpectralGenerator == NULL)
    {
        cout << "can't find the waveSpectralGenerator function file." << endl;
        return -2;
    }
    else
    {
        cout << "The waveSpectralGenerator funtion file has been found!" << endl;
    }
    /*----------------链接完成，开始对各参数进行赋值----------------*
    double H_s = 3.2;
    double T_p = 7.9;
    
    double* Hs=&H_s;
    double* Tp=&T_p;
    double S[201] = {0};
    double Omega[201] = { 0 };
    double S2[200] = { 0 };
    double dm;

    waveSpectralGenerator(Hs, Tp, S, S2, Omega,&dm);

 

    /*--------------------------调用生成波浪--------------------------*




    /*-----------------调用动态链接库（显示链接）---------------------*
    typedef void(_stdcall* WAVEGENERATOR)(double* H_s, double* T_p, double* T,  double* dt, int* N, double* x,  double S[2001], double S2[2000], double Omega[2001], double* wave,double* T_i);
    HINSTANCE hLibrary = LoadLibrary(L"timeDomainSolution_by_fortran.dll");
   
    if (hLibrary == NULL)
    {
        cout << "can't find the wavegenerator dll file" << endl;
        return -1;
    }
    else
    {
        cout << "the dll file has been found!" << endl;
    }
    
    WAVEGENERATOR waveGenerator = (WAVEGENERATOR)GetProcAddress(hLibrary, "WAVEGENERATOR");
    if (waveGenerator == NULL)
    {
        cout << "can't find the waveGenerator function file." << endl;
        return -2;
    }
    else
    {
        cout << "The waveGenerator funtion file has been found!" << endl;
    }
    /*----------------链接完成，开始对各参数进行赋值----------------*
    double TT = 3000;
    double dtt = 0.02;
    int Num_T;
    Num_T = int(TT / dtt);
    double xx = 1;
    double H_s = 3.2;
    double T_p = 7.9;
    
    double* H_s1=&H_s;
    double* T_p1=&T_p;
    double* T=&TT;
    double* dt=&dtt;
    double* x=&xx;
    double SS[2001] = {0};
    double SS2[2000] = {0};
    double Omega2[2001] = {0};
    //double wave[] = { 0 };
    //double T_i[150000] = { 0 };

    double* wave = new double[Num_T];
    double* T_i = new double[Num_T];
   /* for (int f = 0;f < Num_T;f++)
    {
        wave[f] = 0;
        T_i[f] = 0;
    }
    */
    /*
    for (int f = 0;f <= 20;)
    {
        cout << "wave[" << f << "]=" << wave[f] << "T_i[" << f << "]=" << T_i[f] << endl;
        f = f + 1;
    }
    for (int f = 149991;f < Num_T;)
    {
        cout << "wave[" << f << "]=" << wave[f] << "T_i[" << f << "]=" << T_i[f] << endl;
        f = f + 1;
    }
    cout << wave[Num_T - 1] << ',' << *(wave + Num_T - 2) << ',' << *(wave + 1) << ',' << *wave << endl;
    cout << T_i << ',' << T_i[Num_T - 2] << ',' << T_i[1] << ',' << T_i[0] << endl;

    *

    waveGenerator(H_s1, T_p1, T, dt, &Num_T, x, SS, SS2, Omega2, wave,T_i);
   
    cout << SS[36] << ',' << SS2[36] << endl;
    for (int f = 149991;f < Num_T;)
    {
        cout << "wave[" << f << "]=" << wave[f] << "    T_i[" << f << "]=" << T_i[f] << endl;
        f = f + 1;
    }

    /*-------------------------获取文件行数---------------------------*
    /*-----------------调用动态链接库（显示链接）---------------------*
    typedef void(_stdcall* NUMOFLINE)(int *N);
    /* 
    HINSTANCE hLibrary = LoadLibrary(L"timeDomainSolution_by_fortran.dll");

    if (hLibrary == NULL)
    {
        cout << "can't find the dll file" << endl;
        return -1;
    }
    else
    {
        cout << "the dll file has been found!" << endl;
    }
    *

    NUMOFLINE NumOfLine = (NUMOFLINE)GetProcAddress(hLibrary, "NUMOFLINE");
    if (NumOfLine == NULL)
    {
        cout << "can't find the NumOfLine function file." << endl;
        return -2;
    }
    else
    {
        cout << "The NUMOFLINE funtion file has been found!" << endl;
    }
    /*------------------动态链接库调用完成，定义变量--------------------*

    int Num_Line;
    int* Num_L = &Num_Line;

    NumOfLine(&Num_Line);
    
    /*------------------------------------------------------------------*

    /*------------------获取水动力系数----------------------------------*
    /*------------------加载动态链接库----------------------------------*
    typedef void(_stdcall* READDIFFRACTIONDATA)(int* N,double * dat);

    READDIFFRACTIONDATA ReadDiffractionData = (READDIFFRACTIONDATA)GetProcAddress(hLibrary, "READDIFFRACTIONDATA");
    
    if (ReadDiffractionData == NULL)
    {
        cout << "can't find the ReadDiffractionData function file." << endl;
        return -2;
    }
    else
    {
        cout << "The ReadDiffractionData funtion file has been found!" << endl;
    }
    /*-------------------动态链接库加载完成，输入变量-------------------*

    
    double* dat = new double [Num_Line*7];
    
    cout << "flag::1" << endl;

    ReadDiffractionData(&Num_Line, dat);

   
    double** dat1 = new double* [Num_Line];
        for (int i = 0; i < Num_Line; i++)
        {
            dat1[i] = new double[7];
        }
        int ii = 0;
        for (int i = 0;i < Num_Line;i++)
        {
            for (int j = 0;j < 7;j++)
            {
                dat1[i][j] = dat[ii];
                ii++;
            }
            
        }

    cout << "数据读取完毕，开始输出" << endl;
    for (int i = 0;i < Num_Line;i++)
    {
        cout << dat1[i][0] << "     " << dat1[i][1] << "     " << dat1[i][2] << "     " << dat1[i][3] <<"     "<<dat1[i][4]  << "     " << dat1[i][5] << "     " << dat1[i][6] << endl;
    }


   
	start();
	//Start();

	CloseHandle(event1);
	CloseHandle(event2);
	FreeLibrary(h);
    FreeLibrary(hLibrary);

    delete[] wave;
    delete[] T_i;
    delete[] dat;
    for (int i = 0; i < Num_Line; i++) {
        delete[] dat1[i];
    }
    delete[]dat1;

 */
/*---------------------创建波谱-----------------------

double* p = new double[200];
double Hs = 3.2;
double Tp = 7.9;
readmyfile two;
wavespectral one(Hs, Tp, two.F_n,two.freq);
double* s;
s = one.Spectral();
for (int i = 0; i < two.F_n;i++)
{
    cout << s[i] << endl;
}

one.print();
*/
    typedef int (*STOREI)();
    STOREI storeI = (START)GetProcAddress(h, "storeI");//读取水动力系数，并创建波浪谱
    storeI();
    start();

	return 0;
}
