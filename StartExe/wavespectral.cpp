#include<iostream>
#include <fstream>
#include <cstring>
#include<algorithm>
#include"wavespectral.h"
using namespace std;

/*-----------------------data--------------------------------------*/

data::data()
{
	gama = 3.3;
	sigma_a = 0.07;
	sigma_b = 0.09;
	pi = 3.1415926;
	temp = 0.23 + 0.0336 * gama - 0.185 * pow(1.9 + gama,-1);
	beta_j = 0.06238 *pow(temp, -1) * (1.094 - 0.01915 * log(gama));
	C2F = 1 / (2 * pi);

}

void data::print()
{
	cout << beta_j << endl;
}

/*------------------------------wavespectral----------------------------------*/


wavespectral::wavespectral(double Hs, double Tp, int N, double * f)
{
	this->Hs = Hs;
	this->Tp = Tp;
	this->n = N;
	this->f_p = 1 / Tp;
	this->f = f;
}
double * wavespectral::Spectral()
{
	double* s = new double[n];
	double sigma;
	for (int i = 0; i < n;i++)
	{
		if (f[i] <= f_p)
		{
			sigma = sigma_a;
		}
		else if (f[i] > f_p)
		{
			sigma = sigma_b;
		}
		
		s[i] = beta_j * pow(Hs, 2) * pow(Tp, -4) * pow(f[i], -5) * exp(-(5 / 4) * pow(Tp * f[i], -4)) * pow(gama, exp(-pow((f[i] / f_p) - 1, 2) / (2 * pow(sigma, 2))));

	}
	spec = s;
	return spec;
}

void wavespectral::print2()
{
	cout << beta_j<<"   "<<sigma_a <<"		"<< temp << endl ;
}

/*--------------------------------readmyfile------------------------------------*/

readmyfile::readmyfile()
{
	//统计文本有多少行
	char c;
	this->line = 0;
	fstream finl("Spar.3", ios::in);
	if (!finl)
	{
		cerr << "cannot open hydro parameter file!" << endl;
	}
	while (finl.get(c))
	{
		if (c == '\n')
			line++; //行数加1
	}
	finl.close();
	
	//读取数据
	double** I = new double* [line];
	F_n = 0;
	double item = 0;
	for (int i = 0; i < line; i++)
	{
		I[i] = new double[7];
	}
	ifstream myfile("Spar.3");
	if (!myfile.is_open())
	{
		cout << "can not open this file!" << endl;
	}
	for (int i = 0; i < line;i++)
	{
		for (int j = 0;j < 7;j++)
		{
			myfile >> I[i][j];
		}
	}
	//统计有多少个频率，这些频率是多少
	for (int i = 0; i < line; i++)
	{
		if (I[i][0] != item)
		{
			item = I[i][0];
			//cout << "item = " << item;
			F_n++;
		}
	}
	double* Tp_local = new double[F_n];
	double* freq1 = new double[F_n];
	item = 0;
	int j = 0;
	for (int i = 0; i < line;i++)
	{
		if (I[i][0] != item)
		{
			item = I[i][0];
			Tp_local[j] = item;
			freq1[j] = 1 / Tp_local[j];
			j++;
		}
	}

	freq = freq1;
	data = I;


	//统计有多少中入射角
	item = 0;
	Angle_n = 0;
	for (int i = 0; i < line; i++)
	{
		if (I[i][1] != item && Tp_local[0]==I[i][0])
		{
			item = I[i][1];
			//cout << "item = " << item;
			Angle_n++;
		}
	}
	double* Angle1 = new double[Angle_n];
	item = 0;
	j = 0;
	for (int i = 0; i < Angle_n;i++)
	{
		if (I[i][1] != item)
		{
			item = I[i][1];
			Angle1[j] = item;
			j++;
		}
	}
	Angle = Angle1;


}

void readmyfile::print()
{
	for (int i = 0; i < 3;i++)
	{
		cout << data[i][0] << "    " << data[i][1] << "    " << data[i][2] << "    " << data[i][3] << "    " << data[i][4] << "    " << data[i][5] << "    " << data[i][6] << endl;
	}
	cout <<"F_n = "<< F_n << endl;
	cout << "I[1][0]=" << data[1][0] << endl;

}

/*-------------------------force----------------------*/
force::force(double* Angle, int Angle_n)
{
	this->Angle = Angle;
	this->Angle_n = Angle_n;

}

void force::find_data()
{
	double* dm = new double[Angle_n];
	for (int i = 0;i < Angle_n;i++)
	{
		dm[i] = abs(theta - Angle[i]);
	}

	int p = min_element(dm, dm + Angle_n)-dm;
	int q = 0;
	if (theta>Angle[p])
	{
		q = p + 1;
	}
	else if (theta == Angle[p])
	{
		q = p;
	}
	else if (theta < Angle[p])
	{
		q = p - 1;
	}

	double percent;
	percent = (theta - Angle[p]) / (Angle[q] - Angle[p]);



}