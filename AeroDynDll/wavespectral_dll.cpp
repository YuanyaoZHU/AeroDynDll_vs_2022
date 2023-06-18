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
	C2F = 1.0 / (2.0 * pi);
	F2W = 2.0 * pi;
	TwoPi = 2.0 * pi;

}

void data::print()
{
	cout << beta_j << endl;
}

/*------------------------------wavespectral----------------------------------*/


wavespectral::wavespectral(double Hs, double Tp, int N, double * f)//NΪƵ�ʸ�����fΪƵ��
{
	cout << "����wavespctral���ڲ�" << endl;
	this->Hs = Hs;
	this->Tp = Tp;
	this->n = N;
	this->f_p = 1 / Tp;
	this->f = f;
	this->spec = 0;
	if (this->Tp / sqrt(this->Hs) <= 3.6)//���ݡ�Dynamics Modeling and Loads Analysis of an Offshore Floating Wind Turbine���ж��ڲ���٤��ֵ��ȷ������ȡֵ
	{
		this->gama = 5;
	}
	else if (this->Tp / sqrt(this->Hs) > 3.6 && this->Tp / sqrt(this->Hs) <= 5)
	{
		this->gama = exp(5.75 - 1.15 * this->Tp / sqrt(this->Hs));
	}
	else
	{
		this->gama = 1;
	}
	//------------------������������������Լ��벨������------------//
	double* ecc1 = new double[2*N];//ע�⣺���new������ĸ����͸�ֵ�ĸ����Բ��ϻᵼ������ϵͳ�����쳣������system���в�������
	double* ecc2 = new double[2*N];
	double seed1;
	double seed2;
	
	default_random_engine random1;
	default_random_engine random2;
	cout << "������U1���ӣ�" << endl;
	cin >> seed1;
	cout << "������U2���ӣ�" << endl;
	cin >> seed2;
	random1.seed(seed1);
	uniform_real_distribution<double> dis1(0.0, 1.0);
	for (int i = 0;i < 2*N; i++)
	{
		ecc1[i] = dis1(random1);
		cout << ecc1[i]  << "    ";
	}
	cout << endl;
	ec1 = ecc1;
	random2.seed(seed2);
	uniform_real_distribution<double> dis2(0.0, 1.0);
	for (int i = 0;i < 2*N;i++)
	{
		ecc2[i] = dis2(random2);
		cout << ecc2[i]  << "    ";
	}
	cout << endl;
	ec2 = ecc2;
	cout << "����������ɹ���" << endl;
	//�����������ϣ�
	
}
double * wavespectral::Spectral()
{
	//double s[100];
	double* s=new double[n];
	double* s1 = new double[n];
	double* s2 = new double[n];
	double* sw=new double[n];
	double sigma;
	ofstream OutSpec("wavespectrum.txt");
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
		//��Ƶ�ʱ�ʾ���׺�Ƶ�ʱ�ʾ���ײ�2*pi���������������˼��乤��Ӧ�á���148ҳ������ݣ�������õ���Ƶ��Ϊ���������ֵ����FAST�����в��õ��ǽ�Ƶ�ʡ�
		s1[i] = 5.0 / 16.0 * pow(Hs, 2) * Tp * pow(f[i] * Tp, -5);
		s2[i] = exp(-(5 / 4) * pow(Tp * f[i], -4)) * (1 - 0.287 * log(gama)) * pow(gama, exp(-0.5 * (pow(((f[i] * Tp - 1) / sigma), 2))));
		//cout << "���׽����ɹ���" << endl;
		//s[i] = beta_j * pow(Hs, 2) * pow(Tp, -4) * pow(f[i], -5) * exp(-(5 / 4) * pow(Tp * f[i], -4)) * pow(gama, exp(-pow((f[i] / f_p) - 1, 2) / (2 * pow(sigma, 2))));
		s[i] = s1[i] * s2[i];
		sw[i] = s[i] / 2 / pi;

		if (OutSpec.is_open())
		{
			OutSpec << f[i] << "    " << sw[i] << endl;
		}
	}
	spec = sw;//��������������sw������Ϊ���ý�Ƶ�ʷ�ʽ�������㣻��������s������Ϊ����Ƶ�ʽ������㣡
	OutSpec.close();
	cout << "���׽����ɹ���" << endl;
	return spec;
}


void wavespectral::print2()
{
	cout << beta_j<<"   "<<sigma_a <<"		"<< temp << endl ;
}



/*-------------------------------wavespectral2-----------------------------------*/
wavespectral2::wavespectral2(double Hs, double Tp, double DT)
{
	//cout << "����wavespctral2���ڲ�" << endl;
	this->Hs = Hs;
	this->Tp = Tp;
	this->DT = DT;
			
}

double wavespectral2::Spectral(double omega)
{
	
		if (omega <= TwoPi/Tp)//sigma��ȡֵ
		{
			this->sigma = sigma_a;
		}
		else if (omega > TwoPi / Tp)
		{
			this->sigma = sigma_b;
		}

		if (Tp / sqrt(Hs) <= 3.6)//gama��ȡֵ
		{
			gama = 5;
		}
		else if (Tp / sqrt(Hs) > 3.6 && Tp / sqrt(Hs) <= 5)
		{
			gama = exp(5.75 - 1.15 * Tp / sqrt(Hs));
		}
		else if (Tp / sqrt(Hs) > 5)
		{
			gama = 1;
		}

	double C0 = C2F * 5.0 / 16.0 * pow(Hs, 2) * Tp;
	double C1 = pow(omega * Tp / TwoPi, -5);
	double C2 = exp(-(5.0 / 4.0) * pow(Tp * omega/TwoPi, -4));
	double C3 = 1.0 - 0.287 * log(gama);
	double C4 = pow(gama, exp(-0.5 * (pow(((omega * Tp/TwoPi - 1) / sigma), 2))));
	this->S = C0 * C1 * C2 * C3 * C4;
	return S;
}

/*--------------------------------readmyfile------------------------------------*/

readmyfile::readmyfile()
{
	//ͳ���ı��ж�����
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
			line++; //������1
	}
	finl.close();
	
	//��ȡ����
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
	//ͳ���ж��ٸ�Ƶ�ʣ���ЩƵ���Ƕ���
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
	period = Tp_local;
	freq = freq1;
	data = I;


	//ͳ���ж����������
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
	cout << "Angle_n = " << Angle_n << endl;
	double* Angle1 = new double[Angle_n];
	item = 0;
	j = 0;
	for (int i = 0; i < line;i++)
	{
		if (I[i][1] != item && Tp_local[0] == I[i][0])
		{
			item = I[i][1];
			Angle1[j] = item;
			//cout << "Angle1[" << j << "] = " << Angle1[j] << endl;
			j++;
		}
	}
	Angle = Angle1;
	myfile.close();


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

/*--------------------------readrandom---------------------------*/
readrandom::readrandom()
{
	//ͳ���ı��ж�����
	char c;
	this->line = 0;
	fstream finl("random_value.txt", ios::in);
	if (!finl)
	{
		cerr << "cannot open random_value.txt file!" << endl;
	}
	while (finl.get(c))
	{
		if (c == '\n')
			line++; //������1
	}
	finl.close();

	//��ȡ����
	double** I = new double* [line];
	double item = 0;
	for (int i = 0; i < line; i++)
	{
		I[i] = new double[2];
	}
	ifstream myfile("random_value.txt");
	if (!myfile.is_open())
	{
		cout << "can not open random_value.txt file!" << endl;

	}
	for (int i = 0; i < line;i++)
	{
		for (int j = 0;j < 2;j++)
		{
			myfile >> I[i][j];
		}
	}

	data = I;

}

/*-------------------------force----------------------*/
Find::Find(double* Angle, int Angle_n, double theta)
{
	this->Angle = Angle;
	this->Angle_n = Angle_n;
	this->theta = theta;

	double* dm = new double[this->Angle_n];
	for (int i = 0;i < this->Angle_n;i++)
	{
		dm[i] = abs(this->theta - this->Angle[i]);
	}

	p = min_element(dm, dm + this->Angle_n) - dm;
	q = 0;
	if (this->theta > this->Angle[p])
	{
		q = p + 1;
	}
	else if (this->theta == this->Angle[p])
	{
		q = p;
	}
	else if (this->theta < this->Angle[p])
	{
		q = p - 1;
	}
	percent = (this->theta - this->Angle[p]) / (this->Angle[q] - this->Angle[p]);


}

void Find::find_data()
{
	double* dm = new double[Angle_n];
	for (int i = 0;i < Angle_n;i++)
	{
		dm[i] = abs(theta - Angle[i]);
	}

	p = min_element(dm, dm + Angle_n)-dm;
	q = 0;
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
	if (theta == Angle[p])
	{
		percent = 1;
	}
	else
	{
		percent = (theta - Angle[p]) / (Angle[q] - Angle[p]);
	}

}

waveNumber::waveNumber(double h)
{
	this->h = h;
}

double waveNumber::newton(double w)
{
	this->omega = w;
	double k1;
	double k0=10;
	double delta_k;
	double f;
	double f1;
	double i=0;
	//cout << "����newton!" << endl;
	do
	{
		f = F(k0);
		f1 = F1(k0, ec);
		k1 = k0-f/f1;
		delta_k = k1 - k0;
		k0 = k1;
		i++;
		if (i >= 100000)
		{
			cout << "newton���������������µ���" << endl;
			system("pause");
			break;
		}
	} while (abs(delta_k) > ec);
	return k1;
}

double waveNumber::F(double k)
{
	double f;
	f = k * tanh(k * h) - pow(omega,2)/ g;
	return f;
}

double waveNumber::F1(double k, double dk)
{
	double f1;
	f1 = (F(k + dk) - F(k)) / dk;
	return f1;
}

writeOut::writeOut(double time, double motion[6], double RoterSpeed, double GenTorq, double Fairlead_Tension[9], double BaseLoads[6], double TowerTopLoads[6])
{
	fstream in;
	if (time == 0)
	{
		in.open("motion.txt", ios::out);
		if (!in.is_open())
		{
			cout << "���ܴ�motion.txt�ļ�" << endl;
			system("pause");
		}
		in.close();
	}
	in.open("motion.txt", ios::app);
	if (!in.is_open())
	{
		cout << "���ܴ�motion.txt�ļ�" << endl;
		system("pause");
	} 
	in << time << "    " << motion[0] << "    " << motion[1] << "    " << motion[2] << "    " << motion[3] << "    " << motion[4] << "    " << motion[5] << "    " << RoterSpeed << "    " << GenTorq ;
	for(int i=0;i<=8;i++)
	{
		in << "    " << Fairlead_Tension[i];
	}
	for (int i = 0;i < 6;i++)
	{
		in << "    " << BaseLoads[i];
	}
	for (int i = 0;i < 6;i++)
	{
		in << "    " << TowerTopLoads[i];
	}
	in << endl;

	in.close();

}

writeOut_AeroDyn_force::writeOut_AeroDyn_force(double time, double blade_a[3][19], double blade_b[3][19], double blade_c[3][19])
{
	this->time = time;
	for (int j = 0;j < 19;j++)
	{
		for (int i = 0;i < 3;i++)
		{
			this->blade_a[i][j] = blade_a[i][j];
			this->blade_b[i][j] = blade_b[i][j];
			this->blade_c[i][j] = blade_c[i][j];
		}
	}
}

int writeOut_AeroDyn_force::write()
{
	fstream in;
	if (time == 0)
	{
		in.open("AeroDyn_blade_force.txt", ios::out);
		if (!in.is_open())
		{
			cout << "���ܴ�AeroDyn_blade_a_force.txt�ļ�" << endl;
			system("pause");
		}
		in.close();
	}
	in.open("AeroDyn_blade_force.txt", ios::app);
	if (!in.is_open())
	{
		cout << "���ܴ�AeroDyn_blade_force.txt�ļ�" << endl;
		system("pause");
	}
	in << time << "    ";
	for (int j = 0; j <= 18;j++)
	{
		for (int i = 0;i < 3;i++)
		{
			in << blade_a[i][j] << "		";
		}
	}
	for (int j = 0; j <= 18;j++)
	{
		for (int i = 0;i < 3;i++)
		{
			in << blade_b[i][j] << "		";
		}
	}
	for (int j = 0; j <= 18;j++)
	{
		for (int i = 0;i < 3;i++)
		{
			in << blade_c[i][j] << "		";
		}
	}
	in << endl;
	
	return 0;

}

writeOut_Blade_force::writeOut_Blade_force(double time, int num_blade, double blade1_x[19], double blade1_y[19], double blade1_z[19])
{
	this->time = time;
	this->numblade = num_blade;
	for (int i = 0;i < 19;i++)
	{
		this->blade1_x[i] = blade1_x[i];
		this->blade1_y[i] = blade1_y[i];
		this->blade1_z[i] = blade1_z[i];
	}
}
int writeOut_Blade_force::write()
{
	fstream in;
	if (time == 0 && numblade == 1)
	{
		in.open("blade_force1.txt", ios::out);
		if (!in.is_open())
		{
			cout << "���ܴ�blade_force.txt�ļ�" << endl;
			system("pause");
		}
		in.close();
	}
	else if(time == 0 && numblade == 2)
	{
		in.open("blade_force2.txt", ios::out);
		if (!in.is_open())
		{
			cout << "���ܴ�blade_force2.txt�ļ�" << endl;
			system("pause");
		}
		in.close();
	}
	else if (time == 0 && numblade == 3)
	{
		in.open("blade_force3.txt", ios::out);
		if (!in.is_open())
		{
			cout << "���ܴ�blade_force3.txt�ļ�" << endl;
			system("pause");
		}
		in.close();
	}

	if (numblade == 1)
	{
		in.open("blade_force1.txt", ios::app);
		if (!in.is_open())
		{
			cout << "���ܴ�blade_force1.txt�ļ�" << endl;
			system("pause");
		}
	}
	else if (numblade == 2)
	{
		in.open("blade_force2.txt", ios::app);
		if (!in.is_open())
		{
			cout << "���ܴ�blade_force2.txt�ļ�" << endl;
			system("pause");
		}
	}
	else if (numblade == 3)
	{
		in.open("blade_force3.txt", ios::app);
		if (!in.is_open())
		{
			cout << "���ܴ�blade_force3.txt�ļ�" << endl;
			system("pause");
		}
	}
	in << time<<"    ";
	for (int i = 0; i <= 18;i++)
	{
		in << blade1_x[i] << "		"<<blade1_y[i] << "		" << blade1_z[i]<<"		";
	}
	in << endl;

	return 0;
}



morison_velocity::morison_velocity(double time, double motion[3], double WaveTmax,double WaveDT, double waterDepth, double Hs, double Tp, double U1[15000], double U2[15000],double incide_angle)
{
	/// <summary>
	/// ���㲨�������ˮ�ʵ��ٶ�
	/// </summary>
	/// <param name="time">��ǰʱ��</param>
	/// <param name="motion">��ǰ�ڵ�λ��</param>
	/// <param name="WaveTmax">���������</param>
	/// <param name="WaveDT"></param>
	/// <param name="waterDepth">ˮ��</param>
	/// <param name="Hs">���岨��</param>
	/// <param name="Tp">�׷�����</param>
	/// <param name="U1"></param>
	/// <param name="U2"></param>
	/// <param name="incide_angle"></param>
	this->incide_angle = incide_angle;
	for (int i = 0;i < 3;i++)
	{
		this->motion[i] = motion[i];
	}
	
	distance = motion[0] * cos(incide_angle * pi / 180) + motion[1] * sin(incide_angle * pi / 180);
	
	//N = WaveTmax / WaveDT;
	delta_omega = 2* pi / WaveTmax;//��Ƶ�ʵĲ�����ʽ�����Դ�FAST��Hydrodyn�������ļ����ҵ���Ӧ����
	N = 3 * 2 * pi / Tp / delta_omega;//�������3��˼��Ҫȡ3�����׷����ڵ�Ƶ��������м���
	wavespectral2 wave(Hs, Tp, WaveDT);
	waveNumber waveN(waterDepth);//�������������
	for (int i = 0; i < N ; i++)
	{
		omega = (i + 1) * 2 * pi / WaveTmax;
		S = wave.Spectral(omega);
		//cout << "S=" << S << endl;
		k_w = waveN.newton(omega);
		//cout << "k_w = " << k_w << endl;
		AW = sqrt(2 * S * delta_omega);//��һ���Ǹ��ݴ����뺣�󹤳̻����غɵ�18ҳ������д��
		//cout << "AW=" << AW << endl;
		f1 = sqrt(-2 * log(U1[i]));//BOX-Muller���˹�������е�һ��
		//cout << "f1=" << f1 << endl;
		f2 = cos(omega * time + 2 * pi * U2[i] - k_w * distance);//���˵�cos����е�2*pi*U2[i]�ǰ������еĸ�����
		//cout << "f2 = " << f2 << endl;
		f3 = omega * cosh(k_w * (motion[2] + waterDepth)) / sinh(k_w * waterDepth);
		if (isinf(f3))
		{
			f3 = 0;
		}
		f4 = -sin(omega * time + 2 * pi * U2[i] - k_w * distance);
		f5 = omega * sinh(k_w * (motion[2] + waterDepth)) / sinh(k_w * waterDepth);
		if (isinf(f5))
		{
			f5 = 0;
		}
		temp1 = f1 * AW * f2 * f3;
		if (abs(temp1) < 1.0E-6 || isinf(temp1) ||isnan(temp1))//||abs(temp1) >100 || motion[2] > 0
		{
			temp1 = 0;
		}
		//zeta = zeta + f1 * AW * f2 * f3;
		temp2 = f1 * AW * f4 * f5;
		if (abs(temp2) < 1.0E-6 || isinf(temp2) || isnan(temp2))//abs(temp2) >100||motion[2]>0
		{
			temp2 = 0;
		}
		zeta = zeta + temp1;
		zeta1 = zeta1 + temp2;
		
		if (isnan(zeta) || isnan(zeta1) || isinf(zeta) || isinf(zeta1))
		{
			cout << isnan(zeta) << "    " << isnan(zeta1) << "    " << isinf(zeta) << "    " << isinf(zeta1) << endl;
			cout << "zeta or zeta1 = nan!" << endl;
			cout << zeta << endl;
			cout << zeta1 << endl;
			cout << temp1 << endl;
			cout << temp2 << endl;
			system("pause");

			zeta = 0;
			zeta1 = 0;
		}

	}
}

morison_velocity::~morison_velocity()
{
	//cout << "�ͷ�morison_velocity" << endl;
}

double* morison_velocity::Velocity()
{
	if (motion[2] > 0)
	{
		velocity[0] = 0;
		velocity[1] = 0;
		velocity[2] = 0;
	}
	else
	{
		velocity[0] = zeta * cos(incide_angle * pi / 180);
		velocity[1] = zeta * sin(incide_angle * pi / 180);
		velocity[2] = zeta1;
	}
	
	return velocity;
}

int morison_velocity::checkValue(double value)
{
	if (isnan(value)|| isinf(value)||value>3)
	{
		cout << "the number is nan" << endl;
		return -1;
	}
	return 0;
}

void morison_velocity::morison_velocity2(double time, double motion[3], double WaveTmax, double WaveDT, double waterDepth, double Hs, double Tp, double U1[15000], double U2[15000], double incide_angle)
{
	zeta = 0;
	zeta1 = 0;
	this->incide_angle = incide_angle;
	distance = motion[0] * cos(incide_angle * pi / 180) + motion[1] * sin(incide_angle * pi / 180);
	N = WaveTmax / WaveDT;
	delta_omega = pi / WaveTmax;//��Ƶ�ʵĲ�����ʽ�����Դ�FAST��Hydrodyn�������ļ����ҵ���Ӧ����
	wavespectral2 wave(Hs, Tp, WaveDT);
	waveNumber waveN(waterDepth);//�������������

	for (int i = 0; i < 3;i++)
	{
		cout << "motion[" << i << "]=" << motion[i] << endl;
	}
	for (int i = 0; i < N / 2; i++)
	{
		omega = (i + 1) * 2 * pi / WaveTmax;
		S = wave.Spectral(omega);
		cout << "S=" << S << endl;
		k_w = waveN.newton(omega);
		cout << "k_w = " << k_w << endl;
		AW = sqrt(2 * S * delta_omega);
		cout << "AW=" << AW << endl;
		f1 = sqrt(-2 * log(U1[i]));
		cout << "f1=" << f1 << endl;
		f2 = cos(omega * time + 2 * pi * U2[i] - k_w * distance);//���˵�cos��
		cout << "f2 = " << f2 << endl;
		f3 = omega * cosh(k_w * (motion[2] + waterDepth)) / sinh(k_w * waterDepth);
		cout << "f3 = " << f3 << endl;
		cout << "omega=" << omega << "    cosh(k_w * (motion[2] + waterDepth))=" << cosh(k_w * (motion[2] + waterDepth)) << "    sinh(k_w * waterDepth)=" << sinh(k_w * waterDepth) << endl;
		for (int i = 0; i < 3;i++)
		{
			cout << "motion[" << i << "]=" << motion[i] << endl;
		}
		if (isnan(f3))
		{
			f3 = 0;
		}
		f4 = cos(omega * time + 2 * pi * U2[i] - k_w * distance);
		f5 = omega * sinh(k_w * (motion[2] + waterDepth)) / sinh(k_w * waterDepth);
		cout << "f4 = " << f4 << endl;
		cout << "f5 = " << f5 << endl;
		if (isnan(f5))
		{
			f5 = 0;
		}
		temp1 = f1 * AW * f2 * f3;
		cout << "temp1 = " << temp1 << endl;
		if (abs(temp1) < 1.0E-6 || isnan(temp1)||motion[2]>0)
		{
			temp1 = 0;
		}
		//zeta = zeta + f1 * AW * f2 * f3;
		
		temp2 = f1 * AW * f4 * f5;
		cout << "temp2 = " << temp2 << endl;
		if (abs(temp2) < 1.0E-6 || isnan(temp2)||motion[2]>0)
		{
			temp2 = 0;
		}
		zeta = zeta + temp1;
		cout << "zeta = " << zeta << endl;
		zeta1 = zeta1 + temp2;
		cout << "zeta1 = " << zeta1 << endl;
		if (isnan(zeta) || isnan(zeta1)|| isinf(zeta) || isinf(zeta1))
		{
			zeta = 0;
			zeta1 = 0;
		}

	}
}

