#include<iostream>
#include<cmath>
#include<string>
#include<Windows.h>
#include<stdlib.h>
#include<stdio.h>
#include"pch.h"
#include"SharedDll.h"
#include"wavespectral_dll.h"

using namespace std;

//#pragma comment (lib,"../../FloaterZhu/TimeDomainSolution/timeDomainSolution_by_fortran/timeDomainSolution_by_fortran/x64/Release/timeDomainSolution_by_fortran.lib")
//extern "C" _declspec(dllimport) int RETARDA(double* val);

#pragma data_seg("shared")

double theta = NULL;//�������Թ����ڴ�ı������и�ֵ������ñ������ᱻ���빲���ڴ�
double bla[3][19] = { 0 };//ҶƬa����
double blb[3][19] = { 0 };//ҶƬb����
double blc[3][19] = { 0 };//ҶƬc����
double force_tower[3][11] = { 0 };//�����ϵ���
double mo = NULL;
double velocity[10][2500] = { 0 };
double retard[10][150] = { 0 };//���ù������
int flag_Retard = 0;
double simulation_time = NULL;//ģ��ʱ��
double I[22200][7] = { 0 };
int line = NULL;
double Spec[100][2] = { 0 };
double Angle[37] = { 0 };
int Angle_n = NULL;
int F_n = NULL;
double rand_value1[200] = { 0 };
double rand_value2[200] = { 0 };
double t_global=-1.0;
double t_retard = -1.0;
double towerroot[3] = { 0 };//��Ͳ����ԭ���λ��
double towerrootAngle[3] = { 0 };//��Ͳ����ԭ����ƶ��Ƕ�
double towerrootVel[3] = { 0 };//��Ͳ����ԭ����ٶ�
double towerrootAngleVel[3] = { 0 };//��Ͳ����ԭ��Ľ��ٶ�
double Hs = NULL;
double Tp = NULL;
int wave_mod = NULL;
double Tp_1 = NULL;
double Tp_2 = NULL;
double percent_Tp = NULL;
double ragular_wave_phase = NULL;
double incide_angle = NULL;
double time2 = -1.0;
double time3 = -1.0;
double time4 = -1.0;
double WaveTmax = NULL;
double WaveDT = NULL;
double U1[15000] = { 0 };
double U2[15000] = { 0 };
double freq[200] = { 0 };
double waterDepth = NULL;//ˮ��
double watervelocity[3][300] = { 0 };//��¼ˮ�ʵ��ٶ�
int MorisonN = 0;
int MorisonI = 0;
int TF_watervelocity = 1;



#pragma data_seg()
#pragma comment(linker,"/SECTION:shared,RWS")

/*
����������˼·Ϊ������modelica����dll�е�WRITETHETA_GETFORCE��������modelica�е�ҶƬת��д��theta�У�
modelica��ʱ���ڵȴ�״̬����ʱAerodyn��ȡtheta�ǣ����������Ժ������뵽dll�������������У�
��ʱmodelica��ȡ�������ҶƬ���˶�������µ�theta�ǽ�����һ��ѭ��  
*/

extern "C" _declspec(dllexport) int WF_GT1(double input_theta,double t,double *towerRootMotion ,double *towerRootAngle, double *towerRootVel, double *towerRootAngleVel, double *fa, double* fb, double* fc, double* ft, double *m)
{
    cout << "============================================================" << endl;
    theta = input_theta;
    for (int i = 0;i < 3;i++)
    {
        towerroot[i] = towerRootMotion[i];
        towerrootAngle[i] = towerRootAngle[i];
        towerrootVel[i] = towerRootVel[i];
        towerrootAngleVel[i] = towerRootAngleVel[i];
    }
    
    cout << "Finish Write Theta:" << theta << endl;
    cout << "Modelica Time Step:" << t << endl;
    if (t > t_global)  // 2022.8.18�޸ģ�ԭ����if (t != t_global)����ҪĿ����Ϊ�˽��Modelica���㲽�س������⡣
    {
        HANDLE event1 = OpenEvent(EVENT_ALL_ACCESS, FALSE, L"д��Ƕ�");
        if (event1 == NULL)
        {
            cout << "Failed open event1" << endl;
            GetLastError();
            getchar();
            return 0;
        }
        else
        {
            cout << "OpenEvent 'event1' Success!" << endl;
        }

        SetEvent(event1);//�����źţ�����FAST����Ѿ��������
        cout << "Send signal to tell FAST theta have been writen..." << endl;

        HANDLE event2 = OpenEvent(EVENT_ALL_ACCESS, FALSE, L"д������");
        WaitForSingleObject(event2, INFINITE);//�ȴ�FAST����������
        if (event2 == NULL)
        {
            cout << "Failed Open event2" << endl;
            GetLastError();
            getchar();
            return 0;
        }
        else
        {
            cout << "OpenEvent 'event2' Success!" << endl;
        }
        int p = 0;
        int q = 0;
        for (int i = 0;i < 3; i++)
        {
            for (int j = 0;j < 19; j++)
            {

                fa[p] = bla[i][j];
                fb[p] = blb[i][j];
                fc[p] = blc[i][j];
                p++;
            }
            for (int j = 0; j < 11; j++)
            {
                ft[q] = force_tower[i][j];
                q++;
            }
        }

        *m = mo;
        
        cout << "*m=" << *m << endl;
        cout << "read the force success!" << endl;

        CloseHandle(event1);
        CloseHandle(event2);
        t_global = t;
    }
    else if (t <= t_global)
    {
        int p = 0;
        int q = 0;
        for (int i = 0;i < 3; i++)
        {
            for (int j = 0;j < 19; j++)
            {

                fa[p] = bla[i][j];
                fb[p] = blb[i][j];
                fc[p] = blc[i][j];
                p++;
            }
            for (int j = 0; j < 11; j++)
            {
                ft[q] = force_tower[i][j];
                q++;
            }
        }
        *m = mo;

        cout << " The time is same, use the last value to get the force." << endl;
    }
    
    return 0;
}

extern "C" _declspec(dllexport) int GETTHETA(double* theta2, double* towerRootMotion, double* towerRootAngle, double* towerRootVel, double* towerRootAngleVel)
{
    HANDLE event3 = OpenEvent(EVENT_ALL_ACCESS, FALSE, L"д��Ƕ�");
    WaitForSingleObject(event3, INFINITE);
    cout << "���յ���д��Ƕȡ��ź�" << endl;

    *theta2 = theta;
    for (int i = 0;i < 3;i++)
    {
        towerRootMotion[i] = towerroot[i];
        towerRootAngle[i] = towerrootAngle[i];
        towerRootVel[i] = towerrootVel[i];
        towerRootAngleVel[i] = towerrootAngleVel[i];
    }
    
    cout << "����Ƕȳɹ���" << *theta2 << endl;
    CloseHandle(event3);

    return 0;

}

extern "C" _declspec(dllexport) int WRFO1(double *time, double Fa[19][3],double Fb[19][3], double Fc[19][3], double FT[11][3], double *moment)
{
    for (int i = 0;i < 3; i++) 
    {
        for (int j = 0;j < 19; j++) 
        {
            bla[i][j] = Fa[j][i];
            blb[i][j] = Fb[j][i];
            blc[i][j] = Fc[j][i];

            //cout << "bla[" << i << "][" << j << "]=" << bla[i][j] << "          Fa[" << j << "][" << i << "]=" << Fa[j][i] << endl;
            //cout << "blb[" << i << "][" << j << "]=" << blb[i][j] << "          Fb[" << j << "][" << i << "]=" << Fb[j][i] << endl;
            //cout << "blc[" << i << "][" << j << "]=" << blc[i][j] << "          Fc[" << j << "][" << i << "]=" << Fa[j][i] << endl;

        }
        for (int j = 0;j < 11;j++)
        {
            force_tower[i][j] = FT[j][i];
            //cout << "force_tower[" << i << "][" << j << "]=" << force_tower[i][j] << endl;;
        }

    }
    double time1 = *time;
    writeOut_AeroDyn_force writeAeroDynForce(time1, bla, blb, blc);
    writeAeroDynForce.write();
    mo = *moment;
    cout << "mo=" << mo << endl;
    cout << "moment_dll=" << *moment << endl;
    cout << "�Ѿ�д������" << endl;
    HANDLE event4 = OpenEvent(EVENT_ALL_ACCESS, FALSE, L"д������");
    SetEvent(event4);
    CloseHandle(event4);

    return 0;
    
}

extern "C" _declspec(dllexport) int Start()
{
    cout << "����dll��,�س�����������" << endl;
    system("pause");

    //getchar();
    cout << "�������У�" << endl;
    return 0;
}

extern "C" _declspec(dllexport) int retard_cal(double acce,double t, double *val_retard)
{
    
    cout << "����retard_cal�ɹ�" << endl;
    typedef void(_stdcall* RETARDA)(double val[1500]);
    HINSTANCE hLibrary = LoadLibrary(L"timeDomainSolution_by_fortran.dll");

    if (hLibrary == NULL)
    {
        cout << "can't find the dll file" << endl;
        return -1;
    }

    RETARDA retarda = (RETARDA)GetProcAddress(hLibrary, "RETARDA");
    if (retarda == NULL)
    {
        cout << "can't find the function file." << endl;
        return -2;
    }
    double val[1500] = {0};
    int p = 0;
    for (int i = 0;i < 10; i++)
    {
        for (int j = 0;j < 150; j++)
        {
            retard[i][j] = val[p];
            p++;

        }
    }

    retarda(val);
    p = 0;
    for (int i = 0;i < 10; i++)
    {
        for (int j = 0;j < 150; j++)
        {
            retard[i][j] = val[p];

            //cout << "retard[" << i << "][" << j << "]=" << retard[i][j] << endl;

            p++;

        }
    }
    cout << "ִ��RETARDA�ɹ�"<<endl;
    //cout <<" val="<< val[1] <<"   val=" << val[2] << endl;
    
    
    
    //Start();
    FreeLibrary(hLibrary);
    double simulation_t;
    //����ģ��ʱ��
    cout << "����ģ��ʱ�䣺" << endl;
    cin >> simulation_t;
    simulation_time = simulation_t;


    return 0;
    
}

extern "C" _declspec(dllexport) int retard_int_cal(double time,double tp,double acce[6], double* val)
{
    cout << "start using retard integrator" << endl;
    cout << "time=" << time << "    tp=" << tp << "    acce[0]=" << acce[0] <<"  acce[1]="<< acce[1] << "  acce[2]=" << acce[2] << "  acce[3]=" << acce[3] << "  acce[4]=" << acce[4] << "  acce[5]=" << acce[5] << endl;

    /////////////////////////////////////////////////////////
    int NS = (int)(simulation_time / tp);//�ܹ��ж��ټ��㲽
    int N = (int)(time / tp);

    cout << "N=  " << N << endl;
   



    ////////////////////////////////////////////////////////
    double* velocity_t = new double[10];
    //double velocity_t[10] = { 0 };//��������ٶȽ��з���
    velocity_t[0] = acce[2];
    velocity_t[1] = acce[0];
    velocity_t[2] = acce[4];
    velocity_t[3] = acce[0];
    velocity_t[4] = acce[4];
    velocity_t[5] = acce[5];
    velocity_t[6] = acce[1];
    velocity_t[7] = acce[3];
    velocity_t[8] = acce[1];
    velocity_t[9] = acce[3];


    /****************************************������̬����
    double** velo = new double*[10];//
    for (int i = 0;i < 10;i++)
    {
        velo[i] = new double[2500];
    }

 
    /***************************************����̬���鸳ֵ*/
    double velo[10][2500] = { 0 };
    if (N >= 2500)  //�����ֲ�������2500ʱ���ٶ����齫�������������ǰ�ƣ���������ĸ���
    {
        //cout << "�����Ѵ���2500��!" << endl;
        if (time != t_retard)
        {
            for (int i = 0; i < 2499;i++)
            {
                for (int j = 0;j < 10;j++)
                {
                    velocity[j][i] = velocity[j][i + 1];
                }
            }
            t_retard = time;
        }
        for (int i = 0;i < 10;i++)
        {
            velocity[i][2499] = velocity_t[i];
        }
        //cout << "Velocity(2499)=" << velocity[1][2498] << "    Velocity(2500)=" << velocity[1][2499] << endl;
        
    }
    else
    {
        
        for (int i = 0;i < 10;i++)
        {
            velocity[i][N] = velocity_t[i];
        }
        //cout << "Velocity(N-1)=" << velocity[1][N] << "    Velocity(N)=" << velocity[1][2499] << endl;
    }
    for (int i = 0;i < 10;i++)//velo��һ���м䴫�ݱ�������ΪC++�ڵ���fortranʱ�����ù����ڴ�����ݽ��н�����������������һ���м䴫�ݱ���velo��
    {
        for (int j = 0;j < 2500;j++)
        {
            velo[i][j] = velocity[i][j];

        }
    }
    /*----------------------------------------------------*/
    /*-----------------���ö�̬���ӿ⣨��ʾ���ӣ�---------------------*/
    typedef void(_stdcall* RETARDINT)(double *T, double *TP, double retarda[10][150], double velo[10][2500], double val[10]);
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
    RETARDINT retardint = (RETARDINT)GetProcAddress(hLibrary, "RETARDINT");
    if (retardint == NULL)
    {
        cout << "can't find the function file." << endl;
        return -2;
    }
    else
    {
        cout << "The funtion file has been found!" << endl;
    }
    /*----------------������ɣ���ʼ�Ը��������и�ֵ----------------*/
    double *T = &time;
    double *TP = &tp;
    double acc[6] = { 0 };
    for (int i = 0;i < 6; i++)
    {
         acc[i] = acce[i];
       
    }

    double valu[10] = { 0 };

    /*
    double** retarda = new double*[10];
    for (int i = 0;i < 10;i++)
        retarda[i] = new double[150];
    */
    double retarda[10][150] = { 0 };

    for (int i = 0;i < 10; i++)
    {
        for (int j = 0; j < 150; j++)
        {
            retarda[i][j] = retard[i][j];
        }
    }

    retardint(T, TP, retarda, velo, valu);


    for (int i = 0;i < 10;i++)
    {
        val[i] = valu[i];
        //cout << "Integration result  =  " << val[i] << "  ";
    }
    cout <<"-------------------------------------------------"<< endl;

    FreeLibrary(hLibrary);
    delete[] velocity_t;
/*
    for (int i = 0;i < 10;i++)//ɾ����̬����
    {
        delete[] velo[i];
        delete[] retarda[i];
    }
    delete[] velo;
    delete[] retarda;
    delete[] velocity_t;
*/

    return 0;
}

extern "C" _declspec(dllexport) int storeI()//��ȡˮ����ϵ����������������
{
    double pi = 3.1415926;
    cout << "����storeI�ɹ�" << endl;
    //Hs = 3.2;
    //Tp = 7.9;
    //wave_mod = 1;
    cout << "������Hs:" << endl;
    cin >> Hs;
    cout << "������Tp:" << endl;
    cin >> Tp;
    cout << "�����벨����ʽ:(1--Jonswap;2--ragular wave;3--still water)" << endl;
    cin >> wave_mod;
    cout << "������ˮ��waterDepth:" << endl;
    cin >> waterDepth;
    cout << "������Morison��Ԫ����MorisonN:" << endl;
    cin >> MorisonN;
    cout << "�Ƿ����ˮ�ʵ��ٶȣ�(1--����;0--�����㣩" << endl;
    cin >> TF_watervelocity;
    //------------------�Բ�����ʽ���ж��壨1��ʾ����Jonswap����������ˣ�2��ʾ���ɹ��򲨣�---------------------//
    if (wave_mod == 1)
    {
        cout << "�����벨�����ʱ��WaveTmax:" << endl;
        cin >> WaveTmax;
        cout << "�����벨��ģ���ʱ�䲽��WaveDT:" << endl;
        cin >> WaveDT;
        cout << "Jonswap����������" << endl;
        ragular_wave_phase = 0;//Jonswap�����в���Ҫ�Թ�����λ���ж���
    }
    else if (wave_mod == 2)
    {
        cout << "��������򲨵���λ" << endl;
        cin >> ragular_wave_phase;
    }
    //-----------------------------------------------------------//

    readmyfile readfile;//��ȡ.3�ļ���readfile�а�����
    cout << "����readmyfile�ɹ���" << endl;
    for (int i = 0;i < readfile.line;i++)
    {
        for (int j = 0;j < 7;j++)
        {
            I[i][j] = readfile.data[i][j];
            //cout << I[i][j] << "    ";
        }
        //cout << endl;
    }
    for (int i = 0;i < 7;i++)
    {
        cout << I[readfile.line - 1][i];
    }
    cout << endl;
    cout << "�����ڴ�������Ѿ�������ɣ�" << endl;
    
    for (int i = 0;i < readfile.Angle_n;i++)
    {
        Angle[i] = readfile.Angle[i];//��ȡˮ�����ļ��е������
    }
    Angle_n = readfile.Angle_n;//��ȡˮ�����ļ�������ǵ�����
    line = readfile.line;
    F_n = readfile.F_n;

    if (F_n > 200)//�����ȡ��Ƶ�ʸ����Ƿ񳬹��趨������
    {
        cout << "ˮ�����ļ��в���Ƶ�ʸ�������200�������ڴ汣���Ƶ�ʸ�����С��200��" << endl;
        system("pause");
    }
    for (int i = 0;i < F_n;i++)//Ϊˮ�����ļ���Ƶ�ʽ��и�ֵ
    {
        freq[i] = readfile.freq[i];
    }
    cout << "Angle_n = " << Angle_n << "    line = " << line << endl;

    readrandom random_value;//��ȡ���ֵ
    if (random_value.line > 15000)
    {
        cout << "������������������15000�����޷�������ȡ������������ݣ�" << endl;
        cout << "random_value.txt�����������Ϊ��" << random_value.line << endl;

        system("pause");
    }

    for (int i = 0; i < random_value.line;i++)
    {
        U1[i] = random_value.data[i][0];
        U2[i] = random_value.data[i][1];

    }

   
    if (wave_mod == 1)//����Jonswap�����ɲ��˲���
    {
        cout << "�������������ǣ�" << endl;
        cin >> incide_angle;
        
        /*-------------------------------����Spar.3�ļ��ṩ���������ɲ���----------------------------------------------*/

        wavespectral wave(Hs, Tp, readfile.F_n, readfile.freq);//�����˳��:Hs ���岨�ߣ� Tp�׷����ڣ�F_nƵ�ʸ�����freqƵ��ֵ
        cout << "wave�����Ѿ�����" << endl;
        
        double* s = wave.Spectral();

        cout << "����Spectral�������" << endl;
        for (int i = 0;i < readfile.F_n;i++)
        {
            Spec[i][0] = 2*pi*readfile.freq[i];//Spec�����һ��Ϊ����Ƶ�ʣ�{С��������������ǿ��ʹ�õ��ǽ�Ƶ��}���ڶ���Ϊ��Ƶ�ʶ�Ӧ�Ĳ���ֵ
            Spec[i][1] = s[i];
            rand_value1[i] = wave.ec1[i];//�����������������λ���и�ֵ
            rand_value2[i] = wave.ec2[i];//�ڶ�����������и�ֵ

            //cout << Spec[i][0] << "    " << Spec[i][1] << "----------" << rand_value[i] << endl;
        }
        cout << "Spectral�������" << endl;
    }
    else if (wave_mod == 2)//���ɹ��򲨵���ز���
    {
        //double d_o[100];
        double* d_o = new double[F_n];//���ڹ����м�����ӽ���������Ƶ�ʵ�����
        cout << "�������������ǣ�" << endl;
        cin >> incide_angle;
        cout << "��������" << endl;
        for (int i = 0;i < readfile.F_n;i++)
        {
            d_o[i] = readfile.freq[i] - 1 / Tp;
            if (d_o[i] == 0)
            {
                Tp_1 = 1/readfile.freq[i];
                Tp_2 = Tp_1;
                percent_Tp = 1;
                break;
            }
            if (i > 0)
            {
                if (d_o[i] * d_o[i - 1] < 0)
                {
                    Tp_1 = 1 / readfile.freq[i - 1];
                    Tp_2 = 1 / readfile.freq[i];
                    percent_Tp = 1-(readfile.freq[i-1]-1/Tp)/(readfile.freq[i-1]-readfile.freq[i]);
                    //percent_Tp = 1;
                }
            }
        }
        delete[] d_o;
    }
    cout << "StoreI�������" << endl;
    return 0;
}

extern "C" _declspec(dllexport) int excitforce(double forwards_angle, double time, double motion[6], double *Force)//forwards_angleΪ�������ǣ�timeΪʱ�䲽��FΪ6�����ɶȴ�������
{
    cout << "step::1------use excitforce*******************************************************" << endl;
    cout << "forwards_angle = " << forwards_angle << endl;
    cout << "time =" << time << endl;
    double Angle_n_1 = Angle_n;
    double *Angle_1 = new double[Angle_n_1];
    double pi = 3.1415926;
    for (int i = 0;i < Angle_n_1;i++)
    {
        Angle_1[i] = Angle[i];
       // cout << Angle_1[i] << endl;
    }
    double forwards_angle_degree;
    if (wave_mod == 1)
    {
        forwards_angle_degree = forwards_angle * 180 / pi + incide_angle;
    }
    else if (wave_mod == 2)
    {
        forwards_angle_degree = forwards_angle * 180 / pi + incide_angle;
    }
    Find FIND(Angle_1, Angle_n_1,forwards_angle_degree);
    FIND.find_data();
    int p = FIND.p;
    int q = FIND.q;
    double percent = FIND.percent;
    double Angle_p = Angle_1[p];
    double Angle_q = Angle_1[q];
    cout << "p = " << p << "    q = " << q << endl;
    cout << "Angle_p = " << Angle_p << "    Angle_q = " << Angle_q << endl;
    //readmyfile readfile;
    double F[6] = { 0 };//���˼��������
    double A1;
    double A2;
    double B1;
    double B2;
    double A;
    double B;
    double a;
    double sigma1;
    double sigma2;
    double sigmaI=0;
    double sigmaII=0;
    double SIGMA=0;
    double sigma;
    double delta_omega;//��Ƶ�ʲ���
    double zeta=0;//����˲ʱ�߶�
    double C0 = 0;
    double S = 0;
    double omega;
    double AW;//���˷�ֵ
    double f1;//���˷�ֵ�ֲ���
    double f2;//����cos��
    double HA;//���˼��������ݺ�����ֵ
    bool Check[4];//���ڼ���ѭ���Ƿ�����Ҫ��������Ҫ������ֹѭ��
    int k_Tp_1;//��¼Tp_1����ӦƵ�ʵļ�����
    int k_Tp_2;//��¼Tp_2����ӦƵ�ʵļ�����
    int k[4];
    double k_w; //wave number ����
    waveNumber waveN(waterDepth);//�������������
    double distance = motion[0] * cos(incide_angle * pi / 180) + motion[1] * sin(incide_angle * pi / 180);
    //double distance = sqrt(pow(motion[0],2)+pow(motion[1],2));

    /* zhu:���޸��ǽ�������е�PrintOut�ӳ�����
    if (time3 != time)
    {
        writeOut write(time, motion);//���˶������motion.txt�ļ���
        time3 = time;
    }
    */
    for (int i = 0;i < 4;i++)
    {
        Check[i] = FALSE;
    }
    fstream OutFile;
    if (time == 0)//����һ�δ�wave.txtʱ����ɾ��ԭ�е��ļ�����
    {
        OutFile.open("wave.txt", ios::out);
        OutFile.close();
    }

    int N;//ʱ����ɢ��������Ƶ����ɢ���������ʱ����ɢ�������һ�룬��N/2.
    
    //cout << "step 2::-----------finish define parameter" << endl;
    
    //for (int freedom = 0; freedom < 6;freedom++)//�Բ�ͬ���ɶȽ�������
    //{
        //cout << "step 3." << freedom << "!---------------------------------------------------!" << endl;
        /*������γ�����Ҫ��������ˮ�����ļ����ҵ����ʵ�ֵ*/
    A = 0;
    A1 = 0;
    A2 = 0;
    B1 = 0;
    B2 = 0;
    sigma1 = 0;
    sigma2 = 0;
    sigma = 0;
    if (F_n == 0)
    {
        cout << "ȫ�ֱ���F_nδ��ֵ" << endl;
        system("pause");
    }
        
    if (wave_mod == 1)
    {
        /*------------�������ȸ������Բ������۽�������---------------------*/
        N = WaveTmax / WaveDT;
        delta_omega = pi / WaveTmax;//��Ƶ�ʵĲ�����ʽ�����Դ�FAST��Hydrodyn�������ļ����ҵ���Ӧ����
        wavespectral2 wave(Hs, Tp, WaveDT);
        for (int i = 0; i < N / 2; i++)
        {
            omega = (i+1) * 2 * pi / WaveTmax;
            S = wave.Spectral(omega);
            AW = sqrt(2 * S * delta_omega);
            f1 = sqrt(-2 * log(U1[i]));
            f2 = cos(omega * time + 2 * pi * U2[i]);//���˵�cos��
            zeta = zeta + f1 * AW * f2;
        }
        if (time2 != time)//�������
        {
            //cout << "zeta = " << zeta << endl;

            OutFile.open("wave.txt", ios::app);//������������ļ�
            if (!OutFile.is_open())
            {
                cout << "can not open wave.txt file!" << endl;
                system("pause");
            }
            OutFile << time << "    " << zeta << endl;
            zeta = 0;
            OutFile.close();
            time2 = time;
        }
        /*--------------------�����������---------------------*/
        /*--------------------���˼���������-------------------*/
        double* d_o = new double[F_n];
        for (int freedom = 0; freedom < 6;freedom++)//�Բ�ͬ���ɶȽ�������
        {
            //cout << "" << endl;
            for (int i = 0;i < F_n;i++)
            {
                d_o[i] = 0;
            }
            for (int i = 0;i < N / 2;i++)
            {
                omega = (i + 1) * 2 * pi / WaveTmax;
                S = wave.Spectral(omega);
                AW = sqrt(2 * S * delta_omega);
                f1 = sqrt(-2 * log(U1[i]));

                //����k_w
                k_w = waveN.newton(omega);
                //cout << "k_w = " << k_w << endl;
                

                for (int j = 0;j < F_n;j++)
                {
                    d_o[j] = freq[j] - omega / 2 / pi;
                    //cout << "freq[" << j << "]=" << freq[j] << "    omega/2/pi = " << omega / 2 / pi << endl;
                    //cout << "d_o[" << j << "]=" << d_o[j] << endl;
                    if (freq[0] > omega / 2 / pi)//��Ƶ�ʹ�Сʱ������ˮ�������в�û��0Ƶ��ϵ������Ҫ���⴦��
                    {
                        k_Tp_1 = -1;
                        k_Tp_2 = 0;
                        percent_Tp = 1 - omega / 2 / pi / freq[0];
                        break;
                    }

                    if (d_o[j] == 0)
                    {
                        Tp_1 = 1 / freq[j];
                        Tp_2 = Tp_1;
                        percent_Tp = 1;
                        k_Tp_1 = j;
                        k_Tp_2 = j;
                        //cout << "k_Tp_1=" << k_Tp_1 << endl;
                        //cout << "k_Tp_2=" << k_Tp_2 << endl;
                        break;
                    }
                    if (j > 0)
                    {
                        if (d_o[j] * d_o[j - 1] < 0)
                        {
                            Tp_1 = 1 / freq[j - 1];
                            Tp_2 = 1 / freq[j];
                            percent_Tp = 1 - (freq[j - 1] - 1 / Tp) / (freq[j - 1] - freq[j]);
                            k_Tp_1 = j - 1;
                            k_Tp_2 = j;
                            //cout << "k_Tp_1=" << k_Tp_1 << endl;
                            //cout << "k_Tp_2=" << k_Tp_2 << endl;
                            break;
                        }
                    }
                }
                A = 0;
                A1 = 0;
                A2 = 0;
                B1 = 0;
                B2 = 0;
                sigma1 = 0;
                sigma2 = 0;
                sigma = 0;
                if (F_n == 0)
                {
                    cout << "ȫ�ֱ���F_nδ��ֵ" << endl;
                    system("pause");
                }
                //Angle_p = 0;
                //Angle_q = 0;
                    
                    
                //������I�л�ȡA1��A2��B1��B2��ֵ
                //cout << "׼������kֵ����" << endl;
                if (k_Tp_1 == -1)
                {
                    A1 = 0;
                    sigma1 = 0;
                    A2 = 0;
                    sigma2 = 0;


                    k[2] = k_Tp_2 * Angle_n * 6 + p * 6 + freedom;
                    //cout << "k[2] = " << k[2] << endl;
                    B1 = I[k[2]][3];
                    sigmaI = I[k[2]][4];


                    k[3] = k_Tp_2 * Angle_n * 6 + q * 6 + freedom;
                    //cout << "k[3] = " << k[3] << endl;
                    B2 = I[k[3]][3];
                    sigmaII = I[k[3]][4];
                }
                else
                {
                    k[0] = k_Tp_1 * Angle_n * 6 + p * 6 + freedom;
                    //cout << "k[0] = " << k[0] << endl;
                    //cout << "line = " << line << endl;
                    A1 = I[k[0]][3];
                    sigma1 = I[k[0]][4];


                    k[1] = k_Tp_1 * Angle_n * 6 + q * 6 + freedom;
                    // cout << "k[1] = " << k[1] << endl;
                    A2 = I[k[1]][3];
                    sigma2 = I[k[1]][4];


                    k[2] = k_Tp_2 * Angle_n * 6 + p * 6 + freedom;
                    // cout << "k[2] = " << k[2] << endl;
                    B1 = I[k[2]][3];
                    sigmaI = I[k[2]][4];


                    k[3] = k_Tp_2 * Angle_n * 6 + q * 6 + freedom;
                    // cout << "k[3] = " << k[3] << endl;
                    B2 = I[k[3]][3];
                    sigmaII = I[k[3]][4];
                }
                /*

                for (int k = 0; k < line;k++)
                {
                    if (Angle_p == I[k][1] && freedom == (I[k][2] - 1) && Tp_1 == I[k][0])
                    {
                        A1 = I[k][3];
                        sigma1 = I[k][4];
                        Check[0] = TRUE;
                        //cout << "A1[" << freedom << "]=" << A1 << endl;
                        //cout << "sigma1[" << freedom << "]=" << sigma1 << endl;
                        //cout << "Tp_1=" << Tp_1 << endl;
                    }
                    if (Angle_q == I[k][1] && freedom == (I[k][2] - 1) && Tp_1 == I[k][0])
                    {
                        A2 = I[k][3];
                        sigma2 = I[k][4];
                        Check[1] = TRUE;
                        //cout << "A2[" << freedom << "]=" << A2 << endl;
                        //cout << "sigma2[" << freedom << "]=" << sigma2 << endl;
                        //cout << "Tp_1=" << Tp_1 << endl;
                    }
                    if (Angle_p == I[k][1] && freedom == (I[k][2] - 1) && Tp_2 == I[k][0])
                    {
                        B1 = I[k][3];
                        sigmaI = I[k][4];
                        Check[2] = TRUE;
                        //cout << "B1[" << freedom << "]=" << B1 << endl;
                        //cout << "sigmaI[" << freedom << "]=" << sigmaI << endl;
                        //cout << "Tp_2=" << Tp_2 << endl;
                    }
                    if (Angle_q == I[k][1] && freedom == (I[k][2] - 1) && Tp_2 == I[k][0])
                    {
                        B2 = I[k][3];
                        sigmaII = I[k][4];
                        Check[3] = TRUE;
                        //cout << "B2[" << freedom << "]=" << B2 << endl;
                        //cout << "sigmaII[" << freedom << "]=" << sigmaII << endl;
                        //cout << "Tp_2=" << Tp_2 << endl;
                    }
                    if (Check[0] == TRUE && Check[1] == TRUE && Check[2] == TRUE && Check[3] == TRUE)
                    {
                        for (int p = 0;p < 4;p++)
                        {
                            Check[p] = FALSE;
                        }
                        break;//���ĸ�ֵ���ҵ��Ժ����˳�����ѭ��
                    }
                }
                //�˲��ֱ���ע��������û�вο�FAST���ʱ�����ݡ�
                */
                //cout << "------------------------------------------------" << endl;
                A = (1 - percent) * A1 + percent * A2;
                sigma = (1 - percent) * sigma1 + percent * sigma2;//���ݽǶ�ȷ�����ݺ����ĸ�ֵ����λ��
                B = (1 - percent) * B1 + percent * B2;
                SIGMA = (1 - percent) * sigmaI + percent * sigmaII;
                //cout << "percent=" << percent << endl;
                //cout << "percent_Tp = " << percent_Tp << endl;
                if (isnan(A) == true || isnan(B) == true || isnan(percent) == true)
                {
                    cout << "A = " << A << "    B = " << B << "    percent=" << percent << endl;
                    system("pause");
                }

                HA = percent_Tp * A + (1 - percent_Tp) * B;
                if (isnan(percent_Tp) == true)
                {
                    cout << "percent_Tp=" << percent << endl;
                    system("pause");
                }
                sigma = percent_Tp * sigma + (1 - percent_Tp) * SIGMA;
                //cout << "A["<<freedom<<"]=" << A << endl;
                // cout << "sigma[" << freedom << "]=" << sigma << endl;
                //cout << "Hs[" << freedom << "]=" << Hs << endl;
                F[freedom] = F[freedom] + f1 * AW * HA * cos(omega * time + 2 * pi * U2[i] + sigma * pi / 180 - k_w*distance);
                if (isnan(F[freedom]) == true)
                {
                    cout << "fΪnan!" << endl;
                    cout << "f1 = " << f1 << "    AW=" << AW << "    HA=" << HA << endl;
                    cout << "percent_Tp = " << percent_Tp << "    A = " << A << "    B=" << B << "    percent=" << percent << endl;
                    system("pause");
                }
                // cout << F[freedom] << endl;
                //F[freedom] = F[freedom] + Hs / 2 * A * cos(2 * pi / Tp * time + ragular_wave_phase * pi / 180 + sigma * pi / 180);

            }
        }
        delete[] d_o;
    }

            
            /*-------------ԭ�����µĲ����������---------------------
            for (int fre = 0; fre < F_n;fre++)//��ÿһ��Ƶ��Ѱ�Һ��ʵ�����
            {

                for (int i = 0; i < line;i++)
                {
                    if (Angle_p == I[i][1] && freedom == I[i][2] - 1 && Spec[fre][0] == 2 * pi / I[i][0])
                    {
                        A1 = I[i][3];
                        sigma1 = I[i][4];
                    }
                    if (Angle_q == I[i][1] && freedom == I[i][2] - 1 && Spec[fre][0] == 2 * pi / I[i][0])
                    {
                        A2 = I[i][3];
                        sigma2 = I[i][4];
                    }

                }
                //cout << "sigma1 = " << sigma1 << "    sigma2 =" << sigma2 << endl;
                //cout << "percent = " << percent << endl;
                A = percent * A1 + (1 - percent) * A2;
                sigma = percent * sigma1 + (1 - percent) * sigma2;//���ݽǶ�ȷ�����ݺ����ĸ�ֵ����λ��
                bool isnumber = isfinite(sigma);
                if (!isnumber)
                {
                    cout << " the sigma is not a number!" << endl;
                    system("pause");
                    //cout << "sigma = " << sigma << endl;
                }

                if (fre == 0)
                {
                    delta_omega = (Spec[1][0] - Spec[0][0]);

                }
                else
                {
                    delta_omega = (Spec[fre][0] - Spec[fre - 1][0]);
                }
                a = sqrt(2 * Spec[fre][1] * delta_omega);

                if (freedom == 0)
                {
                   // cout << "Spec=" << Spec[fre][1] << endl;
                    
                    zeta = zeta + 1.0 / 2.0 / pi * (sqrt(-2.0 * log(rand_value1[fre])) * sqrt(pi * Spec[fre][1]) * cos(Spec[fre][0] * time + rand_value2[fre] * 2.0 * pi));
                    //cout << "zeta = " << zeta << endl;
                }

                F[freedom] = F[freedom] + a * A * cos(Spec[fre][0] * 2 * pi * time + sigma * pi / 180 + rand_value2[fre] * 2 * pi);
                
                isnumber = isfinite(F[freedom]);
                if (!isnumber)
                {
                    cout << " the F[" << freedom << "] is not a number!" << endl;
                    system("pause");
                    //cout << "F["<<freedom<<"] = " << F[freedom] << endl;
                }
            }
            if (freedom == 0 && time2 !=time)
            {                
                //cout << "zeta = " << zeta << endl;
                
                OutFile.open("wave.txt", ios::app);//������������ļ�
                if (!OutFile.is_open())
                {
                    cout << "can not open wave.txt file!" << endl;
                    system("pause");
                }
                OutFile << time << "    " << zeta << endl;
                zeta = 0;
                OutFile.close();
                time2 = time;
            }
            */

        
    else if (wave_mod == 2)
    {
        for (int freedom = 0; freedom < 6;freedom++)//�Բ�ͬ���ɶȽ�������
        {
            //cout << "step 3." << freedom << "!---------------------------------------------------!" << endl;
            A = 0;
            A1 = 0;
            A2 = 0;
            B1 = 0;
            B2 = 0;
            sigma1 = 0;
            sigma2 = 0;
            sigma = 0;
            if (F_n == 0)
            {
                cout << "ȫ�ֱ���F_nδ��ֵ" << endl;
                system("pause");
            }
            //Angle_p = 0;
            //Angle_q = 0;
            for (int i = 0; i < line;i++)
            {
                if (Angle_p == I[i][1] && freedom == (I[i][2] - 1) && Tp_1 == I[i][0])
                {
                    A1 = I[i][3];
                    sigma1 = I[i][4];
                    //cout << "A1[" << freedom << "]=" << A1 << endl;
                    //cout << "sigma1[" << freedom << "]=" << sigma1 << endl;
                    //cout << "Tp_1=" << Tp_1 << endl;
                }
                if (Angle_q == I[i][1] && freedom == (I[i][2] - 1) && Tp_1 == I[i][0])
                {
                    A2 = I[i][3];
                    sigma2 = I[i][4];
                    //cout << "A2[" << freedom << "]=" << A2 << endl;
                    //cout << "sigma2[" << freedom << "]=" << sigma2 << endl;
                    //cout << "Tp_1=" << Tp_1 << endl;
                }
                if (Angle_p == I[i][1] && freedom == (I[i][2] - 1) && Tp_2 == I[i][0])
                {
                    B1 = I[i][3];
                    sigmaI = I[i][4];
                    //cout << "B1[" << freedom << "]=" << B1 << endl;
                    //cout << "sigmaI[" << freedom << "]=" << sigmaI << endl;
                    //cout << "Tp_2=" << Tp_2 << endl;
                }
                if (Angle_q == I[i][1] && freedom == (I[i][2] - 1) && Tp_2 == I[i][0])
                {
                    B2 = I[i][3];
                    sigmaII = I[i][4];
                    //cout << "B2[" << freedom << "]=" << B2 << endl;
                    //cout << "sigmaII[" << freedom << "]=" << sigmaII << endl;
                    //cout << "Tp_2=" << Tp_2 << endl;
                }
            }
            cout << "------------------------------------------------" << endl;
            A = (1 - percent) * A1 + percent * A2;
            sigma = (1 - percent) * sigma1 + percent * sigma2;//���ݽǶ�ȷ�����ݺ����ĸ�ֵ����λ��
            B = (1 - percent) * B1 + percent * B2;
            SIGMA = (1 - percent) * sigmaI + percent * sigmaII;
            //cout << "percent=" << percent << endl;
            //cout << "percent_Tp = " << percent_Tp << endl;
            A = percent_Tp * A + (1 - percent_Tp) * B;
            sigma = percent_Tp * sigma + (1 - percent_Tp) * SIGMA;
            //cout << "A["<<freedom<<"]=" << A << endl;
            // cout << "sigma[" << freedom << "]=" << sigma << endl;
            //cout << "Hs[" << freedom << "]=" << Hs << endl;
            F[freedom] = Hs / 2 * A * cos(2 * pi / Tp * time + ragular_wave_phase * pi / 180 + sigma * pi / 180);
            //cout <<"F["<<freedom<<"]=  " <<F[freedom] << endl;
            bool isnumber = isfinite(F[freedom]);
            if (!isnumber)
            {
                cout << " the F[" << freedom << "] is not a number!" << endl;
                system("pause");
                //cout << "F["<<freedom<<"] = " << F[freedom] << endl;
            }
        }
    }
    else if (wave_mod == 3)
    {
        for (int i = 0; i < 6;i++)
        {
            F[i] = 0;
        }
    }
               
    
    for (int i = 0;i < 6;i++)
    {
        Force[i] = F[i];
        //cout << Force[i]<<"    ";
    }
    cout << endl;
    //cout << "step5" << endl;
    delete[] Angle_1;
   // cout << "step6" << endl;
    return 0;
}

extern "C" _declspec(dllexport) int fluitVelocity(double time, double motion[3], double* Velocity)
{
    /*
    int N;
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
    double zeta =0;
    double zeta1 = 0;
    double velocity[3];
    double distance;
    double k_w;
    double temp1;
    double temp2;
    //cout << "time = " << time << endl;
    /*
    for (int i=0; i < 3;i++)
    {
        cout << "motion[" << i << "]=" << motion[i] << endl;
    }
    

    
    distance = motion[0] * cos(incide_angle * pi / 180) + motion[1] * sin(incide_angle * pi / 180);


    N = WaveTmax / WaveDT;
    delta_omega = pi / WaveTmax;//��Ƶ�ʵĲ�����ʽ�����Դ�FAST��Hydrodyn�������ļ����ҵ���Ӧ����
    wavespectral2 wave(Hs, Tp, WaveDT);
    waveNumber waveN(waterDepth);//�������������
    for (int i = 0; i < N / 2; i++)
    {
        omega = (i + 1) * 2 * pi / WaveTmax;
        S = wave.Spectral(omega);
        //cout << "S=" << S << endl;
        k_w = waveN.newton(omega);
        //cout << "k_w = " << k_w << endl;
        AW = sqrt(2 * S * delta_omega);
        //cout << "AW=" << AW << endl;
        f1 = sqrt(-2 * log(U1[i]));
        //cout << "f1=" << f1 << endl;
        f2 = cos(omega * time + 2 * pi * U2[i]- k_w*distance);//���˵�cos��
        //cout << "f2 = " << f2 << endl;
        f3 = omega * cosh(k_w * (motion[2] + waterDepth)) / sinh(k_w * waterDepth);
        if (isnan(f3))
        {
            f3 = 0;
        }
        f4 = cos(omega * time + 2 * pi * U2[i] - k_w * distance);
        f5 = omega * sinh(k_w * (motion[2] + waterDepth)) / sinh(k_w * waterDepth);
        if (isnan(f5))
        {
            f5 = 0;
        }
        //cout << "motion+depth=" << motion[2] + waterDepth << endl;
        //cout << "waterdepth = " << waterDepth << endl;
        //cout << "cosh = " << cosh(k_w * (motion[2] + waterDepth)) << endl;
        //cout << "sinh = " << sinh(k_w * waterDepth) << endl;
        //cout << "f3 = " << f3 << endl;
        temp1 = f1 * AW * f2 * f3;
        if (abs(temp1) < 1.0E-6 || isnan(temp1))
        {
            temp1 = 0;
        }
        zeta = zeta + f1 * AW * f2 * f3;
        temp2 = f1 * AW * f4 * f5;
        if (abs(temp2) < 1.0E-6 || isnan(temp2))
        {
            temp2 = 0;
        }
        zeta = zeta + temp1;
        zeta1 = zeta1 + temp2;
        if (isnan(zeta) || isnan(zeta1))
        {
            zeta = 0;
            zeta1 = 0;
        }
        
    }
    
    velocity[0] = zeta * cos(incide_angle*pi/180);
    velocity[1] = zeta * sin(incide_angle*pi/180);
    velocity[2] = zeta1;
    
    for (int i = 0;i < 3;i++)
    {
        if (isnan(velocity[i]))
        {
            cout << "velocity[" << i << "]=" << velocity[i] << endl;
            //system("pause");
            velocity[i] = 0;
            cout << "velocity[" << i << "]=" << velocity[i] << endl;
        }
        //cout << "velocity[" << i << "]=" << velocity[i] << endl;
        Velocity[i] = velocity[i];

    }

    */
    double *velocity;
    if (TF_watervelocity == 1)
    {
        if (time != time4)
        {
        
            morison_velocity morison(time, motion, WaveTmax, WaveDT, waterDepth, Hs, Tp, U1, U2, incide_angle);


            velocity = morison.Velocity();


            for (int i = 0;i < 3;i++)
            {
                Velocity[i] = velocity[i];
                //cout << "Velocity [" << i << "]=  " << Velocity[i] << endl;

            }
            //cout << "motion[" << 2 << "]=" << motion[2] << endl;


            watervelocity[0][MorisonI] = Velocity[0];
            watervelocity[1][MorisonI] = Velocity[1];
            watervelocity[2][MorisonI] = Velocity[2];


            int checkMorisonValue;
            for (int i = 0;i < 3;i++)
            {
                checkMorisonValue = morison.checkValue(Velocity[i]);
                if (checkMorisonValue == -1)
                {
                    morison.morison_velocity2(time, motion, WaveTmax, WaveDT, waterDepth, Hs, Tp, U1, U2, incide_angle);
                }
            }

        }
        else
        {
            Velocity[0] = watervelocity[0][MorisonI];
            Velocity[1] = watervelocity[1][MorisonI];
            Velocity[2] = watervelocity[2][MorisonI];
            /*
            for (int i = 0;i < 3;i++)
            {
                cout << "Velocity [" << i << "]" << Velocity[i] << endl;

            }
            cout << "motion[" << 2 << "]=" << motion[2] << endl;
            */
        }
        if (MorisonI == MorisonN)
        {
            MorisonI = 0;
            time4 = time;
        }
        else
        {
            MorisonI = MorisonI + 1;
        }
    }
    else if (TF_watervelocity == 0)
    {
        Velocity[0] = 0;
        Velocity[1] = 0;
        Velocity[2] = 0;
    }
    
    

    return 0;
}

extern "C" _declspec(dllexport) int PrintOut(double time, double motion[6], double RoterSpeed, double GenTorq, double Fairlead_Tension[9], double BaseLoads[6], double TowerTopLoads[6],double *xx)
{
    if (time3 != time)
    {
        writeOut write(time, motion, RoterSpeed, GenTorq, Fairlead_Tension,BaseLoads,TowerTopLoads);//���˶������motion.txt�ļ���
        time3 = time;
        *xx = 1.0;
    }
    return 0;
}

extern "C" _declspec(dllexport) int PrintBladeForce(double time, int blade_num, double blade_x[19], double blade_y[19], double blade_z[19], int *xx)
{
    writeOut_Blade_force write( time,blade_num, blade_x, blade_y, blade_z);
    write.write();
    *xx = 1;
    return 0;
}

