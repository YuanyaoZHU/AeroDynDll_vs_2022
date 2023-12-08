#pragma once
#ifndef SECOND_ORDER_WAVE_FORCE_H
#define SECOND_ORDER_WAVE_FORCE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include"pch.h"
#include"secondOrderWaveForce.h"

struct WaveData {
    double PER[2];
    double BETA[2];
    int I;
    double MOD;
    double PHS;
    double Re;
    double Im;
    // ���Ը���ʵ����Ҫ��Ӹ���Ĳ��������ֶ�
};

class WAMITDataReader {
public:
    WAMITDataReader(const std::string& filename);

    bool readData(std::vector<WaveData>& waveData);

private:
    std::string filename;
};

class WaveForceCalculator {
public:
    double calculateSecondOrderWaveForce(const std::vector<WaveData>& waveData);
};



#endif // !SECOND_ORDER_WAVE_FORCE_H
