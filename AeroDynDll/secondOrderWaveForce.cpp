#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include"pch.h"
#include"secondOrderWaveForce.h"

WAMITDataReader::WAMITDataReader(const std::string& filename) : filename(filename) {}

bool WAMITDataReader::readData(std::vector<WaveData>& waveData) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return false;
	}

	// 读取每行包含的数据
	double per1, per2, beta1, beta2, mod, phs, re, im;
	int i;
	while (file >> per1 >> per2 >> beta1 >> beta2 >> i >> mod >> phs >> re >> im) {
		WaveData data;
		data.PER[0] = per1;
		data.PER[1] = per2;
		data.BETA[0] = beta1;
		data.BETA[1] = beta2;
		data.I = i;
		data.MOD = mod;
		data.PHS = phs;
		data.Re = re;
		data.Im = im;
		waveData.push_back(data);
	}

	file.close();
	return true;
}

double WaveForceCalculator::calculateSecondOrderWaveForce(const std::vector<WaveData>& waveData) {
	// 整理waveData数据

	// 1. 获得per1和per2相同的数据组
	for (const auto& data : waveData) {


	}
	
	
	// 计算平均漂移力
	
}