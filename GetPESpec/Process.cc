#include "Process.h"
#include "GetSpecParameters.h"
#include <iostream>
#include<iomanip>
#include <stdio.h>
#include <string>
#include <fstream>
#include <stdint.h>
#include <cstdio>
#include "TMath.h"
#include <TSystem.h>
#include <TSystem.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
using namespace std;

TFile * File_In;
TFile * File_Out;
TTree * InputTree;

int TotalEventNumber;
int Run_Number;
int Event_Id;
int Trigger_Id;
int Module_Id;
int Board_Id;
int Channal_Id;
uint16_t Time_Stamp[5];
int Trigger_Interval;
int N_Dimention;
uint16_t Channal_Data[1008];
uint16_t Channal_Index[1008];

TH1F * waveform;
TH1F * average_waveform[112];
TH1F * charge_spectrum[112];
TH1F * PENum[112];
TH1F * Amp[112];
TH1F * PeSpectrum;
TH1F * Spec_3inch;
TH1F * Spec_8inch;
TH1F * Spec_20inch;
TH1F * temp_w[112];
double IntegralMap[BUFFER_LENGTH];

TH1F * Inte = new TH1F("IntegralMap", "IntegralMap", 988, 0, 988);
int ID_now;
int Peak_Start[112];
int Peak_End[112];
int Inte_Start[112];
int Inte_End[112];
double Gain[112];
double charge[112];
int Ch_count = 0;
double PeN;
double PeN_3;
double PeN_8;
double PeN_20;

int map_ch[112];
int map_di[112];
int map_PMT[112];
int map_amp[112];
int useful_flag[112];

int ChId;

int TotalProcess(int argc, char **argv)
{
	int num;
	char Input[100], Output[100], RootFile[100], directory[100], temp[100], Gainfile[100], GainOut[100];
	char n[100];
	ofstream out;
	ofstream gainout;
	num = GetFileName(argc, argv, Input, Output, RootFile, Gainfile, GainOut);
	File_Out = new TFile(RootFile, "recreate");
	if (!num)
		return 0;
	TTree * OutputTree = new TTree("Analyzed_Data", "Charge spectrum");
	for (int i = 0; i < 112; i ++)
	{
		char n[100], m[100];
		sprintf(n, "Charge%d", i + 256);
		sprintf(m, "Charge%d/D", i + 256);
		OutputTree->Branch(n, &charge[i], m);
	}
	Initialize();
	GetCaGain(Gainfile);
	for (int i = 0; i < 112; i ++)
	{
		if (map_amp[i] == 1)
			useful_flag[i] = 1;
		else
			useful_flag[i] = 0;
	}
	useful_flag[67] = 0;
	for (int k = 1; k <= num; k ++)
	{
		if (k < 10)
			sprintf(n, "%s000%d.root", Input, k);
		else if (k < 100)
			sprintf(n, "%s00%d.root", Input, k);
		else if (k < 1000)
			sprintf(n, "%s0%d.root", Input, k);
		else
			sprintf(n, "%s%d.root", Input, k);
		strcpy(directory, INPUT_DIRECTORY);
		sprintf(temp, "%s/%s", directory, n);
		File_In = new TFile(temp);
		File_In->ls();
		if (!GetInputTree())
			return 0;
		ID_now = 0;
		for (int i = 0; i < 10000; i ++)
		{
			InputTree->GetEntry(i);
			GetWaveform();
			GetWaveformSum();
		}
		GetPeakWindow();
		//for (int i = 0; i < 112; i ++)
		//    cout << Peak_Start[i] << "  " << Peak_End[i] << endl;
		PeN = 0;
		PeN_20 = 0;
		PeN_3 = 0;
		PeN_8 = 0;
		for (int i = 0; i < TotalEventNumber; i ++)
		{
			InputTree->GetEntry(i);
			GetWaveform();
			GetIntegralMap();
			FillCharge();
			PENum[ChId - 256]->Fill(charge[ChId - 256] / Gain[ChId - 256] / 1.6);
			Amp[ChId - 256]->Fill(waveform->GetBinContent(waveform->GetMaximumBin()));
			if (useful_flag[ChId - 256])
			{
				if (waveform->GetBinContent(waveform->GetMaximumBin()) < 1000)
				{
					PeN += charge[ChId - 256] / Gain[ChId - 256] / 1.6;
					switch(map_di[ChId - 256])
					{
						case (3):
							PeN_3 += charge[ChId - 256] / Gain[ChId - 256] / 1.6;
							break;
						case (8):
							PeN_8 += charge[ChId - 256] / Gain[ChId - 256] / 1.6;
							break;
						case (20):
							PeN_20 += charge[ChId - 256] / Gain[ChId - 256] / 1.6;
							break;                        
					}
				}                    
				else if (ChId - 254 <= 111 && map_PMT[ChId - 254] == map_PMT[ChId - 256])
				{
					PeN += charge[ChId - 254] / Gain[ChId - 254] / 1.6;
					switch(map_di[ChId - 256])
					{
						case (3):
							PeN_3 += charge[ChId - 254] / Gain[ChId - 254] / 1.6;
							break;
						case (8):
							PeN_8 += charge[ChId - 254] / Gain[ChId - 254] / 1.6;
							break;
						case (20):
							PeN_20 += charge[ChId - 254] / Gain[ChId - 254] / 1.6;
							break;                        
					}
				}        
				else if (ChId - 255 <= 111 && map_PMT[ChId - 255] == map_PMT[ChId - 256])
				{
					PeN += charge[ChId - 255] / Gain[ChId - 255] / 1.6;
					switch(map_di[ChId - 255])
					{
						case (3):
							PeN_3 += charge[ChId - 255] / Gain[ChId - 255] / 1.6;
							break;
						case (8):
							PeN_8 += charge[ChId - 255] / Gain[ChId - 255] / 1.6;
							break;
						case (20):
							PeN_20 += charge[ChId - 255] / Gain[ChId - 255] / 1.6;
							break;                        
					}
				}    
				else
				{
					PeN += charge[ChId - 256] / Gain[ChId - 256] / 1.6;
					switch(map_di[ChId - 256])
					{
						case (3):
							PeN_3 += charge[ChId - 256] / Gain[ChId - 256] / 1.6;
							break;
						case (8):
							PeN_8 += charge[ChId - 256] / Gain[ChId - 256] / 1.6;
							break;
						case (20):
							PeN_20 += charge[ChId - 256] / Gain[ChId - 256] / 1.6;
							break;                        
					}
				}   
			}
			if (Ch_count == 112)
			{
				OutputTree->Fill();
				PeSpectrum->Fill(PeN);
				Spec_3inch->Fill(PeN_3);
				Spec_8inch->Fill(PeN_8);
				Spec_20inch->Fill(PeN_20);
				// cout << PeN << endl;
				// if (PeN_20)
				// 	cout << "PeN:" << PeN << " " << "PeN_20:" << PeN_20 << " " << "PeN_3:" << PeN_3 << " " << "PeN_8" << PeN_8 << " " << endl;
				PeN = 0;
				PeN_20 = 0;
				PeN_3 = 0;
				PeN_8 = 0;
				Ch_count = 0;
			}
		}
	}
	out.open(Output, std::ofstream::out);
	gainout.open(GainOut);
	GetGain();
	File_Out->cd();
	average_waveform[56]->Write();
	average_waveform[40]->Write();

	PeSpectrum->Write();
	Spec_3inch->Write();
	Spec_8inch->Write();
	Spec_20inch->Write();
	waveform->Write();
	Inte->Write();
	OutputTree->Write();
	for (int i = 0; i < 112; i ++) {
		gainout << Gain[i] << endl;
	}
	for (int i = 0; i < 5000; i ++)
		out << PeSpectrum->GetBinContent(i + 1) << "\t";
	out << "\n";
	for (int i = 0; i < 5000; i ++)
		out << Spec_3inch->GetBinContent(i + 1) << "\t";
	out << "\n";
	for (int i = 0; i < 5000; i ++)
		out << Spec_8inch->GetBinContent(i + 1) << "\t";
	out << "\n";
	for (int i = 0; i < 5000; i ++)
		out << Spec_20inch->GetBinContent(i + 1) << "\t";
	out << "\n";
	for (int i = 0; i < 112; i ++)
	{
		charge_spectrum[i]->Write();
		PENum[i]->Write();
		Amp[i]->Write();
		//cout << PENum[i]->GetMean() << endl;
		//out << PENum[i]->GetMean() << endl;
	}
	File_Out->Close();
	cout << "Output to TXT" << endl;
	out.close();
	return 0;
}

int GetCaGain(char *file)
{
	ifstream gin;
	int useless;
	gin.open(file);
	if (gin.fail())
	{
		cout << "error file" << endl;
		gin.close();
		return 0;
	}
	for (int i = 0; i < 112; i ++)
	{
		gin >> useless >> Gain[i];
		//cout << Gain[i] << endl;
	}
	return 1;
}

int Initialize()
{
	char name[20];
	for (int i = 0; i < 112; i ++)
	{
		sprintf(name, "Charge_Spectrum %d", i + 256);
		charge_spectrum[i] = new TH1F(name, name, 10000, -50, 150);
		sprintf(name, "Channal %d", i + 256);
		average_waveform[i] = new TH1F(name, name, 988, 0, 988);
		sprintf(name, "PE number %d", i + 256);
		PENum[i] = new TH1F(name, name, 10000, -50, 5000);
		sprintf(name, "Amplitude spectrum %d", i + 256);
		Amp[i] = new TH1F(name, name, 10000, -5, 1000);
		sprintf(name, "Temp Waveform %d", i + 256);
		temp_w[i] = new TH1F(name, name, 988, 0, 988);
	}
	waveform = new TH1F("waveform", "waveform", 988, 0, 988);
	PeSpectrum = new TH1F("PeSpectrum", "PeSpectrum", 5000, 0, 5000);
	Spec_3inch = new TH1F("3-inch PMT PeSpectrum", "3-inch PMT PeSpectrum", 5000, 0, 5000);
	Spec_8inch = new TH1F("8-inch PMT PeSpectrum", "8-inch PMT PeSpectrum", 5000, 0, 5000);
	Spec_20inch = new TH1F("20-inch PMT PeSpectrum", "20-inch PMT PeSpectrum", 5000, 0, 5000);
	ifstream in;
	in.open("cable_map.txt");
	for (int i = 0; i < 112; i ++)
	{
		in >> map_ch[i];
		in >> map_di[i];
		in >> map_PMT[i];
		in >> map_amp[i];
	}
	return 0;
}
int GetWaveform()
{
	ChId = Module_Id * 16 + Board_Id * 4 + Channal_Id;
	double baseline = 0;
	for (int i = 20; i < 1008; i ++)
		waveform->SetBinContent(i - 19, Channal_Data[i]);
	for (int i = 0; i < 100; i ++)
	{
		baseline += waveform->GetBinContent(i + 1);
	}
	baseline /= 100;
	for (int i = 0; i < BUFFER_LENGTH; i ++)
		waveform->SetBinContent(i + 1, waveform->GetBinContent(i + 1) - baseline);
	return 1;
}
int GetWaveformSum()
{
	if (waveform->GetBinContent(waveform->GetMaximumBin()) > 5)
		for (int i = 0; i < BUFFER_LENGTH; i ++)
			average_waveform[ChId - 256]->SetBinContent(i + 1, average_waveform[ChId - 256]->GetBinContent(i + 1) + waveform->GetBinContent(i + 1));
	return 1;
}
int GetPeakWindow()
{
	for (int i = 0; i < 112; i ++)
	{
		//cout << average_waveform[i]->GetMaximumBin() << "  " << average_waveform[i]->GetBinContent(average_waveform[i]->GetMaximumBin()) << endl;
		Peak_Start[i] = average_waveform[i]->GetMaximumBin() - 100;
		Peak_End[i] = average_waveform[i]->GetMaximumBin() + 150;
	}
	return 1;
}
int GetIntegralMap()
{
	double sum = 0;
	for (int i = 0; i < 988; i ++)
	{
		sum += waveform->GetBinContent(i + 1);
		IntegralMap[i] = sum;
		Inte->SetBinContent(i + 1, sum);
		if (i + 1 >= Peak_Start[ChId - 256] && i + 1 <= Peak_End[ChId - 256])
			temp_w[ChId - 256]->SetBinContent(i + 1, waveform->GetBinContent(i + 1));
	}
	int p = temp_w[ChId - 256]->GetMaximumBin();
	//cout << "Start " << Peak_Start[ChId - 256] << endl;
	//cout << "End " << Peak_End[ChId - 256] << endl;
	//cout << "P " << p << endl;
	//cout << temp_w[ChId - 256]->GetBinContent(temp_w[ChId - 256]->GetMaximumBin()) << endl;
	for (int i = 0; i < 100; i ++)
	{
		if (temp_w[ChId - 256]->GetBinContent(p - i) < 0.1 * temp_w[ChId - 256]->GetBinContent(p))
		{
			Inte_Start[ChId - 256] = p - i;
			break;
		}
	}
	for (int i = 0; i < 150; i ++)
	{
		if (temp_w[ChId - 256]->GetBinContent(p + i) < 0.1 * temp_w[ChId - 256]->GetBinContent(p))
		{
			Inte_End[ChId - 256] = p + i;
			break;
		}
	}
	return 0;
}
int FillCharge()
{
	charge[ChId - 256] = (IntegralMap[Inte_End[ChId - 256]] - IntegralMap[Inte_Start[ChId - 256]]) * (map_amp[ChId - 256] == 1 ? 0.1 : (map_amp[ChId - 256] == 2 ? 0.8 : 3.2)) / 50;
	charge_spectrum[ChId - 256]->Fill(charge[ChId - 256]);
	Ch_count ++;
	return 0;
}
int GetGain()
{
	for (int i = 0; i < 112; i ++)
	{
		// Gain[i] = charge_spectrum[i]->GetMaximumBin() * 0.15 / 1.6;
		Gain[i] = CalibGain(charge_spectrum[i]);
	}
	return 0;
}
int GetFileName(int argc, char **argv, char *input_Filename, char *output_FileName, char *rootfile, char *Gainfile, char* GainOut)
{
	char temp[200];
	int num;
	if (argc == 8)
	{
		strcpy(temp, argv[1]);
		strcat(temp, ".0000");
		strcat(temp, argv[2]);
		strcat(temp, ".daq.RAW._SFO-1._");
		sscanf(argv[3], "%d", &num);
		strcpy(input_Filename, temp);
		// cout << input_Filename << endl;
		strcpy(output_FileName, argv[4]);
		strcpy(rootfile, argv[5]);
		strcpy(Gainfile, argv[6]);
		strcpy(GainOut, argv[7]);
	}
	else
	{
		cout << "Wrong number of parameters!" << endl;
		return 0;
	}
	if (!strstr(output_FileName, ".txt"))
	{
		cout << "The output file is not a txt file." << endl;
		return 0;
	}
	if (!strstr(rootfile, ".root"))
	{
		cout << "The output file is not a txt file." << endl;
		return 0;
	}
	return num;
}

int GetInputTree()
{
	File_In->GetObject("fadc_prototype", InputTree);
	if (!InputTree)
	{
		cout << "Get root data in wave error" << endl;
		return 0;
	}
	TotalEventNumber = InputTree->GetEntries(); 
	cout << "Root event total number is " << TotalEventNumber << endl;
	InputTree->SetBranchAddress("runNumber", &Run_Number);
	InputTree->SetBranchAddress("eventID", &Event_Id);
	InputTree->SetBranchAddress("triggerID", &Trigger_Id);
	InputTree->SetBranchAddress("moduleID", &Module_Id);
	InputTree->SetBranchAddress("bdID", &Board_Id);
	InputTree->SetBranchAddress("chID", &Channal_Id);
	InputTree->SetBranchAddress("timeStamp", Time_Stamp);
	cout << "timeStamp" << endl;
	InputTree->SetBranchAddress("triggerInterval", &Trigger_Interval);
	InputTree->SetBranchAddress("nDimension", &N_Dimention);
	InputTree->SetBranchAddress("chData", Channal_Data);
	cout << "chData" << endl;
	InputTree->SetBranchAddress("chIndex", Channal_Index);
	cout << "chIndex" << endl;

	return 1;
}

