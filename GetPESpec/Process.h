#ifndef _PROCESS_H_
#define _PROCESS_H_

#define INPUT_DIRECTORY "/junofs/prototype/Data_prtJUNO/Raw_Data"
#define OUTPUT_DIRECTORY "/junofs/prototype/Data_prtJUNO/Output"
#define BUFFER_LENGTH 988
int TotalProcess(int argc, char **argv);
int Initialize();
int GetWaveform();
int GetWaveformSum();
int GetGain();
int GetPeakWindow();
int GetIntegralMap();
int FillCharge();
int GetGainNeed();
int GetCableMap();
int GetFileName(int argc, char **argv, char *input_Filename, char *outputFileName, char *rootfile, char *Gainfile, char* GainOut); //Both input file name and output file name
int GetInputTree();
int GetCaGain(char *file);

#endif
