#ifndef GetSpecParameters_h
#define GetSpecParameters_h

#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include "math.h"
#include <TCanvas.h>

using namespace std;

Double_t mygaus(Double_t* x, Double_t* par));
double CalibGain(TH1D* &);

#endif
