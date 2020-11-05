#include "GetSpecParameters.h"

Double_t mygaus(Double_t* x, Double_t* par) {
	Double_t xx = x[0];
	Double_t f = par[0] * exp(- 0.5 * ((xx - par[1]) / par[2]) * ((xx - par[1]) / par[2]));
	return f;
}

double CalibGain(TH1D* &InSpec) {
	TH1D* QOut = new TH1D("ChargeSpecOut", "Charge Spectrum", 10000, -50, 150);
	pos = 0;
	int pedestal = InSpec->GetMaximumBin();
	int peak = InSpec->GetBinContent(pedestal);
	int start = 0;
	for (int i = pedestal; i > 0; i --) {
		if (InSpec->GetBinContent(i) < peak * 0.02) {
			start = i;
			break;
		}
	}
	int HalWidth = pedestal - start;
	int end = pedestal + HalWidth * 0.9;
	TF1* GausFunc = new TF1("gaus", "gaus", start * 25, end * 25);
	InSpec->Fit(GausFunc, "R");
	Double_t pars[3];
	GausFunc->GetParameters(pars);
	const Double_t* PedParsEr = GausFunc->GetParErrors();
	TH1D* test = new TH1D("test", "test", 10000, -50, 150);
	for (int i = 0; i < 10000; i ++) {
		Double_t x[1];
		x[0] = i * 25;
		Double_t Spe = InSpec->GetBinContent(i) - mygaus(x, pars) < 0 ? 0 : InSpec->GetBinContent(i) - mygaus(x, pars);
		QOut->SetBinContent(i, Spe);
		// cout << "Debug: " << InSpec->GetBinContent(i) << "\t" << mygaus(x, pars) << endl;
		test->SetBinContent(i, mygaus(x, pars));
		// cout << "Debug: " << i << "\t" << mygaus(x, pars) << endl;
	}
	int Spe_pos = QOut->GetMaximumBin();
	int Spe_peak = QOut->GetBinContent(Spe_pos);
	if (Spe_pos < pedestal + HalWidth) {
		for (int i = Spe_pos; i > 0; i --) {
			if (QOut->GetBinContent(i) < Spe_peak * 0.1) {
				start = i;
				break;
			}
		}
		HalWidth = Spe_pos - start;
		end = Spe_pos + HalWidth;
		cout << start << "\t" << end << endl;
		TF1* tmp = new TF1("tmp_gaus", "gaus", start * 25, end * 25);
		QOut->Fit(tmp, "R");
		Double_t tmp_pars[3];
		tmp->GetParameters(tmp_pars);
		for (int i = 0; i < 10000; i ++) {
			Double_t x[1];
			x[0] = i * 25;
			Double_t Spe = QOut->GetBinContent(i) - mygaus(x, tmp_pars);
			Spe = Spe < 0 ? 0 : Spe;
			QOut->SetBinContent(i, Spe);
		}
	}
	Spe_pos = QOut->GetMaximumBin();
	Spe_peak = QOut->GetBinContent(Spe_pos);
	for (int i = Spe_pos; i < 10000; i ++) {
		if (QOut->GetBinContent(i) < Spe_peak * 0.2) {
			end = i;
			break;
		}
	}
	HalWidth = end - Spe_pos;
	start = Spe_pos - HalWidth * 0.5;
	TF1* SpeFunc = new TF1("Spe_gaus", "gaus", start * 25, end * 25);
	QOut->Fit(SpeFunc, "R");
	Double_t SpePars[3];
	SpeFunc->GetParameters(SpePars);
	const Double_t* SpeParsEr = SpeFunc->GetParErrors();
	Spe_pos = SpePars[1] / 25;
	TH1D* tmp = new TH1D("tmp", "tmp", Spe_pos - pedestal, 0, (Spe_pos - pedestal) * 25);
	for (int i = pedestal; i < Spe_pos; i ++) {
		tmp->SetBinContent(i - pedestal + 1, InSpec->GetBinContent(i));
	}
	// double valley = tmp->GetBinContent(tmp->GetMinimumBin());
	// double peakV = QOut->GetBinContent(QOut->GetMaximumBin());
	// double peakV = SpePars[0];
	// double peakEr = SpeParsEr[0];
	Double_t valleyEr = tmp->GetBinError(tmp->GetMinimumBin());
	double Gain = (SpePars[1] - pars[1]) * 0.01 / 1.6;
	// double GainEr = (SpeParsEr[1] + PedParsEr[1]) * 0.01 / 1.6;
	// double PV = peakV / valley;
	// double PVEr = peakEr / valley + valleyEr * PV / valley;
	// double Resolution = SpePars[2] / (SpePars[1] - pars[1]);
	// double ResEr = SpeParsEr[2] / (SpePars[1] - pars[1])\
	// 				+ SpeParsEr[1] * SpePars[2] / (SpePars[1] - pars[1]) / (SpeParsEr[1] - pars[1])\
	// 				+ PedParsEr[1] * SpePars[2] / (SpePars[1] - pars[1]) / (SpeParsEr[1] - pars[1]);	
	return Gain;
}
