#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TClonesArray.h"

#include "constants.h"
#include "bandhit.h"
#include "clashit.h"
#include "taghit.h"

// For processing data

using namespace std;

int main(int argc, char ** argv){
	// Set style
	gStyle->SetOptFit(1);
		
	if (argc != 4){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [outfile.root] [inputDatafile] [inputHistogramFile] \n";
		return -1;
	}

	// ToF histogram for background normalization
	TH1D * hToF_bac = new TH1D("hToF_bac","hToF_bac",1000,-25,75);

	TFile* inFile = new TFile(argv[2]);
	TTree* inTree = (TTree*)inFile->Get("tagged");

	TFile * outFile = new TFile(argv[1],"RECREATE");

	ifstream histFile(argv[3]);
	

	int nHist1d;	
	histFile >> nHist1d;
	
	vector<TString> histVar;
	vector<TString> histName;
	vector<int> histBins;
	vector<double> histMin;
	vector<double> histMax;
	vector<TString> histTitle;

	TString hVar, hName, hTitle;
	int hBins;
	double hMin, hMax;		

	for(int i = 0; i < nHist1d; i++) {
		histFile >> hVar >> hName >> hBins >> hMin >> hMax >> hTitle;
		histVar.push_back(hVar);
		histName.push_back(hName);
		histBins.push_back(hBins);
		histMin.push_back(hMin);
		histMax.push_back(hMax);
		histTitle.push_back(hTitle);
	}

	int nVarBin;
	histFile >> nVarBin;

	vector<TString> binVar;
	vector<TString> binName;
	vector<double> binMin;
	vector<double> binMax;
	vector<double> binWidth;
	vector<int> nBins;

	TString bVar, bName;
	double bMin, bMax, bWidth;
	int nb;	

	for(int i = 0; i < nVarBin; i++) {
		histFile >> bVar >> bName >> bMin >> bMax >> bWidth;
		binVar.push_back(bVar);
		binName.push_back(bName);
		binMin.push_back(bMin);
		binMax.push_back(bMax);
		binWidth.push_back(bWidth);
		nBins.push_back((binMax[i] - binMin[i] + binWidth[i]/2.)/binWidth[i]);
	}



	vector<TH1F*> hist1d;

	TString thisHistName;	
	TString drawString;
	TString drawCuts = "";
	TString drawOptions = "GOFF";

	for(int i = 0; i < nHist1d; i++) {

		thisHistName = histName[i];
		hist1d.push_back(new TH1F(thisHistName, "", histBins[i], histMin[i], histMax[i]));

		drawString = Form("%s>>%s", histVar[i].Data(), thisHistName.Data());
		drawCuts = "";
		drawOptions = "GOFF";


		inTree->Draw(drawString, drawCuts, drawOptions);

		for(int j = 0; j < nVarBin; j++) {
			for(int k = 0; k < nBins[j]; k++) {

				thisHistName = Form("%s_%s_bin_%i", histName[i].Data(), binName[j].Data(), k);
				hist1d.push_back(new TH1F(thisHistName, "", histBins[i], histMin[i], histMax[i]));

				drawString = Form("%s>>%s", histVar[i].Data(), thisHistName.Data());
				drawCuts = Form("%s >= %f && %s < %f", binVar[j].Data(), binMin[j] + (k*binWidth[j]), binVar[j].Data(), binMin[j] + ((k+1)*binWidth[j]));
				drawOptions = "GOFF";

				inTree->Draw(drawString, drawCuts, drawOptions);

			}
		}
	}
	


	outFile->cd();
	for(int i = 0; i < hist1d.size(); i++) {
		hist1d[i]->Write();	
	}

//	bacnorm.Write("bacnorm");	
//	hToF_bac->Write();

	outFile->Close();
	return 0;
}

























