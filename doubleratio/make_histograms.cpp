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
#include "TCut.h"

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

	TFile* inFile = new TFile(argv[2]);
	TTree* inTree = (TTree*)inFile->Get("tagged");

	TVector3* bacnorm;
	TH1D* hToF_bac;
	bool includesBG = false;
	if (inFile->GetListOfKeys()->Contains("hToF_bac")){
		bacnorm = (TVector3*)inFile->Get("bacnorm");
		hToF_bac = (TH1D*)inFile->Get("hToF_bac");
		includesBG = true;
	}

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
	vector<TString> histUnit;

	TString hVar, hName, hTitle, hUnit;
	int hBins;
	double hMin, hMax;		

	for(int i = 0; i < nHist1d; i++) {
		histFile >> hVar >> hName >> hBins >> hMin >> hMax >> hTitle >> hUnit;
		histVar.push_back(hVar);
		histName.push_back(hName);
		histBins.push_back(hBins);
		histMin.push_back(hMin);
		histMax.push_back(hMax);
		histTitle.push_back(hTitle);
		if(hUnit == ".") hUnit = "";
		histUnit.push_back(hUnit);
	}

	int nVarBin;
	histFile >> nVarBin;

	vector<TString> binVar;
	vector<TString> binName;
	vector<double> binMin;
	vector<double> binMax;
	vector<double> binWidth;
	vector<int> nBins;
	vector<TString> binTitle;
	vector<TString> binUnit;

	TString bVar, bName, bTitle, bUnit;
	double bMin, bMax, bWidth;
	int nb;	

	for(int i = 0; i < nVarBin; i++) {
		histFile >> bVar >> bName >> bMin >> bMax >> bWidth >> bTitle >> bUnit;
		binVar.push_back(bVar);
		binName.push_back(bName);
		binMin.push_back(bMin);
		binMax.push_back(bMax);
		binWidth.push_back(bWidth);
		binTitle.push_back(bTitle);
		if(bUnit == ".") bUnit = "";
		binUnit.push_back(bUnit);
		nBins.push_back((binMax[i] - binMin[i] + binWidth[i]/2.)/binWidth[i]);
	}



	vector<TH1F*> hist1d;

	TString thisHistName;	
	TString drawString;
	TString drawOptions = "GOFF";
	TCut KIN_CUT = "tag.Xp < 1. && tag.Wp < 2.5";

	int thisHist = -1;

	for(int i = 0; i < nHist1d; i++) {

		thisHistName = histName[i];
		hist1d.push_back(new TH1F(thisHistName, "", histBins[i], histMin[i], histMax[i]));
		thisHist++;

		drawString = Form("%s>>%s", histVar[i].Data(), thisHistName.Data());
		
		inTree->Draw(drawString, KIN_CUT, drawOptions);
		
		hist1d[thisHist]->GetXaxis()->SetTitle(Form("%s %s", histTitle[i].Data(), histUnit[i].Data()));
		hist1d[thisHist]->SetTitle("");

		for(int j = 0; j < nVarBin; j++) {
			for(int k = 0; k < nBins[j]; k++) {
				
				thisHistName = Form("%s_%s_bin_%i", histName[i].Data(), binName[j].Data(), k);
				hist1d.push_back(new TH1F(thisHistName, "", histBins[i], histMin[i], histMax[i]));
				thisHist++;			
	
				drawString = Form("%s>>%s", histVar[i].Data(), thisHistName.Data());
				TCut BIN_CUT = Form("%s >= %f && %s < %f", binVar[j].Data(), binMin[j] + (k*binWidth[j]), binVar[j].Data(), binMin[j] + ((k+1)*binWidth[j]));
				
				inTree->Draw(drawString, KIN_CUT&&BIN_CUT, drawOptions);
				
				hist1d[thisHist]->GetXaxis()->SetTitle(Form("%s %s", histTitle[i].Data(), histUnit[i].Data()));
				hist1d[thisHist]->SetTitle(Form("%.2f < %s < %.2f %s", binMin[j] + (k*binWidth[j]), binTitle[j].Data(), binMin[j] + ((k+1)*binWidth[j]), binUnit[j].Data()));

			}
		}
	}
	


	outFile->cd();
	for(int i = 0; i < hist1d.size(); i++) {
		hist1d[i]->Write();	
	}

	if(includesBG) {
		bacnorm->Write("bacnorm");	
		hToF_bac->Write();
	}

	outFile->Close();
	return 0;
}

























