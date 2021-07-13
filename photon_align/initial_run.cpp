#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"

#include "constants.h"
#include "bandhit.h"
#include "TClonesArray.h"
#include "calib_helper.h"
#include "kinematic_cuts.h"

// For processing data

using namespace std;
shiftsReader shifts;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

int main(int argc, char ** argv){
	// Set style
	gStyle->SetOptFit(1);

	if (argc < 4){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [TDC(0) or FADC(1)] [outputTxtfile] [outputPDFfile] [inputDatafiles]\n";
		return -1;
	}

	int TDCorFADC = atoi(argv[1]);

	// Histograms for every run that sum over all bars
	const int nRuns = argc-4;
	TH1D ** ToF_spec 	= new TH1D*[nRuns];
	TF1 ** ToF_fits		= new TF1*[nRuns];
	TF1 ** ToF_fits_it	= new TF1*[nRuns];
	TCanvas ** cRun		= new TCanvas*[nRuns];
	for(int run = 0; run < nRuns; run++){
		ToF_spec[run] = new TH1D( Form("ToF_spec_%i",run), Form("ToF_spec_%i",run), 800, -15, 25);
	}
	cout << "created histograms\n";

	// Load the bar shifts that were calculated at the initial bar code
	shifts.LoadInitBarFadc	("../include/FADC_pass1v0_initbar.txt");
	shifts.LoadInitBar	("../include/TDC_pass1v0_initbar.txt");
	double * FADC_INITBAR 	= (double*) shifts.getInitBarFadc();
	double * TDC_INITBAR 	= (double*) shifts.getInitBar();
	cout << "loaded shifts\n";


	// Loop over all the files that are given to me to get the best statistics per bar
	for( int i = 4 ; i < argc ; i++ ){
		TFile * inFile = new TFile(argv[i]);
		if (!(inFile->GetListOfKeys()->Contains("calib"))){
			cerr << "File has no entries\n";
			return -2;
		}
		TTree * inTree = (TTree*)inFile->Get("calib");

		//	Event info:
		int Runno		= 0;
		double Ebeam		= 0;
		double gated_charge	= 0;
		double livetime		= 0;
		double starttime	= 0;
		double current		= 0;
		bool goodneutron 	= false;
		int nleadindex	 	= -1;
		// 	Neutron info:
		int nMult		= 0;
		TClonesArray * nHits = new TClonesArray("bandhit");
		// 	Event branches:
		inTree->SetBranchAddress("Runno"		,&Runno			);
		inTree->SetBranchAddress("Ebeam"		,&Ebeam			);
		inTree->SetBranchAddress("gated_charge"		,&gated_charge		);
		inTree->SetBranchAddress("livetime"		,&livetime		);
		inTree->SetBranchAddress("starttime"		,&starttime		);
		inTree->SetBranchAddress("current"		,&current		);
		//	Neutron branches:
		inTree->SetBranchAddress("nMult"		,&nMult			);
		inTree->SetBranchAddress("nHits"		,&nHits			);
		inTree->SetBranchAddress("goodneutron"		,&goodneutron		);
		inTree->SetBranchAddress("nleadindex"		,&nleadindex		);

		// Start working on one of the files, looping over all of the events
		cout << "Working on file: " << argv[i] << "\n";
		for( int ev = 0 ; ev < inTree->GetEntries() ; ev++ ){
			if( ev % 1000000 == 0 ) cout << "\ton event " << ev << "\n";

			// Clear all branches before getting the entry from tree
			Runno		= 0;
			Ebeam		= 0;
			gated_charge	= 0;
			livetime	= 0;
			starttime	= 0;
			current		= 0;
			goodneutron 	= false;
			nleadindex	= -1;
			// 	Neutron info:
			nMult		= 0;
			nHits->Clear();
			
			inTree->GetEntry(ev);

			// Only look for "photon" events
			bandhit * this_photon = (bandhit*) nHits->At(0);
			if( nMult 			!= 1 ) continue;
			if( this_photon->getStatus() 	!= 0 ) continue;
			if( TDCorFADC == 0 )
				if( this_photon->getTof() 	== 0 ) continue;
			if( TDCorFADC == 1 )
				if( this_photon->getTofFadc() 	== 0 ) continue;
			if( this_photon->getEdep()	< 2 ) 	continue;

			int sector 	= this_photon->getSector();
			int layer 	= this_photon->getLayer();
			int component 	= this_photon->getComponent();

			double time 	= 0;
			if( TDCorFADC == 0 )
				time 	= this_photon->getTof()		- TDC_INITBAR[sector*100+layer*10+component];
			if( TDCorFADC == 1 )
				time 	= this_photon->getTofFadc()	- FADC_INITBAR[sector*100+layer*10+component];
			double dL 	= this_photon->getDL().Mag();

			double tof = time - dL/cAir;

			ToF_spec[i-4]->Fill(tof);

		} // end loop over events

		inFile->Close();
	}// end loop over files



	ofstream out_file;
	out_file.open(argv[2]);
	TCanvas * c0 = new TCanvas("c0","c0",900,900);
	TString openname = string(argv[3]) + "(";
	c0 -> Print(openname);
	for(int run = 0; run < nRuns; run++){
		// Create canvas
		cRun[run] = new TCanvas(Form("Run %i",run+1),Form("Run %i",run+1),900,900);
		int thisRun = getRunNumber(argv[run+4]);
		ToF_spec[run]->SetTitle(Form("ToF_spec_%i",thisRun));


		if( ToF_spec[run]->Integral() == 0 ){
			out_file << thisRun <<  " " 
				<< 0 << " " << 0 << " " 
				<< 0 << " " << 0 << " "
				<< 0 << " " << -1 << "\n";
			continue;
		}


		// Get the min and max of the fit based on assuming 0.3ns resolution and the peak position
		double sig_guess = 0.3;
		double max = ToF_spec[run]->GetMaximum();
		double max_pos = ToF_spec[run]->GetXaxis()->GetBinCenter( ToF_spec[run]->GetMaximumBin() );

		double min_fit = max_pos - 5;
		double max_fit = max_pos + 2.*sig_guess;
		ToF_fits[run] = new TF1(Form("ToF_fits_%i",run),"pol0+gaus(1)",min_fit,max_fit);

		// Set parameters of the fit before fitting:
		double background_lvl = 0.;
		for( int i = 1; i < 25; i++){
			background_lvl += ToF_spec[run]->GetBinContent(i);
		}
		background_lvl /= 24;
		// background level:
		ToF_fits[run]->SetParameter(0,background_lvl);
		// constant of gaus:
		ToF_fits[run]->SetParameter(1,max);
		// mean of gaus:
		ToF_fits[run]->SetParameter(2,max_pos);
		// sigma of gaus:
		ToF_fits[run]->SetParameter(3,sig_guess);

		ToF_spec[run]->Fit(ToF_fits[run],"QESR");

		// Now do another iteration of fits with the current parameters only if the parameters are reasonable
		if( ToF_fits[run]->GetParameter(3) < 0.6 ){
			sig_guess = ToF_fits[run]->GetParameter(3);
			max_fit = ToF_fits[run]->GetParameter(2) - 5;
			min_fit = ToF_fits[run]->GetParameter(2) + 1.5*sig_guess;
			ToF_fits_it[run] = new TF1(Form("ToF_fits_%i_it",run),"pol0+gaus(1)",min_fit,max_fit);
			ToF_fits_it[run]->SetLineColor(4);
			ToF_fits_it[run]->SetParameter(0, ToF_fits[run]->GetParameter(0) );
			ToF_fits_it[run]->SetParameter(1, ToF_fits[run]->GetParameter(1) );
			ToF_fits_it[run]->SetParameter(2, ToF_fits[run]->GetParameter(2) );
			ToF_fits_it[run]->SetParameter(3, ToF_fits[run]->GetParameter(3) );
			ToF_spec[run]->Fit(ToF_fits_it[run],"QESR");

			out_file << (thisRun) << " "  
				<< ToF_fits_it[run]->GetParameter(0) << " " << ToF_fits_it[run]->GetParameter(1) << " " 
				<< ToF_fits_it[run]->GetParameter(2) << " " << ToF_fits_it[run]->GetParameter(3) << " "
				<< ToF_spec[run]->Integral() << " " << 1 << "\n";
		}
		else{
			out_file << (thisRun) << " "  
				<< ToF_fits[run]->GetParameter(0) << " " << ToF_fits[run]->GetParameter(1) << " " 
				<< ToF_fits[run]->GetParameter(2) << " " << ToF_fits[run]->GetParameter(3) << " "
				<< ToF_spec[run]->Integral() << " " << 0 << "\n";

		}

		cRun[run]->cd(1);
		ToF_spec[run]->Draw();
		//ToF_fits[run]->Draw("SAME");
		cRun[run]->Update();
		cRun[run]->Modified();
		cRun[run] -> Print(argv[3]);
		//cRun[run]->Write();
	}
	TString endname = string(argv[3]) + ")";
	c0 -> Print(endname);

	out_file.close();



	return 0;
}
