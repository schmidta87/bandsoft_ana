#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TClonesArray.h"
#include "TEventList.h"
#include "TCut.h"
#include "TChain.h"

#include "constants.h"
#include "bandhit.h"
#include "clashit.h"
#include "taghit.h"

// For processing data

using namespace std;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

int main(int argc, char ** argv){
	// Set style
	gStyle->SetOptFit(1);

	if (argc < 3){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [outfile.root] [inputDatafiles]\n";
		return -1;
	}

	// Electron cuts
	const int PID = 11;
	const int charge = -1;
	const double EoP_min = 0.17;
	const double EoP_max = 0.3;
	const double Epcal_min = 0.07;
	const double V_min = 15;
	const double W_min = 15;
	const double vtx_min = -8;
	const double vtx_max =  3;
	const double pE_min = 3;
	const double pE_max = 10.6;
	const double Q2_min = 2;
	const double Q2_max = 10;
	const double W2_min = 2*2;
	TCut ePID 	= Form("eHit->getPID() == %i",PID);
	TCut eCharge 	= Form("eHit->getCharge() == %i",charge); 
	TCut eEoP	= Form("eHit->getEoP() > %f && eHit->getEoP() < %f",EoP_min,EoP_max);
	TCut eEpcal	= Form("eHit->getEpcal() > %f",Epcal_min);
	TCut eVW	= Form("eHit->getV() > %f && eHit->getW() > %f",V_min,W_min);
	TCut eVtx	= Form("eHit->getVtz() > %f && eHit->getVtz() < %f",vtx_min,vtx_max);
	TCut eMom	= Form("eHit->getMomentum() > %f && eHit->getMomentum() < %f",pE_min,pE_max);
	TCut eQ2	= Form("eHit->getQ2() > %f && eHit->getQ2() < %f",Q2_min,Q2_max);
	TCut eW		= Form("eHit->getW2() > %f",W2_min);
	TCut inclusive	= ePID && eCharge && eEoP && eEpcal && eVW && eVtx && eMom && eQ2 && eW;
	TString cut = Form("%s",inclusive.GetTitle());

	// Load input files
	TChain* inTree = new TChain("electrons");
	for( int i = 2 ; i < argc; i++ ){
		cout << "Adding file " << argv[i] << endl;
		inTree->Add(argv[i]);
	}
	// Initialize the input branches
	int Runno		= 0;
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	clashit* eHit = new clashit;
	inTree->SetBranchAddress("Runno"		,&Runno			);
	inTree->SetBranchAddress("Ebeam"		,&Ebeam			);
	inTree->SetBranchAddress("gated_charge"		,&gated_charge		);
	inTree->SetBranchAddress("livetime"		,&livetime		);
	inTree->SetBranchAddress("starttime"		,&starttime		);
	inTree->SetBranchAddress("current"		,&current		);
	// 	Electron branches:
	inTree->SetBranchAddress("eHit"			,&eHit			);

	// Output rootfile
	TFile * outFile = new TFile(argv[1],"RECREATE");
	inTree->LoadTree(0);
	TTree * outTree = inTree->CloneTree(0);

	inTree->Draw(">>goodEvents",cut);
	TEventList * goodEvents = (TEventList*) gDirectory->Get("goodEvents");
	int nEvents = goodEvents->GetN();

	for( int ev = 0 ; ev < nEvents ; ev++ ){
		if( ev % 100000 == 0 ) cout << "\ton event " << ev << "\n";

		// Clear all branches before getting the entry from tree
		Runno 		= 0;
		Ebeam		= 0;
		gated_charge	= 0;
		livetime	= 0;
		starttime 	= 0;
		current		= 0;
		eHit->Clear();

		int entry = goodEvents->GetEntry(ev);
		inTree->GetEntry(entry);

		outTree->Fill();

	} // end loop over events


	outTree->Write();
	outFile->Close();
	return 0;
}
