#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

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
		
	if (argc < 4){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [cutsfile] [outfile.root] [inputDatafiles]\n";
		return -1;
	}

	ifstream cutsFile(argv[1]);
	
	double AdcToMeVee = 2300;
	double minMeVee;
	
	double minQ2; 
	double maxQ2; 
	double minW2;
	double minWp; 
	double maxWp; 
	
	double minCosThNQ;
	double maxCosThNQ;
	double minAs;
	double maxAs;

	double minPn;
	double maxPn;

	cutsFile >> minMeVee;
	cutsFile >> minPn >> maxPn;
	cutsFile >> minQ2 >> maxQ2;
	cutsFile >> minW2;
	cutsFile >> minWp >> maxWp;
	cutsFile >> minCosThNQ >> maxCosThNQ;
	cutsFile >> minAs >> maxAs;
	
	cutsFile.close();

	
	TChain* inTree = new TChain("tagged");
	for(int i = 3; i < argc; i++) {
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
	int nMult		= 0;
	TClonesArray* nHit = new TClonesArray("bandhit");
	clashit* eHit = new clashit;
	TClonesArray* tag  = new TClonesArray("taghit");
	inTree->SetBranchAddress("Runno"		,&Runno			);
	inTree->SetBranchAddress("Ebeam"		,&Ebeam			);
	inTree->SetBranchAddress("gated_charge"		,&gated_charge		);
	inTree->SetBranchAddress("livetime"		,&livetime		);
	inTree->SetBranchAddress("starttime"		,&starttime		);
	inTree->SetBranchAddress("current"		,&current		);
	//	Neutron branches:
	inTree->SetBranchAddress("nMult"		,&nMult			);
	inTree->SetBranchAddress("nHits"		,&nHit			);
	// 	Electron branches:
	inTree->SetBranchAddress("eHit"			,&eHit			);
	//	Tagg branches:
	inTree->SetBranchAddress("tag"			,&tag			);


	// Output rootfile and tree
	TFile * outFile = new TFile(argv[2],"RECREATE");
	TTree * outTree = inTree->CloneTree(0);

	// Start working on one of the files, looping over all of the events
	for( int ev = 0 ; ev < inTree->GetEntries() ; ev++ ){
		if( ev % 100000 == 0 ) cout << "\t Event " << ev << "/" << inTree->GetEntries() << "\n";

		// Clear all branches before getting the entry from tree
		gated_charge	= 0;
		livetime	= 0;
		starttime 	= 0;
		nMult		= 0;

		eHit->Clear();

		inTree->GetEntry(ev);

		// Check neutron multiplicity
		if( nMult != 1 ) continue;

		// Get band and tag hit from clones array
		bandhit* this_nHit = (bandhit*)nHit->At(0);
		taghit* this_tag = (taghit*)tag->At(0);
		
		if( this_nHit->getStatus() != 0 ) continue;
		if( this_nHit->getTofFadc() == 0 ) continue;
		if( this_nHit->getEdep() < minMeVee*AdcToMeVee ) continue;
		
		// Check electron information
		if( eHit->getPID() != 11 ) continue;
		if( eHit->getCharge() != -1 ) continue;
		if( eHit->getEoP() < 0.17 ) continue;
		if( eHit->getEoP() > 0.3 ) continue;
		if( eHit->getEpcal() < 0.07 ) continue;
		if( eHit->getV() < 9 ) continue;
		if( eHit->getW() < 9 ) continue;
		if( eHit->getVtz() < -8 ) continue;
		if( eHit->getVtz() > 3 ) continue;
		if( eHit->getMomentum() < 2. ) continue;
	
		
		if( eHit->getQ2() < minQ2 ) continue;
		if( eHit->getQ2() > maxQ2 ) continue;
		if( cos(this_tag->getThetaNQ()) < minCosThNQ ) continue;
		if( cos(this_tag->getThetaNQ()) > maxCosThNQ ) continue;
		if( eHit->getW2() < minW2 ) continue;

		// Now only look at neutrons in our signal region:
		if( this_nHit->getTofFadc() < 0. ) continue;
		if( this_tag->getMomentumN().Mag() < minPn ) continue;
		if( this_tag->getMomentumN().Mag() > maxPn ) continue;
		if( this_tag->getMomentumN().Mag() != this_tag->getMomentumN().Mag() ) continue; // check if NaN
		if( this_tag->getWp() < minWp ) continue;
		if( this_tag->getWp() > maxWp ) continue;
		if( this_tag->getAs() < minAs ) continue;
		if( this_tag->getAs() > maxAs ) continue;

		outTree->Fill();	

	} // end loop over events

	outFile->cd();
	outTree->Write();

	outFile->Close();
	return 0;
}
