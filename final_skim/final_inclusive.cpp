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
#include "kinematic_cuts.h"
#include "bandhit.h"
#include "clashit.h"
#include "taghit.h"

// For processing data

using namespace std;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

int main(int argc, char ** argv){
	if (argc < 2){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [outfile.root] [inputDatafiles]\n";
		return -1;
	}

	// Conditions for a final accepted event
	TCut ePID 	= Form("eHit->getPID() == %i",					ECUT_PID);
	TCut eCharge 	= Form("eHit->getCharge() == %i",				ECUT_charge); 
	TCut eEoP	= Form("eHit->getEoP() > %f && eHit->getEoP() < %f",		ECUT_EoP_min,		ECUT_EoP_max);
	TCut eEpcal	= Form("eHit->getEpcal() > %f",					ECUT_Epcal_min);
	TCut eVW	= Form("eHit->getV() > %f && eHit->getW() > %f",		ECUT_V_min,		ECUT_W_min);
	TCut eVtx	= Form("eHit->getVtz() > %f && eHit->getVtz() < %f",		ECUT_vtx_min,		ECUT_vtx_max);
	TCut eMom	= Form("eHit->getMomentum() > %f && eHit->getMomentum() < %f",	ECUT_pE_min,		ECUT_pE_max);
	TCut eQ2	= Form("eHit->getQ2() > %f && eHit->getQ2() < %f",		ECUT_Q2_min,		ECUT_Q2_max);
	TCut eW		= Form("eHit->getW2() > %f",					ECUT_W2_min);
	TCut inclusive	= ePID && eCharge && eEoP && eEpcal && eVW && eVtx && eMom && eQ2 && eW;
	//TString cut = Form("%s",inclusive.GetTitle());
	TString cut = Form("%s",inclusive.GetTitle());

	// Load input files
	TChain* inTree = new TChain("electrons");
	for( int i = 2 ; i < argc; i++ ){
		cout << "Adding file " << argv[i] << endl;
		inTree->Add(argv[i]);
	}

	// Create an output rootfile
	TFile * outFile = new TFile(argv[1],"RECREATE");
	inTree->LoadTree(0);
	TTree * outTree = inTree->CloneTree(0);

	// Create the good event list from our cuts defined above
	inTree->Draw(">>goodEvents",cut);
	TEventList * goodEvents = (TEventList*) gDirectory->Get("goodEvents");
	int nEvents = goodEvents->GetN();

	// Loop over all the good events and write the output tree
	for( int ev = 0 ; ev < nEvents ; ev++ ){
		if( ev % 100000 == 0 ) cout << "\ton event " << ev << "\n";

		int entry = goodEvents->GetEntry(ev);
		inTree->GetEntry(entry);

		outTree->Fill();

	} // end loop over events


	// Write and close
	outTree->Write();
	outFile->Close();
	return 0;
}
