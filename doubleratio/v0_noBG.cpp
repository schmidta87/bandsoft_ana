#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
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

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

int main(int argc, char ** argv){
	// Set style
	gStyle->SetOptFit(1);

	if (argc < 2){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [outfile.root] [inputDatafiles]\n";
		return -1;
	}



	const double AdcToMeVee = 2300;
	const double MeVee_cut = 5;
	const double Q2_bin_min = 2;
	const double Q2_bin_max = 10;
	const double CosThetaNQ_bin_min = -1;
	const double CosThetaNQ_bin_max = -0.8;
	const double NMomentum_min = 0.2;
	const double NMomentum_max = 0.65;
	const double Wp_min = 1.8;
	const double Al_min = 1.2;
	const double Al_max = 1.6;
	const double Al_bin_width = 0.05;
	int NAl_bins = (Al_max-Al_min)/Al_bin_width;

	// Output rootfile
	TFile * outFile = new TFile(argv[1],"RECREATE");

	// ToF histogram for background normalization
	TH1D * hXp = new TH1D("hXp","hXp",1000,0,1);
	TH1D * hQ2 = new TH1D("hQ2","hQ2",1000,0,10);
	TH2D * h2XpQ2 = new TH2D("h2XpQ2","h2XpQ2",100,0,10,100,0,1);
	TH2D * h2Q2Wp_hi = new TH2D("h2Q2Wp_hi","h2Q2Wp_hi",100,0,10,100,0,4);
	TH2D * h2Q2Wp_lo = new TH2D("h2Q2Wp_lo","h2Q2Wp_lo",100,0,10,100,0,4);
	TH1D * hWp = new TH1D("hWp","hWp",1000,0,10);
	TH1D * hAs = new TH1D("hAs","hAs",600,1.2,1.6);
	TH1D ** hXpBins = new TH1D*[NAl_bins];
	for( int i = 0 ; i < NAl_bins ; i++){
		hXpBins[i] = new TH1D(Form("hXpBins_%i",i),Form("hXpBins_%i",i),100,0,1);
	}
	TH1D * hAs_hi = new TH1D("hAs_hi","hAs_hi",NAl_bins,Al_min,Al_max);
	TH1D * hAs_lo = new TH1D("hAs_lo","hAs_lo",NAl_bins,Al_min,Al_max);


	// Loop over all the files that are given to me to get the best statistics per bar
	for( int i = 2 ; i < argc ; i++ ){

		TFile * inFile = new TFile(argv[i]);
		TTree * inTree;

		if (inFile->GetListOfKeys()->Contains("T")){
			inTree = (TTree*)inFile->Get("T");
		} else if (inFile->GetListOfKeys()->Contains("tagged")) {
			inTree = (TTree*)inFile->Get("tagged");
		} else {
			cerr << "File has no entries\n";
			return -2;
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



		// Start working on one of the files, looping over all of the events
		cout << "Working on file: " << argv[i] << "\n";
		for( int ev = 0 ; ev < inTree->GetEntries() ; ev++ ){
			if( ev % 1000000 == 0 ) cout << "\ton event " << ev << "\n";

			// Clear all branches before getting the entry from tree
			gated_charge	= 0;
			livetime	= 0;
			starttime 	= 0;
			nMult		= 0;

			// Not necessary after switch to TClonesArray
//			for( int thishit = 0; thishit < maxNeutrons ; thishit++){
//				nHit[thishit].Clear();
//				tag[thishit].Clear();
//			}

			eHit->Clear();

			inTree->GetEntry(ev);

			// Check neutron multiplicity
			if( nMult != 1 ) continue;

			// Get band and tag hit from clones array
			bandhit* this_nHit = (bandhit*)nHit->At(0);
			taghit* this_tag = (taghit*)tag->At(0);
			
			if( this_nHit->getStatus() != 0 ) continue;
			if( this_nHit->getTofFadc() == 0 ) continue;
			if( this_nHit->getEdep() < AdcToMeVee*MeVee_cut ) continue;
			
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

			// Define our 2D bins in Q2 and Theta_nq
			if( eHit->getQ2() < Q2_bin_min ) continue;
			if( eHit->getQ2() > Q2_bin_max ) continue;
			if( cos(this_tag->getThetaNQ()) < CosThetaNQ_bin_min ) continue;
			if( cos(this_tag->getThetaNQ()) > CosThetaNQ_bin_max ) continue;
			if( eHit->getW2() < 2*2 ) continue;

			// Now only look at neutrons in our signal region:
			if( this_tag->getMomentumN().Mag() > NMomentum_max ) continue;
			if( this_tag->getMomentumN().Mag() < NMomentum_min ) continue;
			if( this_tag->getMomentumN().Mag() != this_tag->getMomentumN().Mag() ) continue; // check if NaN
			if( this_tag->getWp() < Wp_min ) continue;
			if( this_tag->getAs() < Al_min ) continue;
			if( this_tag->getAs() > Al_max ) continue;

			// Fill full Xp,As distributions
			if( this_tag->getXp() < 0.35 && this_tag->getXp() > 0.25)
				h2Q2Wp_lo->Fill( eHit->getQ2() , this_tag->getWp() );		
			else if( this_tag->getXp() > 0.5 )
				h2Q2Wp_hi->Fill( eHit->getQ2() , this_tag->getWp() );		
	

			if( this_tag->getWp() > 1.8 && this_tag->getWp() < 3 && eHit->getQ2() > 2 && eHit->getQ2() < 6.5 ){
				hXp->Fill( this_tag->getXp() );
				hAs->Fill( this_tag->getAs() );
				hWp->Fill( this_tag->getWp() );
				hQ2->Fill( eHit->getQ2()  );
				h2XpQ2->Fill( eHit->getQ2() , this_tag->getXp() );

				// Fill Xp distribution for bin in alpha_s
				int binAl = ( this_tag->getAs() - Al_min )/Al_bin_width ;
				hXpBins[binAl]->Fill( this_tag->getXp() );
			
				if( this_tag->getXp() > 0.25 && this_tag->getXp() < 0.35 )
					hAs_lo -> Fill( this_tag->getAs() );
				else if( this_tag->getXp() > 0.5 )
					hAs_hi -> Fill( this_tag->getAs() );
			}


		} // end loop over events

		inFile->Close();
	}// end loop over files

	outFile->cd();
	hXp->Write();
	hAs->Write();
	hQ2->Write();
	h2XpQ2->Write();
	h2Q2Wp_hi->Write();
	h2Q2Wp_lo->Write();
	hWp->Write();
	for( int i = 0 ; i < NAl_bins ; i++){
		hXpBins[i]->Write();
	}
	hAs_lo->Write();
	hAs_hi->Write();



	outFile->Close();
	return 0;
}
