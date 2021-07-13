#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TCanvas.h"

#include "clashit.h"
#include "bandhit.h"
#include "taghit.h"

#include "bin_edges.h"

using std::cerr;
using std::cout;

int main( int argc, char** argv){

	if( argc != 2 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Dat File] \n";
		return -1;
	}

	// Load files
	TFile inFile_Dat(argv[1]);

	// Load TTrees from files
	TTree * inTree_Dat = (TTree*) inFile_Dat.Get("electrons");

	// Get normalization from files for the background tree
	TVector3 * datnorm = (TVector3*) inFile_Dat.Get("bacnorm");
	
	// Set the input branches for data
	clashit * dat_eHit 		= new clashit;
	inTree_Dat->SetBranchAddress("eHit"		,&dat_eHit		);

	// Define the histograms for yield to fill:
	TFile * outFile = new TFile("test_yield_inclusive.root","RECREATE");
	TH1D * dat_xb_lowQ2 = new TH1D("dat_xb_lowQ2","dat_xb_lowQ2",bins_Xb,Xb_min,Xb_max);
	TH1D * dat_xb_highQ2 = new TH1D("dat_xb_highQ2","dat_xb_highQ2",bins_Xb,Xb_min,Xb_max);

	// Loop over the data file:
	for( int event = 0 ; event < inTree_Dat->GetEntries() ; ++event ){
		dat_eHit	->Clear();

		inTree_Dat->GetEntry(event);
		
		double Q2	= dat_eHit->getQ2();
		double Xb	= dat_eHit->getXb();

		if( Q2 >= Q2Bins[0] && Q2 < Q2Bins[1] ){
			dat_xb_lowQ2->Fill( Xb );
		}
		else if(Q2 >= Q2Bins[1] && Q2 < Q2Bins[2] ){
			dat_xb_highQ2->Fill( Xb );
		}
		else continue;
	}


	// Print the background subtracted distributions to a pdf file:
	TCanvas * c 		= new TCanvas("c",		"",800,600);
	c->Divide(1,2);
	outFile->cd();

	c->cd(1);
	dat_xb_lowQ2->Write();
	dat_xb_lowQ2->SetLineWidth(3);
	dat_xb_lowQ2->Scale(1);
	dat_xb_lowQ2->Draw();

	c->cd(2);
	dat_xb_highQ2->Write();
	dat_xb_highQ2->SetLineWidth(3);
	dat_xb_highQ2->Scale(1);
	dat_xb_highQ2->Draw();

	c->Print("rawyield_inclusive.pdf");

	outFile->Close();


	inFile_Dat.Close();
	return 1;
}
