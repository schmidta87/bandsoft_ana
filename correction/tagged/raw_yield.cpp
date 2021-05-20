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

	if( argc != 3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Dat File] [Bac File]\n";
		return -1;
	}

	// Load files
	TFile inFile_Dat(argv[1]);
	TFile inFile_Bac(argv[2]);

	// Load TTrees from files
	TTree * inTree_Dat = (TTree*) inFile_Dat.Get("tagged");
	TTree * inTree_Bac = (TTree*) inFile_Bac.Get("tagged");

	// Get normalization from files for the background tree
	TVector3 * datnorm = (TVector3*) inFile_Dat.Get("bacnorm");
	TVector3 * bacnorm = (TVector3*) inFile_Bac.Get("bacnorm");
	//inTree_Bac->SetWeight( datnorm->X() / bacnorm->X() );
	double bac_weight = datnorm->X() / bacnorm->X();
	
	// Set the input branches for data
	clashit * dat_eHit 		= new clashit;
	TClonesArray* dat_tagged 	= new TClonesArray("taghit");
	int dat_nleadindex		= 0;
	inTree_Dat->SetBranchAddress("tag"		,&dat_tagged		);
	inTree_Dat->SetBranchAddress("eHit"		,&dat_eHit		);
	inTree_Dat->SetBranchAddress("nleadindex"	,&dat_nleadindex	);

	// Set the input branches for background
	clashit * bac_eHit 		= new clashit;
	TClonesArray* bac_tagged 	= new TClonesArray("taghit");
	int bac_nleadindex		= 0;
	inTree_Bac->SetBranchAddress("tag"		,&bac_tagged		);
	inTree_Bac->SetBranchAddress("eHit"		,&bac_eHit		);
	inTree_Bac->SetBranchAddress("nleadindex"	,&bac_nleadindex	);


	// Define the histograms for yield to fill:
	TFile * outFile = new TFile("test_yield.root","RECREATE");
	TH1D ** dat_xb_lowQ2 	= new TH1D*[lowQ2_bins];
	TH1D ** dat_xb_highQ2 	= new TH1D*[highQ2_bins];
	TH1D ** bac_xb_lowQ2 	= new TH1D*[lowQ2_bins];
	TH1D ** bac_xb_highQ2 	= new TH1D*[highQ2_bins];
	for( int i = 0 ; i < lowQ2_bins ; ++i){
		dat_xb_lowQ2[i] = new TH1D(Form("dat_xb_lowQ2_%i",i),Form("dat_xb_lowQ2_%i",i),bins_Xb,Xb_min,Xb_max);
		bac_xb_lowQ2[i] = new TH1D(Form("bac_xb_lowQ2_%i",i),Form("bac_xb_lowQ2_%i",i),bins_Xb,Xb_min,Xb_max);
	}
	for( int i = 0 ; i < highQ2_bins; ++i){
		dat_xb_highQ2[i] = new TH1D(Form("dat_xb_highQ2_%i",i),Form("dat_xb_highQ2_%i",i),bins_Xb,Xb_min,Xb_max);
		bac_xb_highQ2[i] = new TH1D(Form("bac_xb_highQ2_%i",i),Form("bac_xb_highQ2_%i",i),bins_Xb,Xb_min,Xb_max);
	}

	// Loop over the data file:
	for( int event = 0 ; event < inTree_Dat->GetEntries() ; ++event ){
		dat_eHit	->Clear();
		dat_tagged	->Clear();
		dat_nleadindex	= 0.;

		inTree_Dat->GetEntry(event);
		
		// Get the correct tag hit
		taghit * this_tag	= (taghit*)  dat_tagged->At(dat_nleadindex);
		
		double As	= this_tag->getAs();
		double Pt	= this_tag->getPt().Mag();
		double Q2	= dat_eHit->getQ2();
		double Xb	= dat_eHit->getXb();

		if( Q2 >= Q2Bins[0] && Q2 < Q2Bins[1] ){
			for( int this_bin = 0 ; this_bin < lowQ2_bins ; ++this_bin ){
				double As_low  = lowQ2_AsPtBins[this_bin][0];
				double As_high = lowQ2_AsPtBins[this_bin][1];
				double Pt_low  = lowQ2_AsPtBins[this_bin][2];
				double Pt_high = lowQ2_AsPtBins[this_bin][3];
				
				if( As < As_low 	) continue;
				if( As >= As_high	) continue;
				if( Pt < Pt_low		) continue;
				if( Pt >= Pt_high	) continue;

				dat_xb_lowQ2[this_bin]->Fill( Xb );
			}
		}
		else if(Q2 >= Q2Bins[1] && Q2 < Q2Bins[2] ){
			for( int this_bin = 0 ; this_bin < highQ2_bins ; ++this_bin ){
				double As_low  = highQ2_AsPtBins[this_bin][0];
				double As_high = highQ2_AsPtBins[this_bin][1];
				double Pt_low  = highQ2_AsPtBins[this_bin][2];
				double Pt_high = highQ2_AsPtBins[this_bin][3];
				
				if( As < As_low 	) continue;
				if( As >= As_high	) continue;
				if( Pt < Pt_low		) continue;
				if( Pt >= Pt_high	) continue;

				dat_xb_highQ2[this_bin]->Fill( Xb );
			}

		}
		else continue;
	}

	// Loop over the background file:
	for( int event = 0 ; event < inTree_Bac->GetEntries() ; ++event ){
		bac_eHit	->Clear();
		bac_tagged	->Clear();
		bac_nleadindex	= 0.;

		inTree_Bac->GetEntry(event);
		
		// Get the correct tag hit
		taghit * this_tag	= (taghit*)  bac_tagged->At(bac_nleadindex);
		
		double As	= this_tag->getAs();
		double Pt	= this_tag->getPt().Mag();
		double Q2	= bac_eHit->getQ2();
		double Xb	= bac_eHit->getXb();

		if( Q2 >= Q2Bins[0] && Q2 < Q2Bins[1] ){
			for( int this_bin = 0 ; this_bin < lowQ2_bins ; ++this_bin ){
				double As_low  = lowQ2_AsPtBins[this_bin][0];
				double As_high = lowQ2_AsPtBins[this_bin][1];
				double Pt_low  = lowQ2_AsPtBins[this_bin][2];
				double Pt_high = lowQ2_AsPtBins[this_bin][3];
				
				if( As < As_low 	) continue;
				if( As >= As_high	) continue;
				if( Pt < Pt_low		) continue;
				if( Pt >= Pt_high	) continue;

				bac_xb_lowQ2[this_bin]->Fill( Xb , bac_weight );
			}
		}
		else if(Q2 >= Q2Bins[1] && Q2 < Q2Bins[2] ){
			for( int this_bin = 0 ; this_bin < highQ2_bins ; ++this_bin ){
				double As_low  = highQ2_AsPtBins[this_bin][0];
				double As_high = highQ2_AsPtBins[this_bin][1];
				double Pt_low  = highQ2_AsPtBins[this_bin][2];
				double Pt_high = highQ2_AsPtBins[this_bin][3];
				
				if( As < As_low 	) continue;
				if( As >= As_high	) continue;
				if( Pt < Pt_low		) continue;
				if( Pt >= Pt_high	) continue;

				bac_xb_highQ2[this_bin]->Fill( Xb , bac_weight );
			}

		}
		else continue;
	}

	// Print the background subtracted distributions to a pdf file:
	TCanvas * c_lowQ2 	= new TCanvas("c_lowQ2",	"",800,600);
	TCanvas * c_highQ2 	= new TCanvas("c_highQ2",	"",800,600);
	c_lowQ2 ->Divide(5,5);
	c_highQ2->Divide(5,5);
	outFile->cd();

	for( int i = 0 ; i < lowQ2_bins ; ++i){
		dat_xb_lowQ2[i]->Write();
		bac_xb_lowQ2[i]->Write();

		if( dat_xb_lowQ2[i]->Integral() < 1 || bac_xb_lowQ2[i]->Integral() < 1 ) continue;
		dat_xb_lowQ2[i]->Add( bac_xb_lowQ2[i] , -1 );
		c_lowQ2->cd(i+1);

		dat_xb_lowQ2[i]->SetLineWidth(3);
		dat_xb_lowQ2[i]->Draw();
	}
	c_lowQ2->Print("rawyield_lowQ2.pdf");

	for( int i = 0 ; i < highQ2_bins ; ++i){
		dat_xb_highQ2[i]->Write();
		bac_xb_highQ2[i]->Write();

		if( dat_xb_highQ2[i]->Integral() < 1 || bac_xb_highQ2[i]->Integral() < 1 ) continue;
		dat_xb_highQ2[i]->Add( bac_xb_highQ2[i] , -1 );
		c_highQ2->cd(i+1);

		dat_xb_highQ2[i]->SetLineWidth(3);
		dat_xb_highQ2[i]->Draw();
	}
	c_highQ2->Print("rawyield_highQ2.pdf");

	outFile->Close();


	inFile_Dat.Close();
	inFile_Bac.Close();
	return 1;
}
