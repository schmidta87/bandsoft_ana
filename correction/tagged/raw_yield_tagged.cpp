#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TCanvas.h"

#include "genpart.h"
#include "clashit.h"
#include "bandhit.h"
#include "taghit.h"

#include "bin_edges.h"

using std::cerr;
using std::cout;

void fillHist( double Q2, double Pt, double Xb, double As, TH1D**** hist , double weight );
void setError( double * err , TH1D * spb , TH1D * bac );

int main( int argc, char** argv){

	if( argc != 4 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Dat File] [Bac File] [Sim File]\n";
		return -1;
	}

	// Load files
	TFile inFile_Dat(argv[1]);
	TFile inFile_Bac(argv[2]);
	TFile inFile_Sim(argv[3]);

	// Load TTrees from files
	TTree * inTree_Dat = (TTree*) inFile_Dat.Get("tagged");
	TTree * inTree_Bac = (TTree*) inFile_Bac.Get("tagged");
	TTree * inTree_Sim = (TTree*) inFile_Sim.Get("tagged");

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

	// Set the input branches for simulation
	clashit * sim_eHit 		= new clashit;
	TClonesArray* sim_tagged 	= new TClonesArray("taghit");
	int sim_nleadindex		= 0;
	inTree_Sim->SetBranchAddress("tag_smeared"		,&sim_tagged		);
	inTree_Sim->SetBranchAddress("eHit_smeared"		,&sim_eHit		);
	inTree_Sim->SetBranchAddress("nleadindex"	,&sim_nleadindex	);

	// Define the histograms for yield to fill:
	TFile * outFile = new TFile("yield_tagged.root","RECREATE");
	TH1D **** h4_dat_as 	= new TH1D***[bins_Q2];
	TH1D **** h4_bac_as	= new TH1D***[bins_Q2];
	TH1D **** h4_sim_as	= new TH1D***[bins_Q2];

	for( int i = 0 ; i < bins_Q2 ; ++i ){ // bins in Q2
		h4_dat_as[i]	= new TH1D**[bins_Pt];
		h4_bac_as[i]	= new TH1D**[bins_Pt];
		h4_sim_as[i]	= new TH1D**[bins_Pt];
		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			h4_dat_as[i][j] 	= new TH1D*[bins_Xb];
			h4_bac_as[i][j] 	= new TH1D*[bins_Xb];
			h4_sim_as[i][j] 	= new TH1D*[bins_Xb];
			for( int k = 0 ; k < bins_Xb ; ++k ){ // bins in Xb
				h4_dat_as[i][j][k] = new TH1D(Form("dat_yield_as_Q2_%i_Pt_%i_Xb_%i",i,j,k),Form("dat_yield_as_Q2_%i_Pt_%i_Xb_%i",i,j,k),bins_As,As_min,As_max);  // histogram in As
				h4_bac_as[i][j][k] = new TH1D(Form("bac_as_Q2_%i_Pt_%i_Xb_%i",i,j,k),Form("bac_as_Q2_%i_Pt_%i_Xb_%i",i,j,k),bins_As,As_min,As_max);  // histogram in As
				h4_sim_as[i][j][k] = new TH1D(Form("sim_yield_as_Q2_%i_Pt_%i_Xb_%i",i,j,k),Form("sim_yield_as_Q2_%i_Pt_%i_Xb_%i",i,j,k),bins_As,As_min,As_max);  // histogram in As
			}
		}
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

		fillHist( Q2, Pt, Xb, As, h4_dat_as , 1. );
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

		fillHist( Q2, Pt, Xb, As, h4_bac_as , bac_weight );
	}

	// Loop over the simulation file:
	for( int event = 0 ; event < inTree_Sim->GetEntries() ; ++event ){
		sim_eHit	->Clear();
		sim_tagged	->Clear();
		sim_nleadindex	= 0.;

		inTree_Sim->GetEntry(event);
		
		// Get the correct tag hit
		taghit * this_tag	= (taghit*)  sim_tagged->At(sim_nleadindex);
		
		double As	= this_tag->getAs();
		double Pt	= this_tag->getPt().Mag();
		double Q2	= sim_eHit->getQ2();
		double Xb	= sim_eHit->getXb();

		fillHist( Q2, Pt, Xb, As, h4_sim_as , 1. );
	}

	TCanvas ** c_Q2 = new TCanvas*[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		c_Q2[i] = new TCanvas(Form("c_Q2_%i",i),"",1200,1200);
		
		c_Q2[i]->Divide( bins_Pt, bins_Xb );
		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			for( int k = 0 ; k < bins_Xb ; ++k ){ // loop over Xb bins

				if( h4_dat_as[i][j][k]->Integral() < 1 || h4_bac_as[i][j][k]->Integral() < 1 ) continue;
				double * errors = new double[bins_As];
				setError( errors , h4_dat_as[i][j][k] , h4_bac_as[i][j][k] );
				h4_dat_as[i][j][k]->Add( h4_bac_as[i][j][k] , -1 );

				for( int bin = 1 ; bin < h4_dat_as[i][j][k]->GetXaxis()->GetNbins(); ++bin ){
					h4_dat_as[i][j][k]->SetBinError( bin , errors[bin-1] );
				}

				c_Q2[i]->cd( (k*2)+1 + j );

				h4_dat_as[i][j][k]->SetLineWidth(3);
				h4_dat_as[i][j][k]->SetMarkerColor(4);
				h4_dat_as[i][j][k]->SetMinimum(0);
				h4_dat_as[i][j][k]->Draw("*P");

				h4_sim_as[i][j][k]->SetLineWidth(1);
				h4_sim_as[i][j][k]->SetLineColor(2);
				h4_sim_as[i][j][k]->SetMinimum(0);

				h4_dat_as[i][j][k]->Write();
				h4_sim_as[i][j][k]->Write();

				delete[] errors;

			}
		}
		c_Q2[i]->Print(Form("yield_tagged_Q2_%i.pdf",i));
	}
	outFile->Close();

	inFile_Dat.Close();
	inFile_Bac.Close();
	inFile_Sim.Close();
	return 1;
}

void fillHist( double Q2, double Pt, double Xb, double As, TH1D**** hist , double weight ){
	

	// If it's larger than max Q2 or smaller than min Q2, return
	if( Q2 > Q2Bins[bins_Q2] 	) return;
	if( Q2 < Q2Bins[0] 		) return;

	// Need to figure out which bin of Q2, Pt, Xb (then fill as As):
	int this_bin_q2 = -1;
	int this_bin_pt = -1;
	int this_bin_xb = -1;
	
	for( int q2_bin = bins_Q2-1 ; q2_bin >= 0; --q2_bin ){
		if( Q2 > Q2Bins[q2_bin] ) this_bin_q2 = q2_bin;
		if( this_bin_q2 != -1 ) break;
	}

	if( Pt < Pt_min			) return;
	if( Pt > Pt_max			) return;
	if( Xb < Xb_min[this_bin_q2]			) return;
	if( Xb > Xb_min[this_bin_q2]+Xb_step*6		) return;

	this_bin_pt = (int) ((Pt - Pt_min)/Pt_step);
	this_bin_xb = (int) ((Xb - Xb_min[this_bin_q2])/Xb_step);

	// Safety clauses
	if( this_bin_q2 == -1 ){ cerr << "how\n"; return; }
	if( this_bin_pt == -1 ){ cerr << "how\n"; return; }
	if( this_bin_xb == -1 ){ cerr << "how\n"; return; }
	if( this_bin_q2 > bins_Q2-1 ){ cerr << "how\n"; return; }
	if( this_bin_pt > bins_Pt-1 ){ cerr << "how\n"; return; }
	if( this_bin_xb > bins_Xb-1 ){ cerr << "how\n"; return; }


	hist[this_bin_q2][this_bin_pt][this_bin_xb]->Fill( As , weight );

}

void setError( double * err , TH1D * spb , TH1D * bac ){
	
	for( int bin = 1 ; bin < spb->GetXaxis()->GetNbins(); ++bin ){
		double s = spb->GetBinContent(bin);
		double b = bac->GetBinContent(bin);

		double e = sqrt( s + b );
		err[bin-1] = e;
	}


	return;
}
