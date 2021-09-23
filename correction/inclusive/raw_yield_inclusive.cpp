#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TCanvas.h"

#include "clashit.h"

#include "bin_edges.h"

using std::cerr;
using std::cout;

void fillHist( double Q2, double Xb, TH1D** hist, double weight );

int main( int argc, char** argv){

	if( argc != 3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Dat File] [Sim File]\n";
		return -1;
	}

	// Load files
	TFile inFile_Dat(argv[1]);
	TFile inFile_Sim(argv[2]);

	// Load TTrees from files
	TTree * inTree_Dat = (TTree*) inFile_Dat.Get("electrons");
	TTree * inTree_Sim = (TTree*) inFile_Sim.Get("electrons");

	// Set the input branches for data
	clashit * dat_eHit 		= new clashit;
	inTree_Dat->SetBranchAddress("eHit"		,&dat_eHit		);

	// Set the input branches for simulation
	clashit * sim_eHit 		= new clashit;
	inTree_Sim->SetBranchAddress("eHit_smeared"		,&sim_eHit		);

	// Define the histograms for yield to fill:
	TFile * outFile = new TFile("yield_inclusive.root","RECREATE");
	TH1D ** h2_dat_xb 	= new TH1D*[bins_Q2];
	TH1D ** h2_sim_xb	= new TH1D*[bins_Q2];

	for( int i = 0 ; i < bins_Q2 ; ++i ){ // bins in Q2
		double minxb = Xb_min[i];
		double maxxb = minxb + bins_Xb*Xb_step;
		h2_dat_xb[i] = new TH1D(Form("dat_yield_Xb_Q2_%i",i),Form("dat_yield_Xb_Q2_%i",i),bins_Xb,minxb,maxxb);  // histogram in As
		h2_sim_xb[i] = new TH1D(Form("sim_yield_Xb_Q2_%i",i),Form("sim_yield_Xb_Q2_%i",i),bins_Xb,minxb,maxxb);  // histogram in As
	}

	// Loop over the data file:
	for( int event = 0 ; event < inTree_Dat->GetEntries() ; ++event ){
		dat_eHit	->Clear();

		inTree_Dat->GetEntry(event);
		double Q2	= dat_eHit->getQ2();
		double Xb	= dat_eHit->getXb();

		fillHist( Q2, Xb, h2_dat_xb , 243984840./3303854. );
	}

	// Loop over the simulation file:
	for( int event = 0 ; event < inTree_Sim->GetEntries() ; ++event ){
		sim_eHit	->Clear();

		inTree_Sim->GetEntry(event);
		double Q2	= sim_eHit->getQ2();
		double Xb	= sim_eHit->getXb();

		fillHist( Q2, Xb, h2_sim_xb , 1. );
	}

	TCanvas ** c_Q2 = new TCanvas*[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		c_Q2[i] = new TCanvas(Form("c_Q2_%i",i),"",1200,1200);
		
		if( h2_dat_xb[i]->Integral() < 1 || h2_sim_xb[i]->Integral() < 1 ) continue;
		
		h2_dat_xb[i]->SetLineWidth(3);
		h2_dat_xb[i]->SetMarkerColor(4);
		h2_dat_xb[i]->SetMinimum(0);
		h2_dat_xb[i]->Draw("*P");

		h2_sim_xb[i]->SetLineWidth(1);
		h2_sim_xb[i]->SetLineColor(2);
		h2_sim_xb[i]->SetMinimum(0);

		h2_dat_xb[i]->Write();
		h2_sim_xb[i]->Write();
		c_Q2[i]->Print(Form("yield_inclusive_Q2_%i.pdf",i));

	}
	outFile->Close();

	inFile_Dat.Close();
	inFile_Sim.Close();
	return 1;
}

void fillHist( double Q2, double Xb, TH1D** hist, double weight ){

	// If it's larger than max Q2 or smaller than min Q2, return
	if( Q2 > Q2Bins[bins_Q2] 	) return;
	if( Q2 < Q2Bins[0] 		) return;

	// Need to figure out which bin of Q2, Pt, Xb (then fill as As):
	int this_bin_q2 = -1;
	
	for( int q2_bin = bins_Q2-1 ; q2_bin >= 0; --q2_bin ){
		if( Q2 > Q2Bins[q2_bin] ) this_bin_q2 = q2_bin;
		if( this_bin_q2 != -1 ) break;
	}

	if( Xb < Xb_min[this_bin_q2]			) return;
	if( Xb > Xb_min[this_bin_q2]+Xb_step*6		) return;

	// Safety clauses
	if( this_bin_q2 == -1 ){ cerr << "how\n"; return; }
	if( this_bin_q2 > bins_Q2-1 ){ cerr << "how\n"; return; }

	hist[this_bin_q2]->Fill( Xb , weight );

}

