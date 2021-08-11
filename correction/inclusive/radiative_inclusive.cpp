#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"

#include "genpart.h"
#include "clashit.h"

#include "bin_edges.h"

using std::cerr;
using std::cout;
void fillHist( double Q2, double Xb, double weight, TH1D** hist , TH1D** counts );
void setError( double * err , TH1D * ravg , TH1D * cnts );

int main( int argc, char** argv){

	if( argc != 2 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Sim File]\n";
		return -1;
	}

	TFile inFile(argv[1]);
	TTree * inTree = (TTree*) inFile.Get("electrons");

	TClonesArray* gen_particles 	= new TClonesArray("genpart");
	int 	genMult 	= 0;
	double	Ebeam		= 0;
	double	rad_weight	= 0;
	inTree->SetBranchAddress("mcParts"		, &gen_particles		);
	inTree->SetBranchAddress("genMult"		, &genMult			);
	inTree->SetBranchAddress("Ebeam"		, &Ebeam			);
	inTree->SetBranchAddress("weight"		, &rad_weight			);

	TFile * outFile = new TFile("radcorrection_inclusive.root","RECREATE");
	TH1D ** h2_gen_xb 	= new TH1D*[bins_Q2];
	TH1D ** h2_cnt_xb	= new TH1D*[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // bins in Q2
		double minxb = Xb_min[i];
		double maxxb = minxb + bins_Xb*Xb_step;
		h2_gen_xb[i] = new TH1D(Form("radcorr_Xb_Q2_%i",i),Form("radcorr_Xb_Q2_%i",i),bins_Xb,minxb,maxxb); 
		h2_cnt_xb[i] = new TH1D(Form("cnt_Xb_Q2_%i",i),Form("cnt_Xb_Q2_%i",i),bins_Xb,minxb,maxxb); 
	}


	for( int event = 0 ; event < inTree->GetEntries() ; ++event ){
		gen_particles	->Clear();
		genMult 	= 0;
		Ebeam 		= 0;
		rad_weight	= 0;
		//if(event > 1000) break;
		inTree->GetEntry(event);

		// Form As, ThetaNQ, W, X with generated:
		genpart * gen_electron 	= (genpart*) gen_particles->At(0);
		genpart * gen_neutron 	= (genpart*) gen_particles->At(1);
		double	gen_Q2		= gen_electron->getQ2();
		double 	gen_Xb		= gen_electron->getXb();

		// sort the generated values:
		fillHist( gen_Q2, gen_Xb, rad_weight, h2_gen_xb , h2_cnt_xb );

	}

	// Print the background subtracted distributions to a pdf file:
	TCanvas ** c_Q2 = new TCanvas*[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		c_Q2[i] = new TCanvas(Form("c_Q2_%i",i),"",1200,1200);
		

		if( h2_gen_xb[i]->Integral() < 1 || h2_cnt_xb[i]->Integral() < 1 ) continue;
		double * errors = new double[bins_Xb];
		setError( errors , h2_gen_xb[i] , h2_cnt_xb[i] );
		h2_gen_xb[i]->Divide( h2_cnt_xb[i] );

		for( int bin = 1 ; bin < h2_gen_xb[i]->GetXaxis()->GetNbins(); ++bin ){
			h2_gen_xb[i]->SetBinError( bin , errors[bin-1] );
		}


		h2_gen_xb[i]->SetLineWidth(3);
		h2_gen_xb[i]->SetMarkerColor(4);
		h2_gen_xb[i]->SetMaximum(2);
		h2_gen_xb[i]->SetMinimum(0);
		h2_gen_xb[i]->Draw("*P");
		TLine * hline = new TLine(Xb_min[i],1,Xb_min[i]+bins_Xb*Xb_step,1);
		hline->SetLineWidth(2);
		hline->SetLineColor(1);
		hline->SetLineStyle(2);
		hline->Draw("SAME");

		h2_gen_xb[i]->Write();
		//h2_cnt_xb[i]->Write();
		delete[] errors;

		c_Q2[i]->Print(Form("radiative_inclusive_Q2_%i.pdf",i));
	}
	outFile->Close();

	inFile.Close();
	return 1;
}

void fillHist( double Q2, double Xb, double weight, TH1D** hist , TH1D** counts ){
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

	// For radiative weight, we cummulatively add the content in the same bin:
	int fill_bin = hist[this_bin_q2]->GetXaxis()->FindBin( Xb );
	double hist_curr_content = hist[this_bin_q2]->GetBinContent( fill_bin );
	double counts_curr_content = counts[this_bin_q2]->GetBinContent( fill_bin );
	
	hist[this_bin_q2]->SetBinContent( fill_bin, hist_curr_content + weight );
	counts[this_bin_q2]->SetBinContent( fill_bin , counts_curr_content + 1 );

}

void setError( double * err , TH1D * ravg , TH1D * cnts ){
	
	for( int bin = 1 ; bin < ravg->GetXaxis()->GetNbins(); ++bin ){

		double R = ravg->GetBinContent(bin);
		double N = cnts->GetBinContent(bin);

		double e = sqrt( R/pow(N,3./2) );
		cout << R << " " << N << " " << e << "\n";
		if( e!=e ) e = 0;
		err[bin-1] = e;
	}


	return;
}
