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
#include "bandhit.h"
#include "taghit.h"

#include "bin_edges.h"

using std::cerr;
using std::cout;

void fillHist( double Q2, double Pt, double Xp, double As, TH1D**** hist , double weight );
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
	TH1D **** h4_dat_xp 	= new TH1D***[bins_Q2];
	TH1D **** h4_bac_xp	= new TH1D***[bins_Q2];
	TH1D **** h4_sim_xp	= new TH1D***[bins_Q2];

	TH1D **** h4_dat_sim_xp 	= new TH1D***[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // bins in Q2
		h4_dat_xp[i]	= new TH1D**[bins_Pt];
		h4_bac_xp[i]	= new TH1D**[bins_Pt];
		h4_sim_xp[i]	= new TH1D**[bins_Pt];
		h4_dat_sim_xp[i]	= new TH1D**[bins_Pt];
		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			h4_dat_xp[i][j] 	= new TH1D*[bins_As];
			h4_bac_xp[i][j] 	= new TH1D*[bins_As];
			h4_sim_xp[i][j] 	= new TH1D*[bins_As];
			h4_dat_sim_xp[i][j] 	= new TH1D*[bins_As];
			for( int k = 0 ; k < bins_As ; ++k ){ // bins in As


				h4_dat_xp[i][j][k] = new TH1D(Form("dat_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("dat_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in Xp
				h4_bac_xp[i][j][k] = new TH1D(Form("bac_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("bac_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in Xp
				h4_sim_xp[i][j][k] = new TH1D(Form("sim_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("sim_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in Xp
				h4_dat_sim_xp[i][j][k] = new TH1D(Form("dat_sim_ratio_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("dat_sim_ratio_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in Xp
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
		double Xp	= this_tag->getXp_WP();
		double Q2	= dat_eHit->getQ2();
		double Pn	= this_tag->getMomentumN().Mag();
		if( Pn < 0.3 ) continue;

		fillHist( Q2, Pt, Xp, As, h4_dat_xp , 1. );
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
		double Xp	= this_tag->getXp_WP();
		double Q2	= bac_eHit->getQ2();
		double Pn	= this_tag->getMomentumN().Mag();
		if( Pn < 0.3 ) continue;

		fillHist( Q2, Pt, Xp, As, h4_bac_xp , bac_weight );
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
		double Xp	= this_tag->getXp_WP();
		double Q2	= sim_eHit->getQ2();
		double Pn	= this_tag->getMomentumN().Mag();
		if( Pn < 0.3 ) continue;

		fillHist( Q2, Pt, Xp, As, h4_sim_xp , 1. );
	}

	TCanvas ** c_Q2 = new TCanvas*[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		c_Q2[i] = new TCanvas(Form("c_Q2_%i",i),"",1200,1200);
		
		c_Q2[i]->Divide( bins_Pt, bins_As );
		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			for( int k = 0 ; k < bins_As ; ++k ){ // loop over As bins

				c_Q2[i]->cd( (k*2)+1 + j );

				if( h4_dat_xp[i][j][k]->Integral() < 1 || h4_bac_xp[i][j][k]->Integral() < 1 || h4_sim_xp[i][j][k]->Integral() < 1 ) continue;
				cout << i << " " << j << " " << k << "\n";
				double * errors = new double[bins_Xp];
				setError( errors , h4_dat_xp[i][j][k] , h4_bac_xp[i][j][k] );
				h4_dat_xp[i][j][k]->Add( h4_bac_xp[i][j][k] , -1. );

				int ref_bin = h4_dat_xp[i][j][k]->FindBin( Xp_Ref[i] );
				double binref = h4_dat_xp[i][j][k]->GetBinContent( ref_bin ) /
					h4_sim_xp[i][j][k]->GetBinContent( ref_bin );
				double ref_s = h4_sim_xp[i][j][k]->GetBinContent(ref_bin);
				double ref_d = h4_dat_xp[i][j][k]->GetBinContent(ref_bin);
				double ref_f = ref_d / ref_s;
				double ref_err = sqrt( ref_d/(ref_s*ref_s) + ref_d*ref_d/pow(ref_s,3) );

				cout << "\t" << ref_bin << " " << binref << "\n";
				if( binref < 0 || !isfinite(binref) ) continue;


				for( int bin = 1 ; bin < h4_dat_xp[i][j][k]->GetXaxis()->GetNbins(); ++bin ){
					h4_dat_xp[i][j][k]->SetBinError( bin , errors[bin-1] );
	
					double bincontent = h4_dat_xp[i][j][k]->GetBinContent(bin) / 
						h4_sim_xp[i][j][k]->GetBinContent(bin) / binref;			
					cout << "\t\t" << bin << " " << bincontent << "\n";
					if( !isfinite(bincontent) || bincontent < 0 ){
						h4_dat_sim_xp[i][j][k]->SetBinContent(bin, 0);
						h4_dat_sim_xp[i][j][k]->SetBinError(bin , 0);
						continue;
					};
					h4_dat_sim_xp[i][j][k]->SetBinContent(bin, bincontent);
					
					double s = h4_sim_xp[i][j][k]->GetBinContent(bin);
					double d = h4_dat_xp[i][j][k]->GetBinContent(bin);
					double f = d/s;
					double err1 = sqrt( d/(s*s) + d*d/(s*s*s) );

					double ratio_err = sqrt( 
						pow(1./ref_f * err1,2) + pow(f * ref_err,2)/pow(ref_f,4) );

					cout << "\t\t\t" << ref_err << " " << err1 << " " << ratio_err << "\n\n";
					h4_dat_sim_xp[i][j][k]->SetBinError( bin , ratio_err );


				}
				
				if( h4_dat_sim_xp[i][j][k]->Integral() < 1 ) continue;

				TString title = Form("%.1f < Q^{2} < %.1f, %.1f < p_{T} < %.1f, %.2f < a_{S} < %.2f",
						Q2Bins[i],Q2Bins[i+1],
						Pt_min + Pt_step*j,Pt_min + Pt_step*(j+1),
						As_min + As_step*k,As_min + As_step*(k+1)	);
				h4_dat_sim_xp[i][j][k]->SetTitle(title);

				h4_dat_sim_xp[i][j][k]->SetLineWidth(3);
				h4_dat_sim_xp[i][j][k]->SetMarkerColor(4);
				h4_dat_sim_xp[i][j][k]->SetMinimum(0);
				h4_dat_sim_xp[i][j][k]->SetMaximum(2.5);
				h4_dat_sim_xp[i][j][k]->GetXaxis()->SetTitle("x'");
				h4_dat_sim_xp[i][j][k]->Draw("*P");
				for( int bin = 1 ; bin < h4_dat_sim_xp[i][j][k]->GetXaxis()->GetNbins(); ++bin ){
					cout << h4_dat_sim_xp[i][j][k]->GetBinContent(bin) << " " 
						<< h4_dat_sim_xp[i][j][k]->GetBinError(bin) << "\n";
				}

				TLine * hline = new TLine(Xp_min,1,Xp_max,1);
				hline->SetLineWidth(2);
				hline->SetLineColor(1);
				hline->SetLineStyle(2);
				hline->Draw("SAME");
				TLine * vline = new TLine(Xp_Ref[i],0,Xp_Ref[i],2.5);
				vline->SetLineWidth(2);
				vline->SetLineColor(1);
				vline->SetLineStyle(2);
				vline->Draw("SAME");

				//h4_sim_xp[i][j][k]->SetLineWidth(1);
				//h4_sim_xp[i][j][k]->SetLineColor(2);
				//h4_sim_xp[i][j][k]->SetMinimum(0);

				h4_dat_xp[i][j][k]->Write();
				h4_sim_xp[i][j][k]->Write();
				h4_dat_sim_xp[i][j][k]->Write();

				delete[] errors;

			}
		}
		c_Q2[i]->Print(Form("tagged_yield_Q2_%i.pdf",i));
	}
	outFile->Close();

	inFile_Dat.Close();
	inFile_Bac.Close();
	inFile_Sim.Close();
	return 1;
}

void fillHist( double Q2, double Pt, double Xp, double As, TH1D**** hist , double weight ){
	

	// If it's larger than max Q2 or smaller than min Q2, return
	if( Q2 > Q2Bins[bins_Q2] 	) return;
	if( Q2 < Q2Bins[0] 		) return;

	// Need to figure out which bin of Q2, Pt, As (then fill as Xp):
	int this_bin_q2 = -1;
	int this_bin_pt = -1;
	int this_bin_as = -1;
	
	for( int q2_bin = bins_Q2-1 ; q2_bin >= 0; --q2_bin ){
		if( Q2 > Q2Bins[q2_bin] ) this_bin_q2 = q2_bin;
		if( this_bin_q2 != -1 ) break;
	}

	if( Pt < Pt_min			) return;
	if( Pt > Pt_max			) return;
	if( As < As_min		 	) return;
	if( As > As_max		 	) return;

	this_bin_pt = (int) ((Pt - Pt_min)/Pt_step);
	this_bin_as = (int) ((As - As_min)/As_step);

	// Safety clauses
	if( this_bin_q2 == -1 ){ cerr << "how\n"; return; }
	if( this_bin_pt == -1 ){ cerr << "how\n"; return; }
	if( this_bin_as == -1 ){ cerr << "how\n"; return; }
	if( this_bin_q2 > bins_Q2-1 ){ cerr << "how\n"; return; }
	if( this_bin_pt > bins_Pt-1 ){ cerr << "how\n"; return; }
	if( this_bin_as > bins_As-1 ){ cerr << "how\n"; return; }

	hist[this_bin_q2][this_bin_pt][this_bin_as]->Fill( Xp , weight );

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
