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
using std::isfinite;
using std::cout;

void fillHist( double Q2, double Pt, double Xp, double As, TH1D**** hist , double weight );
void setError( double * err , TH1D * rec , TH1D * gen );
void background_subtraction(TH1D* dat, TH1D* bac , double Cscale, double NB_sim , double Sigma_Cscale, double Sigma_NB_sim );
void simulation_weighting(TH1D* sim, double Ndata, double Nsim );
bool bad_bar( bandhit * this_n );

int main( int argc, char** argv){

	if( argc != 3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Sim File 1] [Sim File 2]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] << "\n";

	// Load files
	TFile inFile_Sim1(argv[1]);
	TFile inFile_Sim2(argv[2]);

	// Load TTrees from files
	TTree * inTree_Sim1 = (TTree*) inFile_Sim1.Get("tagged");
	TTree * inTree_Sim2 = (TTree*) inFile_Sim2.Get("tagged");

	// Set the input branches for simulation 1
	clashit * sim1_eHit 		= new clashit;
	TClonesArray* sim1_tagged 	= new TClonesArray("taghit");
	TClonesArray* sim1_bandhi 	= new TClonesArray("bandhit");
	int sim1_nleadindex		= 0;
	inTree_Sim1->SetBranchAddress("tag_smeared"		,&sim1_tagged		);
	inTree_Sim1->SetBranchAddress("nHits"			,&sim1_bandhi		);
	inTree_Sim1->SetBranchAddress("eHit_smeared"		,&sim1_eHit		);
	inTree_Sim1->SetBranchAddress("nleadindex"		,&sim1_nleadindex	);

	// Set the input branches for simulation 2
	clashit * sim2_eHit 		= new clashit;
	TClonesArray* sim2_tagged 	= new TClonesArray("taghit");
	TClonesArray* sim2_bandhi 	= new TClonesArray("bandhit");
	int sim2_nleadindex		= 0;
	inTree_Sim2->SetBranchAddress("tag_smeared"		,&sim2_tagged		);
	inTree_Sim2->SetBranchAddress("nHits"			,&sim2_bandhi		);
	inTree_Sim2->SetBranchAddress("eHit_smeared"		,&sim2_eHit		);
	inTree_Sim2->SetBranchAddress("nleadindex"		,&sim2_nleadindex	);

	// Define the histograms for yield to fill:
	TH1D **** h4_sim1_xp	= new TH1D***[bins_Q2];
	TH1D **** h4_sim2_xp	= new TH1D***[bins_Q2];

	TH1D **** h4_sim1_sim2_xp 	= new TH1D***[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // bins in Q2
		h4_sim1_xp[i]	= new TH1D**[bins_Pt];
		h4_sim2_xp[i]	= new TH1D**[bins_Pt];
		h4_sim1_sim2_xp[i]	= new TH1D**[bins_Pt];
		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			h4_sim1_xp[i][j] 	= new TH1D*[bins_As];
			h4_sim2_xp[i][j] 	= new TH1D*[bins_As];
			h4_sim1_sim2_xp[i][j] 	= new TH1D*[bins_As];
			for( int k = 0 ; k < bins_As ; ++k ){ // bins in As

				h4_sim1_xp[i][j][k] = new TH1D(Form("sim1_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("sim1_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in Xp
				h4_sim2_xp[i][j][k] = new TH1D(Form("sim2_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("sim2_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in Xp
				h4_sim1_sim2_xp[i][j][k] = new TH1D(Form("sim1_sim2_ratio_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("sim1_sim2_ratio_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in Xp
			}
		}
	}


	// Loop over the 1st simulation file:
	for( int event = 0 ; event < inTree_Sim1->GetEntries() ; ++event ){
		sim1_eHit	->Clear();
		sim1_tagged	->Clear();
		sim1_bandhi	->Clear();
		sim1_nleadindex	= 0.;

		inTree_Sim1->GetEntry(event);
		
		// Get the correct tag hit
		taghit * this_tag	= (taghit*)  sim1_tagged->At(sim1_nleadindex);
		bandhit * lead_n	= (bandhit*) sim1_bandhi->At(sim1_nleadindex);
		if( lead_n->getEdep() < 10 ) continue;
		
		double As	= this_tag->getAs();
		double Pt	= this_tag->getPt().Mag();
		double Xp	= this_tag->getXp_WP();
		double Q2	= sim1_eHit->getQ2();
		double Pn	= this_tag->getMomentumN().Mag();
		if( Pn < 0.25 ) continue;

		if( bad_bar(lead_n) ) continue;

		//if( lead_n->getSector()==2 && (lead_n->getComponent() > 3 && lead_n->getComponent() < 8 ) 
		//		&& ( lead_n->getX()>90 || lead_n->getX() < -110 ) ) continue;
		//if( lead_n->getSector()==3 && (lead_n->getComponent() > 0 && lead_n->getComponent() < 3 ) 
		//		&& ( lead_n->getX()>80 || lead_n->getX() < 45 ) ) continue;

		fillHist( Q2, Pt, Xp, As, h4_sim1_xp , 1. );
	}
	// Loop over the 2nd simulation file:
	for( int event = 0 ; event < inTree_Sim2->GetEntries() ; ++event ){
		sim2_eHit	->Clear();
		sim2_tagged	->Clear();
		sim2_bandhi	->Clear();
		sim2_nleadindex	= 0.;

		inTree_Sim2->GetEntry(event);
		
		// Get the correct tag hit
		taghit * this_tag	= (taghit*)  sim2_tagged->At(sim2_nleadindex);
		bandhit * lead_n	= (bandhit*) sim2_bandhi->At(sim2_nleadindex);
		if( lead_n->getEdep() < 10 ) continue;
		
		double As	= this_tag->getAs();
		double Pt	= this_tag->getPt().Mag();
		double Xp	= this_tag->getXp_WP();
		double Q2	= sim2_eHit->getQ2();
		double Pn	= this_tag->getMomentumN().Mag();
		if( Pn < 0.25 ) continue;

		if( bad_bar(lead_n) ) continue;

		//if( lead_n->getSector()==2 && (lead_n->getComponent() > 3 && lead_n->getComponent() < 8 ) 
		//		&& ( lead_n->getX()>90 || lead_n->getX() < -110 ) ) continue;
		//if( lead_n->getSector()==3 && (lead_n->getComponent() > 0 && lead_n->getComponent() < 3 ) 
		//		&& ( lead_n->getX()>80 || lead_n->getX() < 45 ) ) continue;

		fillHist( Q2, Pt, Xp, As, h4_sim2_xp , 1. );
	}

	TCanvas ** c_Q2 = new TCanvas*[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		c_Q2[i] = new TCanvas(Form("c_Q2_%i",i),"",1200,1200);
		
		c_Q2[i]->Divide( bins_Pt, bins_As );
		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			for( int k = 0 ; k < bins_As ; ++k ){ // loop over As bins

				c_Q2[i]->cd( (k*2)+1 + j );

				if( h4_sim1_xp[i][j][k]->Integral() < 1 || h4_sim2_xp[i][j][k]->Integral() < 1 ) continue;

				double * errors = new double[bins_Xp];
				setError( errors , h4_sim1_xp[i][j][k] , h4_sim2_xp[i][j][k] );
				h4_sim1_xp[i][j][k]->Divide( h4_sim2_xp[i][j][k] );
				for( int bin = 1 ; bin < h4_sim1_xp[i][j][k]->GetXaxis()->GetNbins(); ++bin ){
					h4_sim1_xp[i][j][k]->SetBinError( bin , errors[bin-1] );
				}


				int ref_bin = h4_sim1_xp[i][j][k]->FindBin( Xp_Ref[k] );
				double binref = h4_sim1_xp[i][j][k]->GetBinContent( ref_bin );

				h4_sim1_xp[i][j][k]->Scale( 1./binref );
				

				TString title = Form("%.1f < Q^{2} < %.1f, %.1f < p_{T} < %.1f, %.2f < a_{S} < %.2f",
						Q2Bins[i],Q2Bins[i+1],
						Pt_min + Pt_step*j,Pt_min + Pt_step*(j+1),
						As_min + As_step*k,As_min + As_step*(k+1)	);
				h4_sim1_xp[i][j][k]->SetTitle(title);

				h4_sim1_xp[i][j][k]->SetLineWidth(3);
				h4_sim1_xp[i][j][k]->SetMarkerColor(4);
				h4_sim1_xp[i][j][k]->SetMinimum(0.5);
				h4_sim1_xp[i][j][k]->SetMaximum(1.3);
				h4_sim1_xp[i][j][k]->GetXaxis()->SetTitle("x'");
				h4_sim1_xp[i][j][k]->Draw("*P");

				TLine * hline = new TLine(Xp_min,1,Xp_max,1);
				hline->SetLineWidth(2);
				hline->SetLineColor(1);
				hline->SetLineStyle(2);
				hline->Draw("SAME");
				TLine * vline = new TLine(Xp_Ref[k],0.5,Xp_Ref[k],1.3);
				vline->SetLineWidth(2);
				vline->SetLineColor(1);
				vline->SetLineStyle(2);
				vline->Draw("SAME");



			}
		}
		c_Q2[i]->Print(Form("sim1_sim2_yield_Q2_%i.pdf",i));
	}

	inFile_Sim1.Close();
	inFile_Sim2.Close();
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

void setError( double * err , TH1D * rec , TH1D * gen ){
	
	for( int bin = 1 ; bin < rec->GetXaxis()->GetNbins(); ++bin ){

		double r = rec->GetBinContent(bin);
		double g = gen->GetBinContent(bin);

		double e = sqrt( r/(g*g) + r*r/(g*g*g) );
		cout << bin << " " << r << " " << g << " " << r/g << " " << e << "\n";
		if( e!=e ) e = 0;
		err[bin-1] = e;
		cout << err[bin-1] << "\n\n";
	}


	return;
}
void background_subtraction(TH1D* dat, TH1D* bac , double Cscale, double NB_sim , double Sigma_Cscale, double Sigma_NB_sim ){
	
	//dat->Add(bac,-1);
	///*
	for( int bin = 1 ; bin < dat->GetXaxis()->GetNbins(); bin++ ){
		double SpB_bin = dat->GetBinContent(bin);
		double NB_bin = bac->GetBinContent(bin); // before any re-weighting

		double Sigma_NB_bin = sqrt(NB_bin);
		//double Sigma_Cscale already defined above
		//double Sigma_NB_sim already defined above

		double i1 = (Cscale/NB_sim)*Sigma_NB_bin;
		double i2 = (NB_bin/NB_sim)*Sigma_Cscale;
		double i3 = (NB_bin*Cscale)/pow(NB_sim,2)*Sigma_NB_sim;
		double Sigma_B = sqrt(i1*i1 + i2*i2 + i3*i3);

		double Sigma = sqrt( SpB_bin + pow(Sigma_B,2) ); // total uncertainty per bin

		dat->SetBinContent(bin, SpB_bin - NB_bin*Cscale/NB_sim );
		dat->SetBinError(bin, Sigma );
	}
	//*/


	return;
}
void simulation_weighting(TH1D* sim, double Ndata, double Nsim ){
	for( int bin = 1 ; bin < sim->GetXaxis()->GetNbins(); bin++ ){
		double Nsim_bin = sim->GetBinContent(bin); // simulation stats in bin before re-weight

		double Sigma_Nsim_bin = sqrt(Nsim_bin);
		double Sigma_Ndata = sqrt(Ndata);
		double Sigma_Nsim = sqrt(Nsim);
		double i1 = (Ndata/Nsim)*Sigma_Nsim_bin;
		double i2 = (Nsim_bin/Nsim)*Sigma_Ndata;
		double i3 = (Nsim_bin*Ndata)/pow(Nsim,2)*Sigma_Nsim;

		double Sigma = sqrt(i1*i1 + i2*i2 + i3*i3);

		sim->SetBinContent(bin, Nsim_bin * Ndata/Nsim );
		sim->SetBinError(bin, Sigma );
	}


	return;
}

bool bad_bar( bandhit * this_n ){
	
	if(	this_n->getSector() == 1 && this_n->getComponent() == 1 				) return true;
	
	/*
	if(	this_n->getSector() == 3 && this_n->getLayer() == 4 && this_n->getComponent() == 2 	) return true;
	if(	this_n->getSector() == 4 && this_n->getLayer() == 3 && this_n->getComponent() == 2 	) return true;
	if(	this_n->getSector() == 4 && this_n->getLayer() == 2 && this_n->getComponent() == 5 	) return true;
	if(	this_n->getSector() == 1 && this_n->getLayer() == 1 && this_n->getComponent() == 1 	) return true;
	if(	this_n->getSector() == 1 && this_n->getLayer() == 2 && this_n->getComponent() == 1 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 2 && this_n->getComponent() == 5 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 2 && this_n->getComponent() == 7 	) return true;
	if(	this_n->getSector() == 4 && this_n->getLayer() == 2 && this_n->getComponent() == 4 	) return true;
	if(	this_n->getSector() == 4 && this_n->getLayer() == 2 && this_n->getComponent() == 6 	) return true;
	if(	this_n->getSector() == 5 && this_n->getLayer() == 2 && this_n->getComponent() == 1 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 3 && this_n->getComponent() == 6 	) return true;
	if(	this_n->getSector() == 4 && this_n->getLayer() == 3 && this_n->getComponent() == 1 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 4 && this_n->getComponent() == 1 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 4 && this_n->getComponent() == 2 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 4 && this_n->getComponent() == 4 	) return true;
	if(	this_n->getSector() == 3 && this_n->getLayer() == 4 && this_n->getComponent() == 6 	) return true;
	if(	this_n->getSector() == 3 && this_n->getLayer() == 5 && this_n->getComponent() == 3 	) return true;
	*/
	return false;
}
