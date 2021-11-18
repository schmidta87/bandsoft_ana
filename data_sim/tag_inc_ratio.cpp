#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"

#include "genpart.h"
#include "clashit.h"
#include "bandhit.h"
#include "taghit.h"

#include "bin_edges.h"

using std::cerr;
using std::isfinite;
using std::cout;

void fillHist( double Q2, double Pt, double Xp, double As, double Xb, TH1D**** hist , TH1D***** hist_xb , TH1D**** hist_Q2, double weight );
void fillArray( double array[bins_Q2][bins_Pt][bins_As][bins_Xp], TH1D**** hist_Q2, double Q2, double Xb , int FLAG );
void background_subtraction(TH1D* dat, TH1D* bac , double Cscale, double NB_sim , double Sigma_Cscale, double Sigma_NB_sim );
void simulation_weighting(TH1D* sim, double Ndata, double Nsim );
bool bad_bar( bandhit * this_n );

// Bin edges of Xb to use for inclusive for our ratio of tagged to inclusive
double DAT_MIN_XB[bins_Q2][bins_Pt][bins_As][bins_Xp];
double DAT_MAX_XB[bins_Q2][bins_Pt][bins_As][bins_Xp];
double SIM_MIN_XB[bins_Q2][bins_Pt][bins_As][bins_Xp];
double SIM_MAX_XB[bins_Q2][bins_Pt][bins_As][bins_Xp];
double DAT_MIN_Q2[bins_Q2][bins_Pt][bins_As][bins_Xp];
double DAT_MAX_Q2[bins_Q2][bins_Pt][bins_As][bins_Xp];
double SIM_MIN_Q2[bins_Q2][bins_Pt][bins_As][bins_Xp];
double SIM_MAX_Q2[bins_Q2][bins_Pt][bins_As][bins_Xp];

int main( int argc, char** argv){

	if( argc != 6 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Dat File] [Bac File] [Sim File] [Inc Dat File] [Inc Sim File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << "\n";


	// Luminosity weights:
	double Q_INC_DAT = 535885;	// mC
	double Q_TAG_DAT = 3.89811e7;	// mC
	double L_TAG_MC = 3.73246e9; 	// fb^-1
	double L_INC_MC = 6.11738e5;	// fb^-1 

	// Load files
	TFile inFile_Dat(argv[1]);
	TFile inFile_Bac(argv[2]);
	TFile inFile_Sim(argv[3]);
	TFile inFile_Inc_Dat(argv[4]);
	TFile inFile_Inc_Sim(argv[5]);

	// Load TTrees from files
	TTree * inTree_Dat = (TTree*) inFile_Dat.Get("tagged");
	TTree * inTree_Bac = (TTree*) inFile_Bac.Get("tagged");
	TTree * inTree_Sim = (TTree*) inFile_Sim.Get("tagged");
	TTree * inTree_Inc_Dat = (TTree*) inFile_Inc_Dat.Get("electrons");
	TTree * inTree_Inc_Sim = (TTree*) inFile_Inc_Sim.Get("electrons");

	// Background normalization
	TVector3 * datnorm = (TVector3*) inFile_Dat.Get("bacnorm");
	TVector3 * bacnorm = (TVector3*) inFile_Bac.Get("bacnorm");
	// normalization uncertainty of the background:
	double Cscale = datnorm->Z(); // this is the same as bacnorm->X() for the data file
	double Sigma_Cscale = (datnorm->Y() - datnorm->X())/2.;  // this is the before-time and after-time levels
	double NB_sim = bacnorm->X();
	double Sigma_NB_sim = sqrt(NB_sim);
		// sigma_Cscale / Cscale ~ 7%
	double Ndata = 0;
	double Nsim = 0;
	double sim_scaling = 0;
	
	// Set the input branches for data
	clashit * dat_eHit 		= new clashit;
	TClonesArray* dat_tagged 	= new TClonesArray("taghit");
	TClonesArray* dat_bandhi 	= new TClonesArray("bandhit");
	int dat_nleadindex		= 0;
	inTree_Dat->SetBranchAddress("tag"		,&dat_tagged		);
	inTree_Dat->SetBranchAddress("nHits"		,&dat_bandhi		);
	inTree_Dat->SetBranchAddress("eHit"		,&dat_eHit		);
	inTree_Dat->SetBranchAddress("nleadindex"	,&dat_nleadindex	);

	// Set the input branches for background
	clashit * bac_eHit 		= new clashit;
	TClonesArray* bac_tagged 	= new TClonesArray("taghit");
	TClonesArray* bac_bandhi 	= new TClonesArray("bandhit");
	int bac_nleadindex		= 0;
	inTree_Bac->SetBranchAddress("tag"		,&bac_tagged		);
	inTree_Bac->SetBranchAddress("nHits"		,&bac_bandhi		);
	inTree_Bac->SetBranchAddress("eHit"		,&bac_eHit		);
	inTree_Bac->SetBranchAddress("nleadindex"	,&bac_nleadindex	);

	// Set the input branches for simulation
	clashit * sim_eHit 		= new clashit;
	TClonesArray* sim_tagged 	= new TClonesArray("taghit");
	TClonesArray* sim_bandhi 	= new TClonesArray("bandhit");
	int sim_nleadindex		= 0;
	inTree_Sim->SetBranchAddress("tag_smeared"		,&sim_tagged		);
	inTree_Sim->SetBranchAddress("nHits"			,&sim_bandhi		);
	inTree_Sim->SetBranchAddress("eHit_smeared"		,&sim_eHit		);
	inTree_Sim->SetBranchAddress("nleadindex"	,&sim_nleadindex	);

	// Set the input branches for inclusive data
	clashit * inc_dat_eHit 		= new clashit;
	inTree_Inc_Dat->SetBranchAddress("eHit"		,&inc_dat_eHit		);

	// Set the input branches for inclusive simulation
	clashit * inc_sim_eHit 		= new clashit;
	inTree_Inc_Sim->SetBranchAddress("eHit_smeared"	,&inc_sim_eHit		);
	

	// Define the histograms for yield to fill:
	TFile * outFile = new TFile("yield_tagged.root","RECREATE");
	TH1D **** h4_dat_xp 	= new TH1D***[bins_Q2];
	TH1D **** h4_bac_xp	= new TH1D***[bins_Q2];
	TH1D **** h4_sim_xp	= new TH1D***[bins_Q2];
	TH1D ***** h4_dat_xb 	= new TH1D****[bins_Q2];
	TH1D ***** h4_bac_xb	= new TH1D****[bins_Q2];
	TH1D ***** h4_sim_xb	= new TH1D****[bins_Q2];

	TH1D **** h4_dat_Q2	= new TH1D***[bins_Pt];
	TH1D **** h4_bac_Q2	= new TH1D***[bins_Pt];
	TH1D **** h4_sim_Q2	= new TH1D***[bins_Pt];
	TH1D **** h4_inc_dat_Q2	= new TH1D***[bins_Pt];
	TH1D **** h4_inc_sim_Q2	= new TH1D***[bins_Pt];

	TH1D **** h4_dat_sim_xp 	= new TH1D***[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // bins in Q2
		h4_dat_xp[i]	= new TH1D**[bins_Pt];
		h4_bac_xp[i]	= new TH1D**[bins_Pt];
		h4_sim_xp[i]	= new TH1D**[bins_Pt];
		h4_dat_sim_xp[i]	= new TH1D**[bins_Pt];
		h4_dat_xb[i]	= new TH1D***[bins_Pt];
		h4_bac_xb[i]	= new TH1D***[bins_Pt];
		h4_sim_xb[i]	= new TH1D***[bins_Pt];
		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			h4_dat_xp[i][j] 	= new TH1D*[bins_As];
			h4_bac_xp[i][j] 	= new TH1D*[bins_As];
			h4_sim_xp[i][j] 	= new TH1D*[bins_As];
			h4_dat_sim_xp[i][j] 	= new TH1D*[bins_As];
			h4_dat_xb[i][j] 	= new TH1D**[bins_As];
			h4_bac_xb[i][j] 	= new TH1D**[bins_As];
			h4_sim_xb[i][j] 	= new TH1D**[bins_As];
			h4_dat_Q2[j]	= new TH1D**[bins_As];
			h4_bac_Q2[j]	= new TH1D**[bins_As];
			h4_sim_Q2[j]	= new TH1D**[bins_As];
			h4_inc_dat_Q2[j]	= new TH1D**[bins_As];
			h4_inc_sim_Q2[j]	= new TH1D**[bins_As];
			for( int k = 0 ; k < bins_As ; ++k ){ // bins in As


				h4_dat_xp[i][j][k] = new TH1D(Form("dat_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("dat_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in Xp
				h4_bac_xp[i][j][k] = new TH1D(Form("bac_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("bac_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in Xp
				h4_sim_xp[i][j][k] = new TH1D(Form("sim_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("sim_yield_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in Xp
				h4_dat_sim_xp[i][j][k] = new TH1D(Form("dat_sim_ratio_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("dat_sim_ratio_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in Xp

				h4_dat_xb[i][j][k] 	= new TH1D*[bins_Xp];
				h4_bac_xb[i][j][k] 	= new TH1D*[bins_Xp];
				h4_sim_xb[i][j][k] 	= new TH1D*[bins_Xp];
				h4_dat_Q2[j][k]		= new TH1D*[bins_Xp];
				h4_bac_Q2[j][k]		= new TH1D*[bins_Xp];
				h4_sim_Q2[j][k]		= new TH1D*[bins_Xp];
				h4_inc_dat_Q2[j][k]		= new TH1D*[bins_Xp];
				h4_inc_sim_Q2[j][k]		= new TH1D*[bins_Xp];
				for( int m = 0 ; m < bins_Xp ; ++m ){
					h4_dat_xb[i][j][k][m] = new TH1D(Form("dat_yield_xb_Q2_%i_Pt_%i_As_%i_Xp_%i",i,j,k,m),Form("dat_yield_xb_Q2_%i_Pt_%i_As_%i_Xp_%i",i,j,k,m),bins_Xb,Xb_min[i],Xb_min[i] + Xb_step*bins_Xb);  // histogram in Xb
					h4_bac_xb[i][j][k][m] = new TH1D(Form("bac_yield_xb_Q2_%i_Pt_%i_As_%i_Xp_%i",i,j,k,m),Form("bac_yield_xb_Q2_%i_Pt_%i_As_%i_Xp_%i",i,j,k,m),bins_Xb,Xb_min[i],Xb_min[i] + Xb_step*bins_Xb);  // histogram in Xb
					h4_sim_xb[i][j][k][m] = new TH1D(Form("sim_yield_xb_Q2_%i_Pt_%i_As_%i_Xp_%i",i,j,k,m),Form("sim_yield_xb_Q2_%i_Pt_%i_As_%i_Xp_%i",i,j,k,m),bins_Xb,Xb_min[i],Xb_min[i] + Xb_step*bins_Xb);  // histogram in Xb
					h4_dat_Q2[j][k][m] = new TH1D(Form("dat_yield_Q2_Pt_%i_As_%i_Xp_%i",j,k,m),Form("dat_yield_Q2_Pt_%i_As_%i_Xp_%i",j,k,m), 16, 2, 6);
					h4_bac_Q2[j][k][m] = new TH1D(Form("bac_yield_Q2_Pt_%i_As_%i_Xp_%i",j,k,m),Form("bac_yield_Q2_Pt_%i_As_%i_Xp_%i",j,k,m), 16, 2, 6);
					h4_sim_Q2[j][k][m] = new TH1D(Form("sim_yield_Q2_Pt_%i_As_%i_Xp_%i",j,k,m),Form("sim_yield_Q2_Pt_%i_As_%i_Xp_%i",j,k,m), 16, 2, 6);
					h4_inc_dat_Q2[j][k][m] = new TH1D(Form("inc_dat_yield_Q2_Pt_%i_As_%i_Xp_%i",j,k,m),Form("inc_dat_yield_Q2_Pt_%i_As_%i_Xp_%i",j,k,m), 16, 2, 6);
					h4_inc_sim_Q2[j][k][m] = new TH1D(Form("inc_sim_yield_Q2_Pt_%i_As_%i_Xp_%i",j,k,m),Form("inc_sim_yield_Q2_Pt_%i_As_%i_Xp_%i",j,k,m), 16, 2, 6);
				}
			}
		}
	}

	// Loop over the data file:
	for( int event = 0 ; event < inTree_Dat->GetEntries() ; ++event ){
		dat_eHit	->Clear();
		dat_tagged	->Clear();
		dat_bandhi	->Clear();
		dat_nleadindex	= 0.;

		inTree_Dat->GetEntry(event);
		
		// Get the correct tag hit
		taghit * this_tag	= (taghit*)  dat_tagged->At(dat_nleadindex);
		bandhit * lead_n	= (bandhit*) dat_bandhi->At(dat_nleadindex);
		
		double As	= this_tag->getAs();
		double Pt	= this_tag->getPt().Mag();
		double Xp	= this_tag->getXp_WP();
		double Q2	= dat_eHit->getQ2();
		double Pn	= this_tag->getMomentumN().Mag();
		double Xb	= dat_eHit->getXb();
		if( lead_n->getEdep() < 10 ) continue;
		if( Pn < 0.25 ) continue;

		if( bad_bar(lead_n) ) continue;
			
	
		//if( lead_n->getSector()==2 && (lead_n->getComponent() > 3 && lead_n->getComponent() < 8 ) 
		//		&& ( lead_n->getX()>90 || lead_n->getX() < -110 ) ) continue;
		//if( lead_n->getSector()==3 && (lead_n->getComponent() > 0 && lead_n->getComponent() < 3 ) 
		//		&& ( lead_n->getX()>80 || lead_n->getX() < 45 ) ) continue;

		fillHist( Q2, Pt, Xp, As, Xb, h4_dat_xp , h4_dat_xb , h4_dat_Q2 , 1. );
	}

	// Loop over the background file:
	for( int event = 0 ; event < inTree_Bac->GetEntries() ; ++event ){
		bac_eHit	->Clear();
		bac_tagged	->Clear();
		bac_bandhi	->Clear();
		bac_nleadindex	= 0.;

		inTree_Bac->GetEntry(event);
		
		// Get the correct tag hit
		taghit * this_tag	= (taghit*)  bac_tagged->At(bac_nleadindex);
		bandhit * lead_n	= (bandhit*) bac_bandhi->At(bac_nleadindex);
		if( lead_n->getEdep() < 10 ) continue;

		double As	= this_tag->getAs();
		double Pt	= this_tag->getPt().Mag();
		double Xp	= this_tag->getXp_WP();
		double Q2	= bac_eHit->getQ2();
		double Pn	= this_tag->getMomentumN().Mag();
		double Xb	= bac_eHit->getXb();
		if( Pn < 0.25 ) continue;

		if( bad_bar(lead_n) ) continue;

		//if( lead_n->getSector()==2 && (lead_n->getComponent() > 3 && lead_n->getComponent() < 8 ) 
		//		&& ( lead_n->getX()>90 || lead_n->getX() < -110 ) ) continue;
		//if( lead_n->getSector()==3 && (lead_n->getComponent() > 0 && lead_n->getComponent() < 3 ) 
		//		&& ( lead_n->getX()>80 || lead_n->getX() < 45 ) ) continue;

		fillHist( Q2, Pt, Xp, As, Xb, h4_bac_xp , h4_bac_xb , h4_bac_Q2 , 1. );
	}

	// Loop over the simulation file:
	for( int event = 0 ; event < inTree_Sim->GetEntries() ; ++event ){
		sim_eHit	->Clear();
		sim_tagged	->Clear();
		sim_bandhi	->Clear();
		sim_nleadindex	= 0.;

		inTree_Sim->GetEntry(event);
		
		// Get the correct tag hit
		taghit * this_tag	= (taghit*)  sim_tagged->At(sim_nleadindex);
		bandhit * lead_n	= (bandhit*) sim_bandhi->At(sim_nleadindex);
		if( lead_n->getEdep() < 10 ) continue;
		
		double As	= this_tag->getAs();
		double Pt	= this_tag->getPt().Mag();
		double Xp	= this_tag->getXp_WP();
		double Q2	= sim_eHit->getQ2();
		double Pn	= this_tag->getMomentumN().Mag();
		double Xb	= sim_eHit->getXb();
		if( Pn < 0.25 ) continue;

		if( bad_bar(lead_n) ) continue;

		//if( lead_n->getSector()==2 && (lead_n->getComponent() > 3 && lead_n->getComponent() < 8 ) 
		//		&& ( lead_n->getX()>90 || lead_n->getX() < -110 ) ) continue;
		//if( lead_n->getSector()==3 && (lead_n->getComponent() > 0 && lead_n->getComponent() < 3 ) 
		//		&& ( lead_n->getX()>80 || lead_n->getX() < 45 ) ) continue;

		fillHist( Q2, Pt, Xp, As, Xb, h4_sim_xp , h4_sim_xb , h4_sim_Q2 , 1. );
	}


	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			for( int k = 0 ; k < bins_As ; ++k ){ // loop over As bins

				if( h4_dat_xp[i][j][k]->Integral() < 1 || h4_bac_xp[i][j][k]->Integral() < 1 || h4_sim_xp[i][j][k]->Integral() < 1 ) continue;

				// Background subtraction of data and mixed:
				background_subtraction( h4_dat_xp[i][j][k],
					       	h4_bac_xp[i][j][k] , Cscale, NB_sim , Sigma_Cscale, Sigma_NB_sim );
				for( int m = 0 ; m < bins_Xp ; ++m ){
					background_subtraction( h4_dat_xb[i][j][k][m],
							h4_bac_xb[i][j][k][m] , Cscale, NB_sim , Sigma_Cscale, Sigma_NB_sim );
					background_subtraction( h4_dat_Q2[j][k][m],
							h4_bac_Q2[j][k][m] , Cscale, NB_sim , Sigma_Cscale, Sigma_NB_sim );
					
					if( h4_dat_xb[i][j][k][m]->Integral() < 1 || h4_dat_Q2[j][k][m]->Integral() < 1 ){
						DAT_MIN_XB[i][j][k][m] = 0.;
						DAT_MAX_XB[i][j][k][m] = 0.;
						DAT_MIN_Q2[i][j][k][m] = 0.;
						DAT_MAX_Q2[i][j][k][m] = 0.;
						continue;
					}
					if( h4_sim_xb[i][j][k][m]->Integral() < 1 || h4_sim_Q2[j][k][m]->Integral() < 1 ){
						SIM_MIN_XB[i][j][k][m] = 0.;
						SIM_MAX_XB[i][j][k][m] = 0.;
						SIM_MIN_Q2[i][j][k][m] = 0.;
						SIM_MAX_Q2[i][j][k][m] = 0.;
						continue;
					}

					// Create the bins for our inclusive histograms:
					DAT_MIN_XB[i][j][k][m] = h4_dat_xb[i][i][k][m]->GetMean() - h4_dat_xb[i][j][k][m]->GetStdDev()*2.;
					DAT_MAX_XB[i][j][k][m] = h4_dat_xb[i][i][k][m]->GetMean() + h4_dat_xb[i][j][k][m]->GetStdDev()*2.;
					SIM_MIN_XB[i][j][k][m] = h4_sim_xb[i][i][k][m]->GetMean() - h4_sim_xb[i][j][k][m]->GetStdDev()*2.;
					SIM_MAX_XB[i][j][k][m] = h4_sim_xb[i][i][k][m]->GetMean() + h4_sim_xb[i][j][k][m]->GetStdDev()*2.;

					DAT_MIN_Q2[i][j][k][m] = h4_dat_Q2[j][k][m]->GetXaxis()->GetBinCenter( h4_dat_Q2[j][k][m]->FindFirstBinAbove(0) )-0.125 ;
					DAT_MAX_Q2[i][j][k][m] = h4_dat_Q2[j][k][m]->GetXaxis()->GetBinCenter( h4_dat_Q2[j][k][m]->FindLastBinAbove(0)  )+0.125 ;
					SIM_MIN_Q2[i][j][k][m] = h4_sim_Q2[j][k][m]->GetXaxis()->GetBinCenter( h4_sim_Q2[j][k][m]->FindFirstBinAbove(0) )-0.125 ;
					SIM_MAX_Q2[i][j][k][m] = h4_sim_Q2[j][k][m]->GetXaxis()->GetBinCenter( h4_sim_Q2[j][k][m]->FindLastBinAbove(0) 	)+0.125 ;

					h4_dat_xb[i][j][k][m]->Write();
					h4_sim_xb[i][j][k][m]->Write();
					h4_dat_Q2[j][k][m]->Write();
					h4_sim_Q2[j][k][m]->Write();
				}
				

				// For all the bins in the data (after background subtraction), divide the data by simultion
				for( int bin = 1 ; bin < h4_dat_xp[i][j][k]->GetXaxis()->GetNbins(); ++bin ){
					// Set the bin content to be: (data/sim) * luminosity scaling
					double bincontent = h4_dat_xp[i][j][k]->GetBinContent(bin) / 
						h4_sim_xp[i][j][k]->GetBinContent(bin) * (L_TAG_MC / Q_TAG_DAT);	
					if( !isfinite(bincontent) || bincontent < 0 ){
						h4_dat_sim_xp[i][j][k]->SetBinContent(bin, 0);
						h4_dat_sim_xp[i][j][k]->SetBinError(bin , 0);
						continue;
					};
					h4_dat_sim_xp[i][j][k]->SetBinContent(bin, bincontent);
					
					// NEED TO CHECK THIS AND MAKE SURE IT'S CORRECT:
					double s = h4_sim_xp[i][j][k]->GetBinContent(bin) ;
					double d = h4_dat_xp[i][j][k]->GetBinContent(bin) ;
					double a = (L_TAG_MC / Q_TAG_DAT);
					double err = sqrt( a*a*d/(s*s) + a*a*d*d/(s*s*s) );

					h4_dat_sim_xp[i][j][k]->SetBinError( bin , err );
				} // loop over Xp bins
			} // loop over As
		} // loop over Pt
	} // loop over Q2


	// Now figure out the scaling of inclusive we have:
	double INC_COUNTS_DAT[bins_Q2][bins_Pt][bins_As][bins_Xp];
	double INC_COUNTS_SIM[bins_Q2][bins_Pt][bins_As][bins_Xp];

	// Loop over the inclusive data file:
	for( int event = 0 ; event < inTree_Inc_Dat->GetEntries() ; ++event ){
		inc_dat_eHit	->Clear();

		inTree_Inc_Dat->GetEntry(event);

		double Q2	= inc_dat_eHit->getQ2();
		double Xb	= inc_dat_eHit->getXb();

		fillArray( INC_COUNTS_DAT, h4_inc_dat_Q2, Q2, Xb , 1);
	}

	// Loop over the inclusive simulation file:
	for( int event = 0 ; event < inTree_Inc_Sim->GetEntries() ; ++event ){
		inc_sim_eHit	->Clear();

		inTree_Inc_Sim->GetEntry(event);
		
		double Q2	= inc_sim_eHit->getQ2();
		double Xb	= inc_sim_eHit->getXb();

		fillArray( INC_COUNTS_SIM, h4_inc_sim_Q2, Q2, Xb , 0);
	}
	
	TCanvas ** c_Q2 = new TCanvas*[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		c_Q2[i] = new TCanvas(Form("c_Q2_%i",i),"",1200,1200);
		c_Q2[i]->Divide( bins_Pt, bins_As );

		TCanvas ** c_compareQ2_dat = new TCanvas*[bins_Pt];
		TCanvas ** c_compareQ2_sim = new TCanvas*[bins_Pt];
		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			c_compareQ2_dat[j] = new TCanvas(Form("c_compareQ2_dat_%i",j),"",3600,2400);
			c_compareQ2_sim[j] = new TCanvas(Form("c_compareQ2_sim_%i",j),"",3600,2400);
			c_compareQ2_dat[j]->Divide( bins_Xp, bins_As );
			c_compareQ2_sim[j]->Divide( bins_Xp, bins_As );
			for( int k = 0 ; k < bins_As ; ++k ){ // loop over As bins
				cout << "[Q2, Pt, As]: " << i << " " << j << " " << k << "\n";
				for( int m = 0 ; m < bins_Xp ; ++m ){ // loop over Xp bins
					
					double content = h4_dat_sim_xp[i][j][k]->GetBinContent(m+1);
					double newcontent = content * (INC_COUNTS_SIM[i][j][k][m] / INC_COUNTS_DAT[i][j][k][m]) * ( Q_INC_DAT / L_INC_MC );

					double f = content;
					double err1 = h4_dat_sim_xp[i][j][k]->GetBinError(m+1);

					double ref_d = INC_COUNTS_DAT[i][j][k][m];
					double ref_s = INC_COUNTS_SIM[i][j][k][m];
					double ref_a = L_INC_MC / Q_INC_DAT;

					double ref_f = ref_d / ref_s * ref_a;
					double ref_err = sqrt( ref_a*ref_a*ref_d/(ref_s*ref_s) + ref_a*ref_a*ref_d*ref_d/pow(ref_s,3) );

					double ratio_err = sqrt( 
						pow(1./ref_f * err1,2) + pow(f * ref_err,2)/pow(ref_f,4) );

					if( !isfinite(newcontent) || newcontent < 0  || !isfinite(ratio_err) || ratio_err < 0 ){
						h4_dat_sim_xp[i][j][k]->SetBinContent(m+1, 0);
						h4_dat_sim_xp[i][j][k]->SetBinError(m+1 , 0);
						continue;
					};

					h4_dat_sim_xp[i][j][k]->SetBinContent(m+1, newcontent);
					h4_dat_sim_xp[i][j][k]->SetBinError(m+1, ratio_err );

					cout << "\t" 
						<< h4_dat_sim_xp[i][j][k]->GetXaxis()->GetBinCenter(m+1) << " "
						<< h4_dat_sim_xp[i][j][k]->GetBinContent(m+1) << " " 
						<< h4_dat_sim_xp[i][j][k]->GetBinError(m+1) << "\n";

					TString title = Form("%.1f < p_{T} < %.1f, %.2f < a_{S} < %.2f, %.2f < x' < %.2f",
							Pt_min + Pt_step*j,Pt_min + Pt_step*(j+1),
							As_min + As_step*k,As_min + As_step*(k+1),
							Xp_min + Xp_step*m,Xp_min + Xp_step*(m+1)	);
					c_compareQ2_dat[j]->cd( (k*bins_Xp)+1 + m );
					h4_dat_Q2[j][k][m]->SetTitle(title);
					h4_dat_Q2[j][k][m]->SetStats(0);
					h4_dat_Q2[j][k][m]->SetLineWidth(2);
					h4_dat_Q2[j][k][m]->SetMarkerColor(4);
					h4_dat_Q2[j][k][m]->GetXaxis()->SetTitle("Q^{2}");
					h4_dat_Q2[j][k][m]->SetMinimum(0);
					h4_dat_Q2[j][k][m]->Draw("*P");
					h4_inc_dat_Q2[j][k][m]->Scale( h4_dat_Q2[j][k][m]->Integral() /h4_inc_dat_Q2[j][k][m]->Integral() );
					h4_inc_dat_Q2[j][k][m]->SetLineColor(2);
					h4_inc_dat_Q2[j][k][m]->Draw("HIST,SAME");

					TLine * vline1 = new TLine(DAT_MIN_Q2[i][j][k][m],0,DAT_MIN_Q2[i][j][k][m],h4_dat_Q2[j][k][m]->GetMaximum());
					TLine * vline2 = new TLine(DAT_MAX_Q2[i][j][k][m],0,DAT_MAX_Q2[i][j][k][m],h4_dat_Q2[j][k][m]->GetMaximum());
					vline1->SetLineWidth(2);
					vline1->SetLineColor(1);
					vline1->SetLineStyle(2);
					vline1->Draw("SAME");
					vline2->SetLineWidth(2);
					vline2->SetLineColor(1);
					vline2->SetLineStyle(2);
					vline2->Draw("SAME");


					TLegend * legend = new TLegend(0.7,0.8,0.9,0.9);
					legend->AddEntry(h4_dat_Q2[j][k][m],	"Tagged","f");
					legend->AddEntry(h4_inc_dat_Q2[j][k][m],"Inclusive","f");
					legend->Draw("SAME");


					c_compareQ2_sim[j]->cd( (k*bins_Xp)+1 + m );
					h4_sim_Q2[j][k][m]->SetTitle(title);
					h4_sim_Q2[j][k][m]->SetStats(0);
					h4_sim_Q2[j][k][m]->SetLineWidth(2);
					h4_sim_Q2[j][k][m]->SetMarkerColor(4);
					h4_sim_Q2[j][k][m]->GetXaxis()->SetTitle("Q^{2}");
					h4_sim_Q2[j][k][m]->SetMinimum(0);
					h4_sim_Q2[j][k][m]->Draw("*P");
					h4_inc_sim_Q2[j][k][m]->Scale( h4_sim_Q2[j][k][m]->Integral() /h4_inc_sim_Q2[j][k][m]->Integral() );
					h4_inc_sim_Q2[j][k][m]->SetLineColor(2);
					h4_inc_sim_Q2[j][k][m]->Draw("HIST,SAME");

					TLine * vline3 = new TLine(SIM_MIN_Q2[i][j][k][m],0,SIM_MIN_Q2[i][j][k][m],h4_sim_Q2[j][k][m]->GetMaximum());
					TLine * vline4 = new TLine(SIM_MAX_Q2[i][j][k][m],0,SIM_MAX_Q2[i][j][k][m],h4_sim_Q2[j][k][m]->GetMaximum());
					vline3->SetLineWidth(2);
					vline3->SetLineColor(1);
					vline3->SetLineStyle(2);
					vline3->Draw("SAME");
					vline4->SetLineWidth(2);
					vline4->SetLineColor(1);
					vline4->SetLineStyle(2);
					vline4->Draw("SAME");

					legend->Draw("SAME");




					h4_inc_dat_Q2[j][k][m]->Write();
					h4_inc_sim_Q2[j][k][m]->Write();

				} // loop over Xp bins

				c_Q2[i]->cd( (k*2)+1 + j );

				if( h4_dat_sim_xp[i][j][k]->Integral() < 1 ) continue;

				TString title = Form("%.1f < Q^{2} < %.1f, %.1f < p_{T} < %.1f, %.2f < a_{S} < %.2f",
						Q2Bins[i],Q2Bins[i+1],
						Pt_min + Pt_step*j,Pt_min + Pt_step*(j+1),
						As_min + As_step*k,As_min + As_step*(k+1)	);
				h4_dat_sim_xp[i][j][k]->SetTitle(title);

				h4_dat_sim_xp[i][j][k]->SetLineWidth(3);
				h4_dat_sim_xp[i][j][k]->SetMarkerColor(4);
				h4_dat_sim_xp[i][j][k]->SetMinimum(0);
				h4_dat_sim_xp[i][j][k]->SetMaximum(40.0);
				h4_dat_sim_xp[i][j][k]->GetXaxis()->SetTitle("x'");
				h4_dat_sim_xp[i][j][k]->Draw("*P");

				TLine * hline = new TLine(Xp_min,1,Xp_max,1);
				hline->SetLineWidth(2);
				hline->SetLineColor(1);
				hline->SetLineStyle(2);
				hline->Draw("SAME");

				h4_dat_sim_xp[i][j][k]->Write();

			} // loop over As
			c_compareQ2_dat[j]->Print(Form("compareQ2_dat_Pt_%i.pdf",j));
			c_compareQ2_sim[j]->Print(Form("compareQ2_sim_Pt_%i.pdf",j));
		} // loop over Pt
		c_Q2[i]->Print(Form("tagged_yield_Q2_%i.pdf",i));
	} // loop over Q2
				
	/*

			}
		}
	}
	*/
	outFile->Close();

	inFile_Dat.Close();
	inFile_Bac.Close();
	inFile_Sim.Close();
	return 1;
}

void fillHist( double Q2, double Pt, double Xp, double As, double Xb, TH1D**** hist , TH1D***** hist_xb , TH1D**** hist_Q2, double weight ){
	

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

	int this_bin_xp = (int) ((Xp - Xp_min)/Xp_step);
	hist_xb[this_bin_q2][this_bin_pt][this_bin_as][this_bin_xp]->Fill( Xb , weight );
	hist_Q2[this_bin_pt][this_bin_as][this_bin_xp]->Fill( Q2, weight );

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

void fillArray( double array[bins_Q2][bins_Pt][bins_As][bins_Xp], TH1D**** hist_Q2, double Q2, double Xb , int FLAG ){

	// If it's larger than max Q2 or smaller than min Q2, return
	if( Q2 > Q2Bins[bins_Q2] 	) return;
	if( Q2 < Q2Bins[0] 		) return;

	// Get bin of Q2:
	int this_bin_q2 = -1;
	for( int q2_bin = bins_Q2-1 ; q2_bin >= 0; --q2_bin ){
		if( Q2 > Q2Bins[q2_bin] ) this_bin_q2 = q2_bin;
		if( this_bin_q2 != -1 ) break;
	}
	if( this_bin_q2 == -1 ){ cerr << "how\n"; return; }

	// Figure out bin of Xb needed:
	for( int j = 0 ; j < bins_Pt ; ++j ){
		for( int k = 0 ; k < bins_As ; ++k ){
			for( int m = 0 ; m < bins_Xp ; ++m ){
				double min = 0.;
				double max = 0.;
				double min_q2 = 0.;
				double max_q2 = 0.;
				if( FLAG == 1 ){
					min = DAT_MIN_XB[this_bin_q2][j][k][m];
					max = DAT_MAX_XB[this_bin_q2][j][k][m];
					min_q2 = DAT_MIN_Q2[this_bin_q2][j][k][m];
					max_q2 = DAT_MAX_Q2[this_bin_q2][j][k][m];
				}
				else if( FLAG == 0 ){
					min = SIM_MIN_XB[this_bin_q2][j][k][m];
					max = SIM_MAX_XB[this_bin_q2][j][k][m];
					min_q2 = SIM_MIN_Q2[this_bin_q2][j][k][m];
					max_q2 = SIM_MAX_Q2[this_bin_q2][j][k][m];
				}

				if( Xb < max && Xb > min && Q2 < max_q2 && Q2 > min_q2 ){
					array[this_bin_q2][j][k][m] += 1;
					hist_Q2[j][k][m]->Fill(Q2);
				}

			}
		}
	}

	return;
}
