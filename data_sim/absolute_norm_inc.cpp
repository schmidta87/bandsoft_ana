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

#include "bin_edges.h"

using std::cerr;
using std::isfinite;
using std::cout;

const int inc_bins_Xb = 6;
const double inc_min_Xb = 0.1;
const double inc_max_Xb = 0.7;
const double inc_step_Xb = (inc_max_Xb - inc_min_Xb)/inc_bins_Xb;

const int inc_bins_Q2 = 13;
const double inc_min_Q2 = 2;
const double inc_max_Q2 = 8.5;
const double inc_step_Q2 = (inc_max_Q2 - inc_min_Q2)/inc_bins_Q2;

const int inc_bins_theta = 7;
const double inc_min_theta = 8.5;
const double inc_max_theta = 29.5;
const double inc_step_theta = (inc_max_theta - inc_min_theta)/inc_bins_theta;

const int inc_bins_phi = 72;
const double inc_min_phi = -180.;
const double inc_max_phi = 180.;
const double inc_step_phi = (inc_max_phi - inc_min_phi)/inc_bins_phi;

const int inc_bins_pe = 9;
const double inc_min_pe = 3.;
const double inc_max_pe = 7.5;
const double inc_step_pe = (inc_max_pe - inc_min_pe)/inc_bins_pe;

void fillHist( TH1D** hist_xb, TH1D** hist_Q2, TH1D*** hist_pe, TH1D*** hist_theta, TH1D*** hist_phi, double Xb, double Q2, double Pe, double Theta, double Phi , double weight );
void drawHists( TH1D* data, TH1D* sim, TString title, TString xtitle );
void drawHistsRatio( TH1D* data, TH1D* sim, TString title, TString xtitle );

int main( int argc, char** argv){

	if( argc != 3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Inc Dat File] [Inc Sim File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] << "\n";


	// Luminosity weights:
	double L_INC_MC = 6.11738e5 / 1E6;		// nb^-1 -> fb^-1 

	double Q_INC_DAT = 535885;		// nC
	double LD2_den = 0.1644; 		// g/cm^3
	double target_L = 5.;			// cm  (-5 to -1)
	double mDeut = 3.3435837724E-24;	// g
	double Coul = 1.60218E-10;		// nC
	double cm2_to_fb = 1E-39;		// 1/cm^2 -> 1/fb
	double num_electron = Q_INC_DAT / Coul;				// # electrons
	double target_density = LD2_den * target_L / mDeut;		// # nucleons / cm^2
	double L_INC_DAT = num_electron * target_density * cm2_to_fb;	// fb^-1

	// Load files
	TFile inFile_Inc_Dat(argv[1]);
	TFile inFile_Inc_Sim(argv[2]);

	// Load TTrees from files
	TTree * inTree_Inc_Dat = (TTree*) inFile_Inc_Dat.Get("electrons");
	TTree * inTree_Inc_Sim = (TTree*) inFile_Inc_Sim.Get("electrons");

	// Set the input branches for inclusive data
	clashit * inc_dat_eHit 		= new clashit;
	inTree_Inc_Dat->SetBranchAddress("eHit"		,&inc_dat_eHit		);

	// Set the input branches for inclusive simulation
	clashit * inc_sim_eHit 		= new clashit;
	inTree_Inc_Sim->SetBranchAddress("eHit_smeared"	,&inc_sim_eHit		);
	

	// Define the histograms for yield to fill:
	// 	bin in Q2(x)		--
	// 	bin in x(Q2)		--
	// 	bin in p(theta,phi)	--
	//	bin in theta(p,phi)	--
	//	bin in phi(theta,p)	--
	TFile * outFile = new TFile("inc_comparison.root","RECREATE");

	TH1D ** h2_dat_xb 	= new TH1D*[inc_bins_Q2]; // plot of Xb in bins of Q2
	TH1D ** h2_sim_xb	= new TH1D*[inc_bins_Q2];
	for( int i = 0 ; i < inc_bins_Q2 ; ++i ){ 
		h2_dat_xb[i] = new TH1D(Form("h2_dat_xb_Q2_%i",i),Form("h2_dat_xb_Q2_%i",i),inc_bins_Xb,inc_min_Xb,inc_max_Xb);
		h2_sim_xb[i] = new TH1D(Form("h2_sim_xb_Q2_%i",i),Form("h2_sim_xb_Q2_%i",i),inc_bins_Xb,inc_min_Xb,inc_max_Xb);
	}
	TH1D ** h2_dat_Q2	= new TH1D*[inc_bins_Xb]; // plot of Q2 in bins of Xb
	TH1D ** h2_sim_Q2	= new TH1D*[inc_bins_Xb];
	for( int j = 0 ; j < inc_bins_Xb ; ++j ){
		h2_dat_Q2[j] = new TH1D(Form("h2_dat_Q2_xb_%i",j),Form("h2_dat_Q2_xb_%i",j),inc_bins_Q2,inc_min_Q2,inc_max_Q2);
		h2_sim_Q2[j] = new TH1D(Form("h2_sim_Q2_xb_%i",j),Form("h2_sim_Q2_xb_%i",j),inc_bins_Q2,inc_min_Q2,inc_max_Q2);
	}


	TH1D *** h3_dat_pe 	= new TH1D**[inc_bins_theta]; // plot of pe in bins of theta, phi
	TH1D *** h3_dat_phi	= new TH1D**[inc_bins_theta]; // plot of phi in bins of theta, pe
	TH1D *** h3_sim_pe 	= new TH1D**[inc_bins_theta]; // plot of pe in bins of theta, phi
	TH1D *** h3_sim_phi	= new TH1D**[inc_bins_theta]; // plot of phi in bins of theta, pe
	for( int i = 0 ; i < inc_bins_theta ; ++i ){
		h3_dat_pe[i]	= new TH1D*[inc_bins_phi];
		h3_sim_pe[i]	= new TH1D*[inc_bins_phi];
		for( int j = 0 ; j < inc_bins_phi ; ++j ){
			h3_dat_pe[i][j] = new TH1D(Form("h3_dat_pe_theta_%i_phi_%i",i,j),Form("h3_dat_pe_theta_%i_phi_%i",i,j),inc_bins_pe,inc_min_pe,inc_max_pe);
			h3_sim_pe[i][j] = new TH1D(Form("h3_sim_pe_theta_%i_phi_%i",i,j),Form("h3_sim_pe_theta_%i_phi_%i",i,j),inc_bins_pe,inc_min_pe,inc_max_pe);
		}
		h3_dat_phi[i]	= new TH1D*[inc_bins_pe];
		h3_sim_phi[i]	= new TH1D*[inc_bins_pe];
		for( int j = 0 ; j < inc_bins_pe ; ++j ){
			h3_dat_phi[i][j] = new TH1D(Form("h3_dat_phi_theta_%i_pe_%i",i,j),Form("h3_dat_phi_theta_%i_pe_%i",i,j),inc_bins_phi,inc_min_phi,inc_max_phi);
			h3_sim_phi[i][j] = new TH1D(Form("h3_sim_phi_theta_%i_pe_%i",i,j),Form("h3_sim_phi_theta_%i_pe_%i",i,j),inc_bins_phi,inc_min_phi,inc_max_phi);
		}
	}
	TH1D *** h3_dat_theta	= new TH1D**[inc_bins_pe]; // plot of theta in bins of pe, phi
	TH1D *** h3_sim_theta	= new TH1D**[inc_bins_pe]; // plot of theta in bins of pe, phi
	for( int i = 0 ; i < inc_bins_pe ; ++i ){
		h3_dat_theta[i]	= new TH1D*[inc_bins_phi];
		h3_sim_theta[i]	= new TH1D*[inc_bins_phi];
		for( int j = 0 ; j < inc_bins_phi ; ++j ){
			h3_dat_theta[i][j] = new TH1D(Form("h3_dat_theta_pe_%i_phi_%i",i,j),Form("h3_dat_theta_pe_%i_phi_%i",i,j),inc_bins_theta,inc_min_theta,inc_max_theta);
			h3_sim_theta[i][j] = new TH1D(Form("h3_sim_theta_pe_%i_phi_%i",i,j),Form("h3_sim_theta_pe_%i_phi_%i",i,j),inc_bins_theta,inc_min_theta,inc_max_theta);
		}
	}

	// Loop over the inclusive data file:
	for( int event = 0 ; event < inTree_Inc_Dat->GetEntries() ; ++event ){
		inc_dat_eHit	->Clear();

		inTree_Inc_Dat->GetEntry(event);

		double Q2	= inc_dat_eHit->getQ2();
		double Xb	= inc_dat_eHit->getXb();
		double pe	= inc_dat_eHit->getMomentum();
		double theta	= inc_dat_eHit->getTheta() * 180./M_PI;
		double phi	= inc_dat_eHit->getPhi() * 180./M_PI;
	
		fillHist( h2_dat_xb, h2_dat_Q2, h3_dat_pe, h3_dat_theta, h3_dat_phi , Xb, Q2, pe, theta, phi , 1./L_INC_DAT);

	}

	// Loop over the inclusive simulation file:
	for( int event = 0 ; event < inTree_Inc_Sim->GetEntries() ; ++event ){
		inc_sim_eHit	->Clear();

		inTree_Inc_Sim->GetEntry(event);
		
		double Q2	= inc_sim_eHit->getQ2();
		double Xb	= inc_sim_eHit->getXb();
		double pe	= inc_sim_eHit->getMomentum();
		double theta	= inc_sim_eHit->getTheta() * 180./M_PI;
		double phi	= inc_sim_eHit->getPhi() * 180./M_PI;

		fillHist( h2_sim_xb, h2_sim_Q2, h3_sim_pe, h3_sim_theta, h3_sim_phi , Xb, Q2, pe, theta, phi , 1./L_INC_MC);
	}
	
	/*
	TCanvas ** c_Q2 = new TCanvas*[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		for( int m = 0 ; m < bins_Xp ; ++m ){ // loop over Xp bins

		} // loop over Xp bins
	} // loop over Q2
	*/
	outFile->cd();	

	TCanvas ** compare_xb = new TCanvas*[inc_bins_Q2];
	for( int i = 0 ; i < inc_bins_Q2 ; ++i ){ 
		compare_xb[i] = new TCanvas(Form("compare_xb_%i",i),Form("compare_xb_%i",i),1200,1200);
		compare_xb[i]->Divide(1,2);
		compare_xb[i]->cd(1);
		
		TString title = Form("%.2f < Q^{2} < %.2f", inc_min_Q2 + inc_step_Q2*i, inc_min_Q2 + inc_step_Q2*(i+1) );
		drawHists( h2_dat_xb[i], h2_sim_xb[i], title, "x_{B}" ); 

		compare_xb[i]->cd(2);
		drawHistsRatio( h2_dat_xb[i], h2_sim_xb[i], title, "x_{B}" ); 

		compare_xb[i]->Print(Form("compare_xb_Q2_%i.pdf",i));
		h2_dat_xb[i]->Write();
		h2_sim_xb[i]->Write();
	}

	TCanvas ** compare_q2 = new TCanvas*[inc_bins_Xb];
	for( int j = 0 ; j < inc_bins_Xb ; ++j ){
		compare_q2[j] = new TCanvas(Form("compare_q2_%i",j),Form("compare_q2_%i",j),1200,1200);
		compare_q2[j]->Divide(1,2);
		compare_q2[j]->cd(1);

		TString title = Form("%.2f < x_B < %.2f", inc_min_Xb + inc_step_Xb*j, inc_min_Xb + inc_step_Xb*(j+1) );
		drawHists( h2_dat_Q2[j], h2_sim_Q2[j], title, "Q^{2} [GeV]^{2}" ); 

		compare_q2[j]->cd(2);
		drawHistsRatio( h2_dat_Q2[j], h2_sim_Q2[j], title, "Q^{2} [GeV]^{2}" ); 

		compare_q2[j]->Print(Form("compare_Q2_xb_%i.pdf",j));
		h2_dat_Q2[j]->Write();
		h2_sim_Q2[j]->Write();
	}


	TCanvas *** compare_pe = new TCanvas**[inc_bins_theta];
	TCanvas *** compare_phi = new TCanvas**[inc_bins_theta];
	for( int i = 0 ; i < inc_bins_theta ; ++i ){
		compare_pe[i] = new TCanvas*[inc_bins_phi];
		for( int j = 0 ; j < inc_bins_phi ; ++j ){
			compare_pe[i][j] = new TCanvas(Form("compare_pe_theta_%i_phi_%i",i,j),Form("compare_pe_theta_%i_phi_%i",i,j),1200,1200);
			compare_pe[i][j]->Divide(1,2);
			TString title = Form("%.2f < Theta_{e} < %.2f , %.2f < Phi_{e} < %.2f", inc_min_theta + inc_step_theta*i , inc_min_theta + inc_step_theta*(i+1),
												inc_min_phi + inc_step_phi*j , inc_min_phi + inc_step_phi*(j+1) );

			// Make pdf plots of pe comparisons
			compare_pe[i][j]->cd(1);
			drawHists( h3_dat_pe[i][j] , h3_sim_pe[i][j], title, "p_{e} [GeV/c]");
			compare_pe[i][j]->cd(2);
			drawHistsRatio( h3_dat_pe[i][j] , h3_sim_pe[i][j], title, "p_{e} [GeV/c]");
			compare_pe[i][j]->Print(Form("compare_pe_theta_%i_phi_%i.pdf",i,j));

			h3_dat_pe[i][j]->Write();
			h3_sim_pe[i][j]->Write();
		}
		compare_phi[i] = new TCanvas*[inc_bins_pe];
		for( int j = 0 ; j < inc_bins_pe ; ++j ){
			compare_phi[i][j] = new TCanvas(Form("compare_phi_theta_%i_pe_%i",i,j),Form("compare_phi_theta_%i_pe_%i",i,j),1200,1200);
			compare_phi[i][j]->Divide(1,2);
			TString title = Form("%.2f < Theta_{e} < %.2f , %.2f < P_{e} < %.2f", inc_min_theta + inc_step_theta*i , inc_min_theta + inc_step_theta*(i+1),
												inc_min_pe + inc_step_pe*j , inc_min_pe + inc_step_pe*(j+1) );

			// Make pdf plots of phi comparisons
			compare_phi[i][j]->cd(1);
			drawHists( h3_dat_phi[i][j] , h3_sim_phi[i][j], title, "Phi_{e} [deg.]");
			compare_phi[i][j]->cd(2);
			drawHistsRatio( h3_dat_phi[i][j] , h3_sim_phi[i][j], title, "Phi [deg.]");
			compare_phi[i][j]->Print(Form("compare_phi_theta_%i_pe_%i.pdf",i,j));

			h3_dat_phi[i][j]->Write();
			h3_sim_phi[i][j]->Write();
		}
	}


	TCanvas *** compare_theta = new TCanvas**[inc_bins_pe];
	for( int i = 0 ; i < inc_bins_pe ; ++i ){
		compare_theta[i] = new TCanvas*[inc_bins_phi];
		for( int j = 0 ; j < inc_bins_phi ; ++j ){
			compare_theta[i][j] = new TCanvas(Form("compare_theta_pe_%i_phi_%i",i,j),Form("compare_theta_pe_%i_phi_%i",i,j),1200,1200);
			compare_theta[i][j]->Divide(1,2);
			TString title = Form("%.2f < P_{e} < %.2f , %.2f < Phi_{e} < %.2f", inc_min_pe + inc_step_pe*i , inc_min_pe + inc_step_pe*(i+1),
												inc_min_phi + inc_step_phi*j , inc_min_phi + inc_step_phi*(j+1) );

			// Make pdf plots of pe comparisons
			compare_theta[i][j]->cd(1);
			drawHists( h3_dat_theta[i][j] , h3_sim_theta[i][j], title, "Theta_{e} [deg.]");
			compare_theta[i][j]->cd(2);
			drawHistsRatio( h3_dat_theta[i][j] , h3_sim_theta[i][j], title, "Theta_{e} [deg.]");
			compare_theta[i][j]->Print(Form("compare_theta_pe_%i_phi_%i.pdf",i,j));

			h3_dat_theta[i][j]->Write();
			h3_sim_theta[i][j]->Write();
		}
	}


	outFile->Close();

	inFile_Inc_Dat.Close();
	inFile_Inc_Sim.Close();
	return 1;
}

void fillHist( TH1D** hist_xb, TH1D** hist_Q2, TH1D*** hist_pe, TH1D*** hist_theta, TH1D*** hist_phi, double Xb, double Q2, double Pe, double Theta, double Phi , double weight ){

	if( Q2 > inc_max_Q2		) return;
	if( Q2 < inc_min_Q2		) return;
	if( Xb > inc_max_Xb		) return;
	if( Xb < inc_min_Xb		) return;
	if( Pe > inc_max_pe		) return;
	if( Pe < inc_min_pe		) return;
	if( Theta > inc_max_theta	) return;
	if( Theta < inc_min_theta	) return;
	if( Phi > inc_max_phi		) return;
	if( Phi < inc_min_phi		) return;

	int this_bin_q2 = -1;
	int this_bin_xb = -1;
	int this_bin_pe = -1;
	int this_bin_theta = -1;
	int this_bin_phi = -1;

	
	this_bin_q2 = (int) ((Q2 - inc_min_Q2)/inc_step_Q2);
	this_bin_xb = (int) ((Xb - inc_min_Xb)/inc_step_Xb);
	this_bin_pe = (int) ((Pe - inc_min_pe)/inc_step_pe);
	this_bin_theta = (int) ((Theta - inc_min_theta)/inc_step_theta);
	this_bin_phi = (int) ((Phi - inc_min_phi)/inc_step_phi);

	hist_Q2[this_bin_xb]->Fill( Q2, weight );
	hist_xb[this_bin_q2]->Fill( Xb, weight );
	hist_pe[this_bin_theta][this_bin_phi]->Fill( Pe, weight );
	hist_theta[this_bin_pe][this_bin_phi]->Fill( Theta, weight );
	hist_phi[this_bin_theta][this_bin_pe]->Fill( Phi, weight );

	return;
}

void drawHists( TH1D* data, TH1D* sim, TString title, TString xtitle ){

	double max = 0;
	if( sim->GetMaximum() > max ) max = sim->GetMaximum();
	if( data->GetMaximum() > max) max = data->GetMaximum();

	data->SetTitle(title);
	data->SetStats(0);
	data->SetLineWidth(2);
	data->SetMarkerColor(4);
	data->GetXaxis()->SetTitle(xtitle);
	data->SetMinimum(0);
	data->SetMaximum(max*1.05);
	data->Draw("*P");
	//sim->Scale( h4_dat_Q2[j][k][m]->Integral() /h4_inc_dat_Q2[j][k][m]->Integral() );
	sim->SetLineColor(2);
	sim->Draw("HIST,SAME");

	return;
}

void drawHistsRatio( TH1D* data, TH1D* sim, TString title, TString xtitle ){

	TH1D * data_copy = (TH1D*) data->Clone();
	TH1D * sim_copy = (TH1D*) sim->Clone();

	data_copy->SetTitle(title);
	data_copy->SetStats(0);
	data_copy->SetLineWidth(2);
	data_copy->SetMarkerColor(4);
	data_copy->GetXaxis()->SetTitle(xtitle);
	data_copy->SetMinimum(0);

	data_copy->Divide(sim_copy);
	data_copy->Draw("ep");
	TLine* line = new TLine(data_copy->GetXaxis()->GetBinCenter(1), 1., data_copy->GetXaxis()->GetBinCenter(data_copy->GetXaxis()->GetNbins()), 1.);
	line->SetLineWidth(3);
	line->SetLineColor(2);
	line->Draw("same");

	data_copy->GetYaxis()->SetRangeUser(0,1.1);
	
	data_copy->GetXaxis()->SetTitle(xtitle);
	data_copy->GetYaxis()->SetTitle("Data/Sim");

	return;
}
