#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"

#include "genpart.h"
#include "clashit.h"

#include "bin_edges.h"

using std::cerr;
using std::isfinite;
using std::cout;

const int inc_bins_theta = 7;
const double inc_min_theta = 8.5;
const double inc_max_theta = 29.5;
const double inc_step_theta = (inc_max_theta - inc_min_theta)/inc_bins_theta;

const int inc_bins_phi = 36;
const double inc_min_phi = -180.;
const double inc_max_phi = 180.;
const double inc_step_phi = (inc_max_phi - inc_min_phi)/inc_bins_phi;

const int inc_bins_pe = 9;
const double inc_min_pe = 3.;
const double inc_max_pe = 7.5;
const double inc_step_pe = (inc_max_pe - inc_min_pe)/inc_bins_pe;

void fillHist( TH1D**** hist_u, TH1D**** hist_v, TH1D**** hist_w, TH1D*** hist_phi, double Pe, double Theta, double Phi, int sector, double u, double v, double w, double weight );
void drawHists( TH1D* data, TH1D* sim, TString title, TString xtitle );
void drawHistsRatio( TH1D* data, TH1D* sim, TString title, TString xtitle );

int main( int argc, char** argv){

	if( argc != 3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Inc Dat File] [Inc Sim File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] << "\n";



	// Load files
	TFile inFile_Inc_Dat(argv[1]);
	TFile inFile_Inc_Sim(argv[2]);

	// Load TTrees from files
	TTree * inTree_Inc_Dat = (TTree*) inFile_Inc_Dat.Get("electrons");
	TTree * inTree_Inc_Sim = (TTree*) inFile_Inc_Sim.Get("electrons");

	// Set the input branches for inclusive data
	clashit * inc_dat_eHit 		= new clashit;
	double livetime = -1;
	double gated_charge = -1;
	inTree_Inc_Dat->SetBranchAddress("eHit"		,&inc_dat_eHit		);
	inTree_Inc_Dat->SetBranchAddress("livetime"	,&livetime		);
	inTree_Inc_Dat->SetBranchAddress("gated_charge"	,&gated_charge		);

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

	TH1D **** h3_dat_tof = new TH1D***[6];
	TH1D **** h3_sim_tof = new TH1D***[6];
	for( int i = 0 ; i < 6 ; ++i ){
		h3_dat_tof[i] = new TH1D**[3];
		h3_sim_tof[i] = new TH1D**[3];
		for( int j = 0 ; j < 2 ; ++j ){
			h3_dat_tof[i][j] = new TH1D*[inc_bins_pe];
			h3_sim_tof[i][j] = new TH1D*[inc_bins_pe];
			for( int k = 0 ; k < inc_bins_pe; ++k ){
				h3_dat_tof[i][j][k] = new TH1D(Form("h3_dat_tof_sec_%i_lay_%i_pe_%i",i,j,k),Form("h3_dat_tof_sec_%i_lay_%i_pe_%i",i,j,k),50,-0.5,49.5);
				h3_sim_tof[i][j][k] = new TH1D(Form("h3_sim_tof_sec_%i_lay_%i_pe_%i",i,j,k),Form("h3_sim_tof_sec_%i_lay_%i_pe_%i",i,j,k),50,-0.5,49.5);
			}
		}
	}

	TH2D *** h3_dat_dc 	= new TH2D**[3]; 		// 3 regions
	TH2D *** h3_sim_dc 	= new TH2D**[3]; 		// 3 regions
	for( int i = 0 ; i < 3 ; ++i ){
		h3_dat_dc[i] = new TH2D*[inc_bins_phi]; 	// phi bins
		h3_sim_dc[i] = new TH2D*[inc_bins_phi]; 	// phi bins
		for( int j = 0 ; j < inc_bins_phi ; ++j ){
			h3_dat_dc[i][j] = new TH2D(Form("h3_dat_dc_region_%i_phi_%i",i,j),Form("h3_dat_dc_region_%i_phi_%i",i,j),
					inc_bins_theta*10,inc_min_theta,inc_max_theta,inc_bins_pe*10,inc_min_pe,inc_max_pe);
			h3_sim_dc[i][j] = new TH2D(Form("h3_sim_dc_region_%i_phi_%i",i,j),Form("h3_sim_dc_region_%i_phi_%i",i,j),
					inc_bins_theta*10,inc_min_theta,inc_max_theta,inc_bins_pe*10,inc_min_pe,inc_max_pe);
		}
	}

	TH1D **** h3_dat_pcal_u	= new TH1D***[inc_bins_theta];
	TH1D **** h3_dat_pcal_v	= new TH1D***[inc_bins_theta];
	TH1D **** h3_dat_pcal_w	= new TH1D***[inc_bins_theta];
	TH1D **** h3_sim_pcal_u	= new TH1D***[inc_bins_theta];
	TH1D **** h3_sim_pcal_v	= new TH1D***[inc_bins_theta];
	TH1D **** h3_sim_pcal_w	= new TH1D***[inc_bins_theta];
	for( int i = 0 ; i < inc_bins_theta ; ++i ){
		h3_dat_pcal_u[i] 	= new TH1D**[inc_bins_pe];
		h3_dat_pcal_v[i] 	= new TH1D**[inc_bins_pe];
		h3_dat_pcal_w[i] 	= new TH1D**[inc_bins_pe];
		h3_sim_pcal_u[i] 	= new TH1D**[inc_bins_pe];
		h3_sim_pcal_v[i] 	= new TH1D**[inc_bins_pe];
		h3_sim_pcal_w[i] 	= new TH1D**[inc_bins_pe];
		for( int j = 0 ; j < inc_bins_pe ; ++j ){
			h3_dat_pcal_u[i][j] 	= new TH1D*[6];
			h3_dat_pcal_v[i][j] 	= new TH1D*[6];
			h3_dat_pcal_w[i][j] 	= new TH1D*[6];
			h3_sim_pcal_u[i][j] 	= new TH1D*[6];
			h3_sim_pcal_v[i][j] 	= new TH1D*[6];
			h3_sim_pcal_w[i][j] 	= new TH1D*[6];
			for( int k = 0 ; k < 6 ; ++k ){
				h3_dat_pcal_u[i][j][k] 	= new TH1D(Form("h3_dat_pcal_u_theta_%i_pe_%i_sector_%i",i,j,k),Form("h3_dat_pcal_u_theta_%i_pe_%i_sector_%i",i,j,k),100,0,450);
				h3_dat_pcal_v[i][j][k] 	= new TH1D(Form("h3_dat_pcal_v_theta_%i_pe_%i_sector_%i",i,j,k),Form("h3_dat_pcal_v_theta_%i_pe_%i_sector_%i",i,j,k),100,0,450);
				h3_dat_pcal_w[i][j][k] 	= new TH1D(Form("h3_dat_pcal_w_theta_%i_pe_%i_sector_%i",i,j,k),Form("h3_dat_pcal_w_theta_%i_pe_%i_sector_%i",i,j,k),100,0,450);
				h3_sim_pcal_u[i][j][k] 	= new TH1D(Form("h3_sim_pcal_u_theta_%i_pe_%i_sector_%i",i,j,k),Form("h3_sim_pcal_u_theta_%i_pe_%i_sector_%i",i,j,k),100,0,450);
				h3_sim_pcal_v[i][j][k] 	= new TH1D(Form("h3_sim_pcal_v_theta_%i_pe_%i_sector_%i",i,j,k),Form("h3_sim_pcal_v_theta_%i_pe_%i_sector_%i",i,j,k),100,0,450);
				h3_sim_pcal_w[i][j][k] 	= new TH1D(Form("h3_sim_pcal_w_theta_%i_pe_%i_sector_%i",i,j,k),Form("h3_sim_pcal_w_theta_%i_pe_%i_sector_%i",i,j,k),100,0,450);
			}
		}
	}

	TH1D *** h3_sim_phi	= new TH1D**[inc_bins_theta]; // plot of phi in bins of theta, pe
	TH1D *** h3_dat_phi	= new TH1D**[inc_bins_theta]; // plot of phi in bins of theta, pe
	for( int i = 0 ; i < inc_bins_theta ; ++i ){
		h3_dat_phi[i]	= new TH1D*[inc_bins_pe];
		h3_sim_phi[i]	= new TH1D*[inc_bins_pe];
		for( int j = 0 ; j < inc_bins_pe ; ++j ){
			h3_dat_phi[i][j] = new TH1D(Form("h3_dat_phi_theta_%i_pe_%i",i,j),Form("h3_dat_phi_theta_%i_pe_%i",i,j),inc_bins_phi*2,inc_min_phi,inc_max_phi);
			h3_sim_phi[i][j] = new TH1D(Form("h3_sim_phi_theta_%i_pe_%i",i,j),Form("h3_sim_phi_theta_%i_pe_%i",i,j),inc_bins_phi*2,inc_min_phi,inc_max_phi);
		}
	}

	// Luminosity weights:
	double L_INC_MC = 6.11738e5 / 1E6;		// nb^-1 -> fb^-1 
	//double L_INC_MC_OSG_TEST = 3.077e5 / 1E6;
	double L_INC_MC_OSG_TEST =  6.11738e5 / 1E6 * (10000./25000.);

	inTree_Inc_Dat->GetEntry(inTree_Inc_Dat->GetEntries() - 1 );
	double Q_INC_DAT = gated_charge;		// nC
	double LD2_den = 0.1644; 		// g/cm^3
	double target_L = 5.;			// cm  (-5 to -1)
	double mDeut = 3.3435837724E-24;	// g
	double Coul = 1.60218E-10;		// nC
	double cm2_to_fb = 1E-39;		// 1/cm^2 -> 1/fb
	double num_electron = Q_INC_DAT / Coul;				// # electrons
	double target_density = LD2_den * target_L / mDeut;		// # nucleons / cm^2
	double L_INC_DAT = num_electron * target_density * cm2_to_fb;	// fb^-1

	// Loop over the inclusive data file:
	for( int event = 0 ; event < inTree_Inc_Dat->GetEntries() ; ++event ){
		inc_dat_eHit	->Clear();
		livetime = -1;

		inTree_Inc_Dat->GetEntry(event);
		if( livetime < 0 ) continue;

		double pe	= inc_dat_eHit->getMomentum();
		double theta	= inc_dat_eHit->getTheta() * 180./M_PI;
		double phi	= inc_dat_eHit->getPhi() * 180./M_PI;
		int sector	= inc_dat_eHit->getSector();
		double u	= inc_dat_eHit->getU();
		double v	= inc_dat_eHit->getV();
		double w	= inc_dat_eHit->getW();

		std::vector<int>	scint_sec = inc_dat_eHit->getScint_sector();
		std::vector<int>	scint_lay = inc_dat_eHit->getScint_layer();
		std::vector<int>	scint_com = inc_dat_eHit->getScint_component();
		int this_bin_pe = (int) ((pe - inc_min_pe)/inc_step_pe);
		for( int scint_hit = 0 ; scint_hit < scint_sec.size() ; ++scint_hit ){
			if( pe > inc_max_pe		) continue;
			if( pe < inc_min_pe		) continue;
			h3_dat_tof[scint_sec[scint_hit]-1][scint_lay[scint_hit]-1][this_bin_pe]->Fill( scint_com[scint_hit] ,  1./(L_INC_DAT*livetime));
		}

		double dc_x1 = inc_dat_eHit->getDC_x1();
		double dc_x2 = inc_dat_eHit->getDC_x2();
		double dc_x3 = inc_dat_eHit->getDC_x3();
		double dc_y1 = inc_dat_eHit->getDC_y1();
		double dc_y2 = inc_dat_eHit->getDC_y2();
		double dc_y3 = inc_dat_eHit->getDC_y3();
		double phi_1 = atan2( dc_y1 , dc_x1 ) * 180./M_PI;
		double phi_2 = atan2( dc_y2 , dc_x2 ) * 180./M_PI;
		double phi_3 = atan2( dc_y3 , dc_x3 ) * 180./M_PI;
		int this_bin_phi;

		if( phi_1 > inc_min_phi && phi_1 < inc_max_phi ){
			this_bin_phi = (int) ((phi_1 - inc_min_phi)/inc_step_phi);
			h3_dat_dc[0][this_bin_phi] -> Fill( theta , pe , 1./(L_INC_DAT*livetime));
		}
		if( phi_2 > inc_min_phi && phi_2 < inc_max_phi ){
			this_bin_phi = (int) ((phi_2 - inc_min_phi)/inc_step_phi);
			h3_dat_dc[1][this_bin_phi] -> Fill( theta , pe , 1./(L_INC_DAT*livetime));
		}
		if( phi_3 > inc_min_phi && phi_3 < inc_max_phi ){
			this_bin_phi = (int) ((phi_3 - inc_min_phi)/inc_step_phi);
			h3_dat_dc[2][this_bin_phi] -> Fill( theta , pe , 1./(L_INC_DAT*livetime));
		}
	
		fillHist( h3_dat_pcal_u, h3_dat_pcal_v, h3_dat_pcal_w , h3_dat_phi ,  pe, theta, phi, sector, u, v, w, 1./(L_INC_DAT*livetime));

	}

	// Loop over the inclusive simulation file:
	for( int event = 0 ; event < inTree_Inc_Sim->GetEntries() ; ++event ){
		inc_sim_eHit	->Clear();

		inTree_Inc_Sim->GetEntry(event);
		
		double pe	= inc_sim_eHit->getMomentum();
		double theta	= inc_sim_eHit->getTheta() * 180./M_PI;
		double phi	= inc_sim_eHit->getPhi() * 180./M_PI;
		int sector	= inc_sim_eHit->getSector();
		double u	= inc_sim_eHit->getU();
		double v	= inc_sim_eHit->getV();
		double w	= inc_sim_eHit->getW();

		std::vector<int>	scint_sec = inc_sim_eHit->getScint_sector();
		std::vector<int>	scint_lay = inc_sim_eHit->getScint_layer();
		std::vector<int>	scint_com = inc_sim_eHit->getScint_component();
		int this_bin_pe = (int) ((pe - inc_min_pe)/inc_step_pe);
		for( int scint_hit = 0 ; scint_hit < scint_sec.size() ; ++scint_hit ){
			if( pe > inc_max_pe		) continue;
			if( pe < inc_min_pe		) continue;
			h3_sim_tof[scint_sec[scint_hit]-1][scint_lay[scint_hit]-1][this_bin_pe]->Fill( scint_com[scint_hit] , 1./(L_INC_MC_OSG_TEST));
		}

		double dc_x1 = inc_sim_eHit->getDC_x1();
		double dc_x2 = inc_sim_eHit->getDC_x2();
		double dc_x3 = inc_sim_eHit->getDC_x3();
		double dc_y1 = inc_sim_eHit->getDC_y1();
		double dc_y2 = inc_sim_eHit->getDC_y2();
		double dc_y3 = inc_sim_eHit->getDC_y3();
		double phi_1 = atan2( dc_y1 , dc_x1 ) * 180./M_PI;
		double phi_2 = atan2( dc_y2 , dc_x2 ) * 180./M_PI;
		double phi_3 = atan2( dc_y3 , dc_x3 ) * 180./M_PI;
		int this_bin_phi;

		if( phi_1 > inc_min_phi && phi_1 < inc_max_phi ){
			this_bin_phi = (int) ((phi_1 - inc_min_phi)/inc_step_phi);
			h3_sim_dc[0][this_bin_phi] -> Fill( theta , pe , 1./(L_INC_MC_OSG_TEST));
		}
		if( phi_2 > inc_min_phi && phi_2 < inc_max_phi ){
			this_bin_phi = (int) ((phi_2 - inc_min_phi)/inc_step_phi);
			h3_sim_dc[1][this_bin_phi] -> Fill( theta , pe , 1./(L_INC_MC_OSG_TEST));
		}
		if( phi_3 > inc_min_phi && phi_3 < inc_max_phi ){
			this_bin_phi = (int) ((phi_3 - inc_min_phi)/inc_step_phi);
			h3_sim_dc[2][this_bin_phi] -> Fill( theta , pe , 1./(L_INC_MC_OSG_TEST));
		}
	
		fillHist( h3_sim_pcal_u, h3_sim_pcal_v, h3_sim_pcal_w , h3_sim_phi , pe, theta, phi, sector, u, v, w, 1./(L_INC_MC_OSG_TEST));
	}
	
	outFile->cd();	

	TCanvas *** compare_phi = new TCanvas**[inc_bins_theta];
	for( int i = 0 ; i < inc_bins_theta ; ++i ){
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

	TCanvas **** compare_tof = new TCanvas***[6];
	for( int i = 0 ; i < 6 ; ++i ){
		compare_tof[i] = new TCanvas**[3];
		for( int j = 0 ; j < 2 ; ++j ){
			compare_tof[i][j] = new TCanvas*[inc_bins_pe];
			for( int k = 0 ; k < inc_bins_pe; ++k){

				compare_tof[i][j][k] = new TCanvas(Form("compare_tof_sec_%i_layer_%i_pe_%i",i,j,k),Form("compare_tof_sec_%i_layer_%i_pe_%i",i,j,k),1200,1200);
				compare_tof[i][j][k] -> Divide(1,2);
				TString title = Form("Sector %i , Layer %i , %.2f < P_{e} < %.2f", i+1, j+1, inc_min_pe + inc_step_pe*k , inc_min_pe + inc_step_pe*(k+1));
				compare_tof[i][j][k] -> cd(1);
				drawHists( h3_dat_tof[i][j][k] , h3_sim_tof[i][j][k], title, "ToF Bar ID");
				compare_tof[i][j][k] -> cd(2);
				drawHistsRatio( h3_dat_tof[i][j][k] , h3_sim_tof[i][j][k], title, "ToF Bar ID");
				compare_tof[i][j][k]->Print(Form("compare_tof_sector_%i_layer_%i_pe_%i.pdf",i+1,j+1,k));
			}
		}
	}

	TCanvas *** compare_dc = new TCanvas**[3];
	for( int i = 0 ; i < 3 ; ++i){
		compare_dc[i] = new TCanvas*[inc_bins_phi];
		for( int j = 0 ; j < inc_bins_phi ; ++j ){
			compare_dc[i][j] = new TCanvas(Form("compare_dc_region_%i_phi_%i",i,j),Form("compare_dc_region_%i_phi_%i",i,j),1200,1200);
			compare_dc[i][j]->Divide(2,2);

			TString title = Form("%.2f < Region %i Phi_{e} < %.2f", inc_min_phi + inc_step_phi*j, (i+1) , inc_min_phi + inc_step_phi*(j+1));

			compare_dc[i][j]->cd(1);
			h3_dat_dc[i][j]->SetTitle(title);
			h3_dat_dc[i][j]->SetStats(0);
			h3_dat_dc[i][j]->GetXaxis()->SetTitle("Theta_{e} [deg.]");
			h3_dat_dc[i][j]->GetYaxis()->SetTitle("p_{e} [GeV/c]");
			h3_dat_dc[i][j]->Draw("colz");

			compare_dc[i][j]->cd(2);
			h3_sim_dc[i][j]->SetTitle(title);
			h3_sim_dc[i][j]->SetStats(0);
			h3_sim_dc[i][j]->GetXaxis()->SetTitle("Theta_{e} [deg.]");
			h3_sim_dc[i][j]->GetYaxis()->SetTitle("p_{e} [GeV/c]");
			h3_sim_dc[i][j]->Draw("colz");

			compare_dc[i][j]->cd(3);
			TH2D * data_copy = (TH2D*) h3_dat_dc[i][j]->Clone();
			TH2D * sim_copy = (TH2D*) h3_sim_dc[i][j]->Clone();
			data_copy->Divide(sim_copy);
			title = title + " , Data/Sim Ratio";
			data_copy->SetTitle(title);
			data_copy->SetStats(0);
			data_copy->SetMaximum(3);
			data_copy->SetMinimum(0.01);
			data_copy->Draw("colz");
			
			compare_dc[i][j]->Print(Form("compare_dc_region_%i_phi_%i.pdf",(i+1),j));


		}
	}

	TCanvas **** compare_u = new TCanvas***[inc_bins_theta];
	TCanvas **** compare_v = new TCanvas***[inc_bins_theta];
	TCanvas **** compare_w = new TCanvas***[inc_bins_theta];
	for( int i = 0 ; i < inc_bins_theta ; ++i ){
		compare_u[i] = new TCanvas**[inc_bins_pe];
		compare_v[i] = new TCanvas**[inc_bins_pe];
		compare_w[i] = new TCanvas**[inc_bins_pe];
		for( int j = 0 ; j < inc_bins_pe ; ++j ){
			compare_u[i][j] = new TCanvas*[6];
			compare_v[i][j] = new TCanvas*[6];
			compare_w[i][j] = new TCanvas*[6];
			for( int k = 0 ; k < 6 ; ++k ){
				compare_u[i][j][k] = new TCanvas(Form("compare_u_theta_%i_pe_%i_sector_%i",i,j,k),Form("compare_u_theta_%i_pe_%i_sector_%i",i,j,k),1200,1200);
				compare_v[i][j][k] = new TCanvas(Form("compare_v_theta_%i_pe_%i_sector_%i",i,j,k),Form("compare_v_theta_%i_pe_%i_sector_%i",i,j,k),1200,1200);
				compare_w[i][j][k] = new TCanvas(Form("compare_w_theta_%i_pe_%i_sector_%i",i,j,k),Form("compare_w_theta_%i_pe_%i_sector_%i",i,j,k),1200,1200);

				compare_u[i][j][k]->Divide(1,2);
				compare_v[i][j][k]->Divide(1,2);
				compare_w[i][j][k]->Divide(1,2);
				TString title = Form("Sector %i , %.2f < Theta_{e} < %.2f , %.2f < P_{e} < %.2f", (k+1), inc_min_theta + inc_step_theta*i , inc_min_theta + inc_step_theta*(i+1),
													inc_min_pe + inc_step_pe*j , inc_min_pe + inc_step_pe*(j+1) );

				// Make pdf plots of phi comparisons
				compare_u[i][j][k]->cd(1);
				drawHists( h3_dat_pcal_u[i][j][k] , h3_sim_pcal_u[i][j][k], title, "PCal_{u} [cm]");
				compare_u[i][j][k]->cd(2);
				drawHistsRatio( h3_dat_pcal_u[i][j][k] , h3_sim_pcal_u[i][j][k], title, "PCal_{u} [cm]");
				compare_u[i][j][k]->Print(Form("compare_u_theta_%i_pe_%i_sector_%i.pdf",i,j,(k+1)));

				// Make pdf plots of phi comparisons
				compare_v[i][j][k]->cd(1);
				drawHists( h3_dat_pcal_v[i][j][k] , h3_sim_pcal_v[i][j][k], title, "PCal_{v} [cm]");
				compare_v[i][j][k]->cd(2);
				drawHistsRatio( h3_dat_pcal_v[i][j][k] , h3_sim_pcal_v[i][j][k], title, "PCal_{v} [cm]");
				compare_v[i][j][k]->Print(Form("compare_v_theta_%i_pe_%i_sector_%i.pdf",i,j,(k+1)));

				// Make pdf plots of phi comparisons
				compare_w[i][j][k]->cd(1);
				drawHists( h3_dat_pcal_w[i][j][k] , h3_sim_pcal_w[i][j][k], title, "PCal_{w} [cm]");
				compare_w[i][j][k]->cd(2);
				drawHistsRatio( h3_dat_pcal_w[i][j][k] , h3_sim_pcal_w[i][j][k], title, "PCal_{w} [cm]");
				compare_w[i][j][k]->Print(Form("compare_w_theta_%i_pe_%i_sector_%i.pdf",i,j,(k+1)));
			}

		}
	}

	outFile->Close();

	inFile_Inc_Dat.Close();
	inFile_Inc_Sim.Close();
	return 1;
}

void fillHist( TH1D**** hist_u, TH1D**** hist_v, TH1D**** hist_w, TH1D*** hist_phi, double Pe, double Theta, double Phi, int sector, double u, double v, double w, double weight ){

	if( Pe > inc_max_pe		) return;
	if( Pe < inc_min_pe		) return;
	if( Theta > inc_max_theta	) return;
	if( Theta < inc_min_theta	) return;
	if( Phi > inc_max_phi		) return;
	if( Phi < inc_min_phi		) return;

	int this_bin_pe = -1;
	int this_bin_theta = -1;
	int this_bin_phi = -1;

	
	this_bin_pe = (int) ((Pe - inc_min_pe)/inc_step_pe);
	this_bin_theta = (int) ((Theta - inc_min_theta)/inc_step_theta);
	this_bin_phi = (int) ((Phi - inc_min_phi)/inc_step_phi);

	hist_u[this_bin_theta][this_bin_pe][sector-1]->Fill( u, weight );
	hist_v[this_bin_theta][this_bin_pe][sector-1]->Fill( v, weight );
	hist_w[this_bin_theta][this_bin_pe][sector-1]->Fill( w, weight );

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

	data_copy->GetYaxis()->SetRangeUser(0,1);
	
	data_copy->GetXaxis()->SetTitle(xtitle);
	data_copy->GetYaxis()->SetTitle("Data/Sim");

	return;
}
