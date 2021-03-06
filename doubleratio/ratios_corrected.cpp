#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"

#include "../include/bin_edges.h"


void ratios_corrected(TString path){


	//Create Eq12 and then Eq14 results from BAND writeup

	//Uncertainty Array for corrections coeffficents. Set all to 0 at the moment
	double coeff_errors_tagged[bins_As];
	double coeff_errors_inclusive[bins_Xb];
	for (int Asbin = 0; Asbin < bins_As; Asbin++ ) {
		coeff_errors_tagged[Asbin] = 0;
	}
	for (int Xbbin = 0; Xbbin < bins_Xb; Xbbin++ ) {
		coeff_errors_inclusive[Xbbin] = 0;
	}

	//Files for Tagged and Inclusive data and simulation
	TFile * inFileDatTag = new TFile(path+"/yield_tagged.root");
	TFile * inFileDatInc = new TFile(path+"/yield_inclusive.root");
	//TFile * inFileSimInc = new TFile(inFileSimIncName);

	//Files for Bin Migration correction
	TFile * inFileCorrBMTag = new TFile(path+"/migcorrection_tagged.root");
	TFile * inFileCorrBMInc = new TFile(path+"/migcorrection_inclusive.root");

	//Files for Radiative corrections
	TFile * inFileCorrRadTag = new TFile(path+"/radcorrection_tagged.root");
	TFile * inFileCorrRadInc = new TFile(path+"/radcorrection_inclusive.root");

	//Files for Acceptance/Phase Space corrections
	TFile * inFileCorrAccTag = new TFile(path+"/acceptance_tagged.root");
	TFile * inFileCorrAccInc = new TFile(path+"/acceptance_inclusive.root");

	//Print file content for checking
	//	inFileDatTag->ls();
	//	inFileDatInc->ls();
	//	inFileCorrBMTag->ls();
	//	inFileCorrBMInc->ls();
	//	inFileCorrRadTag->ls();
	//	inFileCorrRadInc->ls();
	//	inFileCorrAccTag->ls();
	//	inFileCorrAccInc->ls();

	//Histos for Yields
	TH1D **** h4_data_tagged_yield = new TH1D***[bins_Q2];
	TH1D **** h4_simulation_tagged_yield = new TH1D***[bins_Q2];
	TH1D ** h2_data_inclusive_yield = new TH1D*[bins_Q2];
	TH1D ** h2_simulation_inclusive_yield = new TH1D*[bins_Q2];
	//Histos for Bin Migration Corrections
	TH1D **** h4_tagged_binmigration = new TH1D***[bins_Q2];
	TH1D ** h2_inclusive_binmigration = new TH1D*[bins_Q2];
	//Histos for Radiative Corrections
	TH1D **** h4_tagged_radiative = new TH1D***[bins_Q2];
	TH1D ** h2_inclusive_radiative = new TH1D*[bins_Q2];
	//Histos for Acceptance Corrections
	TH1D **** h4_tagged_acceptance = new TH1D***[bins_Q2];
	TH1D ** h2_inclusive_acceptance = new TH1D*[bins_Q2];

	//Histos for Corrected Data
	TH1D **** h4_data_corrected_tagged = new TH1D***[bins_Q2];
	TH1D **** h4_sim_corrected_tagged = new TH1D***[bins_Q2];
	TH1D ** h2_data_corrected_inclusive = new TH1D*[bins_Q2];
	TH1D ** h2_sim_corrected_inclusive = new TH1D*[bins_Q2];

	for( int i = 0 ; i < bins_Q2 ; ++i ){ // bins in Q2
		//Getting Yields
		h2_data_inclusive_yield[i] = (TH1D*) inFileDatInc->Get(Form("dat_yield_Xb_Q2_%i",i)); //x-axis is xB
		h2_simulation_inclusive_yield[i] = (TH1D*) inFileDatInc->Get(Form("sim_yield_Xb_Q2_%i",i)); //x-axis is xB
		h4_data_tagged_yield[i]	= new TH1D**[bins_Pt];
		h4_simulation_tagged_yield[i]	= new TH1D**[bins_Pt];
		h2_data_corrected_inclusive[i] = (TH1D*) h2_data_inclusive_yield[i]->Clone(); //x-axis is xB
		h2_sim_corrected_inclusive[i] = (TH1D*) h2_simulation_inclusive_yield[i]->Clone(); //x-axis is xB
		h2_data_corrected_inclusive[i]->SetTitle(Form("dat_corrected_Xb_Q2_%i",i));
		h2_sim_corrected_inclusive[i]->SetTitle(Form("sim_corrected_Xb_Q2_%i",i));
		h2_data_corrected_inclusive[i]->SetName(Form("dat_corrected_Xb_Q2_%i",i));
		h2_sim_corrected_inclusive[i]->SetName(Form("sim_corrected_Xb_Q2_%i",i));
		h4_data_corrected_tagged[i]	= new TH1D**[bins_Pt];
		h4_sim_corrected_tagged[i]	= new TH1D**[bins_Pt];
		//Getting Bin migration
		h2_inclusive_binmigration[i] = (TH1D*) inFileCorrBMInc->Get(Form("migcorr_Xb_Q2_%i",i)); //x-axis is xB
		h4_tagged_binmigration[i]	= new TH1D**[bins_Pt];
		//Getting Radiative corrections
		h2_inclusive_radiative[i] = (TH1D*) inFileCorrRadInc->Get(Form("radcorr_Xb_Q2_%i",i)); //x-axis is xB
		h4_tagged_radiative[i]	= new TH1D**[bins_Pt];
		//Getting Acceptance corrections
		h2_inclusive_acceptance[i] = (TH1D*) inFileCorrAccInc->Get(Form("acceptance_inclusive_xb_Q2_%i",i)); //x-axis is xB
		h4_tagged_acceptance[i]	= new TH1D**[bins_Pt];

		//Set Errors for coeffficents to 0
		h2_inclusive_binmigration[i]->SetError(coeff_errors_inclusive);
		h2_inclusive_radiative[i]->SetError(coeff_errors_inclusive);
		h2_inclusive_acceptance[i]->SetError(coeff_errors_inclusive);

		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			h4_data_tagged_yield[i][j] 	= new TH1D*[bins_Xb];
			h4_simulation_tagged_yield[i][j] 	= new TH1D*[bins_Xb];
			h4_data_corrected_tagged[i][j]	= new TH1D*[bins_Xb];
			h4_sim_corrected_tagged[i][j]	= new TH1D*[bins_Xb];
			h4_tagged_binmigration[i][j]	= new TH1D*[bins_Xb];
			h4_tagged_radiative[i][j]	= new TH1D*[bins_Xb];
			h4_tagged_acceptance[i][j]	= new TH1D*[bins_Xb];

			for( int k = 0 ; k < bins_Xb ; ++k ){ // bins in Xb
				//Getting Yields
				h4_data_tagged_yield[i][j][k] = (TH1D*) inFileDatTag->Get(Form("dat_yield_as_Q2_%i_Pt_%i_Xb_%i",i,j,k)); //x axis is alphaS
				h4_simulation_tagged_yield[i][j][k] = (TH1D*) inFileDatTag->Get(Form("sim_yield_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));  //x axis is alphaS
				h4_data_corrected_tagged[i][j][k] = (TH1D*) h4_data_tagged_yield[i][j][k]->Clone();
				h4_sim_corrected_tagged[i][j][k] = (TH1D*) h4_simulation_tagged_yield[i][j][k]->Clone();
				h4_data_corrected_tagged[i][j][k]->SetTitle(Form("dat_corrected_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_sim_corrected_tagged[i][j][k]->SetTitle(Form("sim_corrected_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_data_corrected_tagged[i][j][k]->SetName(Form("dat_corrected_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_sim_corrected_tagged[i][j][k]->SetName(Form("sim_corrected_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				//Getting Corrections
				h4_tagged_binmigration[i][j][k] = (TH1D*) inFileCorrBMTag->Get(Form("migcorr_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));  //x axis is alphaS
				h4_tagged_radiative[i][j][k] = (TH1D*) inFileCorrRadTag->Get(Form("radcorr_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));  //x axis is alphaS
				h4_tagged_acceptance[i][j][k] = (TH1D*) inFileCorrAccTag->Get(Form("acceptance_tagged_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));  //x axis is alphaS

				//Set Errors for coeffficents to 0
				h4_tagged_binmigration[i][j][k]->SetError(coeff_errors_tagged);
				h4_tagged_radiative[i][j][k]->SetError(coeff_errors_tagged);
				h4_tagged_acceptance[i][j][k]->SetError(coeff_errors_tagged);

				//just for debug
				//			cout << " for i j k " << i << " " << j << " " << k << " " <<  h4_tagged_binmigration[i][j][k]->GetBinContent(2) << " " <<  h4_tagged_acceptance[i][j][k]->GetBinContent(2) << " "
				//			<<  h4_tagged_radiative[i][j][k]->GetBinContent(2) << " " << endl;
			}
		}
		//DEBUG:
		//	cout << "for i : " << h2_data_inclusive_yield[i]->GetNbinsX() << endl;

	}

	//Now we have all data and correction factors loaded
	//Create new histograms for data and inclusive divided by all corrections factors, hence creating nominator and denominator separately for Eq 12.

	for( int i = 0 ; i < bins_Q2 ; ++i ){ // bins in Q2
		//Divide by Bin Migration Correction
		h2_data_corrected_inclusive[i]->Divide(h2_inclusive_binmigration[i]);
		h2_sim_corrected_inclusive[i]->Divide(h2_inclusive_binmigration[i]);
		//Divide by Radiative Correction
		h2_data_corrected_inclusive[i]->Divide(h2_inclusive_radiative[i]);
		h2_sim_corrected_inclusive[i]->Divide(h2_inclusive_radiative[i]);
		//Divide by Phase Space Correction
		h2_data_corrected_inclusive[i]->Divide(h2_inclusive_acceptance[i]);
		h2_sim_corrected_inclusive[i]->Divide(h2_inclusive_acceptance[i]);

		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			for( int k = 0 ; k < bins_Xb ; ++k ){ // bins in Xb
				//Divide by Bin Migration Correction
				h4_data_corrected_tagged[i][j][k]->Divide(h4_tagged_binmigration[i][j][k]);
				h4_sim_corrected_tagged[i][j][k]->Divide(h4_tagged_binmigration[i][j][k]);
				//Divide by Radiative Correction
				h4_data_corrected_tagged[i][j][k]->Divide(h4_tagged_radiative[i][j][k]);
				h4_sim_corrected_tagged[i][j][k]->Divide(h4_tagged_radiative[i][j][k]);
				//Divide by Phase Space Correction
				h4_data_corrected_tagged[i][j][k]->Divide(h4_tagged_acceptance[i][j][k]);
				h4_sim_corrected_tagged[i][j][k]->Divide(h4_tagged_acceptance[i][j][k]);

			}
		}
	}
	//Histos for Eq 12
	TH1D **** h4_data_singleratio = new TH1D***[bins_Q2];
	TH1D **** h4_sim_singleratio = new TH1D***[bins_Q2];
	TCanvas * c_data_corrected[bins_Q2][bins_Pt][bins_Xb];

	//Show intermediate histograms for nominator and denominator. Create single ratio histograms
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins

		h4_data_singleratio[i]	= new TH1D**[bins_Pt];
		h4_sim_singleratio[i]	= new TH1D**[bins_Pt];

		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			h4_data_singleratio[i][j]	= new TH1D*[bins_Xb];
			h4_sim_singleratio[i][j]	= new TH1D*[bins_Xb];
			for( int k = 0 ; k < bins_Xb ; ++k ){ // loop over Xb bina

				c_data_corrected[i][j][k]= new TCanvas(Form("test_%i_%i_%i",i,j,k),"",1200,1200);	

				h4_data_singleratio[i][j][k]= (TH1D*) h4_data_corrected_tagged[i][j][k]->Clone();
				h4_sim_singleratio[i][j][k] = (TH1D*) h4_sim_corrected_tagged[i][j][k]->Clone();
				h4_data_singleratio[i][j][k]->SetTitle(Form("dat_singleratio_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_sim_singleratio[i][j][k]->SetTitle(Form("sim_singleratio_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_data_singleratio[i][j][k]->SetName(Form("dat_singleratio_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_sim_singleratio[i][j][k]->SetName(Form("sim_singleratio_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));

				c_data_corrected[i][j][k]->cd();
				h4_data_corrected_tagged[i][j][k]->SetLineWidth(3);
				h4_data_corrected_tagged[i][j][k]->SetMarkerColor(4);
				//		h4_data_corrected_tagged[i][j][k]->SetMaximum(2);
				h4_data_corrected_tagged[i][j][k]->SetMinimum(0);
				h4_data_corrected_tagged[i][j][k]->Draw("*P");
				h4_data_corrected_tagged[i][j][k]->GetXaxis()->SetTitle("a_{S}");

				c_data_corrected[i][j][k]->Print(Form("data_tagged_corrected_Q2_%i_Pt_%i_Xb_%i.pdf",i,j,k));

			}
		}
	}
	//Canvas for inclusive data corrected
	TCanvas *c_Q2_data_inc_corrected = new TCanvas("c_Q2_data_inc_corrected","",1200,1200);
	c_Q2_data_inc_corrected->Divide(1,bins_Q2);
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		c_Q2_data_inc_corrected->cd(i+1);
		h2_data_corrected_inclusive[i]->SetLineWidth(3);
		h2_data_corrected_inclusive[i]->SetMarkerColor(4);
		//		h2_data_corrected_inclusive[i]->SetMaximum(2);
		h2_data_corrected_inclusive[i]->SetMinimum(0);
		h2_data_corrected_inclusive[i]->Draw("*P");
		h2_data_corrected_inclusive[i]->GetXaxis()->SetTitle("x_{B}");
	}
	c_Q2_data_inc_corrected->Print("data_inclusive_corrected.pdf");

	//Now use corrected histograms to fill single ratios from Equation 12
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bin
		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			for( int k = 0 ; k < bins_Xb ; ++k ){ // bins in Xb
				if (h2_data_corrected_inclusive[i]->GetBinContent(k+1)!=0) {
					h4_data_singleratio[i][j][k]->Scale(1./h2_data_corrected_inclusive[i]->GetBinContent(k+1)); //x-axis is alphaS
					h4_sim_singleratio[i][j][k]->Scale(1./h2_data_corrected_inclusive[i]->GetBinContent(k+1));
				}
				else {
					h4_data_singleratio[i][j][k]->Scale(0);
					h4_sim_singleratio[i][j][k]->Scale(0);
				}
			}
		}
	}

	//Histos for Eq 14
	TH1D **** h4_data_doubleratio = new TH1D***[bins_Q2];
	TH1D **** h4_sim_doubleratio = new TH1D***[bins_Q2];

	//Create Canvas for single ratio plots and Clone data for later double ratio plot
	TCanvas * c_data_singleratio[bins_Q2][bins_Pt][bins_Xb];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins

		h4_data_doubleratio[i]	= new TH1D**[bins_Pt];
		h4_sim_doubleratio[i]	= new TH1D**[bins_Pt];

		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			h4_data_doubleratio[i][j]	= new TH1D*[bins_Xb];
			h4_sim_doubleratio[i][j]	= new TH1D*[bins_Xb];
			for( int k = 0 ; k < bins_Xb ; ++k ){ // loop over Xb bina

				c_data_singleratio[i][j][k] = new TCanvas(Form("c_singleratio_%i_%i_%i",i,j,k),"",1200,1200);

				h4_data_doubleratio[i][j][k]= (TH1D*) h4_data_singleratio[i][j][k]->Clone();
				h4_sim_doubleratio[i][j][k] = (TH1D*) h4_sim_singleratio[i][j][k]->Clone();
				h4_data_doubleratio[i][j][k]->SetTitle(Form("dat_doubleratio_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_sim_doubleratio[i][j][k]->SetTitle(Form("sim_doubleratio_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_data_doubleratio[i][j][k]->SetName(Form("dat_doubleratio_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_sim_doubleratio[i][j][k]->SetName(Form("sim_doubleratio_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_data_doubleratio[i][j][k]->Sumw2();
				h4_sim_doubleratio[i][j][k]->Sumw2();

				c_data_singleratio[i][j][k]->cd();
				h4_data_singleratio[i][j][k]->SetLineWidth(3);
				h4_data_singleratio[i][j][k]->SetMarkerColor(4);
				//	h4_data_singleratio[i][j][k]->SetMaximum(2);
				h4_data_singleratio[i][j][k]->SetMinimum(0);
				h4_sim_singleratio[i][j][k]->SetLineWidth(3);
				h4_sim_singleratio[i][j][k]->SetMarkerColor(2);
				//	h4_sim_singleratio[i][j][k]->Draw("*P");
				h4_data_singleratio[i][j][k]->Draw("*P");
				h4_data_singleratio[i][j][k]->GetXaxis()->SetTitle("a_{S}");

				c_data_singleratio[i][j][k]->Print(Form("data_tagged_singleratio_Q2_%i_Pt_%i_Xb_%i.pdf",i,j,k));
			}
		}
	}

	//Now make double ratio by dividing each single ratio histogram by the first xB bin hence h4_data_singleratio[i][j][0]
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bin
		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			for( int k = 0 ; k < bins_Xb ; ++k ){ // loop over Xb bina
				h4_data_doubleratio[i][j][k]->Divide(	h4_data_singleratio[i][j][0]);
				h4_sim_doubleratio[i][j][k]->Divide(	h4_sim_singleratio[i][j][0]);
				//cout << " for i j " << i << " " << j << " "  << " " <<  1./h2_data_corrected_inclusive[i]->GetBinContent(k+1) << endl;
			}
		}
	}

	//Histos for data to simulation quadro-ratio of doble ratios :)
	TH1D **** h4_data_sim_doubleratios = new TH1D***[bins_Q2];
	//Canvas to plot double ratios overlayed data and sim
	TCanvas * c_data_doubleratio[bins_Q2][bins_Pt][bins_Xb];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins

		h4_data_sim_doubleratios[i] = new TH1D**[bins_Pt];
		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			h4_data_sim_doubleratios[i][j] = new TH1D*[bins_Xb];
			for( int k = 0 ; k < bins_Xb ; ++k ){ // loop over Xb bina

				c_data_doubleratio[i][j][k] = new TCanvas(Form("c_doubleratio_%i_%i_%i",i,j,k),"",1200,1200);

				h4_data_sim_doubleratios[i][j][k] =(TH1D*) h4_data_doubleratio[i][j][k]->Clone();
				h4_data_sim_doubleratios[i][j][k]->SetTitle(Form("data_sim_quadratio_as_Q2_%i_Pt_%i",i,j));
				h4_data_sim_doubleratios[i][j][k]->SetName(Form("data_sim_quadratio_as_Q2_%i_Pt_%i",i,j));

				c_data_doubleratio[i][j][k]->cd();
				h4_data_doubleratio[i][j][k]->SetLineWidth(3);
				h4_data_doubleratio[i][j][k]->SetMarkerColor(4);
				h4_data_doubleratio[i][j][k]->SetMaximum(2);
				h4_data_doubleratio[i][j][k]->SetMinimum(0);
				h4_sim_doubleratio[i][j][k]->SetMaximum(2);
				h4_sim_doubleratio[i][j][k]->SetMinimum(0);
				h4_sim_doubleratio[i][j][k]->SetLineWidth(3);
				h4_sim_doubleratio[i][j][k]->SetMarkerColor(2);
				h4_sim_doubleratio[i][j][k]->Draw("E1");
				h4_data_doubleratio[i][j][k]->Draw("SAMEP");
				h4_sim_doubleratio[i][j][k]->GetXaxis()->SetTitle("a_{S}");



				c_data_doubleratio[i][j][k]->Print(Form("data_tagged_doubleratio_Q2_%i_Pt_%i_Xb_%i.pdf",i,j,k));

			}
		}
	}

	//Now make data to sim double ratio by dividing each double ratio histogram. h4_data_sim_doubleratios[i][j][k] has cloned data from data double ratio
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bin
		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			for( int k = 0 ; k < bins_Xb ; ++k ){ // loop over Xb bina
				h4_data_sim_doubleratios[i][j][k]->Divide(h4_sim_doubleratio[i][j][k]);
			}
		}
	}

	//Canvas to plot quad ratios of data/sim doubleratios
	TCanvas * c_data_quadratio[bins_Q2][bins_Pt][bins_Xb];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins

		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			for( int k = 0 ; k < bins_Xb ; ++k ){ // loop over Xb bina

				c_data_quadratio[i][j][k] = new TCanvas(Form("c_quadratio_%i_%i_%i",i,j,k),"",1200,1200);


				c_data_quadratio[i][j][k]->cd();
				h4_data_sim_doubleratios[i][j][k]->SetLineWidth(3);
				h4_data_sim_doubleratios[i][j][k]->SetMarkerColor(4);
				h4_data_sim_doubleratios[i][j][k]->SetMaximum(2);
				h4_data_sim_doubleratios[i][j][k]->SetMinimum(0);
				h4_data_sim_doubleratios[i][j][k]->Draw("E1");
				h4_data_sim_doubleratios[i][j][k]->GetXaxis()->SetTitle("a_{S}");

				c_data_quadratio[i][j][k]->Print(Form("data_tagged_quadratio_Q2_%i_Pt_%i_Xb_%i.pdf",i,j,k));
			}
		}
	}

	TH1D **** h4_data_doubleratios_fixedAs = new TH1D***[bins_Q2];
	TH1D **** h4_sim_doubleratios_fixedAs = new TH1D***[bins_Q2];
	TCanvas * c_data_doubleratios_fixedAs[bins_Q2][bins_Pt][bins_As];
	// Now I want to take the double ratios and loop over Q2, Pt, and take bins of fixed alpha and plot as a function of xB
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		h4_data_doubleratios_fixedAs[i] = new TH1D**[bins_Pt];
		h4_sim_doubleratios_fixedAs[i] = new TH1D**[bins_Pt];
		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			h4_data_doubleratios_fixedAs[i][j] = new TH1D*[bins_As];
			h4_sim_doubleratios_fixedAs[i][j] = new TH1D*[bins_As];
			
			for( int k = 0 ; k < bins_As ; ++k ){ // for bins in As

				c_data_quadratio[i][j][k] = new TCanvas(Form("c_doubleratios_fixedAs_%i_%i_%i",i,j,k),"",1200,1200);
				h4_data_doubleratios_fixedAs[i][j][k] = new TH1D(Form("data_doubleratios_fixedAs_%i_%i_%i",i,j,k),"",bins_Xp,Xp_min,Xp_min+bins_Xp*Xp_step);
				h4_sim_doubleratios_fixedAs[i][j][k] = new TH1D(Form("sim_doubleratios_fixedAs_%i_%i_%i",i,j,k),"",bins_Xp,Xp_min,Xp_min+bins_Xp*Xp_step);

				// Now we need to grab each As bin from all the Xb plots and store it here:
				for( int m = 0 ; m < bins_Xb ; ++m ){ // loop over Xb bina
					double this_as_bin_content = h4_data_doubleratio[i][j][m] -> GetBinContent( k+1 ); // grabs the As bin

					double this_as_bin_content_sim = h4_sim_doubleratio[i][j][m] -> GetBinContent( k+1 ); // grabs the As bin
					h4_sim_doubleratios_fixedAs[i][j][k]->SetBinContent( m+1 , this_as_bin_content_sim );
					if( this_as_bin_content == 0 || this_as_bin_content_sim == 0 ) continue;

					double this_xb = Xb_min[i] + m*Xb_step;
					double this_as = As_min + k*As_step;
					double this_pt = Pt_min + j*Pt_step;
					double this_xp = this_xb / (2-this_as);
					int fill_bin = h4_data_doubleratios_fixedAs[i][j][k]->FindBin( this_xp );

					h4_data_doubleratios_fixedAs[i][j][k]->SetBinContent( fill_bin , this_as_bin_content/this_as_bin_content_sim );

				}
				c_data_quadratio[i][j][k]->cd();
				h4_data_doubleratios_fixedAs[i][j][k]->SetLineWidth(3);
				h4_data_doubleratios_fixedAs[i][j][k]->SetMarkerColor(4);
				h4_data_doubleratios_fixedAs[i][j][k]->SetMaximum(2);
				h4_data_doubleratios_fixedAs[i][j][k]->SetMinimum(0);
				h4_data_doubleratios_fixedAs[i][j][k]->Draw("E1");
				h4_data_doubleratios_fixedAs[i][j][k]->GetXaxis()->SetTitle("x'");

				//h4_sim_doubleratios_fixedAs[i][j][k]->SetLineWidth(3);
				//h4_sim_doubleratios_fixedAs[i][j][k]->SetLineColor(2);
				//h4_sim_doubleratios_fixedAs[i][j][k]->Draw("SAME,E1");

				c_data_quadratio[i][j][k]->Print(Form("data_tagged_doubleratios_fixedAs_Q2_%i_Pt_%i_As_%i.pdf",i,j,k));
			}
		}
	}



	TFile *storehistos = new TFile("ratio_histograms.root","RECREATE");
	storehistos->cd();

	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		h2_data_corrected_inclusive[i]->Write();//Denominator Eq12
		h2_sim_corrected_inclusive[i]->Write();//Denominator Eq12

		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			for( int k = 0 ; k < bins_Xb ; ++k ){ // loop over Xb bin
				h4_data_corrected_tagged[i][j][k]->Write(); //Nominator Eq12
				h4_sim_corrected_tagged[i][j][k]->Write(); //Nominator Eq12
				h4_data_singleratio[i][j][k]->Write(); //Eq 12 result
				h4_sim_singleratio[i][j][k]->Write(); //Eq 12 result
				h4_data_doubleratio[i][j][k]->Write(); //Eq 15 result
				h4_sim_doubleratio[i][j][k]->Write(); //Eq 15 result
				h4_data_sim_doubleratios[i][j][k]->Write(); //Data to sim double ratios
			}
		}
	}
	storehistos->Write();
	storehistos->Close();

}
