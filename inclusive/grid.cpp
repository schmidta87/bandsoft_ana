#include <iostream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TClonesArray.h"

#include "constants.h"
#include "clashit.h"
#include "genpart.h"

TTree* readTree(TFile* inFile);
void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
using namespace std;

int main(int argc, char ** argv){

	if (argc<4){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [SimFile] [DataFile] [OutputFile]\n";
		return -1;
	}

	// Electron cuts
	const int PID = 11;
	const int charge = -1;
	const double EoP_min = 0.17;
	const double EoP_max = 0.3;
	const double Epcal_min = 0.07;
	const double V_min = 15;
	const double W_min = 15;
	const double vtx_min = -8;
	const double vtx_max =  3;
	const double pE_min = 0;
	const double pE_max = 10.2;
	const double Q2_min = 2;
	const double Q2_max = 10;
	const double W2_min = 2*2;
	TCut ePID 	= Form("eHit->getPID() == %i",PID);
	TCut eCharge 	= Form("eHit->getCharge() == %i",charge); 
	TCut eEoP	= Form("eHit->getEoP() > %f && eHit->getEoP() < %f",EoP_min,EoP_max);
	TCut eEpcal	= Form("eHit->getEpcal() > %f",Epcal_min);
	TCut eVW	= Form("eHit->getV() > %f && eHit->getW() > %f",V_min,W_min);
	TCut eVtx	= Form("eHit->getVtz() > %f && eHit->getVtz() < %f",vtx_min,vtx_max);
	TCut eMom	= Form("eHit->getMomentum() > %f && eHit->getMomentum() < %f",pE_min,pE_max);
	TCut eQ2	= Form("eHit->getQ2() > %f && eHit->getQ2() < %f",Q2_min,Q2_max);
	TCut eW		= Form("eHit->getW2() > %f",W2_min);
	TCut inclusive	= ePID && eCharge && eEoP && eEpcal && eVW && eVtx && eQ2 && eW;
	//TCut inclusive	= ePID && eCharge && eMom && eQ2 && eW;
	TString cut = Form("weight*%s",inclusive.GetTitle());
	// MC cuts
	TCut mcMom	= Form("mcParts[0]->getMomentum() > %f && mcParts[0]->getMomentum() < %f",pE_min,pE_max);
	TCut mcQ2	= Form("mcParts[0]->getQ2() > %f && mcParts[0]->getQ2() < %f",Q2_min,Q2_max);
	TCut mcW	= Form("mcParts[0]->getW2() > %f",W2_min);
	TCut mc_inclusive = mcQ2 && mcW;
	TString mc_cut = Form("weight*%s",mc_inclusive.GetTitle());

	TFile * inSimFile = new TFile(argv[1]);
	TTree * inSimTree = (TTree*) inSimFile->Get("electrons");

	TFile * outFile = new TFile(argv[3],"RECREATE");

	//int q2_bins = 8;
	int q2_bins = 1;
	double w_bins[3] = {2,3,4.5};
	TH1D ** h1gen_w_q2s = new TH1D*[q2_bins];
	TH1D ** h1rec_w_q2s = new TH1D*[q2_bins];
	for( int i = 0 ; i < q2_bins ; i++ ){
		//h1gen_w_q2s[i] = new TH1D(Form("h1gen_w_q2s_%i",i),Form("h1gen_w_q2s_%i",i),10,2,4.5);
		//h1rec_w_q2s[i] = new TH1D(Form("h1rec_w_q2s_%i",i),Form("h1rec_w_q2s_%i",i),10,2,4.5);
		h1gen_w_q2s[i] = new TH1D(Form("h1gen_w_q2s_%i",i),Form("h1gen_w_q2s_%i",i),2,w_bins);
		h1rec_w_q2s[i] = new TH1D(Form("h1rec_w_q2s_%i",i),Form("h1rec_w_q2s_%i",i),2,w_bins);
	}
	int pE_bins = 8;
	TH1D ** h1gen_thetaE_pEs = new TH1D*[pE_bins];
	TH1D ** h1rec_thetaE_pEs = new TH1D*[pE_bins];
	for( int i = 0 ; i < q2_bins ; i++ ){
		h1gen_thetaE_pEs[i] = new TH1D(Form("h1gen_thetaE_pEs_%i",i),Form("h1gen_thetaE_pEs_%i",i),40,0,40);
		h1rec_thetaE_pEs[i] = new TH1D(Form("h1rec_thetaE_pEs_%i",i),Form("h1rec_thetaE_pEs_%i",i),40,0,40);
	}
	
	TCanvas * c_w_q2 = new TCanvas("c_w_q2");
	TCanvas * c_w_q2_rat = new TCanvas("c_w_q2_rat");
	c_w_q2->Divide(2,4);
	c_w_q2_rat->Divide(2,4);
	for( int i = 0 ; i < q2_bins ; i++ ){
		//double Q2_bin_min = Q2_min + 1*i;
		//double Q2_bin_max = Q2_min + 1*(i+1);
		//TCut Q2_refined = Form("eHit->getQ2() > %f && eHit->getQ2() < %f ",Q2_bin_min,Q2_bin_max);
		//TCut new_cut = inclusive && Q2_refined;
		//cut = Form("weight*%s",new_cut.GetTitle());
		c_w_q2->cd(i+1);
		inSimTree->Draw(Form("sqrt(eHit->getW2()) >> h1rec_w_q2s_%i",i)			,cut);
		//h1rec_w_q2s[i]->SetTitle(Form("%i < Q2 < %i",(int)Q2_bin_min,(int)Q2_bin_max));
		h1rec_w_q2s[i]->SetTitle(Form("Q2 > %f",Q2_min));

		//TCut mc_Q2_refined = Form("mcParts[0]->getQ2() > %f && mcParts[0]->getQ2() < %f",Q2_bin_min,Q2_bin_max);
		//TCut mc_new_cut = mc_inclusive && mc_Q2_refined;
		//mc_cut = Form("weight*%s",mc_new_cut.GetTitle());
		inSimTree->Draw(Form("sqrt(mcParts[0]->getW2()) >> h1gen_w_q2s_%i",i)		,mc_cut);

		label1D(h1rec_w_q2s[i],h1gen_w_q2s[i],"W [GeV]","Counts");


		c_w_q2_rat->cd(i+1);
		label1D_ratio(h1rec_w_q2s[i],h1gen_w_q2s[i],"W [GeV]","Rec/Gen");


	}
	c_w_q2->Update();

	
	TCanvas * c_thetaE_pE = new TCanvas("c_thetaE_pE");
	TCanvas * c_thetaE_pE_rat = new TCanvas("c_thetaE_pE_rat");
	c_thetaE_pE->Divide(2,4);
	c_thetaE_pE_rat->Divide(2,4);
	for( int i = 0 ; i < q2_bins ; i++ ){
		double pE_bin_min = pE_min + 1*i;
		double pE_bin_max = pE_min + 1*(i+1);
		TCut pE_refined = Form("eHit->getMomentum() > %f && eHit->getMomentum() < %f ",pE_bin_min,pE_bin_max);
		TCut new_cut = inclusive && pE_refined;
		cut = Form("weight*%s",new_cut.GetTitle());
		c_thetaE_pE->cd(i+1);
		inSimTree->Draw(Form("eHit->getTheta()*180./TMath::Pi() >> h1rec_thetaE_pEs_%i",i)	,cut);
		h1rec_thetaE_pEs[i]->SetTitle(Form("%i < |p_{e}| < %i",(int)pE_bin_min,(int)pE_bin_max));

		TCut mc_pE_refined = Form("mcParts[0]->getMomentum() > %f && mcParts[0]->getMomentum() < %f",pE_bin_min,pE_bin_max);
		TCut mc_new_cut = mc_inclusive && mc_pE_refined;
		mc_cut = Form("weight*%s",mc_new_cut.GetTitle());
		inSimTree->Draw(Form("mcParts[0]->getTheta()*180./TMath::Pi() >> h1gen_thetaE_pEs_%i",i)	,mc_cut);

		label1D(h1rec_thetaE_pEs[i],h1gen_thetaE_pEs[i],"Theta e [Deg.]","Counts");


		c_thetaE_pE_rat->cd(i+1);
		label1D_ratio(h1rec_thetaE_pEs[i],h1gen_thetaE_pEs[i],"Theta e [Deg.]","Counts");


	}
	c_thetaE_pE->Update();
	
	TCanvas * c0 = new TCanvas("c0","c0",700,900);
	c0 -> Modified();
	c0 -> Update(); 
	c0 -> Print(Form("corrections_inclusive.pdf(" ));
	c_w_q2		->Print(Form("corrections_inclusive.pdf(" ));
	c_w_q2_rat	->Print(Form("corrections_inclusive.pdf(" ));
	c_thetaE_pE	->Print(Form("corrections_inclusive.pdf(" ));
	c_thetaE_pE_rat	->Print(Form("corrections_inclusive.pdf(" ));
	c0 -> Print(Form("corrections_inclusive.pdf)" ));
	


	outFile->cd();
	outFile->Close();

	return 1;
}

void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel){
	data->SetLineColor(1);
	data->SetLineWidth(3);
	data->SetStats(0);

	sim->SetLineColor(9);
	sim->SetLineWidth(3);
	sim->SetStats(0);
	//sim->Scale(data->Integral() / sim->Integral() );

	data->Draw("h,p");
	sim->Draw("hist,same");

	double max1 = data->GetMaximum()*1.1;
	double max2 = sim->GetMaximum()*1.1;
	data->GetYaxis()->SetRangeUser(0,max(max1,max2));
	
	data->GetXaxis()->SetTitle(xlabel);
	data->GetYaxis()->SetTitle(ylabel);

	TLegend * legend = new TLegend(0.7,0.8,0.9,0.9);
	//legend->AddEntry(data,"Radiation On","f");
	//legend->AddEntry(sim,"Radiation Off","f");
	legend->AddEntry(data,"Rec","f");
	legend->AddEntry(sim,"Gen","f");
	legend->Draw();

	return;
}
void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel){
	TH1D * data_copy = (TH1D*) data->Clone();
	data_copy->Sumw2();
	TH1D * sim_copy = (TH1D*) sim->Clone();
	
	data_copy->SetLineColor(1);
	data_copy->SetLineWidth(3);
	data_copy->SetStats(0);

	sim_copy->SetLineColor(9);
	sim_copy->SetLineWidth(3);
	sim_copy->SetStats(0);
	//sim_copy->Scale(data_copy->Integral() / sim_copy->Integral() );

	data_copy->Divide(sim_copy);
	data_copy->Draw("ep");
	TLine* line = new TLine(data_copy->GetXaxis()->GetBinCenter(1), 1., data_copy->GetXaxis()->GetBinCenter(data_copy->GetXaxis()->GetNbins()), 1.);
	line->SetLineWidth(3);
	line->SetLineColor(2);
	line->Draw("same");

	double max1 = data_copy->GetMaximum()*1.1;
	data_copy->GetYaxis()->SetRangeUser(0.,1.1);
	
	data_copy->GetXaxis()->SetTitle(xlabel);
	data_copy->GetYaxis()->SetTitle(ylabel);

	return;
}
