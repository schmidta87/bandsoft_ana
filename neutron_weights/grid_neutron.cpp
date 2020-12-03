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
	const double pE_max = 10.6;
	const double Q2_min = 2;
	const double Q2_max = 10;
	const double W2_min = 2*2;
	const double gen_pN_min = 0.3;
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
	TString cut = Form("%s",inclusive.GetTitle());
	// MC cuts
	//TCut mcMom	= Form("mcParts[0]->getMomentum() > %f && mcParts[0]->getMomentum() < %f",pE_min,pE_max);
	//TCut mcQ2	= Form("mcParts[0]->getQ2() > %f && mcParts[0]->getQ2() < %f",Q2_min,Q2_max);
	//TCut mcW	= Form("mcParts[0]->getW2() > %f",W2_min);
	//TCut mc_inclusive = mcQ2 && mcW;
	//TString mc_cut = Form("weight*%s",mc_inclusive.GetTitle());
	// Neutron cuts
	const int MeVee_cut = 5;
	const double pN_min = 0.3;
	TCut nMult	= Form("nMult == 1");
	TCut nStatus	= Form("nHits[0]->getStatus() == 0 ");
	TCut nEdep	= Form("nHits[0]->getEdep()/1E4 > %i",MeVee_cut);
	TCut nTof	= Form("nHits[0]->getTofFadc() != 0");
	TCut nMom	= Form("tag[0]->getMomentumN().Mag() > %f",pN_min);
	TCut tagged	= nMult && nStatus && nEdep && nTof && nMom;
	TString tag_cut = Form("%s",tagged.GetTitle());


	TFile * inSimFile = new TFile(argv[1]);
	TTree * inSimTree = (TTree*) inSimFile->Get("tagged");

	TFile * outFile = new TFile(argv[3],"RECREATE");

	double as_bins[3] = {1.3,1.4,1.6};
	int n_as_bins = 2;
	TH1D ** h1tag_as_w = new TH1D*[n_as_bins];
	TH1D ** h1inc_as_w = new TH1D*[n_as_bins];
	double w_bins[3] = {2,3,4.5};
	for( int i = 0 ; i < n_as_bins ; i++ ){
		h1tag_as_w[i] = new TH1D(Form("h1tag_as_w_%i",i),Form("h1tag_as_w_%i",i),2,w_bins);
		h1inc_as_w[i] = new TH1D(Form("h1inc_as_w_%i",i),Form("h1inc_as_w_%i",i),2,w_bins);
	}
	
	// AlphaS calculation by hand until we add it to our mcgen tree:
	// (sqrt(mcParts[1].getMomentum()*mcParts[1].getMomentum() + 0.939*0.939) 
	// 	- mcParts[1].getMomentum()*( 
	// 	cos(mcParts[1].getTheta())*cos(mcParts[0].getTheta()) + 
	// 	sin(mcParts[1].getTheta())*cos(mcParts[1].getPhi())*sin(mcParts[0].getTheta())*cos(mcParts[0].getPhi()) + 
	// 	sin(mcParts[1].getTheta())*sin(mcParts[1].getPhi())*sin(mcParts[0].getTheta())*sin(mcParts[0].getPhi()) ) 
	// )/0.939

	TCanvas * c_as_w = new TCanvas("c_as_w");
	TCanvas * c_as_w_rat = new TCanvas("c_as_w_rat");
	c_as_w->Divide(2,4);
	c_as_w_rat->Divide(2,4);
	for( int i = 0 ; i < n_as_bins ; i++ ){
		double As_bin_min = as_bins[i];
		double As_bin_max = as_bins[i+1];

		TCut As_refined = Form("(sqrt(mcParts[1].getMomentum()*mcParts[1].getMomentum() + 0.939*0.939) - mcParts[1].getMomentum()*( cos(mcParts[1].getTheta())*cos(mcParts[0].getTheta()) + sin(mcParts[1].getTheta())*cos(mcParts[1].getPhi())*sin(mcParts[0].getTheta())*cos(mcParts[0].getPhi()) + sin(mcParts[1].getTheta())*sin(mcParts[1].getPhi())*sin(mcParts[0].getTheta())*sin(mcParts[0].getPhi()) ) )/0.939 > %f && (sqrt(mcParts[1].getMomentum()*mcParts[1].getMomentum() + 0.939*0.939) - mcParts[1].getMomentum()*( cos(mcParts[1].getTheta())*cos(mcParts[0].getTheta()) + sin(mcParts[1].getTheta())*cos(mcParts[1].getPhi())*sin(mcParts[0].getTheta())*cos(mcParts[0].getPhi()) + sin(mcParts[1].getTheta())*sin(mcParts[1].getPhi())*sin(mcParts[0].getTheta())*sin(mcParts[0].getPhi()) ) )/0.939 < %f ",As_bin_min,As_bin_max);
		TCut new_cut = inclusive && As_refined;
		cut = Form("%s",new_cut.GetTitle());
		c_as_w->cd(i+1);
		inSimTree->Draw(Form("sqrt(eHit->getW2()) >> h1inc_as_w_%i",i)			,cut);
		h1inc_as_w[i]->SetTitle(Form("%f < alpha_{S} < %f",As_bin_min,As_bin_max));

		TCut As_refined_rec = Form("tag[0].getAs() > %f && tag[0].getAs() < %f",As_bin_min,As_bin_max);
		TCut tag_new_cut = inclusive && tagged && As_refined_rec;
		tag_cut = Form("%s",tag_new_cut.GetTitle());
		inSimTree->Draw(Form("sqrt(eHit->getW2()) >> h1tag_as_w_%i",i)			,tag_cut);

		label1D(h1inc_as_w[i],h1tag_as_w[i],"W [GeV]","Counts");


		c_as_w_rat->cd(i+1);
		label1D_ratio(h1tag_as_w[i],h1inc_as_w[i],"W [GeV]","Tag/Inclusive");


	}
	c_as_w->Update();

	
	
	TCanvas * c0 = new TCanvas("c0","c0",700,900);
	c0 -> Modified();
	c0 -> Update(); 
	c0 -> Print(Form("corrections_inclusive.pdf(" ));
	c_as_w		->Print(Form("corrections_inclusive.pdf(" ));
	c_as_w_rat	->Print(Form("corrections_inclusive.pdf(" ));
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
	legend->AddEntry(data,"Inclusive","f");
	legend->AddEntry(sim,"Tagged","f");
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
