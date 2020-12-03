
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

#include "constants.h"
#include "clashit.h"

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
	const double pE_min = 3;
	const double pE_max = 10.6;
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
	TCut inclusive	= ePID && eCharge && eEoP && eEpcal && eVW && eVtx && eMom && eQ2 && eW;
	//TCut inclusive	= ePID && eCharge && eMom && eQ2 && eW;
	TString cut = Form("weight*%s",inclusive.GetTitle());

	TFile * inSimFile = new TFile(argv[1]);
	TFile * inDatFile = new TFile(argv[2]);
	TTree * inSimTree = (TTree*)	readTree(inSimFile);
	TTree * inDatTree = (TTree*)	readTree(inDatFile);
	
	// Define histograms for the output
	TFile * outFile = new TFile(argv[3],"RECREATE");
	TH1D * h1dat_p_e	= new TH1D("h1dat_p_e",		"",250,0,10);
	TH1D * h1sim_p_e	= new TH1D("h1sim_p_e",		"",250,0,10);

	TH1D * h1dat_px_e	= new TH1D("h1dat_px_e",		"",250,-5,5);
	TH1D * h1sim_px_e	= new TH1D("h1sim_px_e",		"",250,-5,5);
	TH1D * h1dat_py_e	= new TH1D("h1dat_py_e",		"",250,-5,5);
	TH1D * h1sim_py_e	= new TH1D("h1sim_py_e",		"",250,-5,5);
	TH1D * h1dat_pz_e	= new TH1D("h1dat_pz_e",		"",250,0,10);
	TH1D * h1sim_pz_e	= new TH1D("h1sim_pz_e",		"",250,0,10);

	TH1D * h1dat_th_e	= new TH1D("h1dat_th_e",	"",250,0,40);
	TH1D * h1sim_th_e	= new TH1D("h1sim_th_e",	"",250,0,40);

	TH1D * h1dat_ph_e	= new TH1D("h1dat_ph_e",	"",250,-200,200);
	TH1D * h1sim_ph_e	= new TH1D("h1sim_ph_e",	"",250,-200,200);

	TH1D * h1dat_q		= new TH1D("h1dat_q",		"",250,0,10);
	TH1D * h1sim_q		= new TH1D("h1sim_q",		"",250,0,10);

	TH1D * h1dat_th_q	= new TH1D("h1dat_th_q",	"",250,0,50);
	TH1D * h1sim_th_q	= new TH1D("h1sim_th_q",	"",250,0,50);

	TH1D * h1dat_ph_q	= new TH1D("h1dat_ph_q",	"",250,-200,200);
	TH1D * h1sim_ph_q	= new TH1D("h1sim_ph_q",	"",250,-200,200);

	TH1D * h1dat_nu		= new TH1D("h1dat_nu",		"",250,0,10);
	TH1D * h1sim_nu		= new TH1D("h1sim_nu",		"",250,0,10);

	TH1D * h1dat_Q2		= new TH1D("h1dat_Q2",		"",250,0,10);
	TH1D * h1sim_Q2		= new TH1D("h1sim_Q2",		"",250,0,10);

	TH1D * h1dat_xB		= new TH1D("h1dat_xB",		"",250,0,1);
	TH1D * h1sim_xB		= new TH1D("h1sim_xB",		"",250,0,1);

	TH1D * h1dat_W		= new TH1D("h1dat_W",		"",250,0,5);
	TH1D * h1sim_W		= new TH1D("h1sim_W",		"",250,0,5);



	TCanvas * c_p_e = new TCanvas("c_p_e");
	c_p_e->Divide(1,2);
	c_p_e->cd(1);
	inDatTree->Draw("eHit->getMomentum() >> h1dat_p_e"								,cut);
	inSimTree->Draw("eHit->getMomentum() >> h1sim_p_e"								,cut);
	label1D(	h1dat_p_e,	h1sim_p_e,	"|p_e| [GeV]",	"Counts"	);
	c_p_e->cd(2);
	label1D_ratio(	h1dat_p_e,      h1sim_p_e,      "|p_e| [GeV]",  "Data/Sim"	);
	c_p_e->Update();

	TCanvas * c_p_comp = new TCanvas("c_p_comp");
	c_p_comp->Divide(2,2);
	c_p_comp->cd(1);
	inDatTree->Draw("eHit->getMomentum()*sin(eHit->getTheta())*cos(eHit->getPhi()) >> h1dat_px_e"                                    	,cut);
        inSimTree->Draw("eHit->getMomentum()*sin(eHit->getTheta())*cos(eHit->getPhi()) >> h1sim_px_e"                                    	,cut);
	label1D(        h1dat_px_e,	h1sim_px_e,	"p_e X [GeV]",	"Counts"	);
	c_p_comp->cd(2);
	inDatTree->Draw("eHit->getMomentum()*sin(eHit->getTheta())*sin(eHit->getPhi()) >> h1dat_py_e"                                    	,cut);
        inSimTree->Draw("eHit->getMomentum()*sin(eHit->getTheta())*sin(eHit->getPhi()) >> h1sim_py_e"                                    	,cut);
	label1D(        h1dat_py_e,	h1sim_py_e,	"p_e Y [GeV]",	"Counts"	);
	c_p_comp->cd(3);
	inDatTree->Draw("eHit->getMomentum()*cos(eHit->getTheta()) >> h1dat_pz_e"                                               	,cut);
        inSimTree->Draw("eHit->getMomentum()*cos(eHit->getTheta()) >> h1sim_pz_e"                                               	,cut);
	label1D(        h1dat_pz_e,	h1sim_pz_e,	"p_e Z [GeV]",	"Counts"	);
	c_p_comp->cd(2);
	c_p_comp->Update();
	

	TCanvas * c_th_e = new TCanvas("c_th_e");
	c_th_e->Divide(1,2);
	c_th_e->cd(1);
	inDatTree->Draw("eHit->getTheta()*180./TMath::Pi() >> h1dat_th_e"					,cut);
	inSimTree->Draw("eHit->getTheta()*180./TMath::Pi() >> h1sim_th_e"					,cut);
	label1D(        h1dat_th_e,	h1sim_th_e,	"Theta e [deg.]",	"Counts"	);
	c_th_e->cd(2);
	label1D_ratio(	h1dat_th_e,     h1sim_th_e,     "Theta e [deg.]",       "Data/Sim"	);
	c_th_e->Update();	

	TCanvas * c_ph_e = new TCanvas("c_ph_e");
	c_ph_e->Divide(1,2);
	c_ph_e->cd(1);
	inDatTree->Draw("eHit->getPhi()*180./TMath::Pi() >> h1dat_ph_e"						,cut);
	inSimTree->Draw("eHit->getPhi()*180./TMath::Pi() >> h1sim_ph_e"						,cut);
	label1D(        h1dat_ph_e,	h1sim_ph_e,	"Phi e [deg.]",	"Counts"	);
	c_ph_e->cd(2);
	label1D_ratio(	h1dat_ph_e,     h1sim_ph_e,     "Phi e [deg.]", "Data/Sim"	);
	c_ph_e->Update();

	TCanvas * c_q = new TCanvas("c_q");
	c_q->Divide(1,2);
	c_q->cd(1);
	inDatTree->Draw("eHit->getQ() >> h1dat_q"									,cut);
	inSimTree->Draw("eHit->getQ() >> h1sim_q"									,cut);
	label1D(        h1dat_q,	h1sim_q,	"|q| [GeV]",	"Counts"	);
	c_q->cd(2);
	label1D_ratio(	h1dat_q,        h1sim_q,        "|q| [GeV]",    "Data/Sim"	);
	c_q->Update();

	TCanvas * c_th_q = new TCanvas("c_th_q");
	c_th_q->Divide(1,2);
	c_th_q->cd(1);
	inDatTree->Draw("eHit->getThetaQ()*180./TMath::Pi() >> h1dat_th_q"					,cut);
	inSimTree->Draw("eHit->getThetaQ()*180./TMath::Pi() >> h1sim_th_q"					,cut);
	label1D(        h1dat_th_q,	h1sim_th_q,	"Theta q [deg.]",	"Counts"	);
	c_th_q->cd(2);
	label1D_ratio(	h1dat_th_q,     h1sim_th_q,     "Theta q [deg.]",       "Data/Sim"	);
	c_th_q->Update();

	TCanvas * c_ph_q = new TCanvas("c_ph_q");
	c_ph_q->Divide(1,2);
	c_ph_q->cd(1);
	inDatTree->Draw("eHit->getPhiQ()*180./TMath::Pi() >> h1dat_ph_q"						,cut);
	inSimTree->Draw("eHit->getPhiQ()*180./TMath::Pi() >> h1sim_ph_q"						,cut);
	label1D(        h1dat_ph_q,	h1sim_ph_q,	"Phi q [deg.]",	"Counts"	);
	c_ph_q->cd(2);
	label1D_ratio(	h1dat_ph_q,     h1sim_ph_q,     "Phi q [deg.]", "Data/Sim"	);
	c_ph_q->Update();

	TCanvas * c_nu = new TCanvas("c_nu");
	c_nu->Divide(1,2);
	c_nu->cd(1);
	inDatTree->Draw("eHit->getOmega() >> h1dat_nu"								,cut);
	inSimTree->Draw("eHit->getOmega() >> h1sim_nu"								,cut);
	label1D(        h1dat_nu,	h1sim_nu,	"Nu [GeV]",	"Counts"	);
	c_nu->cd(2);
	label1D_ratio(	h1dat_nu,       h1sim_nu,       "Nu [GeV]",     "Data/Sim"	);
	c_nu->Update();
	
	TCanvas * c_Q2 = new TCanvas("c_Q2");
	c_Q2->Divide(1,2);
	c_Q2->cd(1);
	inDatTree->Draw("eHit->getQ2() >> h1dat_Q2"								,cut);
	inSimTree->Draw("eHit->getQ2() >> h1sim_Q2"								,cut);
	label1D(        h1dat_Q2,	h1sim_Q2,	"Q2 [GeV^2]",	"Counts"	);
	c_Q2->cd(2);
	label1D_ratio(	h1dat_Q2,       h1sim_Q2,       "Q2 [GeV^2]",   "Data/Sim"	);
	c_Q2->Update();

	TCanvas * c_xB = new TCanvas("c_xB");
	c_xB->Divide(1,2);
	c_xB->cd(1);
	inDatTree->Draw("eHit->getXb() >> h1dat_xB"								,cut);
	inSimTree->Draw("eHit->getXb() >> h1sim_xB"								,cut);
	label1D(        h1dat_xB,	h1sim_xB,	"xB",	"Counts"	);
	c_xB->cd(2);
	label1D_ratio(	h1dat_xB,       h1sim_xB,       "xB",   "Data/Sim"	);
	c_xB->Update();

	TCanvas * c_W = new TCanvas("c_W");
	c_W->Divide(1,2);
	c_W->cd(1);
	inDatTree->Draw("sqrt(eHit->getW2()) >> h1dat_W"								,cut);
	inSimTree->Draw("sqrt(eHit->getW2()) >> h1sim_W"								,cut);
	label1D(        h1dat_W,	h1sim_W,	"W [GeV]",	"Counts"	);
	c_W->cd(2);
	label1D_ratio(	h1dat_W,        h1sim_W,        "W [GeV]",      "Data/Sim"	);
	c_W->Update();

	TCanvas * c0 = new TCanvas("c0","c0",700,900);
	c0 -> Modified();
	c0 -> Update(); 
	c0 -> Print(Form("comparison_test.pdf(" ));
	c_p_e	->Print(Form("comparison_test.pdf(" ));
	c_p_comp->Print(Form("comparison_test.pdf(" ));
	c_th_e	->Print(Form("comparison_test.pdf(" ));
	c_ph_e	->Print(Form("comparison_test.pdf(" ));
	c_q	->Print(Form("comparison_test.pdf(" ));
	c_th_q	->Print(Form("comparison_test.pdf(" ));
	c_ph_q	->Print(Form("comparison_test.pdf(" ));
	c_nu	->Print(Form("comparison_test.pdf(" ));
	c_Q2	->Print(Form("comparison_test.pdf(" ));
	c_xB	->Print(Form("comparison_test.pdf(" ));
	c_W	->Print(Form("comparison_test.pdf(" ));
	c0 -> Print(Form("comparison_test.pdf)" ));
	
	

	outFile->cd();
	c_p_e->Write();
	c_p_comp->Write();
	c_th_e->Write();
	c_ph_e->Write();
	c_q->Write();
	c_th_q->Write();
	c_ph_q->Write();
	c_nu->Write();
	c_Q2->Write();
	c_xB->Write();
	c_W->Write();
	outFile->Close();

	return 1;
}

TTree* readTree(TFile* inFile){
	TTree * inTree = (TTree*)inFile->Get("electrons");
	// Initialize the input branches
	int Runno		= 0;
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	double weight		= 0;
	int nMult		= 0;
	clashit* eHit = new clashit;
	inTree->SetBranchAddress("Runno"		,&Runno			);
	inTree->SetBranchAddress("Ebeam"		,&Ebeam			);
	inTree->SetBranchAddress("gated_charge"		,&gated_charge		);
	inTree->SetBranchAddress("livetime"		,&livetime		);
	inTree->SetBranchAddress("starttime"		,&starttime		);
	inTree->SetBranchAddress("current"		,&current		);
	inTree->SetBranchAddress("weight"		,&weight		);
	// 	Electron branches:
	inTree->SetBranchAddress("eHit"			,&eHit			);
	
	return inTree;
}

void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel){
	data->SetLineColor(1);
	data->SetLineWidth(3);
	data->SetStats(0);

	sim->SetLineColor(9);
	sim->SetLineWidth(3);
	sim->SetStats(0);
	sim->Scale(data->Integral() / sim->Integral() );

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
	legend->AddEntry(data,"Data","f");
	legend->AddEntry(sim,"Sim","f");
	legend->Draw();

	return;
}
void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel){
	TH1D * data_copy = (TH1D*) data->Clone();
	TH1D * sim_copy = (TH1D*) sim->Clone();
	
	data_copy->SetLineColor(1);
	data_copy->SetLineWidth(3);
	data_copy->SetStats(0);

	sim_copy->SetLineColor(9);
	sim_copy->SetLineWidth(3);
	sim_copy->SetStats(0);
	sim_copy->Scale(data_copy->Integral() / sim_copy->Integral() );

	data_copy->Divide(sim_copy);
	data_copy->Draw("ep");
	TLine* line = new TLine(data_copy->GetXaxis()->GetBinCenter(1), 1., data_copy->GetXaxis()->GetBinCenter(data_copy->GetXaxis()->GetNbins()), 1.);
	line->SetLineWidth(3);
	line->SetLineColor(2);
	line->Draw("same");

	double max1 = data_copy->GetMaximum()*1.1;
	data_copy->GetYaxis()->SetRangeUser(0.5,1.5);
	
	data_copy->GetXaxis()->SetTitle(xlabel);
	data_copy->GetYaxis()->SetTitle(ylabel);

	return;
}
