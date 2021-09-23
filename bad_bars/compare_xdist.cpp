#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>
using namespace std;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
void background_subtraction(TH1D* dat, TH1D* bac , double Cscale, double NB_sim , double Sigma_Cscale, double Sigma_NB_sim );
void simulation_weighting(TH1D* sim, double Ndata, double Nsim );
void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
int main( int argc, char** argv ){

	if( argc != 4 ){
		return -1;
	}

	// Files
	TFile * dat = new TFile(argv[1]);
	TFile * bac = new TFile(argv[2]);
	TFile * sim = new TFile(argv[3]);

	// Trees
	TTree * datTree	= (TTree*) dat->Get("tagged");
	TTree * bacTree	= (TTree*) bac->Get("tagged");
	TTree * simTree	= (TTree*) sim->Get("tagged");

	// Background normalization
	TVector3 * datnorm = (TVector3*) dat->Get("bacnorm");
	TVector3 * bacnorm = (TVector3*) bac->Get("bacnorm");
	// normalization uncertainty of the background:
	double Cscale = datnorm->Z(); // this is the same as bacnorm->X() for the data file
	double Sigma_Cscale = (datnorm->Y() - datnorm->X())/2.;  // this is the before-time and after-time levels
	double NB_sim = bacnorm->X();
	double Sigma_NB_sim = sqrt(NB_sim);
		// sigma_Cscale / Cscale ~ 7%
	double Ndata = 0;
	double Nsim = 0;
	double sim_scaling = 0;

	// Add a TCut on momentum:
	TCut mom_cut = "tag[nleadindex]->getMomentumN().Mag() > 0.3 && tag[nleadindex]->getMomentumN().Theta() < 168.5*TMath::Pi()/180.";

	// Let's get the absolute normalization for the simulation using some pN plot:
	TH1D * pn_dat = new TH1D("pn_dat","pn_dat",100,0,2);
	TH1D * pn_bac = new TH1D("pn_bac","pn_bac",100,0,2);
	TH1D * pn_sim = new TH1D("pn_sim","pn_sim",100,0,2);
	datTree->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_dat" , mom_cut );
	bacTree->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_bac" , mom_cut );
	simTree->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_sim" , mom_cut );
	background_subtraction( pn_dat , pn_bac , Cscale, NB_sim, Sigma_Cscale, Sigma_NB_sim );
	Ndata = pn_dat->Integral(); 		// total data we have
	Nsim = pn_sim->Integral(); 		// total simulation we have
	sim_scaling = Ndata/Nsim;		// scale of the simulation bin
	// Simulation weighting
	simulation_weighting( pn_sim, Ndata, Nsim );
	double full_simscale = sim_scaling;
	
	// Histogram of all our chi2 values:
	TH1D * chi2_vals = new TH1D("chi2_vals","",40,0,20);

	// Now we want to loop over all bars and create histograms to get some sort of "quality" statistic:
	TH1D **** xn_dat = new TH1D***[5];
	TH1D **** xn_bac = new TH1D***[5];
	TH1D **** xn_sim = new TH1D***[5];
	TCanvas * cTotal = new TCanvas("cTotal","cTotal",800,600);
	cTotal->cd();
	label1D(pn_dat, pn_sim, "p_{n} [GeV/c]", "Counts");
	cTotal->Print("xn_distributions_10pt2.pdf(");
	for( int layer = 1 ; layer < 6 ; layer++ ){
		xn_dat[layer-1] = new TH1D**[5];
		xn_bac[layer-1] = new TH1D**[5];
		xn_sim[layer-1] = new TH1D**[5];
		for( int sector = 1 ; sector < 6 ; sector++ ){
			xn_dat[layer-1][sector-1] = new TH1D*[slc[layer-1][sector-1]];
			xn_bac[layer-1][sector-1] = new TH1D*[slc[layer-1][sector-1]];
			xn_sim[layer-1][sector-1] = new TH1D*[slc[layer-1][sector-1]];
			for( int component = 1 ; component < slc[layer-1][sector-1]+1 ; component++ ){
				xn_dat[layer-1][sector-1][component-1] = 
					new TH1D(Form("xn_dat_%i_%i_%i",layer,sector,component),Form("xn_dat_%i_%i_%i",layer,sector,component),18,-90,90);
				xn_dat[layer-1][sector-1][component-1]->Sumw2();
				xn_bac[layer-1][sector-1][component-1] = 
					new TH1D(Form("xn_bac_%i_%i_%i",layer,sector,component),Form("xn_bac_%i_%i_%i",layer,sector,component),18,-90,90);
				xn_sim[layer-1][sector-1][component-1] = 
					new TH1D(Form("xn_sim_%i_%i_%i",layer,sector,component),Form("xn_sim_%i_%i_%i",layer,sector,component),18,-90,90);
		
				// Cut for this specific bar:
				TCut bar_cut = Form("nHits[nleadindex]->getLayer()==%i && nHits[nleadindex]->getSector()==%i && nHits[nleadindex]->getComponent()==%i",
							layer,sector,component);

				TCanvas * c = new TCanvas(Form("c_%i_%i_%i",layer,sector,component),"",800,600);
				datTree->Draw(Form("nHits[nleadindex]->getDL().X() >> xn_dat_%i_%i_%i",layer,sector,component), mom_cut && bar_cut );
				bacTree->Draw(Form("nHits[nleadindex]->getDL().X() >> xn_bac_%i_%i_%i",layer,sector,component), mom_cut && bar_cut );
				simTree->Draw(Form("nHits[nleadindex]->getDL().X() >> xn_sim_%i_%i_%i",layer,sector,component), mom_cut && bar_cut );
				if( 	xn_dat[layer-1][sector-1][component-1]->Integral() == 0 ||
					xn_bac[layer-1][sector-1][component-1]->Integral() == 0 ||
					xn_sim[layer-1][sector-1][component-1]->Integral() == 0 ) continue;


				
				background_subtraction( xn_dat[layer-1][sector-1][component-1] , 
						xn_bac[layer-1][sector-1][component-1] , Cscale, NB_sim, Sigma_Cscale, Sigma_NB_sim );
				xn_sim[layer-1][sector-1][component-1]->Scale( full_simscale );

				double rescaling = xn_dat[layer-1][sector-1][component-1]->Integral() / xn_sim[layer-1][sector-1][component-1]->Integral();
				//xn_sim[layer-1][sector-1][component-1]->Scale( rescaling );


				// Calculate some sort of chi2 with data vs simulation:
				double chi2 = 0;
				int npts = 0;
				for( int bin = 1 ; bin < xn_dat[layer-1][sector-1][component-1]->GetNbinsX()+1 ; bin++){
					double datval = xn_dat[layer-1][sector-1][component-1]->GetBinContent(bin);
					double simval = xn_sim[layer-1][sector-1][component-1]->GetBinContent(bin);
					double daterr = xn_dat[layer-1][sector-1][component-1]->GetBinError(bin);
					double toterr = sqrt( daterr*daterr + simval );
					double this_chi2 = pow( (datval-simval)/toterr ,2 );
					cout << "\tbar details: " << datval << " " << simval << " " << daterr << " " << this_chi2 << " " << xn_dat[layer-1][sector-1][component-1]->GetBinCenter(bin) << "\n";
					cout << flush;
					if( datval == 0 || simval == 0 || daterr == 0 || (daterr!=daterr) ) continue;
					chi2 += this_chi2;
					npts += 1;
				}
				chi2 /= npts;
				chi2_vals->Fill( chi2 );

				// Set the title for this histogram
				int barID = sector*100 + layer*10 + component;
				xn_sim[layer-1][sector-1][component-1]->SetTitle(Form("Bar = %i , Chi = %f , C_{new} = %f",barID,chi2,rescaling));
				
				cout << "plotting bar: " << layer << " " << sector << " " << component << " " << chi2 << "\n" << flush;
				// Draw and format:
				xn_sim[layer-1][sector-1][component-1]->SetLineColor(2);
				xn_sim[layer-1][sector-1][component-1]->SetLineWidth(1);
				xn_sim[layer-1][sector-1][component-1]->SetStats(0);
				if( xn_sim[layer-1][sector-1][component-1]->GetMaximum() > xn_dat[layer-1][sector-1][component-1]->GetMaximum() ){ 
					xn_sim[layer-1][sector-1][component-1]->SetMaximum( xn_sim[layer-1][sector-1][component-1]->GetMaximum()*1.1 ); 
				}
				else{
					xn_sim[layer-1][sector-1][component-1]->SetMaximum( xn_dat[layer-1][sector-1][component-1]->GetMaximum()*1.1 ); 
				}
				xn_sim[layer-1][sector-1][component-1]->Draw("HIST");

				xn_dat[layer-1][sector-1][component-1]->SetLineColor(4);
				xn_dat[layer-1][sector-1][component-1]->SetMarkerColor(4);
				xn_dat[layer-1][sector-1][component-1]->SetMarkerStyle(8);
				xn_dat[layer-1][sector-1][component-1]->SetMarkerSize(1);
				xn_dat[layer-1][sector-1][component-1]->SetStats(0);
				xn_dat[layer-1][sector-1][component-1]->Draw("P,SAME");

				xn_sim[layer-1][sector-1][component-1]->SetMinimum(-1);

				c->Print("xn_distributions_10pt2.pdf");
				c->Close();
			}
		}
	}

	TCanvas * c_chi2 = new TCanvas("c_chi2","c_chi2",800,600);
	c_chi2->cd();
	chi2_vals->Draw("HIST");
	c_chi2->Print("xn_distributions_10pt2.pdf");
	c_chi2->Close();
	cTotal->Print("xn_distributions_10pt2.pdf)");
	cTotal->Close();

	dat->Close();
	bac->Close();
	sim->Close();
	return 0;
}
void background_subtraction(TH1D* dat, TH1D* bac , double Cscale, double NB_sim , double Sigma_Cscale, double Sigma_NB_sim ){
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
void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel){
	data->SetLineColor(4);
	data->SetMarkerColor(4);
	data->SetMarkerStyle(8);
	data->SetMarkerSize(1);
	data->SetStats(0);

	sim->SetLineColor(2);
	sim->SetLineWidth(1);
	sim->SetStats(0);
	sim->SetTitle("");
	data->SetTitle("");

	sim->Draw("hist");
	data->Draw("p,same");

	for( int bin = 1 ; bin < data->GetXaxis()->GetNbins() ; ++bin ){
		cerr << "Data: " << data->GetBinCenter(bin) << " " << data->GetBinContent(bin) << " " << data->GetBinError(bin) << "\n";
		cerr << "Sim: " << sim->GetBinCenter(bin) << " " << sim->GetBinContent(bin) << " " << sim->GetBinError(bin) << "\n";
	}

	double max1 = data->GetMaximum()*1.1;
	double max2 = sim->GetMaximum()*1.1;
	sim->GetYaxis()->SetRangeUser(0,max(max1,max2));
	
	sim->GetXaxis()->SetTitle(xlabel);
	sim->GetYaxis()->SetTitle(ylabel);

	TLegend * legend = new TLegend(0.7,0.8,0.9,0.9);
	//legend->AddEntry(data,"Radiation On","f");
	//legend->AddEntry(sim,"Radiation Off","f");
	legend->AddEntry(data,"Data","f");
	legend->AddEntry(sim,"Sim","f");
	legend->Draw();

	return;
}
