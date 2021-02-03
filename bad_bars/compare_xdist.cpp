#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TCut.h"
#include "TCanvas.h"

#include <iostream>
using namespace std;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
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
	bacTree->SetWeight( datnorm->X() / bacnorm->X() );

	// Add a TCut on momentum:
	TCut mom_cut = "tag[nleadindex]->getMomentumN().Mag() > 0.3";

	// Let's get the absolute normalization for the simulation using some pN plot:
	TH1D * pn_dat = new TH1D("pn_dat","pn_dat",100,0,1);
	TH1D * pn_bac = new TH1D("pn_bac","pn_bac",100,0,1);
	TH1D * pn_sim = new TH1D("pn_sim","pn_sim",100,0,1);
	datTree->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_dat" , mom_cut );
	bacTree->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_bac" , mom_cut );
	simTree->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_sim" , mom_cut );
	pn_dat->Add(pn_bac,-1);
	double full_simscale = pn_dat->Integral() / pn_sim->Integral();
	
	// Histogram of all our chi2 values:
	TH1D * chi2_vals = new TH1D("chi2_vals","",40,0,20);

	// Now we want to loop over all bars and create histograms to get some sort of "quality" statistic:
	TH1D **** xn_dat = new TH1D***[5];
	TH1D **** xn_bac = new TH1D***[5];
	TH1D **** xn_sim = new TH1D***[5];
	TCanvas * cTotal = new TCanvas("cTotal","cTotal",800,600);
	cTotal->Print("xn_distributions.pdf(");
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
					new TH1D(Form("xn_dat_%i_%i_%i",layer,sector,component),Form("xn_dat_%i_%i_%i",layer,sector,component),30,-120,120);
				xn_dat[layer-1][sector-1][component-1]->Sumw2();
				xn_bac[layer-1][sector-1][component-1] = 
					new TH1D(Form("xn_bac_%i_%i_%i",layer,sector,component),Form("xn_bac_%i_%i_%i",layer,sector,component),30,-120,120);
				xn_sim[layer-1][sector-1][component-1] = 
					new TH1D(Form("xn_sim_%i_%i_%i",layer,sector,component),Form("xn_sim_%i_%i_%i",layer,sector,component),30,-120,120);
		
				// Cut for this specific bar:
				TCut bar_cut = Form("nHits[nleadindex]->getLayer()==%i && nHits[nleadindex]->getSector()==%i && nHits[nleadindex]->getComponent()==%i",
							layer,sector,component);

				TCanvas * c = new TCanvas("c","",800,600);
				datTree->Draw(Form("nHits[nleadindex]->getDL().X() >> xn_dat_%i_%i_%i",layer,sector,component), mom_cut && bar_cut );
				bacTree->Draw(Form("nHits[nleadindex]->getDL().X() >> xn_bac_%i_%i_%i",layer,sector,component), mom_cut && bar_cut );
				simTree->Draw(Form("nHits[nleadindex]->getDL().X() >> xn_sim_%i_%i_%i",layer,sector,component), mom_cut && bar_cut );

				xn_dat[layer-1][sector-1][component-1]->Add(xn_bac[layer-1][sector-1][component-1],-1);
				xn_sim[layer-1][sector-1][component-1]->Scale( full_simscale );

				double rescaling = xn_dat[layer-1][sector-1][component-1]->Integral() / xn_sim[layer-1][sector-1][component-1]->Integral();

				// Calculate some sort of chi2 with data vs simulation:
				double chi2 = 0;
				int npts = 0;
				for( int bin = 1 ; bin < xn_dat[layer-1][sector-1][component-1]->GetNbinsX() ; bin++){
					double datval = xn_dat[layer-1][sector-1][component-1]->GetBinContent(bin);
					double simval = xn_sim[layer-1][sector-1][component-1]->GetBinContent(bin);
					double daterr = xn_dat[layer-1][sector-1][component-1]->GetBinError(bin);
					if( datval == 0 || simval == 0 || daterr == 0 || (daterr!=daterr) ) continue;
					chi2 += pow( (datval-simval)/daterr ,2);
					npts += 1;
				}
				chi2 /= npts;
				chi2_vals->Fill( chi2 );

				// Set the title for this histogram
				int barID = sector*100 + layer*10 + component;
				xn_sim[layer-1][sector-1][component-1]->SetTitle(Form("Bar = %i , Chi = %f , C_{new} = %f",barID,chi2,rescaling));
				
				cout << "\tworking on bar: " << layer << " " << sector << " " << component << " " << chi2 << "\n";
				// Draw and format:
				xn_sim[layer-1][sector-1][component-1]->SetLineColor(2);
				xn_sim[layer-1][sector-1][component-1]->SetLineWidth(1);
				xn_sim[layer-1][sector-1][component-1]->SetStats(0);
				xn_sim[layer-1][sector-1][component-1]->Draw("HIST");

				xn_dat[layer-1][sector-1][component-1]->SetLineColor(4);
				xn_dat[layer-1][sector-1][component-1]->SetMarkerColor(4);
				xn_dat[layer-1][sector-1][component-1]->SetMarkerStyle(8);
				xn_dat[layer-1][sector-1][component-1]->SetMarkerSize(1);
				xn_dat[layer-1][sector-1][component-1]->SetStats(0);
				xn_dat[layer-1][sector-1][component-1]->Draw("P,SAME");

				c->Print("xn_distributions.pdf");
				c->Close();

			}
		}
	}

	TCanvas * c_chi2 = new TCanvas("c_chi2","c_chi2",800,600);
	c_chi2->cd();
	chi2_vals->Draw("HIST");
	c_chi2->Print("xn_distributions.pdf");
	c_chi2->Close();

	cTotal->Print("xn_distributions.pdf)");
	cTotal->Close();

	dat->Close();
	bac->Close();
	sim->Close();
	return 0;
}
