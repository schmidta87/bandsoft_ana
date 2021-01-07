void compare_twomethods(TString inDat1, TString inBac1, TString inDat2, TString inBac2 ){

	// Define some function used
	void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);

	// Load TFiles
	TFile * inFileDat1 = new TFile(inDat1);
	TFile * inFileBac1 = new TFile(inBac1);
	TFile * inFileDat2 = new TFile(inDat2);
	TFile * inFileBac2 = new TFile(inBac2);

	// Get the TTrees
	TTree * inTreeDat1 = (TTree*) inFileDat1->Get("tagged");
	TTree * inTreeBac1 = (TTree*) inFileBac1->Get("tagged");
	TTree * inTreeDat2 = (TTree*) inFileDat2->Get("tagged");
	TTree * inTreeBac2 = (TTree*) inFileBac2->Get("tagged");

	// Get and set the background normalization
	TVector3 * datnorm1 = (TVector3*)inFileDat1->Get("bacnorm");
	TVector3 * bacnorm1 = (TVector3*)inFileBac1->Get("bacnorm");
	inTreeBac1->SetWeight( datnorm1->X() / bacnorm1->X() );
	TVector3 * datnorm2 = (TVector3*)inFileDat2->Get("bacnorm");
	TVector3 * bacnorm2 = (TVector3*)inFileBac2->Get("bacnorm");
	inTreeBac2->SetWeight( datnorm2->X() / bacnorm2->X() );

	// Define histograms we want to plot:
	TH1D * pn_dat1 = new TH1D("pn_dat1","pn_dat1",40,0.2,0.6);
	TH1D * pn_bac1 = new TH1D("pn_bac1","pn_bac1",40,0.2,0.6);
	TH1D * pn_dat2 = new TH1D("pn_dat2","pn_dat2",40,0.2,0.6);
	TH1D * pn_bac2 = new TH1D("pn_bac2","pn_bac2",40,0.2,0.6);
	
	// Draw the full pn distribution
	TCanvas * c_pn = new TCanvas("c_pn","",800,600);
	//c_pn->Divide(2,2);
	c_pn->cd(1);
	inTreeDat1->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_dat1","tag[nleadindex]->getMomentumN().Mag() > 0.3");
	inTreeBac1->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_bac1","tag[nleadindex]->getMomentumN().Mag() > 0.3");
	inTreeDat2->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_dat2","tag[nleadindex]->getMomentumN().Mag() > 0.3");
	inTreeBac2->Draw("tag[nleadindex]->getMomentumN().Mag() >> pn_bac2","tag[nleadindex]->getMomentumN().Mag() > 0.3");


	pn_dat1->Add( pn_bac1 , -1 );
	pn_dat2->Add( pn_bac2 , -1 );

	label1D(pn_dat1,pn_dat2,"|p_{n}| [GeV/c]","Counts");

	c_pn->SaveAs("full_pn.pdf");



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

	sim->Draw("hist");
	data->Draw("p,same");

	double max1 = data->GetMaximum()*1.1;
	double max2 = sim->GetMaximum()*1.1;
	sim->GetYaxis()->SetRangeUser(0,max(max1,max2));
	
	sim->GetXaxis()->SetTitle(xlabel);
	sim->GetYaxis()->SetTitle(ylabel);

	TLegend * legend = new TLegend(0.7,0.8,0.9,0.9);
	legend->AddEntry(data,"Data New Method","f");
	legend->AddEntry(sim,"Data Old Method","f");
	legend->Draw();

	return;
}

