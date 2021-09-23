void aspT(TString inDat, TString inBac, TString inSim){

	// Define some function used
	void label2D(TH2D* data, TH2D* sim, TString xlabel, TString ylabel);
	void label2D_ratio(TH2D* data, TH2D* sim, TString xlabel, TString ylabel, double ymin , double ymax );

	// Load TFiles
	TFile * inFileDat = new TFile(inDat);
	TFile * inFileBac = new TFile(inBac);
	TFile * inFileSim = new TFile(inSim);

	// Get the TTrees
	TTree * inTreeDat = (TTree*) inFileDat->Get("tagged");
	TTree * inTreeBac = (TTree*) inFileBac->Get("tagged");
	TTree * inTreeSim = (TTree*) inFileSim->Get("tagged");

	// Get and set the background normalization
	TVector3 * datnorm = (TVector3*)inFileDat->Get("bacnorm");
	TVector3 * bacnorm = (TVector3*)inFileBac->Get("bacnorm");
	inTreeBac->SetWeight( datnorm->Z() / bacnorm->X() );

	// Define histograms we want to plot:
	TH2D * aspt_dat = new TH2D(Form("aspt_dat"),"",30,1.35,1.65,20,0.,0.1);
	TH2D * aspt_bac = new TH2D(Form("aspt_bac"),"",30,1.35,1.65,20,0.,0.1);
	TH2D * aspt_sim = new TH2D(Form("aspt_sim"),"",30,1.35,1.65,20,0.,0.1);

	// Draw the full aspt distribution
	double sim_scaling = 0;
	TCut pTcut = "tag[nleadindex]->getPt().Mag() < 0.1";

	TCanvas * c_aspt = new TCanvas("c_aspt","",800,600);

	c_aspt->cd(1);
	inTreeDat->Draw(Form("tag[nleadindex]->getPt().Mag() : tag[nleadindex]->getAs() >> aspt_dat"),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut );
	inTreeBac->Draw(Form("tag[nleadindex]->getPt().Mag() : tag[nleadindex]->getAs() >> aspt_bac"),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut );
	inTreeSim->Draw(Form("tag[nleadindex]->getPt().Mag() : tag[nleadindex]->getAs() >> aspt_sim"),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut );

	// Background subraction
	aspt_dat->Add(aspt_bac,-1);

	c_aspt->SetRightMargin(0.09);
	c_aspt->SetLeftMargin(0.15);
	c_aspt->SetBottomMargin(0.15);
	
	aspt_dat->SetTitle("");
	label2D(aspt_dat,aspt_sim,"#alpha_{S}","p_{T} [GeV/c]");
	c_aspt->SaveAs("as_vs_pt.pdf");



	return;
}

void label2D(TH2D* data, TH2D* sim, TString xlabel, TString ylabel){
	data->SetLineColor(4);
	data->SetMarkerColor(4);
	data->SetMarkerStyle(8);
	data->SetMarkerSize(1);
	data->SetStats(0);
	data->GetXaxis()->SetTickLength(0.05);
	data->GetYaxis()->SetTickLength(0.05);
	data->SetMinimum(0);

	sim->SetLineColor(2);
	sim->SetLineWidth(1);
	sim->SetStats(0);


	//sim->Draw("hist");
	data->Draw("hist,col");
	data->GetXaxis()->SetTitle(xlabel);
	data->GetYaxis()->SetTitle(ylabel);
	data->SetTitleSize(0.07,"xyz");

	//double max1 = data->GetMaximum()*1.1;
	//double max2 = sim->GetMaximum()*1.1;
	//sim->GetYaxis()->SetRangeUser(0,max(max1,max2));
	

	//TLegend * legend = new TLegend(0.7,0.8,0.9,0.9);
	//legend->AddEntry(data,"Radiation On","f");
	//legend->AddEntry(sim,"Radiation Off","f");
	//legend->AddEntry(data,"Data","f");
	//legend->AddEntry(sim,"Sim","f");
	//legend->Draw();

	return;
}

void label2D_ratio(TH2D* data, TH2D* sim, TString xlabel, TString ylabel, double ymin , double ymax ){
	gStyle->SetOptFit(1);
	
	TH1D * data_copy = (TH1D*) data->Clone();
	TH1D * sim_copy = (TH1D*) sim->Clone();
	
	data_copy->SetLineColor(1);
	data_copy->SetLineWidth(3);
	//data_copy->SetStats(0);

	sim_copy->SetLineColor(9);
	sim_copy->SetLineWidth(3);
	//sim_copy->SetStats(0);
	//sim_copy->Scale(data_copy->Integral() / sim_copy->Integral() );

	data_copy->Divide(sim_copy);
	data_copy->SetTitle( sim_copy->GetTitle() );
	data_copy->Draw("ep");
	TLine* line = new TLine(data_copy->GetXaxis()->GetBinCenter(1), 1., data_copy->GetXaxis()->GetBinCenter(data_copy->GetXaxis()->GetNbins()), 1.);
	line->SetLineWidth(3);
	line->SetLineColor(2);
	line->Draw("same");

	data_copy->GetYaxis()->SetRangeUser(ymin,ymax);
	
	data_copy->GetXaxis()->SetTitle(xlabel);
	data_copy->GetYaxis()->SetTitle(ylabel);

	return;
}
