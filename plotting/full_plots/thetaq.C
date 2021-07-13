void thetaq(TString inDat, TString inBac, TString inSim){

	// Define some function used
	void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
	void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel, double ymin , double ymax );

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
	inTreeBac->SetWeight( datnorm->X() / bacnorm->X() );

	// Define histograms we want to plot:
	TH1D ** thetaq_dat = new TH1D*[3];
	TH1D ** thetaq_bac = new TH1D*[3];
	TH1D ** thetaq_sim = new TH1D*[3];
	for(int i = 0 ; i < 3 ; i++){
		thetaq_dat[i] = new TH1D(Form("thetaq_dat_%i",i),"",40,0,20);
		thetaq_bac[i] = new TH1D(Form("thetaq_bac_%i",i),"",40,0,20);
		thetaq_sim[i] = new TH1D(Form("thetaq_sim_%i",i),"",40,0,20);
	}

	// Draw the full thetaq distribution
	TCanvas * c_thetaq = new TCanvas("c_thetaq","",800,600);
	double sim_scaling = 0;
	c_thetaq->Divide(3,2);
	for( int i = 0 ; i < 3 ; i++){
		TCut pTcut = "";
		TString pTtitle = "Full pT";
		if( i == 1 ){
			pTcut = "tag[nleadindex]->getPt().Mag() < 0.1";
			pTtitle = "Low pT";
		}
		if( i == 2 ){
			pTcut = "tag[nleadindex]->getPt().Mag() >= 0.1";
			pTtitle = "High pT";
		}

		c_thetaq->cd(i+1);
		inTreeDat->Draw(Form("tag[nleadindex]->getMomentumQ().Theta()*180./TMath::Pi() >> thetaq_dat_%i",i),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut);
		inTreeBac->Draw(Form("tag[nleadindex]->getMomentumQ().Theta()*180./TMath::Pi() >> thetaq_bac_%i",i),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut);
		inTreeSim->Draw(Form("tag[nleadindex]->getMomentumQ().Theta()*180./TMath::Pi() >> thetaq_sim_%i",i),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut);

		// Background subraction
		thetaq_dat[i]->Add(thetaq_bac[i],-1);

		// Simulation scaling only from no pT cut distribution (i.e. from full distribution)
		double full_simnorm = (double)thetaq_dat[0]->Integral() / thetaq_sim[0]->Integral();
		if( i == 0 ) sim_scaling = full_simnorm;
		thetaq_sim[i]->Scale( sim_scaling );
		
		
		thetaq_sim[i]->SetTitle(pTtitle+Form(", C_{sim} = %f, ",sim_scaling));
		label1D(thetaq_dat[i],thetaq_sim[i],"Theta_{q} [Deg.]","Counts");

		c_thetaq->cd(4+i);
		label1D_ratio(thetaq_dat[i],thetaq_sim[i],"Theta_{q} [Deg.]","Data/Sim",0,2);
	}


	c_thetaq->SaveAs("full_thetaq.pdf");

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
	//legend->AddEntry(data,"Radiation On","f");
	//legend->AddEntry(sim,"Radiation Off","f");
	legend->AddEntry(data,"Data","f");
	legend->AddEntry(sim,"Sim","f");
	legend->Draw();

	return;
}

void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel, double ymin , double ymax ){
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
