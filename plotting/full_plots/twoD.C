void twoD(TString inDat, TString inBac, TString inSim){

	TCut pNcut = "tag[nleadindex]->getMomentumN().Mag() < 0.32 && tag[nleadindex]->getMomentumN().Mag() > 0.25 && !(nHits[nleadindex]->getSector()==1 && nHits[nleadindex]->getComponent()==1) && nHits[nleadindex]->getEdep()>10";
	TCut pNcut_sim = "tag_smeared[nleadindex]->getMomentumN().Mag() < 0.32 && tag_smeared[nleadindex]->getMomentumN().Mag() > 0.25 && !(nHits[nleadindex]->getSector()==1 && nHits[nleadindex]->getComponent()==1) && nHits[nleadindex]->getEdep()>10";
	TCut pTcut[3] = {"tag[nleadindex]->getPt().Mag() < 0.2","tag[nleadindex]->getPt().Mag()< 0.1","tag[nleadindex]->getPt().Mag()>0.1 && tag[nleadindex]->getPt().Mag()<0.2"};
	TCut pTcut_sim[3] = {"tag_smeared[nleadindex]->getPt().Mag() < 0.2",
	     			"tag_smeared[nleadindex]->getPt().Mag()< 0.1",
				"tag_smeared[nleadindex]->getPt().Mag()>0.1 && tag_smeared[nleadindex]->getPt().Mag()<0.2"};

	// Define some function used
	void label2D(TH2D* hist, TString xlabel, TString ylabel);
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
	inTreeBac->SetWeight( datnorm->X() / bacnorm->X() );

	// Define histograms we want to plot:
	TH2D ** twoD_dat = new TH2D*[3];
	TH2D ** twoD_bac = new TH2D*[3];
	TH2D ** twoD_sim = new TH2D*[3];
	for(int i = 0 ; i < 3 ; i++){
		twoD_dat[i] = new TH2D(Form("twoD_dat_%i",i),"",30,-120,120,30,-50,100);
		twoD_bac[i] = new TH2D(Form("twoD_bac_%i",i),"",30,-120,120,30,-50,100);
		twoD_sim[i] = new TH2D(Form("twoD_sim_%i",i),"",30,-120,120,30,-50,100);
	}

	// Draw the full twoD distribution
	TCanvas * c_twoD = new TCanvas("c_twoD","",800,600);
	double sim_scaling = 0;
	c_twoD->Divide(3,3);
	for( int i = 0 ; i < 3 ; i++){
		TString pTtitle = "Full pT";
		if( i == 1 ){
			pTtitle = "Low pT";
		}
		if( i == 2 ){
			pTtitle = "High pT";
		}

		c_twoD->cd(i+1);
		inTreeDat->Draw(Form("nHits[nleadindex]->getDL().Y() : nHits[nleadindex]->getDL().X() >> twoD_dat_%i",i),pNcut && pTcut[i] );
		inTreeBac->Draw(Form("nHits[nleadindex]->getDL().Y() : nHits[nleadindex]->getDL().X() >> twoD_bac_%i",i),pNcut && pTcut[i] );
		inTreeSim->Draw(Form("nHits[nleadindex]->getDL().Y() : nHits[nleadindex]->getDL().X() >> twoD_sim_%i",i),pNcut_sim && pTcut_sim[i] );

		// Background subraction
		twoD_dat[i]->Add(twoD_bac[i],-1);

		// Simulation scaling only from no pT cut distribution (i.e. from full distribution)
		double full_simnorm = (double)twoD_dat[i]->Integral() / twoD_sim[i]->Integral();
		if( i == 0 ) sim_scaling = full_simnorm;
		twoD_sim[i]->Scale( full_simnorm );
		
		
		c_twoD->cd(i+1);
		twoD_dat[i]->SetTitle(pTtitle+Form(", Data, C_{sim} = %f, ",full_simnorm));
		label2D(twoD_dat[i],"x [cm]","y [cm]");

		c_twoD->cd(i+4);
		twoD_sim[i]->SetTitle(pTtitle+Form(", Sim, C_{sim} = %f, ",full_simnorm));
		label2D(twoD_sim[i],"x [cm]","y [cm]");

		c_twoD->cd(7+i);
		label2D_ratio(twoD_dat[i],twoD_sim[i],"x [cm]","y [cm]",0,2);
	}


	c_twoD->SaveAs("full_twoD.pdf");

	return;
}

void label2D(TH2D* hist, TString xlabel, TString ylabel){

	hist->SetStats(0);

	hist->Draw("colz");

	
	hist->GetXaxis()->SetTitle(xlabel);
	hist->GetYaxis()->SetTitle(ylabel);

	return;
}

void label2D_ratio(TH2D* data, TH2D* sim, TString xlabel, TString ylabel, double ymin , double ymax ){
	gStyle->SetOptFit(1);
	
	TH2D * data_copy = (TH2D*) data->Clone();
	TH2D * sim_copy = (TH2D*) sim->Clone();
	
	data_copy->SetLineColor(1);
	data_copy->SetLineWidth(3);
	//data_copy->SetStats(0);

	sim_copy->SetLineColor(9);
	sim_copy->SetLineWidth(3);
	//sim_copy->SetStats(0);
	//sim_copy->Scale(data_copy->Integral() / sim_copy->Integral() );

	data_copy->Divide(sim_copy);
	data_copy->SetTitle( "Data/Sim Ratio" );

	data_copy->SetMaximum(3);
	data_copy->SetMinimum(-1);
	data_copy->Draw("colz");
	//TLine* line = new TLine(data_copy->GetXaxis()->GetBinCenter(1), 1., data_copy->GetXaxis()->GetBinCenter(data_copy->GetXaxis()->GetNbins()), 1.);
	//line->SetLineWidth(3);
	//line->SetLineColor(2);
	//line->Draw("same");

	//data_copy->GetYaxis()->SetRangeUser(ymin,ymax);
	
	data_copy->GetXaxis()->SetTitle(xlabel);
	data_copy->GetYaxis()->SetTitle(ylabel);

	return;
}
