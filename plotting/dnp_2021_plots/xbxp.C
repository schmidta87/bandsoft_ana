void xbxp(TString inDat, TString inBac, TString inSim){
	const int 	bins_As 		= 3;
	const double 	As_step 		= 0.1;
	const double 	As_min			= 1.3;
	const double 	As_max 			= As_min + As_step*bins_As; 


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
	inTreeBac->SetWeight( datnorm->X() / bacnorm->X() );

	// Define histograms we want to plot:
	TH2D ** xbxp_dat = new TH2D*[3];
	TH2D ** xbxp_bac = new TH2D*[3];
	TH2D ** xbxp_sim = new TH2D*[3];
	for(int i = 0 ; i < 3 ; i++){
		xbxp_dat[i] = new TH2D(Form("xbxp_dat_%i",i),"",40,0.1,0.5,70,0.1,0.8);
		xbxp_bac[i] = new TH2D(Form("xbxp_bac_%i",i),"",40,0.1,0.5,70,0.1,0.8);
		xbxp_sim[i] = new TH2D(Form("xbxp_sim_%i",i),"",40,0.1,0.5,70,0.1,0.8);
	}

	// Draw the full xbxp distribution
	double sim_scaling = 0;
	TCut pTcut = "tag[nleadindex]->getPt().Mag() < 0.1";
	for( int i = 0 ; i < bins_As ; i++){
		TCanvas * c_xbxp = new TCanvas("c_xbxp","",800,600);
		double this_min_as = As_min + i*(As_max - As_min)/bins_As;
		double this_max_as = As_min + (i+1)*(As_max - As_min)/bins_As;
		TCut this_as_cut = Form("tag[nleadindex]->getAs() > %f && tag[nleadindex]->getAs() < %f",this_min_as,this_max_as);

		c_xbxp->cd(1);
		inTreeDat->Draw(Form("tag[nleadindex]->getXp() : eHit->getXb() >> xbxp_dat_%i",i),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut && this_as_cut);
		inTreeBac->Draw(Form("tag[nleadindex]->getXp() : eHit->getXb() >> xbxp_bac_%i",i),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut && this_as_cut);
		inTreeSim->Draw(Form("tag[nleadindex]->getXp() : eHit->getXb() >> xbxp_sim_%i",i),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut && this_as_cut);

		// Background subraction
		xbxp_dat[i]->Add(xbxp_bac[i],-1);

		// Simulation scaling only from no pT cut distribution (i.e. from full distribution)
		//double full_simnorm = (double)xbxp_dat[0]->Integral() / xbxp_sim[0]->Integral();
		//if( i == 0 ) sim_scaling = full_simnorm;
		//xbxp_sim[i]->Scale( sim_scaling );
		
		
		//xbxp_sim[i]->SetTitle(pTtitle+Form(", C_{sim} = %f, ",sim_scaling));
		c_xbxp->SetRightMargin(0.09);
		c_xbxp->SetLeftMargin(0.15);
		c_xbxp->SetBottomMargin(0.15);
		
		TString current_title = Form("%.2f < #alpha_{S} < %.2f",this_min_as,this_max_as);
		xbxp_dat[i]->SetTitle(current_title);
		label2D(xbxp_dat[i],xbxp_sim[i],"x_{B}","x'");
		c_xbxp->SaveAs(Form("xb_vs_xp_asBin%i.pdf",i));

		//c_xbxp->cd(4+i);
		//label2D_ratio(xbxp_dat[i],xbxp_sim[i],"x'","Data/Sim",0,2);
	}



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
