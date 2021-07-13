void virt(TString inDat, TString inBac, TString inSim){

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
	TH1D ** virt_dat = new TH1D*[3];
	TH1D ** virt_bac = new TH1D*[3];
	TH1D ** virt_sim = new TH1D*[3];
	for(int i = 0 ; i < 3 ; i++){
		virt_dat[i] = new TH1D(Form("virt_dat_%i",i),"",25,-0.6,-0.1);
		virt_bac[i] = new TH1D(Form("virt_bac_%i",i),"",25,-0.6,-0.1);
		virt_sim[i] = new TH1D(Form("virt_sim_%i",i),"",25,-0.6,-0.1);
	}

	// Draw the full virt distribution
	TCanvas * c_virt = new TCanvas("c_virt","",800,600);
	double sim_scaling = 0;
	c_virt->Divide(3,2);
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

		c_virt->cd(i+1);
		inTreeDat->Draw(Form("( pow(1.8756 - sqrt( pow(tag[nleadindex]->getMomentumN().Mag(),2) + pow(0.9396 ,2) ) ,2) - pow( tag[nleadindex]->getMomentumN().Mag() , 2 ) - pow(0.93827,2) )/pow(0.93827,2) >> virt_dat_%i",i),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut);
		inTreeBac->Draw(Form("( pow(1.8756 - sqrt( pow(tag[nleadindex]->getMomentumN().Mag(),2) + pow(0.9396 ,2) ) ,2) - pow( tag[nleadindex]->getMomentumN().Mag() , 2 ) - pow(0.93827,2) )/pow(0.93827,2) >> virt_bac_%i",i),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut);
		inTreeSim->Draw(Form("( pow(1.8756 - sqrt( pow(tag[nleadindex]->getMomentumN().Mag(),2) + pow(0.9396 ,2) ) ,2) - pow( tag[nleadindex]->getMomentumN().Mag() , 2 ) - pow(0.93827,2) )/pow(0.93827,2) >> virt_sim_%i",i),"tag[nleadindex]->getMomentumN().Mag() > 0.3" && pTcut);

		// Background subraction
		virt_dat[i]->Add(virt_bac[i],-1);

		// Simulation scaling only from no pT cut distribution (i.e. from full distribution)
		double full_simnorm = (double)virt_dat[0]->Integral() / virt_sim[0]->Integral();
		if( i == 0 ) sim_scaling = full_simnorm;
		virt_sim[i]->Scale( sim_scaling );
		
		
		virt_sim[i]->SetTitle(pTtitle+Form(", C_{sim} = %f, ",sim_scaling));
		label1D(virt_dat[i],virt_sim[i],"Virtuality^{2}","Counts");

		c_virt->cd(4+i);
		label1D_ratio(virt_dat[i],virt_sim[i],"Virtuality^{2}","Data/Sim",0,2);
	}


	c_virt->SaveAs("full_virt.pdf");

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
					
	TLegend * legend = new TLegend(0.1,0.8,0.3,0.9);
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
