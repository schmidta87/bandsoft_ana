void xB_Q2_inc(TString inDat, TString inSim){

	// Define some function used
	void label2D(TH2D* data, TString xlabel, TString ylabel);

	// Load TFiles
	TFile * inFileDat = new TFile(inDat);
	TFile * inFileSim = new TFile(inSim);

	// Get the TTrees
	TTree * inTreeDat = (TTree*) inFileDat->Get("electrons");
	TTree * inTreeSim = (TTree*) inFileSim->Get("electrons");

	// Define histograms we want to plot:
	TH2D * xB_Q2_dat = new TH2D("xB_Q2_dat","",150,0.1,0.7,150,2,8);

	// Draw the full Q2 distribution
	TCanvas * c_xB_Q2 = new TCanvas("c_xB_Q2","",800,600);

	c_xB_Q2->cd(1);
	inTreeDat->Draw("eHit->getQ2() : eHit->getXb() >> xB_Q2_dat");

	xB_Q2_dat->SetTitle("Inclusive Kinematics");
	label2D(xB_Q2_dat,"x_{B}","Q^{2}");

	c_xB_Q2->SaveAs("full_xB_Q2-inc.pdf");

	return;
}

void label2D(TH2D* data, TString xlabel, TString ylabel){
	//data->SetLineColor(4);
	//data->SetMarkerColor(4);
	//data->SetMarkerStyle(8);
	//data->SetMarkerSize(1);
	data->SetStats(0);

	data->Draw("colz");

	
	data->GetXaxis()->SetTitle(xlabel);
	data->GetYaxis()->SetTitle(ylabel);


	return;
}
