void Q2_W2_inc(TString inDat, TString inSim){

	// Define some function used
	void label2D(TH2D* data, TString xlabel, TString ylabel);

	// Load TFiles
	TFile * inFileDat = new TFile(inDat);
	TFile * inFileSim = new TFile(inSim);

	// Get the TTrees
	TTree * inTreeDat = (TTree*) inFileDat->Get("electrons");
	TTree * inTreeSim = (TTree*) inFileSim->Get("electrons");

	// Define histograms we want to plot:
	TH2D * Q2_W2_dat = new TH2D("Q2_W2_dat","",150,2,8,90,2,3.8);

	// Draw the full Q2 distribution
	TCanvas * c_Q2_W2 = new TCanvas("c_Q2_W2","",800,600);

	c_Q2_W2->cd(1);
	inTreeDat->Draw("sqrt(eHit->getW2()) : eHit->getQ2() >> Q2_W2_dat");

	Q2_W2_dat->SetTitle("Inclusive Kinematics");
	label2D(Q2_W2_dat,"Q^{2}","W");

	c_Q2_W2->SaveAs("full_Q2_W2-inc.pdf");

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
