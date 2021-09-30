void tof_raw(TString inputFileNameDat, TString inputFileNameMix){

	
	TFile inFileDat = TFile(inputFileNameDat);
	TFile inFileMix = TFile(inputFileNameMix);

	TH1D *ToF_dat = (TH1D*) inFileDat.Get("hToF_full");
	TH1D *ToF_mix = (TH1D*) inFileMix.Get("hToF_full");

	TVector3 *datnorm = (TVector3*) inFileDat.Get("bacnorm");
	TVector3 *bacnorm = (TVector3*) inFileMix.Get("bacnorm");

	ToF_mix->Scale( datnorm->Y() / bacnorm->X() );


	TCanvas * c = new TCanvas("c","",800,600);
	c->cd();
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.15);

	ToF_dat->Rebin(16);
	ToF_mix->Rebin(16);

	ToF_dat->SetMinimum(0);
	ToF_mix->SetMinimum(0);

	ToF_dat->GetXaxis()->SetLabelSize(0.05);
	ToF_dat->GetYaxis()->SetLabelSize(0.05);


	ToF_dat->SetTitle("");
	ToF_dat->GetXaxis()->SetTitle("ToF [ns]");
	ToF_dat->GetXaxis()->SetTitleSize(0.65);
	ToF_dat->GetYaxis()->SetTitle("Counts");
	ToF_dat->GetYaxis()->SetTitleSize(0.65);

	ToF_dat->GetXaxis()->SetRangeUser(12,100);
	ToF_mix->GetXaxis()->SetRangeUser(12,100);

	ToF_dat->SetLineWidth(3);
	ToF_mix->SetLineWidth(3);
	ToF_mix->SetLineColor(2);

	ToF_dat->SetStats(0);
	ToF_mix->SetStats(0);

	ToF_dat->Draw();
	ToF_mix->Draw("H,SAME");

	for(int i = 1 ; i <= ToF_dat->GetXaxis()->GetNbins(); ++i ){
		if( ToF_dat->GetXaxis()->GetBinCenter(i) < 12 ) continue;
		if( ToF_dat->GetXaxis()->GetBinCenter(i) > 100 ) continue;
		cout << ToF_dat->GetXaxis()->GetBinCenter(i) << " " << ToF_dat->GetBinContent(i) << " " << ToF_mix->GetBinContent(i) << "\n";
	}

	c->SaveAs("ToFHist.pdf");


	inFileDat.Close();
	inFileMix.Close();
	return;
}
