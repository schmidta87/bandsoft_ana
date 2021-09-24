#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TVector3.h"
#include "TH1.h"
#include "TChain.h"

using namespace std;
int main(int argc, char ** argv){

	if (argc < 3){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [opt 0 (MC) 1 (dat) 2 (bac) ] [inputDatafiles]\n";
		return -1;
	}
	int MC_DATA_OPT = atoi(argv[1]);
	TFile * output = NULL;
	if( MC_DATA_OPT == 0 ) output = new TFile("SIM-CHAINED.root","RECREATE");
	if( MC_DATA_OPT == 1 ) output = new TFile("DAT-CHAINED.root","RECREATE");
	if( MC_DATA_OPT == 2 ) output = new TFile("MIX-CHAINED.root","RECREATE");
	
	TChain * tagged = new TChain("tagged");
	TH1D * hToF_full = new TH1D("hToF_full","hToF_full",8000,-100,300);
	TH1D * hThetaNQ = new TH1D("hThetaNQ","hThetaNQ",8000,140,180);
	TH1D * hThetaN 	= new TH1D("hThetaN","hThetaN",8000,140,180);
	TH1D * hPe = new TH1D("hPe","hPe",4000,3,7);
	

	double total_normalization_x = 0;
	double total_normalization_y = 0;
	double total_normalization_z = 0;
	for( int fi = 2 ; fi < argc ; fi++){
		cout << argv[fi] << "\n";
		// Add file to chain:
		tagged->Add(argv[fi]);

		// Open file to grab histo and norm vector:
		TFile * thisFile = new TFile(argv[fi]);

		TH1D * thisToF = (TH1D*)thisFile->Get("hToF_full");
		TH1D * thisThetaNQ = (TH1D*)thisFile->Get("hThetaNQ");
		TH1D * thisThetaN = (TH1D*)thisFile->Get("hThetaN");
		TH1D * thisPe = (TH1D*)thisFile->Get("hPe");
		
		hToF_full->Add( thisToF , 1 );
		hThetaNQ->Add( thisThetaNQ , 1);
		hThetaN->Add( thisThetaN, 1 );
		hPe->Add( thisPe , 1 );
		
		if( MC_DATA_OPT > 0 ){
			TVector3 * thisNorm = (TVector3*)thisFile->Get("bacnorm");
			total_normalization_x	+=	thisNorm->X();
			total_normalization_y	+=	thisNorm->Y();
			total_normalization_z	+=	thisNorm->Z();
			delete thisNorm;
		}	

		delete thisToF;
		delete thisThetaNQ;
		delete thisThetaN;
		delete thisPe;
		thisFile->Close();
	}

	output->cd();
	tagged->CloneTree(-1,"fast");
	
	// Set the background normalization
	TVector3 * bacnorm = new TVector3();
	bacnorm->SetXYZ( total_normalization_x, total_normalization_y, total_normalization_z );
	bacnorm->Write("bacnorm");

	hToF_full->Write("hToF_full");
	hThetaNQ->Write("hThetaNQ");
	hThetaN->Write("hThetaN");
	hPe->Write("hPe");

	output->Write();
	output->Close();

	return 0;
}
