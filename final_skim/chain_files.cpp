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

	if (argc < 2){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [inputDatafiles]\n";
		return -1;
	}

	TFile * output = new TFile("CHAINED.root","RECREATE");
	TChain * tagged = new TChain("tagged");
	TH1D * hToF_full = new TH1D("hToF_full","hToF_full",8000,-100,300);
	

	
	double total_normalization_x = 0;
	double total_normalization_y = 0;
	double total_normalization_z = 0;
	for( int fi = 1 ; fi < argc ; fi++){
		cout << argv[fi] << "\n";
		// Add file to chain:
		tagged->Add(argv[fi]);

		// Open file to grab histo and norm vector:
		TFile * thisFile = new TFile(argv[fi]);

		TH1D * thisToF = (TH1D*)thisFile->Get("hToF_full");
		hToF_full->Add( thisToF , 1 );

		TVector3 * thisNorm = (TVector3*)thisFile->Get("bacnorm");
		total_normalization_x	+=	thisNorm->X();
		total_normalization_y	+=	thisNorm->Y();
		total_normalization_z	+=	thisNorm->Z();

		delete thisNorm;
		delete thisToF;
		thisFile->Close();
	}

	output->cd();
	tagged->CloneTree(-1,"fast");
	
	// Set the background normalization
	TVector3 * bacnorm = new TVector3();
	bacnorm->SetXYZ( total_normalization_x, total_normalization_y, total_normalization_z );
	bacnorm->Write("bacnorm");

	hToF_full->Write("hToF_full");

	output->Write();
	output->Close();

	return 0;
}
