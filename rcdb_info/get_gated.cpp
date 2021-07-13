#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH1.h"
#include "TClonesArray.h"


#include "RCDB/Connection.h"

#include "constants.h"

using namespace std;


int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 2 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [inputFiles] \n\n";
		cerr << "\t\t[inputFile] = ____.root ____.root ____.root ...\n\n";
		return -1;
	}

	double total_charge = 0;
	// Load input file
	for( int i = 1 ; i < argc ; i++ ){
		
		TFile * inFile = new TFile(argv[i]);
		TTree * inTree = (TTree*)inFile->Get("tagged");

		double gated_charge = 0;
		inTree->SetBranchAddress("gated_charge",	&gated_charge);


		int nEvents = inTree->GetEntries();
		inTree->GetEntry( nEvents-1 );
		cout << argv[i] << " " << gated_charge << "\n";

		total_charge += gated_charge;


	}// end loop over files
	cout << total_charge << "\n";

	return 0;
}



