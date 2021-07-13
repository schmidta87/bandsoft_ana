#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TVector3.h"

using namespace std;
int main(int argc, char ** argv){

	if (argc < 2){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [inputDatafiles]\n";
		return -1;
	}
	
	double total_normalization = 0;
	for( int fi = 1 ; fi < argc ; fi++){
		TFile * thisFile = new TFile(argv[fi]);

		TVector3 * thisNorm = (TVector3*)thisFile->Get("bacnorm");
		total_normalization	+=	thisNorm->X();

		delete thisNorm;
		thisFile->Close();
	}
	cout.precision(10);
	cout << total_normalization << "\n";
	


	return 0;
}
