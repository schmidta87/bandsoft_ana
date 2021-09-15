#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TVector3.h"

using namespace std;
int main(int argc, char ** argv){

	if (argc < 1){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [inputDatafiles]\n";
		return -1;
	}
	
	double total_normalization_x = 0;
	double total_normalization_y = 0;
	double total_normalization_z = 0;
	for( int fi = 1 ; fi < argc ; fi++){
		cout << argv[fi] << "\n";
		TFile * thisFile = new TFile(argv[fi]);

		TVector3 * thisNorm = (TVector3*)thisFile->Get("bacnorm");
		total_normalization_x	+=	thisNorm->X();
		total_normalization_y	+=	thisNorm->Y();
		total_normalization_z	+=	thisNorm->Z();

		delete thisNorm;
		thisFile->Close();
	}
	cout.precision(10);
	cout << total_normalization_x << "\n";
	cout << total_normalization_y << "\n";
	cout << total_normalization_z << "\n";
	


	return 0;
}
