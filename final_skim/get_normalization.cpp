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
			<< "\t./code [opt] [inputDatafiles]\n"
			<< "\t\t [opt==0: data, getZ()]\n"
			<< "\t\t [opt==1: bkgd, getX()]\n";
		return -1;
	}
	int opt = atoi(argv[1]);
	
	double total_normalization = 0;
	for( int fi = 2 ; fi < argc ; fi++){
		cout << argv[fi] << "\n";
		TFile * thisFile = new TFile(argv[fi]);

		TVector3 * thisNorm = (TVector3*)thisFile->Get("bacnorm");
		if( opt == 0)
			total_normalization	+=	thisNorm->Z();
		else if(opt == 1)
			total_normalization	+=	thisNorm->X();

		delete thisNorm;
		thisFile->Close();
	}
	cout.precision(10);
	cout << total_normalization << "\n";
	


	return 0;
}
