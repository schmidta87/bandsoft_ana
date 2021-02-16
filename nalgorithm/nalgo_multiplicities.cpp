#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "bandhit.h"
#include "kinematic_cuts.h"

using namespace std;
int main( int argc, char** argv){
	if( argc != 3 ){
		cerr << "Wrong number of arguments. Please instead use:\n";
		cerr << "\t./code [InputSkim] [OutputRootFile]\n";
		return 0;
	}

	// Load input file:
	TFile * inFile = new TFile(argv[1]);
	TTree * inTree = (TTree*) inFile->Get("tagged");
	// 	Neutron info:
	int 		nMult		= 0;
	TClonesArray * 	nHits 		= new TClonesArray("bandhit");
	bool 		goodneutron 	= false;
	int 		nleadindex	= -1;
	//	Neutron branches:
	inTree->SetBranchAddress("nMult"		,&nMult			);
	inTree->SetBranchAddress("nHits"		,&nHits			);
	inTree->SetBranchAddress("goodneutron"		,&goodneutron		);
	inTree->SetBranchAddress("nleadindex"		,&nleadindex		);

	// Create output file:
	TFile * outFile = new TFile(argv[2],"RECREATE");
	TTree * outTree = new TTree("nalgo","Investigating BAND neutron algo");
	int		new_nMult 	= 0;
	int		layer		= 0;
	double		xDiff		= 0;
	double		yDiff		= 0;
	double		zDiff		= 0;
	double		tDiff		= 0;
	outTree->Branch("new_nMult"			,&new_nMult	);
	outTree->Branch("layer"				,&layer		);
	outTree->Branch("xDiff"				,&xDiff		);
	outTree->Branch("yDiff"				,&yDiff		);
	outTree->Branch("zDiff"				,&zDiff		);
	outTree->Branch("tDiff"				,&tDiff		);


	int veto = 0;
	for( int event = 0; event < inTree->GetEntries(); ++event ){
		nMult = -1;
		nHits->Clear();
		goodneutron = false;
		nleadindex = -1;
		// New info:
		new_nMult = 0;
		layer = -1;
		xDiff = -666;
		yDiff = -666;
		zDiff = -666;
		tDiff = -666;

		inTree->GetEntry(event);

		int closest_idx = -1;
		double closest_dL = 1e10;
		for( int ihit = 0 ; ihit < nHits->GetEntriesFast() ; ++ihit ){
			bandhit * thishit = (bandhit*) nHits->At(ihit);
			if( thishit->getEdep() > DataAdcToMeVee*2 ){
				++new_nMult;
				
				if( thishit->getDL().Mag() < closest_dL ){
					closest_idx = ihit;
					closest_dL = thishit->getDL().Mag();
				}
			}

			
		}
		if( new_nMult == 0 ) continue;

		// Now for these events, roughly:
		// 	83.4% are 1-hit events
		// 	11.1% are 2-hit events -- up to this covers ~95% of all events
		// 	2.5%  are 3-hit events -- up to this covers ~97% of all events
		// 	1.2%  are 4-hit events -- up to this covers ~98% of all events
		// 	<1%   are 5-hit events
		
		if( new_nMult == 2 ){
			int layer1 = ((bandhit*)nHits->At(0) )->getLayer();
			int layer2 = ((bandhit*)nHits->At(1) )->getLayer();
			if( layer1 == 6 && layer2 == 6 ) layer=0;
			else if( layer1 == 6 || layer2 == 6) layer=1;
			else{ 
				layer=2; 
				xDiff = fabs( ((bandhit*)nHits->At(0))->getX() - ((bandhit*)nHits->At(1))->getX() );
				yDiff = fabs( ((bandhit*)nHits->At(0))->getY() - ((bandhit*)nHits->At(1))->getY() );
				zDiff = fabs( ((bandhit*)nHits->At(0))->getZ() - ((bandhit*)nHits->At(1))->getZ() );
				tDiff = fabs( ((bandhit*)nHits->At(0))->getTof() - ((bandhit*)nHits->At(1))->getTof() );

			}
		}

		outTree->Fill();

	}

	// Cleanup
	inFile->Close();
	outTree->Write();
	outFile->Close();
	return 1;
}
