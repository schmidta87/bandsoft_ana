#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unordered_map>
#include <vector>

#include "bandhit.h"
#include "clashit.h"
#include "taghit.h"

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TChain.h"
#include "TVectorT.h"
#include "TRandom3.h"
#include "constants.h"
#include "kinematic_cuts.h"
#include "TClonesArray.h"


using namespace std;

//clashit (const clashit &input) {x = p2.x; y = p2.y; }


int main(int argc, char ** argv){
	if (argc != 4){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./event_mixing [outputRootfile] [inputElectronSkim] [inputTaggedSkim]\n";
		return -1;
	}

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("tagged","BAND Neutrons and CLAS Electrons");
	//	Event info:
	int Runno		= 0;
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	bool goodneutron = false;
	int nleadindex = -1;
	double weight		= 0;
	//	MC info:
	int genMult		= 0;
	TClonesArray * mcParts = new TClonesArray("genpart");
	TClonesArray &saveMC = *mcParts;
	// 	Neutron info:
	int nMult		= 0;
	TClonesArray * nHits = new TClonesArray("bandhit");
	TClonesArray &saveHit = *nHits;
	//	Electron info:
	clashit eHit;
	//	Tagged info:
	TClonesArray * tags = new TClonesArray("taghit");
	TClonesArray &saveTags = *tags;
	// 	Event branches:
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("current"	,&current		);
	outTree->Branch("weight"	,&weight		);
	//	Neutron branches:
	outTree->Branch("nMult"		,&nMult			);
	outTree->Branch("nHits"		,&nHits			);
	//Branches to store if good Neutron event and leadindex
	outTree->Branch("goodneutron"		,&goodneutron	);
	outTree->Branch("nleadindex"		,&nleadindex			);
	//	Electron branches:
	outTree->Branch("eHit"		,&eHit			);
	//	Tagged branches:
	outTree->Branch("tag"		,&tags			);
	//	MC branches:
	outTree->Branch("genMult"	,&genMult		);
	outTree->Branch("mcParts"	,&mcParts		);

	// Setup rand for picking particles
	TRandom3 * myRand = new TRandom3(0);


	// Define 4D binned ToF plots and bin edges and kin cuts from kinematic_cuts.h file
	const double max_CosTheta_nq 	= cos(NCUT_THETANQ_min);	// bit non-intuitive but min angle = 150 which is a "max" costhetanq of -0.8
	const double min_CosTheta_nq 	= cos(NCUT_THETANQ_max);	// and then max angle = 180 which is a "min" costhetanq of -1

	// Define background and signal edges from from kinematic_cuts.h file
	const double bkgrd_min 		= NCUT_BACK_Tof_min;
	const double bkgrd_max 		= NCUT_BACK_Tof_max;
	const double signal_min 	= NCUT_Tof_min;
	const double signal_max 	= NCUT_Tof_max;

	// Loop over the neutron and electron skim files given
	TFile * inFile_e = new TFile(argv[2]);
	TTree * inTree_e = (TTree*)inFile_e->Get("electrons");
		// combine all the neutron files
	TChain* inTree_n = new TChain("tagged");
	for( int i = 3 ; i < argc; i++ ){
		cout << "Adding file " << argv[i] << endl;
		inTree_n->Add(argv[i]);
	}

	clashit * input_eHit 		= new clashit;
	TClonesArray* input_nHit 	= new TClonesArray("bandhit");
	TClonesArray* input_tag 	= new TClonesArray("taghit");
	int input_nleadindex		= -1;
	bool input_ngoodneutron		= false;
	double input_Ebeam		= 0;

	inTree_e->SetBranchAddress("eHit",		&input_eHit);
	inTree_e->SetBranchAddress("Ebeam",		&input_Ebeam);
	inTree_n->SetBranchAddress("nHits",		&input_nHit);
	inTree_n->SetBranchAddress("tag",		&input_tag);
	inTree_n->SetBranchAddress("nleadindex",	&input_nleadindex);
	inTree_n->SetBranchAddress("goodneutron",	&input_ngoodneutron);



	// Read entire tree into memory due to issue with jumping around in tree as we read it..
	vector<clashit> electron_list;
	vector<bandhit> neutron_list;
	vector<double>  Ebeam_list;
	int neutrons_saved = 0;
	cout << "Saving neutrons to mem...\n";
	for( int neutron = 0 ; neutron < inTree_n->GetEntries() ; neutron++ ){
		input_nHit->Clear();
		input_tag->Clear();
		input_nleadindex = -1;
		input_ngoodneutron = false;

		inTree_n->GetEntry(neutron);
	
		// Require a good neutron event:
		if(	input_ngoodneutron !=	NCUT_goodneutron	) continue;
		if(	input_nleadindex ==	NCUT_leadindex		) continue;
		
		// Grab our neutron:
		bandhit * lead_n = (bandhit*) input_nHit->At(input_nleadindex);
		bandhit copy_n;
		copy_n = *lead_n;

		// Check that the neutron we have is in our background region, in the CosThetaNQ bin, and in the Q2 bin
		if( lead_n->getTof() < NCUT_BACK_Tof_min || lead_n->getTof() > NCUT_BACK_Tof_max ) continue;

		// And then double check that it's above an Edep cut and has good status:
		if(	lead_n->getEdep() < NCUT_Edep 			) continue;


		// If it passes all these cuts then just push it back 
		//  -- 	the number of events we have in this list is the number of
		// 	background events in our background region
		neutron_list.push_back( copy_n );
	}

	cout << "Saving electrons to mem...\n";
	for( int electron = 0 ; electron < inTree_e->GetEntries() ; electron++ ){
		input_eHit->Clear();
		input_Ebeam = 0;

		inTree_e->GetEntry(electron);

		// Require a good electron event:
		if(	input_eHit->getPID() 						!= ECUT_PID 	) 	continue;
		if(	input_eHit->getCharge()						!= ECUT_charge 	)	continue;
		if(	input_eHit->getEoP() < ECUT_EoP_min || input_eHit->getEoP() > ECUT_EoP_max	) 	continue;
		if(	input_eHit->getEpcal() 						< ECUT_Epcal_min)	continue;
		if(	input_eHit->getV() < ECUT_V_min 	|| input_eHit->getW() < ECUT_W_min	) 	continue;
		if(	input_eHit->getVtz() < ECUT_vtx_min	|| input_eHit->getVtz() > ECUT_vtx_max	) 	continue;
		if(	input_eHit->getMomentum() < ECUT_pE_min || input_eHit->getMomentum() > ECUT_pE_max ) 	continue;
		if(	input_eHit->getQ2() < ECUT_Q2_min	|| input_eHit->getQ2() > ECUT_Q2_max	)	continue;
		if(	input_eHit->getW2() < ECUT_W2_min						)	continue;

		// now just get it and save it
		clashit copy_e;
		copy_e = *input_eHit;
		electron_list.push_back( copy_e );
		Ebeam_list.push_back( input_Ebeam );
	}

	// Loop over all neutron events
	for( int neutron = 0 ; neutron < neutron_list.size() ; neutron++ ){
		// Loop over 100 random electrons for each neutron:
		if( neutron % 1000 == 0 ) cout << "working on neutron " << neutron << "\n";
		for( int i_electron = 0 ; i_electron < 100; ++i_electron ){
			
			
			int nSignalBunches = (signal_max - signal_min)/BEAM_BUNCH;
			// Grab a random electron 
			int electron = myRand->Rndm() * electron_list.size();
			// Read in the stored variables
			double this_Ebeam 	= Ebeam_list	[electron];
			double p_e		= electron_list	[electron].getMomentum();
			double theta_e		= electron_list	[electron].getTheta();
			double phi_e		= electron_list	[electron].getPhi();
			double vtz 		= electron_list	[electron].getVtz();
			double theta_n		= neutron_list	[neutron].getDL().Theta();
			double phi_n		= neutron_list	[neutron].getDL().Phi();
			double dL		= neutron_list	[neutron].getDL().Mag();
			double bkg_Tof		= neutron_list	[neutron].getTof();

			
			// Figure out what bunch this background is in
			int this_bBunch		= (bkg_Tof - bkgrd_min )/BEAM_BUNCH;

			// Create all electron quantites only once
			TVector3 beamVec(0,0,this_Ebeam);
			TVector3 eVec;	eVec.SetMagThetaPhi(p_e,theta_e,phi_e);
			TVector3 qVec;	qVec = beamVec - eVec;
			double q 	= electron_list[electron].getQ();
			double theta_q  = electron_list[electron].getThetaQ();
			double phi_q 	= electron_list[electron].getPhiQ();
			double nu 	= electron_list[electron].getOmega();
			double Q2 	= electron_list[electron].getQ2();
			double xB	= electron_list[electron].getXb();
			double W2	= electron_list[electron].getW2();
				// Calculate the nq angles
			TVector3 norm_scatter = qVec.Cross( beamVec );
			norm_scatter 	= norm_scatter.Unit();

			// Loop over all signal bunches to push our background into
			for( int thisBunch = 0 ; thisBunch < (signal_max-signal_min)/BEAM_BUNCH ; thisBunch++ ){	

				// Shift ToF by current background bunch time and desired signal bunch time
				double bunch_shift	= (signal_min + BEAM_BUNCH*thisBunch) - (this_bBunch*BEAM_BUNCH + bkgrd_min);

				double Tof 	= bkg_Tof + bunch_shift;
				double TofpM	= Tof / (dL/100.);
				double beta 	= (1./TofpM) * (1./cAir) * (100./1);
				double p_n 	= mN / sqrt( 1./pow(beta,2) - 1. );
				
				TVector3 nVec;	nVec.SetMagThetaPhi(p_n,theta_n,phi_n);
				double E_n 	= sqrt( mN*mN + p_n*p_n );

				TVector3 norm_reaction = qVec.Cross( nVec );
				norm_reaction 	= norm_reaction.Unit();
				double phi_nq 	= norm_scatter.Angle( norm_reaction );
				double theta_nq = nVec.Angle( qVec );
				double CosTheta_nq = cos(theta_nq);
				TVector3 direction = norm_scatter.Cross(norm_reaction);
				if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
				}
				else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
					phi_nq *= (-1);
				}

				double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_n) - 2.*nu*E_n + 2.*q*p_n*cos(theta_nq);
				double Wp = sqrt(W_primeSq);
				double As = (E_n - p_n*CosTheta_nq)/mN;

				// Different definitions of x'
				// Bonus definition is the default
				double Xp = eHit.getQ2()/(2.*( eHit.getOmega()*(mD-E_n) + p_n*eHit.getQ()*CosTheta_nq));
				// W' definition
				double Xp_WP  = eHit.getQ2()/(W_primeSq - mN*mN + eHit.getQ2());
				// Bjorken definition
				double Xp_Bj  = eHit.getXb()/(2. - As);
				// PRC definition
				double Ei = mD - E_n;
				double ps_plus = mD/2. * As;
				double virt = (Ei*Ei - p_n*p_n - mN*mN)/(mN*mN);
				double p_plus = mD - ps_plus;
				double q_plus = eHit.getOmega() - eHit.getQ();
				double tP = virt * mN * mN;
				double Xp_PRC = (eHit.getQ2() - (q_plus/p_plus)*tP)/(W_primeSq - mN*mN + eHit.getQ2() - (q_plus/p_plus)*tP);

				TVector3 Pt;
				TVector3 pN_par_q = nVec.Dot(qVec) / (qVec.Mag2()) * qVec;
				Pt = nVec - pN_par_q;


				// Clear ehit branches
				eHit.Clear();
				Ebeam 		= 0;
				gated_charge 	= 0;
				livetime 	= 0;
				starttime 	= 0;
				current 	= 0;
				weight 		= 1;
				// Clear nhit branches
				nHits->Clear();
				nMult 		= 1;
				goodneutron 	= NCUT_goodneutron;
				nleadindex 	= 0;
				// Clear tag branches
				tags->Clear();

				// Save electron info
				Ebeam = this_Ebeam;
				eHit = electron_list[electron];

				// Save neutron info with new ToF
				neutron_list[neutron].setTof(Tof);
				neutron_list[neutron].setTofFadc(Tof);
				new(saveHit[0]) bandhit;
				saveHit[0] = &neutron_list[neutron];

				// Save tag info with new vectors
				taghit new_tag;
				new_tag.setMomentumE	(eVec 		);
				new_tag.setMomentumN	(nVec		);
				new_tag.setMomentumQ	(qVec		);
				new_tag.setMomentumB	(beamVec	);

				new_tag.setPhiNQ	(phi_nq		);
				new_tag.setThetaNQ	(theta_nq	);
				new_tag.setWp		(Wp		);
				new_tag.setAs		(As		);
				new_tag.setPt		(Pt		);
				new_tag.setXp		(Xp		);
				new_tag.setXp_WP	(Xp_WP		);
				new_tag.setXp_Bj	(Xp_Bj		);
				new_tag.setXp_PRC	(Xp_PRC		);
				new(saveTags[0]) taghit;
				saveTags[0] = &new_tag;


				outTree->Fill();

			} // end loop over bunches
		} // end loop over 100 electrons
	} // end loop over neutrons
	
	// Write output file:
	outFile->cd();
	outTree->Write();
	outFile->Close();


	return 0;
}
