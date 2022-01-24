#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"

#include "genpart.h"
#include "clashit.h"
#include "bandhit.h"
#include "particles.h"

#include "kinematic_cuts.h"

using namespace std;
bool pointsToBand(double theta,double phi,double z_m);



int main(int argc, char** argv) {

	if( argc < 4 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [MC (0) or Data (1)] [outputFile] [inputFiles]... \n\n";
		return -1;
	}

	int MC_DATA_OPT = atoi(argv[1]);

	TChain inChain("qe");
	for(int k = 3; k < argc; k++){
		inChain.Add(argv[k]);
	}

	// Set the input branches
	clashit *	eHit 		= new clashit;
	TClonesArray* 	nHits 		= new TClonesArray("bandhit");
	TClonesArray*	charged_parts	= new TClonesArray("particles");
	int nleadindex			= 0;
	int goodneutron			= 0;
	double Ebeam 			= 0.;
	int charged_mult		= 0;
	double weight			= 0.;
	if( MC_DATA_OPT ){ // 1 == DATA
		inChain.SetBranchAddress("eHit"		,&eHit	);
	}
	else{ // 0 == SIM
		inChain.SetBranchAddress("eHit_smeared"	,&eHit	);
	}
	inChain.SetBranchAddress("nHits"		,&nHits		);
	inChain.SetBranchAddress("nleadindex"		,&nleadindex	);
	inChain.SetBranchAddress("goodneutron"		,&goodneutron	);
	inChain.SetBranchAddress("Ebeam"		,&Ebeam		);
	inChain.SetBranchAddress("charged_mult"		,&charged_mult	);
	inChain.SetBranchAddress("charged_particles"	,&charged_parts	);
	inChain.SetBranchAddress("weight"		,&weight 	);

	// Create output tree
	TFile * outFile = new TFile(argv[2],"RECREATE");
	TTree * outTree = new TTree("qe","CLAS & BAND QE physics");
	TVector3 beamVec, eVec, qVec, pVec, missVec, nVec;
	double Emiss, Pmiss, Mmiss, Thetamiss, ThetaPQ, Edep_n, Tof_n, out_weight;
	bool ptb;
	bool valid_n;
	outTree->Branch("beamVec"		,&beamVec	);
	outTree->Branch("eVec"			,&eVec		);
	outTree->Branch("qVec"			,&qVec		);
	outTree->Branch("pVec"			,&pVec		);
	outTree->Branch("ThetaPQ"		,&ThetaPQ	);
	outTree->Branch("missVec"		,&missVec	);
	outTree->Branch("Emiss"			,&Emiss		);
	outTree->Branch("Pmiss"			,&Pmiss		);
	outTree->Branch("Mmiss"			,&Mmiss		);
	outTree->Branch("Thetamiss"		,&Thetamiss	);
	outTree->Branch("ptb"			,&ptb		);
	outTree->Branch("nVec"			,&nVec		);
	outTree->Branch("valid_n"		,&valid_n	);
	outTree->Branch("Edep_n"		,&Edep_n	);
	outTree->Branch("Tof_n"			,&Tof_n		);
	outTree->Branch("weight"		,&out_weight	);

	std::vector<TH1D*> hists;
	TH1D * h1_e_v		= new TH1D("h1_e_v"	,"h1_e_v"	,100,0,250);
	TH1D * h1_e_w		= new TH1D("h1_e_w"	,"h1_e_w"	,100,0,250);
	TH1D * h1_e_eop		= new TH1D("h1_e_eop"	,"h1_e_eop"	,100,0,0.5);
	TH1D * h1_e_epcal	= new TH1D("h1_e_epcal"	,"h1_e_epcal"	,100,0,0.5);
	TH1D * h1_e_p		= new TH1D("h1_e_p"	,"h1_e_p"	,100,0,2);
	TH1D * h1_e_vtz		= new TH1D("h1_e_vtz"	,"h1_e_vtz"	,100,-10,0);
	TH1D * h1_p_chi2	= new TH1D("h1_p_chi2"	,"h1_p_chi2"	,100,-10,10);
	TH1D * h1_p_theta	= new TH1D("h1_p_theta"	,"h1_p_theta"	,100,0,180);
	TH1D * h1_p_vtz		= new TH1D("h1_p_vtz"	,"h1_p_vtz"	,100,-10,10);
	TH1D * h1_p_pmult	= new TH1D("h1_p_pmult"	,"h1_p_pmult"	,100,0,10);
	TH1D * h1_pmiss		= new TH1D("h1_pmiss"	,"h1_pmiss"	,100,0,1);
	TH1D * h1_thetamiss	= new TH1D("h1_thetamiss","h1_thetamiss",100,0,180);
	hists.push_back(	h1_e_v		);
        hists.push_back(	h1_e_w		);
        hists.push_back(	h1_e_eop	);
        hists.push_back(	h1_e_epcal	);
        hists.push_back(	h1_e_p		);
        hists.push_back(	h1_e_vtz	);
        hists.push_back(	h1_p_chi2	);
        hists.push_back(	h1_p_theta	);
        hists.push_back(	h1_p_vtz	);
        hists.push_back(	h1_p_pmult	);
        hists.push_back(	h1_pmiss	);
        hists.push_back(	h1_thetamiss	);




	int eep_counter = 0;
	int eepn_counter = 0;
	int fin = inChain.GetEntries();
	for(int event = 0; event < fin; event++){
		eHit		->Clear();
		nHits		->Clear();
		charged_parts	->Clear();
		nleadindex 	= -1;
		goodneutron	= 0.;
		Ebeam 		= 0.;
		charged_mult 	= 0;
		weight 		= 0.;
		beamVec.Clear();
		eVec.Clear();
		qVec.Clear();
		pVec.Clear();
		ThetaPQ = 0;
		missVec.Clear();
		Emiss = 0;
		Pmiss = 0;
		Mmiss = 0;
		Thetamiss = 0;
		ptb = false;
		nVec.Clear();
		nVec.SetMagThetaPhi(0,0,0);
		valid_n = false;
		Edep_n = 0;
		Tof_n =0;
		out_weight = 0;
		inChain.GetEntry(event);

		//Display completed  
		if((event%100000) == 0){
			cerr << (event*100.)/fin <<"% complete \n";
		}

		out_weight = weight;

		////////////////////////////////////////////////
		// electron PID cuts
		////////////////////////////////////////////////
		h1_e_v		->Fill(eHit->getV()		,weight	);
		h1_e_w		->Fill(eHit->getW()		,weight	);
		h1_e_eop	->Fill(eHit->getEoP()		,weight	);
		h1_e_epcal	->Fill(eHit->getEpcal()		,weight	);
		h1_e_p		->Fill(eHit->getMomentum()	,weight	);
		h1_e_vtz	->Fill(eHit->getVtz()		,weight	);

		if(eHit->getV() < ECUT_V_min)		continue; 
		if(eHit->getW() < ECUT_W_min)		continue; 
		if(eHit->getEoP() < ECUT_EoP_min)	continue; 
		if(eHit->getEoP() > ECUT_EoP_max)	continue; 
		if(eHit->getEpcal() < ECUT_Epcal_min)	continue;
		if(eHit->getMomentum() < 1)		continue; 
		if(eHit->getMomentum() > Ebeam)		continue; 
		if(eHit->getVtz() < ECUT_vtx_min)	continue; 
		if(eHit->getVtz() > ECUT_vtx_max)	continue; 

		beamVec.SetXYZ( 0,0,Ebeam);
		eVec.SetMagThetaPhi( eHit->getMomentum() , eHit->getTheta() , eHit->getPhi() );
		qVec.SetMagThetaPhi( eHit->getQ() , eHit->getThetaQ() , eHit->getPhiQ() );

		////////////////////////////////////////////////
		// proton PID cuts
		////////////////////////////////////////////////

		int numProtons = 0;
		int indexP = -1;
		for(int thisPart = 0; thisPart < charged_mult; thisPart++){
			particles * this_particle = (particles*) charged_parts->At(thisPart);

			int PID = this_particle->getPID();
			int charge = this_particle->getCharge();
			int status = this_particle->getStatus();
			// Ignore any protons from the central detector
			if( PID == 2212 && charge == 1 && status > 1000 && status < 4000){
				
				// Do a chi2 cut on PID:
				double chi2pid = this_particle->getChi2();
				h1_p_chi2	->Fill(chi2pid					,weight );
				h1_p_theta	->Fill(this_particle->getTheta()* 180./M_PI 	,weight );
				h1_p_vtz	->Fill(this_particle->getVtz()			,weight );
				if( chi2pid < -4 || chi2pid > 8 ) 			continue;	// TODO we can refine this
				if( this_particle->getTheta() * 180./M_PI > 35	)	continue;
	
				// Do a ThetaPQ cut on the lead
				pVec.SetMagThetaPhi( 	this_particle->getMomentum(),
							this_particle->getTheta(),
							this_particle->getPhi()			);
				ThetaPQ = pVec.Angle( qVec );
				//if( ThetaPQ * 180./M_PI < 25 ) 	continue;	// TODO can come back to this later

				indexP = thisPart;
				++numProtons;
			}

		}
		h1_p_pmult	->Fill(numProtons	,weight );
		if( numProtons != 1 ) continue;
		particles * proton = (particles*) charged_parts->At(indexP);
		pVec.SetMagThetaPhi( 	proton->getMomentum(),proton->getTheta(),proton->getPhi()  );
		ThetaPQ = pVec.Angle( qVec );
		double Ep = sqrt( mP*mP + pVec.Mag2() );
		double Ee = sqrt( mE*mE + eVec.Mag2() );

		missVec = beamVec - eVec - pVec;
		Emiss = Ebeam + mD - Ee - Ep;
		Pmiss = missVec.Mag();
		Mmiss = sqrt((Emiss * Emiss) - missVec.Mag2() );
		Thetamiss = missVec.Theta();
		ptb = pointsToBand( missVec.Theta() , missVec.Phi() , -3 );


		////////////////////////////////////////////////
		// (e,ep) event selection cuts
		////////////////////////////////////////////////
		h1_pmiss	->Fill(Pmiss			,weight );
		h1_thetamiss	->Fill(Thetamiss* 180./M_PI 	,weight );
		if( Pmiss < 0.25 )			continue;
		if( Thetamiss * 180./M_PI < 80 )	continue;

		if( ptb && (Mmiss > 0.9 && Mmiss < 1.1)) ++eep_counter;
		////////////////////////////////////////////////
		// neutron PID
		////////////////////////////////////////////////
		if( goodneutron == NCUT_goodneutron && nleadindex != NCUT_leadindex ){

			bandhit * this_neutron = (bandhit *) nHits->At(nleadindex);
			//if( this_neutron->getEdep () > NCUT_Edep ){
				//bool fiducial1 = !(this_neutron->getSector()==1 && (this_neutron->getX() >  70 || this_neutron->getX() < -70));
				//bool fiducial2 = !(this_neutron->getSector()==2 && (this_neutron->getX() >  90 || this_neutron->getX() < -90));
				//bool fiducial3 = !(this_neutron->getSector()==3 && (this_neutron->getX() >  90 || this_neutron->getX() <  60));
				//bool fiducial4 = !(this_neutron->getSector()==4 && (this_neutron->getX() < -90 || this_neutron->getX() > -60));
				//bool fiducial5 = !(this_neutron->getSector()==5 && (this_neutron->getX() >  90 || this_neutron->getX() < -90));
				//bool fiducial = fiducial1 || fiducial2 || fiducial3 || fiducial4 || fiducial5;
				//if( fiducial == false )					continue;
				bool thetaFiducial = this_neutron->getDL().Theta() > 168.5*M_PI/180.;
				//if( thetaFiducial == false ) 				continue;
				
				nVec.SetMagThetaPhi( 	this_neutron->getMomentumN().Mag() ,
							this_neutron->getDL().Theta(),
							this_neutron->getDL().Phi()		);
				if( ptb && (Mmiss > 0.9 && Mmiss < 1.1)) ++eepn_counter;
				valid_n = true;
				Edep_n = this_neutron->getEdep();
				Tof_n = this_neutron->getTof();
			//}

		}
		outTree->Fill();
		

		/*
		
		////////////////////////////////////////////////
		//(e,e'n) Cuts
		////////////////////////////////////////////////
		if( goodneutron != NCUT_goodneutron)			continue;
		if( nleadindex == NCUT_leadindex )			continue;
		bandhit * this_neutron = (bandhit *) nHits->At(nleadindex);
		if( this_neutron->getEdep () < NCUT_Edep )		continue;
		bool fiducial1 = !(this_neutron->getSector()==1 && (this_neutron->getX() >  70 || this_neutron->getX() < -70));
		bool fiducial2 = !(this_neutron->getSector()==2 && (this_neutron->getX() >  90 || this_neutron->getX() < -90));
		bool fiducial3 = !(this_neutron->getSector()==3 && (this_neutron->getX() >  90 || this_neutron->getX() <  60));
		bool fiducial4 = !(this_neutron->getSector()==4 && (this_neutron->getX() < -90 || this_neutron->getX() > -60));
		bool fiducial5 = !(this_neutron->getSector()==5 && (this_neutron->getX() >  90 || this_neutron->getX() < -90));
		bool fiducial = fiducial1 || fiducial2 || fiducial3 || fiducial4 || fiducial5;
		if( fiducial == false )					continue;
		bool thetaFiducial = this_neutron->getDL().Theta() > 168.5*M_PI/180.;
		if( thetaFiducial == false ) 				continue;
		*/
		
		// PID is good in FTOF
		// 	ThetaPQ cut on lead
		// PID not good in CD but should improve
		// Chi2PID cut not great in CD 
		//
		// look at how many protons pass the criteria - thetaPQ included
		// 	- ToF cut instead of Chi2PID ( not a big difference in FD over CD )
		// 	- 
	}
	cout << eep_counter << " " << eepn_counter << "\n";


	//inFile->Close();
	outFile->cd();
	for( int h = 0 ; h < hists.size() ; ++h ) hists[h]->Write();
	outTree->Write();
	outFile->Close();
	cout<<"Finished making file: "<<argv[1]<<"\n";

	return 0;
}

bool pointsToBand(double theta,double phi,double z_m){
	//double z = z_m*100; // from m to cm
	double z = z_m;

	// Numbers taken from band/src/main/java/org/jlab/rec/band/constants/Parameters.java
	double thickness  = 7.2;                                // thickness of each bar (cm)
	double layerGap[] = {7.94, 7.62, 7.94, 7.62, 7.3};      // gap between center of neighbouring layers (cm), 1-2, 2-3, 3-4, 4-5, 5-6

	// Numbers taken from clas-band-calib/bin/src/org/clas/fcmon/band/BANDConstants.java
	double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};

	// Distance from ideal target to upstream end of BAND
	// (from BAND survey report, 02/18/2019)
	double zUpst = (-302.69-302.69-302.57-302.64)/4.; // [cm]

	// Distance from ideal target to downstream end of layer 5
	double zDown = (zUpst + 5*thickness) - z_m;

	double rho   = zDown/cos(theta);
	double xDown = rho*sin(theta)*cos(phi);
	double yDown = rho*sin(theta)*sin(phi);

	double globalX = (-240.5-240.5+241.0+243.7)/4.; // [cm] --> Not using this yet (need to make sure we have the right coordinate system)
	double globalY = (-211.0+228.1-210.6+228.1)/4.; // [cm]

	// Sector boundaries
	double topSec1  = globalY + 13*thickness;
	double topSec2  = globalY + 10*thickness;
	double topSec34 = globalY +  3*thickness;
	double topSec5  = globalY -  3*thickness;
	double downSec5 = globalY -  5*thickness;

	if( yDown >= topSec1 || yDown <= downSec5 ) return 0;

	if(             (yDown < topSec1  && yDown >= topSec2  && fabs(xDown) < bandlen[0]/2. )||
			(yDown < topSec2  && yDown >= topSec34 && fabs(xDown) < bandlen[1]/2. )||
			(yDown < topSec34 && yDown >= topSec5  && fabs(xDown) < bandlen[1]/2. && fabs(xDown) > bandlen[1]/2.-bandlen[2])||
			(yDown < topSec5  && yDown >= downSec5 && fabs(xDown) < bandlen[4]/2. )
	  ){
		return 1;
	}
	return 0;
}

