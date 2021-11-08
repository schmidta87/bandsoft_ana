#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TEventList.h"
#include "TCut.h"
#include "TChain.h"
#include "TVector3.h"

#include "constants.h"
#include "kinematic_cuts.h"
#include "bandhit.h"
#include "clashit.h"
#include "taghit.h"
#include "genpart.h"

using namespace std;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
void getMcTag( TClonesArray * const gen , const double Ebeam, taghit * gen_taghit );
bool pointsToBand(double theta,double phi,double z_cm);
bool checkTagHit(taghit * this_taghit);

int main(int argc, char ** argv){
	if (argc < 2){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [sim gen file] [sim reco file] \n\n";
		return -1;
	}

	TFile * inFileGen 	= new TFile(argv[1]);
	TTree * inTreeGen 	= (TTree*) inFileGen->Get("tagged");
	double Ebeam = 0.;
	clashit * eHit = NULL;
	TClonesArray * gen 	= new TClonesArray("genpart");
	TClonesArray * nHits	= new TClonesArray("bandhit");
	int nMult = 0;
	bool goodneutron = false;
	int nleadindex = 0;
	TClonesArray * tags	= new TClonesArray("taghit");
	inTreeGen->SetBranchAddress("eHit_smeared",	&eHit	);
	inTreeGen->SetBranchAddress("Ebeam",	&Ebeam	);
	inTreeGen->SetBranchAddress("mcParts",	&gen	);
	inTreeGen->SetBranchAddress("nMult",	&nMult	);
	inTreeGen->SetBranchAddress("nHits",	&nHits	);
	inTreeGen->SetBranchAddress("goodneutron",	&goodneutron);
	inTreeGen->SetBranchAddress("nleadindex",	&nleadindex);
	inTreeGen->SetBranchAddress("tag",		&tags);

	//TFile * inFileRec = new TFile(argv[2]);
	//TTree * inTreeRec = (TTree*) inFileRec->Get("tagged");
	//TClonesArray * rec	= new TClonesArray("tag");
	//inTreeRec->SetBranchAddress("tag",	rec	);

	int GoodElectron = 0;
	int GoodGeneratedTag = 0;
	int GoodNeutronReco = 0;
	int GoodTaggedReco = 0;

	for( int gen_event = 0 ; gen_event < inTreeGen->GetEntries() ; ++gen_event ){
		eHit->Clear();
		Ebeam = 0.;
		gen->Clear();
		nHits->Clear();
		nMult = 0;
		goodneutron = false;
		nleadindex = 0;
		tags->Clear();
		
		inTreeGen->GetEntry(gen_event);
		//cout << "Event " << gen_event << "\n";


		//cout << eHit->getPID() << "\n";
		//cout << eHit->getCharge() << "\n";
		//cout << eHit->getEoP() << "\n";
		//cout << eHit->getEpcal() << "\n";
		//cout << eHit->getV() << "\n";
		//cout << eHit->getVtz() << "\n";
		//cout << eHit->getMomentum() << "\n";
		//cout << eHit->getQ2() << "\n";
		//cout << eHit->getW2() << "\n";
		// Check if eHit exists and passes our final cut selection
		if( eHit->getPID() != ECUT_PID ) 							continue; // PID check
		if( eHit->getCharge() != ECUT_charge)							continue; // charge check
		if( eHit->getEoP() < ECUT_EoP_min || eHit->getEoP() > ECUT_EoP_max ) 			continue; // E/p check
		if( eHit->getEpcal() < ECUT_Epcal_min )							continue; // Epcal check
		if( eHit->getV() < ECUT_V_min )								continue; // Pcal V check
		if( eHit->getVtz() < ECUT_vtx_min || eHit->getVtz() > ECUT_vtx_max ) 			continue; // vtx check
		if( eHit->getMomentum() < ECUT_pE_min || eHit->getMomentum() > ECUT_pE_max )		continue; // pE check
		if( eHit->getQ2() < ECUT_Q2_min || eHit->getQ2() > ECUT_Q2_max )			continue; // Q2 check
		if( eHit->getW2() < ECUT_W2_min )							continue; // W2 check
		//cout << "\te reconstructed\n";
		++GoodElectron;

		// Get the GENERATED tag hit info:
		taghit gen_taghit;
		getMcTag( gen, Ebeam, &gen_taghit );
	
		//cout << pointsToBand( gen_taghit.getMomentumN().Theta() , gen_taghit.getMomentumN().Phi() , 0 )  << "\n";
		//cout << gen_taghit.getMomentumN().Theta() << "\n";
		//cout << gen_taghit.getThetaNQ() << "\n";
		//cout << gen_taghit.getMomentumN().Mag() << "\n";
		//cout << gen_taghit.getWp()  << "\n";
		//cout << gen_taghit.getAs()  << "\n";
		//cout << nMult << "\n";
		// Ask if it has the correct tag kinematics for our final cut selection:
		if( checkTagHit( &gen_taghit ) == false ) continue;
		//cout << "\tgood tag generated\n";
		++GoodGeneratedTag;
		
		// Now check if our neutron was reconstructed:
		if( nMult == 0 ) 				continue; // nMult check
		if( goodneutron != NCUT_goodneutron ) 		continue; // goodneutron cut
		if( nleadindex == NCUT_leadindex ) 		continue; // nleadindex cut
		bandhit * lead_n = (bandhit*) nHits->At(nleadindex);
		if( lead_n->getEdep() < 10 )			continue; // Edep cut
		++GoodNeutronReco;


		// Lastly, see if the event had a good reconstructed tag hit:
		taghit * rec_taghit = (taghit*) tags->At(nleadindex);
		if( checkTagHit( rec_taghit ) == false ) continue;
		++GoodTaggedReco;

		
		//cout << "\tnHits not empty\n";
	
	}
	cout << GoodElectron << " " << GoodGeneratedTag << " " << GoodNeutronReco << " " << GoodTaggedReco << "\n";


	inFileGen->Close();
	//inFileRec->Close();
	return 0;
}

void getMcTag( TClonesArray * const gen , const double Ebeam, taghit * gen_taghit ){
	TVector3 nVec;
	TVector3 eVec;
	TVector3 beamVec(0,0,Ebeam);
	for( int genidx = 0 ; genidx < gen->GetEntriesFast() ; ++genidx ){
		genpart * this_gen_hit = (genpart*) gen->At(genidx);
		int PID = this_gen_hit->getPID();

		if( PID == 2112 ){
			nVec.SetMagThetaPhi( this_gen_hit->getMomentum() , this_gen_hit->getTheta() , this_gen_hit->getPhi() );
		}
		else if( PID == 11 ){
			eVec.SetMagThetaPhi( this_gen_hit->getMomentum() , this_gen_hit->getTheta() , this_gen_hit->getPhi() );
		}
		else{
			continue;
		}
	}
	TVector3 qVec;	qVec = beamVec - eVec;
	double q 	= qVec.Mag();
	double theta_q  = qVec.Theta();
	double phi_q 	= qVec.Phi();
	double nu 	= Ebeam - sqrt( pow(eVec.Mag(),2) + mE*mE );
	double Q2 	= qVec.Mag()*qVec.Mag() - pow(nu,2);
	double xB	= Q2/(2.*mP*nu);
	double W2	= mP*mP - Q2 + 2*nu*mP;
	TVector3 norm_scatter = qVec.Cross( beamVec );
	norm_scatter 	= norm_scatter.Unit();

	// With neutron and electron vectors set, let's calculate everything else:
	double p_n = nVec.Mag();
	double E_n = sqrt( mN*mN + p_n*p_n );

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
	double Xp = Q2/(2.*( nu*(mD-E_n) + p_n*q*CosTheta_nq));
	// W' definition
	double Xp_WP  = Q2/(W_primeSq - mN*mN + Q2);
	// Bjorken definition
	double Xp_Bj  = xB/(2. - As);
	// PRC definition
	double Ei = mD - E_n;
	double ps_plus = mD/2. * As;
	double virt = (Ei*Ei - p_n*p_n - mN*mN)/(mN*mN);
	double p_plus = mD - ps_plus;
	double q_plus = nu - q;
	double tP = virt * mN * mN;
	double Xp_PRC = (Q2 - (q_plus/p_plus)*tP)/(W_primeSq - mN*mN + Q2 - (q_plus/p_plus)*tP);

	TVector3 Pt;
	TVector3 pN_par_q = nVec.Dot(qVec) / (qVec.Mag2()) * qVec;
	Pt = nVec - pN_par_q;

	gen_taghit->setMomentumE	(eVec 		);
	gen_taghit->setMomentumN	(nVec		);
	gen_taghit->setMomentumQ	(qVec		);
	gen_taghit->setMomentumB	(beamVec	);

	gen_taghit->setPhiNQ	(phi_nq		);
	gen_taghit->setThetaNQ	(theta_nq	);
	gen_taghit->setWp	(Wp		);
	gen_taghit->setAs	(As		);
	gen_taghit->setPt	(Pt		);
	gen_taghit->setXp	(Xp		);
	gen_taghit->setXp_WP	(Xp_WP		);
	gen_taghit->setXp_Bj	(Xp_Bj		);
	gen_taghit->setXp_PRC	(Xp_PRC		);
}

bool pointsToBand(double theta,double phi,double z_cm){

  double inset = 10;//inset distance from edges of BAND in [cm]
  double z = z_cm;

  // Numbers taken from band/src/main/java/org/jlab/rec/band/constants/Parameters.java  
  double thickness  = 7.3;                                // thickness of each bar (cm)
  double layerGap[] = {7.94, 7.62, 7.94, 7.62, 7.3};      // gap between center of neig\
							  hbouring layers (cm), 1-2, 2-3, 3-4, 4-5, 5-6                                                 

  // Numbers taken from clas-band-calib/bin/src/org/clas/fcmon/band/BANDConstants.java  
  double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};

  // Distance from ideal target to upstream end of BAND                                 
  // (from BAND survey report, 02/18/2019)                                              
  double zUpst = (-302.69-302.69-302.57-302.64)/4.; // [cm]                             

  // Distance from ideal target to downstream end of layer 5                            
  double zDown = (zUpst + 5*thickness) - z_cm;

  double rho   = zDown/cos(theta);
  double xDown = rho*sin(theta)*cos(phi);
  double yDown = rho*sin(theta)*sin(phi);

  double globalX = (-240.5-240.5+241.0+243.7)/4./10.; // [cm] --> Not using this yet (n\
						      eed to make sure we have the right coordinate system)                                         
  double globalY = (-211.0+228.1-210.6+228.1)/4./10.; // [cm]

  // Sector boundaries                                                                  
  double topSec1  = globalY + 12*thickness; // ****** THIS WILL CUT OUT THE VERY TOP OF THE DETECTOR, I.E. BARS 1 IN SECTOR 1 FOR ALL LAYERS
  double topSec2  = globalY + 10*thickness;
  double topSec34 = globalY +  3*thickness;
  double topSec5  = globalY -  3*thickness;
  double downSec5 = globalY -  5*thickness;

  if( yDown >= topSec1 || yDown <= downSec5 ) return 0;

  if(
     (yDown < topSec1  && yDown >= topSec2  && fabs(xDown) < bandlen[0]/2. - inset) ||
     ( (yDown < topSec2  && yDown >= topSec34 && fabs(xDown) < bandlen[1]/2. - inset) &&   !(yDown >= topSec2 - inset && yDown < topSec2  && fabs(xDown) > (bandlen[0]/2.  - inset) ) && !(yDown >= topSec34 && yDown < topSec34 + inset && abs(xDown) < bandlen[1]/2 - bandlen[2] + inset) ) ||
     (yDown < topSec34 && yDown >= topSec5  && fabs(xDown) < bandlen[1]/2.- inset && fabs(xDown) > (bandlen[1]/2.-bandlen[2] + inset)) ||
     ( (yDown < topSec5  && yDown >= downSec5 && fabs(xDown) < bandlen[4]/2. - inset) && !(yDown >= topSec5 -inset && yDown < topSec5 + inset && abs(xDown) < bandlen[1]/2 - bandlen[2] + inset) )

      ){
    return 1;
    }
  return 0;
}


bool checkTagHit(taghit * this_taghit){
	if( pointsToBand( this_taghit->getMomentumN().Theta() , this_taghit->getMomentumN().Phi() , 0 ) == 0 )		return false; // 10cm fiducial cut AND cutting out the top bars
	if( this_taghit->getMomentumN().Theta() 	> 168.5*M_PI/180. 	) 					return false;
	if( this_taghit->getThetaNQ() < NCUT_THETANQ_min 	|| this_taghit->getThetaNQ() > NCUT_THETANQ_max	)	return false; // thetaNQ cut
	if( this_taghit->getMomentumN().Mag() < NCUT_Pn_min || this_taghit->getMomentumN().Mag() > NCUT_Pn_max )	return false; // pN cut
	if( this_taghit->getMomentumN().Mag() != this_taghit->getMomentumN().Mag()	)				return false; // pN not a number cut
	if( this_taghit->getWp() < NCUT_Wp_min || this_taghit->getWp() > NCUT_Wp_max )					return false; // Wp cut
	if( this_taghit->getAs() < NCUT_As_min || this_taghit->getAs() > NCUT_As_max )					return false; // As cut
	if( this_taghit->getMomentumN().Mag() < 0.25 )									return false; // Additional pN min cut

	return true;

}
