#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"

#include "genpart.h"
#include "clashit.h"
#include "bandhit.h"
#include "taghit.h"

#include "bin_edges.h"

using std::cerr;
using std::cout;
void calculate_AsPt( double & As, double & Pt , double & Xp_WP , const double Ebeam, genpart * const electron, genpart * const neutron );
void calculate_AsPt_recN( double & As, double & Pt , double & Xp_WP , const double Ebeam, genpart * const electron, taghit * const neutron );
void calculate_AsPt_recBoth( double & As, double & Pt , double & Xp_WP , const double Ebeam, clashit * const electron, taghit * const neutron );
void calculate_AsPt_recE( double & As, double & Pt , double & Xp_WP , const double Ebeam, clashit * const electron, genpart * const neutron );
void fillHist( double Q2, double Pt, double Xp, double As, TH1D**** hist , double weight );
void setError( double * err , TH1D * rec , TH1D * gen );
bool bad_bar( bandhit * this_n );

int main( int argc, char** argv){

	if( argc != 2 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Sim File]\n";
		return -1;
	}
	
	TFile inFile(argv[1]);
	TTree * inTree = (TTree*) inFile.Get("tagged");

	TClonesArray* gen_particles 	= new TClonesArray("genpart");
	clashit * rec_electron 		= new clashit;
	TClonesArray* rec_tagged 	= new TClonesArray("taghit");
	TClonesArray* rec_neutrons 	= new TClonesArray("bandhit");
	int 	genMult 	= 0;
	int	nleadindex 	= 0;
	double	Ebeam		= 0;
	inTree->SetBranchAddress("mcParts"		, &gen_particles		);
	inTree->SetBranchAddress("genMult"		, &genMult			);
	inTree->SetBranchAddress("eHit"			, &rec_electron			);
	inTree->SetBranchAddress("tag"			, &rec_tagged			);
	inTree->SetBranchAddress("nHits"		, &rec_neutrons			);
	inTree->SetBranchAddress("nleadindex"		, &nleadindex			);
	inTree->SetBranchAddress("Ebeam"		, &Ebeam			);
	
	TFile * outFile = new TFile("migcorrection_tagged.root","RECREATE");
	TH1D **** h4_gen_as 	= new TH1D***[bins_Q2];
	TH1D **** h4_rec_as	= new TH1D***[bins_Q2];
	TH2D **** h4_twoD_cutonrec	= new TH2D***[bins_Q2];
	TH2D **** h4_twoD_cutongen	= new TH2D***[bins_Q2];

	for( int i = 0 ; i < bins_Q2 ; ++i ){ // bins in Q2
		h4_gen_as[i]	= new TH1D**[bins_Pt];
		h4_rec_as[i]	= new TH1D**[bins_Pt];
		h4_twoD_cutonrec[i]	= new TH2D**[bins_Pt];
		h4_twoD_cutongen[i]	= new TH2D**[bins_Pt];
		for( int j = 0 ; j < bins_Pt ; ++j ){ // bins in Pt
			h4_gen_as[i][j] 	= new TH1D*[bins_As];
			h4_rec_as[i][j] 	= new TH1D*[bins_As];
			h4_twoD_cutonrec[i][j]	= new TH2D*[bins_As];
			h4_twoD_cutongen[i][j]	= new TH2D*[bins_As];
			for( int k = 0 ; k < bins_As ; ++k ){ // bins in As
				h4_gen_as[i][j][k] = new TH1D(Form("gen_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("gen_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in As
				h4_rec_as[i][j][k] = new TH1D(Form("migcorr_xp_Q2_%i_Pt_%i_As_%i",i,j,k),Form("migcorr_xp_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max);  // histogram in As

				h4_twoD_cutonrec[i][j][k] 	= new TH2D(Form("mig_cutonrec_Q2_%i_Pt_%i_As_%i",i,j,k),Form("mig_cutonrec_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max,bins_Xp,Xp_min,Xp_max);
				h4_twoD_cutongen[i][j][k] 	= new TH2D(Form("mig_cutongen_Q2_%i_Pt_%i_As_%i",i,j,k),Form("mig_cutongen_Q2_%i_Pt_%i_As_%i",i,j,k),bins_Xp,Xp_min,Xp_max,bins_Xp,Xp_min,Xp_max);

			}
		}
	}
	

	for( int event = 0 ; event < inTree->GetEntries() ; ++event ){
		gen_particles	->Clear();
		rec_electron	->Clear();
		rec_tagged	->Clear();
		rec_neutrons	->Clear();
		genMult 	= 0;
		nleadindex	= 0;
		Ebeam 		= 0;
		//if(event > 1000) break;
		inTree->GetEntry(event);

		// Get the reconstructed variables
		taghit * rec_tag	= (taghit*)  rec_tagged->At(nleadindex);
		double rec_Pn		= rec_tag->getMomentumN().Mag();

		bandhit * lead_n	= (bandhit*)	rec_neutrons->At(nleadindex);
		if( rec_Pn < 0.25 ) 		continue;
		if( bad_bar(lead_n) ) 		continue;
		if( lead_n->getEdep() < 10 ) 	continue;


		double rec_As		= rec_tag->getAs();
		double rec_Pt 		= rec_tag->getPt().Mag();
		double rec_Xp_WP	= rec_tag->getXp_WP();
		double rec_Q2		= rec_electron->getQ2();

		// Form As, ThetaNQ, W, X with generated:
		genpart * gen_electron 	= (genpart*) gen_particles->At(0);
		genpart * gen_neutron 	= (genpart*) gen_particles->At(1);
		double	gen_Q2		= gen_electron->getQ2();

		double gen_As, gen_Pt, gen_Xp_WP;
		calculate_AsPt( gen_As, gen_Pt, gen_Xp_WP, Ebeam, gen_electron, gen_neutron );

		// Use generated electron but reconstructed neutron:
		double genErecN_As, genErecN_Pt, genErecN_Xp_WP;
		double genErecN_Q2 = gen_Q2;
		calculate_AsPt_recN( genErecN_As, genErecN_Pt, genErecN_Xp_WP, Ebeam, gen_electron, rec_tag );

		// Use reconstructed electron but generated neutron:
		double recEgenN_As, recEgenN_Pt, recEgenN_Xp_WP;
		double recEgenN_Q2 = rec_Q2;
		calculate_AsPt_recE( recEgenN_As, recEgenN_Pt, recEgenN_Xp_WP, Ebeam, rec_electron , gen_neutron );

		// Use reconstructed electron AND reconstructed neutron to do crosscheck:
		double recErecN_As, recErecN_Pt, recErecN_Xp_WP;
		double recErecN_Q2 = rec_Q2;
		calculate_AsPt_recBoth( recErecN_As, recErecN_Pt, recErecN_Xp_WP, Ebeam, rec_electron, rec_tag );
		//cout << recErecN_As - rec_As << " " << recErecN_Pt - rec_Pt << " " << recErecN_Xp_WP - rec_Xp << "\n";

		// sort the generated values:
		fillHist( gen_Q2, gen_Pt, gen_Xp_WP, gen_As, h4_gen_as , 1.);

		// sort the reconstructed values:
		fillHist( rec_Q2, rec_Pt, rec_Xp_WP, rec_As, h4_rec_as , 1. );

		int gen_bin_pt = (int) ((gen_Pt - Pt_min)/Pt_step);
		int gen_bin_as = (int) ((gen_As - As_min)/As_step);
		int rec_bin_pt = (int) ((rec_Pt - Pt_min)/Pt_step);
		int rec_bin_as = (int) ((rec_As - As_min)/As_step);

		int gen_bin_q2 = -1;
		int rec_bin_q2 = -1;
		for( int q2_bin = bins_Q2-1 ; q2_bin >= 0; --q2_bin ){
			if( gen_Q2 > Q2Bins[q2_bin] ) gen_bin_q2 = q2_bin;
			if( rec_Q2 > Q2Bins[q2_bin] ) rec_bin_q2 = q2_bin;
		}


		if( 	gen_bin_pt < bins_Pt && gen_bin_pt >= 0 &&
			gen_bin_as < bins_As && gen_bin_as >= 0 &&
			gen_bin_q2 < bins_Q2 && gen_bin_q2 >= 0 ){
				h4_twoD_cutongen[gen_bin_q2][gen_bin_pt][gen_bin_as]->Fill( gen_Xp_WP , rec_Xp_WP );
		}
		if( 	rec_bin_pt < bins_Pt && rec_bin_pt >= 0 &&
			rec_bin_as < bins_As && rec_bin_as >= 0 &&
			rec_bin_q2 < bins_Q2 && rec_bin_q2 >= 0 ){
				h4_twoD_cutonrec[rec_bin_q2][rec_bin_pt][rec_bin_as]->Fill( gen_Xp_WP , rec_Xp_WP );
		}

		
	}

	TCanvas ** c_Q2 = new TCanvas*[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){ // loop over Q2 bins
		c_Q2[i] = new TCanvas(Form("c_Q2_%i",i),"",1200,1200);
		
		c_Q2[i]->Divide( bins_Pt, bins_As );
		for( int j = 0 ; j < bins_Pt ; ++j){ // loop over Pt bins
			for( int k = 0 ; k < bins_As ; ++k ){ // loop over Xb bins


				if( h4_gen_as[i][j][k]->Integral() < 1 || h4_rec_as[i][j][k]->Integral() < 1 ) continue;
				double * errors = new double[bins_Xp];
				setError( errors , h4_rec_as[i][j][k] , h4_gen_as[i][j][k] );
				h4_rec_as[i][j][k]->Divide( h4_gen_as[i][j][k] );
				for( int bin = 1 ; bin <= h4_rec_as[i][j][k]->GetXaxis()->GetNbins(); ++bin ){
					h4_rec_as[i][j][k]->SetBinError( bin , errors[bin-1] );
				}

				int ref_bin = h4_rec_as[i][j][k]->FindBin( Xp_Ref[k] );
				double binref = h4_rec_as[i][j][k]->GetBinContent( ref_bin );

				h4_rec_as[i][j][k]->Scale( 1./binref );



				c_Q2[i]->cd( (k*2)+1 + j );

				gStyle->SetPaintTextFormat("4.1f m");
				h4_rec_as[i][j][k]->SetLineWidth(3);
				h4_rec_as[i][j][k]->SetMarkerColor(4);
				h4_rec_as[i][j][k]->SetMaximum(2);
				h4_rec_as[i][j][k]->SetMinimum(0);
				h4_rec_as[i][j][k]->Draw("*P");
				TLine * hline = new TLine(Xp_min,1,Xp_max,1);
				hline->SetLineWidth(2);
				hline->SetLineColor(1);
				hline->SetLineStyle(2);
				hline->Draw("SAME");

				//h4_gen_as[i][j][k]->Write();
				h4_rec_as[i][j][k]->Write();

				/*
				for( int biny = 1 ; biny <= h4_twoD_cutonrec[i][j][k]->GetYaxis()->GetNbins(); ++biny ){
					double rowsum = 0.0;
					for( int binx = 1 ; binx <= h4_twoD_cutonrec[i][j][k]->GetXaxis()->GetNbins(); ++binx ){
						rowsum += h4_twoD_cutonrec[i][j][k]->GetBinContent(binx,biny);
					}
					for( int binx = 1 ; binx <= h4_twoD_cutonrec[i][j][k]->GetXaxis()->GetNbins(); ++binx ){
						double thiscontent = h4_twoD_cutonrec[i][j][k]->GetBinContent(binx,biny);
						if( thiscontent <= 0 || rowsum <= 0 ) continue;
						h4_twoD_cutonrec[i][j][k]->SetBinContent(binx,biny,thiscontent/rowsum);
					}
				}
				*/
				/*
				for( int binx = 1 ; binx <= h4_twoD_cutonrec[i][j][k]->GetXaxis()->GetNbins(); ++binx ){
					double colsum = 0.0;
					for( int biny = 1 ; biny <= h4_twoD_cutonrec[i][j][k]->GetYaxis()->GetNbins(); ++biny ){
						colsum += h4_twoD_cutonrec[i][j][k]->GetBinContent(binx,biny);
					}
					for( int biny = 1 ; biny <= h4_twoD_cutonrec[i][j][k]->GetYaxis()->GetNbins(); ++biny ){
						double thiscontent = h4_twoD_cutonrec[i][j][k]->GetBinContent(binx,biny);
						if( thiscontent <= 0 || colsum <= 0 ) continue;
						h4_twoD_cutonrec[i][j][k]->SetBinContent(binx,biny,thiscontent/colsum);
					}
				}
				*/

				h4_twoD_cutongen[i][j][k]->GetXaxis()->SetTitle("x'_{gen}");
				h4_twoD_cutongen[i][j][k]->GetYaxis()->SetTitle("x'_{rec}");
				h4_twoD_cutongen[i][j][k]->SetMarkerSize(1.8);
				h4_twoD_cutonrec[i][j][k]->GetXaxis()->SetTitle("x'_{gen}");
				h4_twoD_cutonrec[i][j][k]->GetYaxis()->SetTitle("x'_{rec}");
				h4_twoD_cutonrec[i][j][k]->SetMarkerSize(1.8);
				h4_twoD_cutongen[i][j][k]->Write();
				h4_twoD_cutonrec[i][j][k]->Write();

				delete[] errors;

			}
		}
		c_Q2[i]->Print(Form("binmigration_tagged_Q2_%i.pdf",i));
	}
	outFile->Close();

	inFile.Close();
	return 1;
}


void calculate_AsPt( double & As, double & Pt , double & Xp_WP , const double Ebeam, genpart * const electron, genpart * const neutron ){

	TVector3 	beamVec(0,0,Ebeam);
	TVector3	eVec; 
	eVec.SetMagThetaPhi( electron->getMomentum(), electron->getTheta(), electron->getPhi() );
	TVector3	qVec; qVec = beamVec - eVec;
	TVector3	nVec; 
	nVec.SetMagThetaPhi( neutron->getMomentum(), neutron->getTheta(), neutron->getPhi() );
	nVec.Unit();


	TVector3 norm_scatter = qVec.Cross( beamVec );
	norm_scatter 	= norm_scatter.Unit();

	TVector3 norm_reaction = qVec.Cross( nVec );
	norm_reaction 	= norm_reaction.Unit();

	double phi_nq 		= norm_scatter.Angle( norm_reaction );
	double theta_nq 	= nVec.Angle( qVec );
	double gen_CosThetaNQ 	= cos(theta_nq);

	TVector3 direction = norm_scatter.Cross(norm_reaction);
	if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
	}
	else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
		phi_nq *= (-1);
	}

	nVec.SetMagThetaPhi( neutron->getMomentum(), neutron->getTheta(), neutron->getPhi() );
	double p_n = neutron->getMomentum();

	double E_n 	= sqrt( mN*mN + p_n*p_n );
	double W_primeSq = mD*mD - electron->getQ2() + mN*mN + 2.*mD*(electron->getOmega()-E_n) - 2.*electron->getOmega()*E_n + 2.*electron->getQ()*p_n*cos(theta_nq);
	double Wp = sqrt(W_primeSq);
	double Xp = electron->getQ2()/(2.*( electron->getOmega()*(mD-E_n) + p_n*electron->getQ()*gen_CosThetaNQ));
	Xp_WP = electron->getQ2()/(W_primeSq - mN*mN + electron->getQ2());

	TVector3 Pt_vec;
	TVector3 pN_par_q = nVec.Dot(qVec) / (qVec.Mag2()) * qVec;

	As = (E_n - p_n*gen_CosThetaNQ)/mN;
	Pt_vec = nVec - pN_par_q;
	Pt = Pt_vec.Mag();
}
void calculate_AsPt_recN( double & As, double & Pt , double & Xp_WP , const double Ebeam, genpart * const electron, taghit * const neutron ){
	TVector3 	beamVec(0,0,Ebeam);
	TVector3	eVec; 
	eVec.SetMagThetaPhi( electron->getMomentum(), electron->getTheta(), electron->getPhi() );
	TVector3	qVec; qVec = beamVec - eVec;
	TVector3	nVec; 
	nVec = (TVector3) neutron->getMomentumN();
	nVec.Unit();


	TVector3 norm_scatter = qVec.Cross( beamVec );
	norm_scatter 	= norm_scatter.Unit();

	TVector3 norm_reaction = qVec.Cross( nVec );
	norm_reaction 	= norm_reaction.Unit();

	double phi_nq 		= norm_scatter.Angle( norm_reaction );
	double theta_nq 	= nVec.Angle( qVec );
	double gen_CosThetaNQ 	= cos(theta_nq);

	TVector3 direction = norm_scatter.Cross(norm_reaction);
	if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
	}
	else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
		phi_nq *= (-1);
	}

	nVec = (TVector3) neutron->getMomentumN();
	double p_n = neutron->getMomentumN().Mag();

	double E_n 	= sqrt( mN*mN + p_n*p_n );
	double W_primeSq = mD*mD - electron->getQ2() + mN*mN + 2.*mD*(electron->getOmega()-E_n) - 2.*electron->getOmega()*E_n + 2.*electron->getQ()*p_n*cos(theta_nq);
	double Wp = sqrt(W_primeSq);
	double Xp = electron->getQ2()/(2.*( electron->getOmega()*(mD-E_n) + p_n*electron->getQ()*gen_CosThetaNQ));
	Xp_WP = electron->getQ2()/(W_primeSq - mN*mN + electron->getQ2());

	TVector3 Pt_vec;
	TVector3 pN_par_q = nVec.Dot(qVec) / (qVec.Mag2()) * qVec;

	As = (E_n - p_n*gen_CosThetaNQ)/mN;
	Pt_vec = nVec - pN_par_q;
	Pt = Pt_vec.Mag();

}

void calculate_AsPt_recE( double & As, double & Pt , double & Xp_WP , const double Ebeam, clashit * const electron, genpart * const neutron ){

	TVector3 	beamVec(0,0,Ebeam);
	TVector3	eVec; 
	eVec.SetMagThetaPhi( electron->getMomentum(), electron->getTheta(), electron->getPhi() );
	TVector3	qVec; qVec = beamVec - eVec;
	TVector3	nVec; 
	nVec.SetMagThetaPhi( neutron->getMomentum(), neutron->getTheta(), neutron->getPhi() );
	nVec.Unit();


	TVector3 norm_scatter = qVec.Cross( beamVec );
	norm_scatter 	= norm_scatter.Unit();

	TVector3 norm_reaction = qVec.Cross( nVec );
	norm_reaction 	= norm_reaction.Unit();

	double phi_nq 		= norm_scatter.Angle( norm_reaction );
	double theta_nq 	= nVec.Angle( qVec );
	double gen_CosThetaNQ 	= cos(theta_nq);

	TVector3 direction = norm_scatter.Cross(norm_reaction);
	if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
	}
	else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
		phi_nq *= (-1);
	}

	nVec.SetMagThetaPhi( neutron->getMomentum(), neutron->getTheta(), neutron->getPhi() );
	double p_n = neutron->getMomentum();

	double E_n 	= sqrt( mN*mN + p_n*p_n );
	double W_primeSq = mD*mD - electron->getQ2() + mN*mN + 2.*mD*(electron->getOmega()-E_n) - 2.*electron->getOmega()*E_n + 2.*electron->getQ()*p_n*cos(theta_nq);
	double Wp = sqrt(W_primeSq);
	double Xp = electron->getQ2()/(2.*( electron->getOmega()*(mD-E_n) + p_n*electron->getQ()*gen_CosThetaNQ));
	Xp_WP = electron->getQ2()/(W_primeSq - mN*mN + electron->getQ2());

	TVector3 Pt_vec;
	TVector3 pN_par_q = nVec.Dot(qVec) / (qVec.Mag2()) * qVec;

	As = (E_n - p_n*gen_CosThetaNQ)/mN;
	Pt_vec = nVec - pN_par_q;
	Pt = Pt_vec.Mag();
}
void calculate_AsPt_recBoth( double & As, double & Pt , double & Xp_WP , const double Ebeam, clashit * const electron, taghit * const neutron ){
	TVector3 	beamVec(0,0,Ebeam);
	TVector3	eVec; 
	eVec.SetMagThetaPhi( electron->getMomentum(), electron->getTheta(), electron->getPhi() );
	TVector3	qVec; qVec = beamVec - eVec;
	TVector3	nVec; 
	nVec = (TVector3) neutron->getMomentumN();
	nVec.Unit();


	TVector3 norm_scatter = qVec.Cross( beamVec );
	norm_scatter 	= norm_scatter.Unit();

	TVector3 norm_reaction = qVec.Cross( nVec );
	norm_reaction 	= norm_reaction.Unit();

	double phi_nq 		= norm_scatter.Angle( norm_reaction );
	double theta_nq 	= nVec.Angle( qVec );
	double gen_CosThetaNQ 	= cos(theta_nq);

	TVector3 direction = norm_scatter.Cross(norm_reaction);
	if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
	}
	else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
		phi_nq *= (-1);
	}

	nVec = (TVector3) neutron->getMomentumN();
	double p_n = neutron->getMomentumN().Mag();

	double E_n 	= sqrt( mN*mN + p_n*p_n );
	double W_primeSq = mD*mD - electron->getQ2() + mN*mN + 2.*mD*(electron->getOmega()-E_n) - 2.*electron->getOmega()*E_n + 2.*electron->getQ()*p_n*cos(theta_nq);
	double Wp = sqrt(W_primeSq);
	double Xp = electron->getQ2()/(2.*( electron->getOmega()*(mD-E_n) + p_n*electron->getQ()*gen_CosThetaNQ));
	Xp_WP = electron->getQ2()/(W_primeSq - mN*mN + electron->getQ2());

	TVector3 Pt_vec;
	TVector3 pN_par_q = nVec.Dot(qVec) / (qVec.Mag2()) * qVec;

	As = (E_n - p_n*gen_CosThetaNQ)/mN;
	Pt_vec = nVec - pN_par_q;
	Pt = Pt_vec.Mag();

}

void fillHist( double Q2, double Pt, double Xp, double As, TH1D**** hist , double weight ){
	

	// If it's larger than max Q2 or smaller than min Q2, return
	if( Q2 > Q2Bins[bins_Q2] 	) return;
	if( Q2 < Q2Bins[0] 		) return;

	// Need to figure out which bin of Q2, Pt, As (then fill as Xp):
	int this_bin_q2 = -1;
	int this_bin_pt = -1;
	int this_bin_as = -1;
	
	for( int q2_bin = bins_Q2-1 ; q2_bin >= 0; --q2_bin ){
		if( Q2 > Q2Bins[q2_bin] ) this_bin_q2 = q2_bin;
		if( this_bin_q2 != -1 ) break;
	}

	if( Pt < Pt_min			) return;
	if( Pt > Pt_max			) return;
	if( As < As_min		 	) return;
	if( As > As_max		 	) return;

	this_bin_pt = (int) ((Pt - Pt_min)/Pt_step);
	this_bin_as = (int) ((As - As_min)/As_step);

	// Safety clauses
	if( this_bin_q2 == -1 ){ cerr << "how\n"; return; }
	if( this_bin_pt == -1 ){ cerr << "how\n"; return; }
	if( this_bin_as == -1 ){ cerr << "how\n"; return; }
	if( this_bin_q2 > bins_Q2-1 ){ cerr << "how\n"; return; }
	if( this_bin_pt > bins_Pt-1 ){ cerr << "how\n"; return; }
	if( this_bin_as > bins_As-1 ){ cerr << "how\n"; return; }

	hist[this_bin_q2][this_bin_pt][this_bin_as]->Fill( Xp , weight );

}
void setError( double * err , TH1D * rec , TH1D * gen ){
	
	for( int bin = 1 ; bin < rec->GetXaxis()->GetNbins(); ++bin ){

		double r = rec->GetBinContent(bin);
		double g = gen->GetBinContent(bin);

		double e = sqrt( r/(g*g) + r*r/(g*g*g) );
		cout << bin << " " << r << " " << g << " " << r/g << " " << e << "\n";
		if( e!=e ) e = 0;
		err[bin-1] = e;
		cout << err[bin-1] << "\n\n";
	}


	return;
}

bool bad_bar( bandhit * this_n ){
	
	if(	this_n->getSector() == 1 && this_n->getComponent() == 1 				) return true;
	
	/*
	if(	this_n->getSector() == 3 && this_n->getLayer() == 4 && this_n->getComponent() == 2 	) return true;
	if(	this_n->getSector() == 4 && this_n->getLayer() == 3 && this_n->getComponent() == 2 	) return true;
	if(	this_n->getSector() == 4 && this_n->getLayer() == 2 && this_n->getComponent() == 5 	) return true;
	if(	this_n->getSector() == 1 && this_n->getLayer() == 1 && this_n->getComponent() == 1 	) return true;
	if(	this_n->getSector() == 1 && this_n->getLayer() == 2 && this_n->getComponent() == 1 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 2 && this_n->getComponent() == 5 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 2 && this_n->getComponent() == 7 	) return true;
	if(	this_n->getSector() == 4 && this_n->getLayer() == 2 && this_n->getComponent() == 4 	) return true;
	if(	this_n->getSector() == 4 && this_n->getLayer() == 2 && this_n->getComponent() == 6 	) return true;
	if(	this_n->getSector() == 5 && this_n->getLayer() == 2 && this_n->getComponent() == 1 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 3 && this_n->getComponent() == 6 	) return true;
	if(	this_n->getSector() == 4 && this_n->getLayer() == 3 && this_n->getComponent() == 1 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 4 && this_n->getComponent() == 1 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 4 && this_n->getComponent() == 2 	) return true;
	if(	this_n->getSector() == 2 && this_n->getLayer() == 4 && this_n->getComponent() == 4 	) return true;
	if(	this_n->getSector() == 3 && this_n->getLayer() == 4 && this_n->getComponent() == 6 	) return true;
	if(	this_n->getSector() == 3 && this_n->getLayer() == 5 && this_n->getComponent() == 3 	) return true;
	*/
	return false;
}
