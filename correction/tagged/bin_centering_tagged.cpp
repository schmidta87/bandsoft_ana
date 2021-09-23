// Take in:
// 	- simulated yield (per bin)
// 	- bin-migration correction (per bin)
// 	- radiative correction (per bin)
//	- acceptance correction (per bin)
//	- cross-section for average point (per bin)
// Output:
// 	- bin-centering correction

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"

#include "../../include/bin_edges.h"

void loadDiffCS( double mem[bins_Q2][bins_Pt][bins_Xb][bins_As] ){

	for( int i = 0 ; i < bins_Q2 ; ++i ){
		for( int j = 0 ; j < bins_Pt ; ++j ){
			for( int k = 0 ; k < bins_Xb ; ++k ){
				for( int m = 0 ; m < bins_As ; ++m ){
					mem[i][j][k][m] = 0.;
				}
			}
		}
	}

	std::string line;
	std::ifstream f;
	f.open("../point_cs/sim_bin_avgpoints_pointcs.txt");
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			int i_binq2, i_binpt, i_binxb, i_binas;
			double i_q2, i_pt, i_xb, i_as, diffcs;
			i_binq2 = i_binpt = i_binxb = i_binas = 
				i_q2 = i_pt = i_xb = i_as = diffcs = 0.;
			ss >> i_binq2 >> i_binpt >> i_binxb >> i_binas >>
				i_q2 >> i_pt >> i_xb >> i_as >> diffcs;
			
			mem[i_binq2][i_binpt][i_binxb][i_binas] = diffcs;
		}
	}
	f.close();


	
	return;
}

void bin_centering_tagged(TString path){

	double TaggedDiffCS[bins_Q2][bins_Pt][bins_Xb][bins_As];
	void loadDiffCS( double mem[bins_Q2][bins_Pt][bins_Xb][bins_As] );
	loadDiffCS( TaggedDiffCS );

	//Files for Tagged yield simulation
	TFile * inFileYieldTag = new TFile(path+"/yield_tagged.root");

	//Files for Bin Migration correction
	TFile * inFileCorrBMTag = new TFile(path+"/migcorrection_tagged.root");

	//Files for Radiative corrections
	TFile * inFileCorrRadTag = new TFile(path+"/radcorrection_tagged.root");

	//Files for Acceptance/Phase Space corrections
	TFile * inFileCorrAccTag = new TFile(path+"/acceptance_tagged.root");

	//Histos for copying
	TH1D **** h4_simulation_tagged_yield 	= new TH1D***[bins_Q2];
	TH1D **** h4_tagged_binmigration 	= new TH1D***[bins_Q2];
	TH1D **** h4_tagged_radiative 		= new TH1D***[bins_Q2];
	TH1D **** h4_tagged_acceptance 		= new TH1D***[bins_Q2];
	for( int i = 0 ; i < bins_Q2 ; ++i ){
		h4_simulation_tagged_yield 	[i]= new TH1D**[bins_Pt];
		h4_tagged_binmigration 		[i]= new TH1D**[bins_Pt];
		h4_tagged_radiative 		[i]= new TH1D**[bins_Pt];
		h4_tagged_acceptance 		[i]= new TH1D**[bins_Pt];
		for( int j = 0 ; j < bins_Pt ; ++j ){
			h4_simulation_tagged_yield 	[i][j]= new TH1D*[bins_Xb];
			h4_tagged_binmigration 		[i][j]= new TH1D*[bins_Xb];
			h4_tagged_radiative 		[i][j]= new TH1D*[bins_Xb];
			h4_tagged_acceptance 		[i][j]= new TH1D*[bins_Xb];
			for( int k = 0 ; k < bins_Xb ; ++k ){
				h4_simulation_tagged_yield 	[i][j][k]=(TH1D*) inFileYieldTag	->Get(Form("sim_yield_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_tagged_binmigration 		[i][j][k]=(TH1D*) inFileCorrBMTag	->Get(Form("migcorr_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_tagged_radiative 		[i][j][k]=(TH1D*) inFileCorrRadTag	->Get(Form("radcorr_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_tagged_acceptance 		[i][j][k]=(TH1D*) inFileCorrAccTag	->Get(Form("acceptance_tagged_as_Q2_%i_Pt_%i_Xb_%i",i,j,k));
				h4_simulation_tagged_yield[i][j][k]->Divide( h4_tagged_binmigration[i][j][k]	);
				h4_simulation_tagged_yield[i][j][k]->Divide( h4_tagged_radiative[i][j][k]	);
				h4_simulation_tagged_yield[i][j][k]->Divide( h4_tagged_acceptance[i][j][k]	);

				for( int m = 0 ; m < bins_As ; ++m ){
					double curr_bin_content = h4_simulation_tagged_yield[i][j][k]->GetBinContent( m+1 ) / (3.94707e+09/(1E-06));
					double diff_cs = TaggedDiffCS[i][j][k][m];
					if( diff_cs == 0 || diff_cs != diff_cs ) continue;
					cout << i << " " << j << " " << k << " " << m << " " << curr_bin_content << " " << diff_cs << "\n";
					h4_simulation_tagged_yield[i][j][k]->SetBinContent( m+1 , curr_bin_content / diff_cs );
					h4_simulation_tagged_yield[i][j][k]->SetBinError( m+1 , 0 );
					
				}

				if( k > 0 ){
					// Now we can take ratios of these bin-centering corrections for x / reference-x
					h4_simulation_tagged_yield[i][j][k]->Divide( h4_simulation_tagged_yield[i][j][0] );

					TCanvas * c1 = new TCanvas(Form("c_bincentering_Q2_%i_Pt_%i_Xb_%i",i,j,k),"",1200,1200);
					c1->cd();
					h4_simulation_tagged_yield[i][j][k]->Draw("*P");
					h4_simulation_tagged_yield[i][j][k]->GetXaxis()->SetTitle("a_{S}");

					c1->Print(Form("bincentering_Q2_%i_Pt_%i_Xb_%i.pdf",i,j,k));
				}
			}
		}
	}

	return;
}

