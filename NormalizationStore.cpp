#include <TMinuit.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>      
#include "TComplex.h" 
#include "TMath.h"  
#include "TLorentzVector.h"
#include <map>
#include "cstdlib"
#include <vector>
#include <utility>
#include <time.h> 
#include "TH1F.h"
#include "TH2F.h"
#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
#include <TRandom3.h>
#include "TVirtualFitter.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
//#include "omp.h"

using namespace std; 
  
int num_reso ; //including LASS/Kpi S-wave
int allnum_reso ;

map<string, double> amplmap ;
map<string, double> phasemap ;
map<string, double> KMamplmap ;
map<string, double> KMphasemap ;
map<string, double> errKMamplmap ;
map<string, double> errKMphasemap ;
map<string, double> kmatAR ;
map<string, double> kmatPH ;
map<string, double> lassparammap ;
map<string, double> erramplmap ;
map<string, double> errphasemap ;

const double mD0 = 1.86483;
const double mKs = 0.49761;
const double mPi = 0.13957;

double prob_dcs1=0., prob_dcs2=0., prob_dcs3=0., prob_dcs4=0.;

int fixflag_ar[50], fixflag_ph[50];
double effi[200][200];
int datanumber_ks, datanumber_kl;
int taggednormnumber_ks[3], taggednormnumber_kl[3] ;
Double_t normnumber_ks=0., normnumber_kl=0.;

TComplex value0(0.,0.);
TComplex value1(0.,0.);

TComplex Z[20], L[5], F[5], I[5][5];
TComplex Z_dcs[20];

double mod2L_ks[5][3], mod2F_ks[5][3];
double mod2L_kl[5][3], mod2F_kl[5][3];

double probZ_ks[20][3], probZ_dcs_ks[20][3];
double reinterZ_ks[20][20][3], iminterZ_ks[20][20][3], reinterKM0_ks[5][5][3], iminterKM0_ks[5][5][3], reinterKM1_ks[5][5][3], iminterKM1_ks[5][5][3], reinterI_ks[5][5][3], iminterI_ks[5][5][3] ;
double reinterZ_cfdcs_ks[20][20][3], iminterZ_cfdcs_ks[20][20][3];
double reinterZ_dcs_ks[20][20][3], iminterZ_dcs_ks[20][20][3];
double  reinterZval0_ks[20][5][3], iminterZval0_ks[20][5][3], reinterZval1_ks[20][5][3], iminterZval1_ks[20][5][3] ;
double  reinterval0Z_ks[20][5][3], iminterval0Z_ks[20][5][3], reinterval1Z_ks[20][5][3], iminterval1Z_ks[20][5][3] ;
double  reinterZval0_dcs_ks[20][5][3], iminterZval0_dcs_ks[20][5][3], reinterZval1_dcs_ks[20][5][3], iminterZval1_dcs_ks[20][5][3] ;

double probZ_kl[20][3], probZ_dcs_kl[20][3];
double reinterZ_kl[20][20][3], iminterZ_kl[20][20][3], reinterKM0_kl[5][5][3], iminterKM0_kl[5][5][3], reinterKM1_kl[5][5][3], iminterKM1_kl[5][5][3], reinterI_kl[5][5][3], iminterI_kl[5][5][3] ;
double reinterZ_cfdcs_kl[20][20][3], iminterZ_cfdcs_kl[20][20][3];
double reinterZ_dcs_kl[20][20][3], iminterZ_dcs_kl[20][20][3];
double  reinterZval0_kl[20][5][3], iminterZval0_kl[20][5][3], reinterZval1_kl[20][5][3], iminterZval1_kl[20][5][3] ;
double  reinterval0Z_kl[20][5][3], iminterval0Z_kl[20][5][3], reinterval1Z_kl[20][5][3], iminterval1Z_kl[20][5][3] ;
double  reinterZval0_dcs_kl[20][5][3], iminterZval0_dcs_kl[20][5][3], reinterZval1_dcs_kl[20][5][3], iminterZval1_dcs_kl[20][5][3] ;

std::vector<TComplex> u1j;
std::vector<TComplex> u1jdcs;
std::vector<TComplex> U1j;

// Define the complex coupling constants
double g[5][5]; // g[Physical pole]Decay channel]
int TAG_ks[50000], TAG_kl[50000];

Double_t CMplus2_ks[200000], CMminus2_ks[200000];
std::map<string, Double_t> CMulfacbw_re_ks[200000];
std::map<string, Double_t> CMulfacbw_im_ks[200000];
std::map<string, Double_t> CLASS_contribution_re_ks[200000];
std::map<string, Double_t> CLASS_contribution_im_ks[200000];

std::vector <double> CU1j_re_ks ;
std::vector <double> CU1j_im_ks ; 
std::vector <double> CBB_ks[40000] ;

TComplex A0_ks[50000];
TComplex A1_ks[50000];
TComplex A2_ks[50000];
TComplex A3_ks[50000];
TComplex A4_ks[50000];

TComplex U1jCal_ks[50000][5];

Double_t CMplus2_kl[50000], CMminus2_kl[50000];
std::map<string, Double_t> CMulfacbw_re_kl[50000];
std::map<string, Double_t> CMulfacbw_im_kl[50000];
std::map<string, Double_t> CLASS_contribution_re_kl[50000];
std::map<string, Double_t> CLASS_contribution_im_kl[50000];

std::vector <double> CU1j_re_kl ;
std::vector <double> CU1j_im_kl ; 
std::vector <double> CBB_kl[50000] ;

TComplex A0_kl[50000];
TComplex A1_kl[50000];
TComplex A2_kl[50000];
TComplex A3_kl[50000];
TComplex A4_kl[50000];

TComplex U1jCal_kl[50000][5];

string resonance[30], state[30];
double fixedpar[50];

TH2D *histo = new TH2D("histo", "histo", 100, 0.3, 3.0, 100, 0.3, 3.0);
TH1F *histu1j = new TH1F("histu1j", "", 100, -2.0, 2.0);

double normprobcal_ks[3], normprobcal_kl[3], normprobcal_ksbkg[3];

TH1D *pion = new TH1D("pion gamma","",10,0.,1.) ;
Double_t pi_track_gamma[10], pi_pid_gamma[10] ;
double out[4], out_ks[4];
Double_t piplus_trk_gamma, piminus_trk_gamma, piplus_pid_gamma, piminus_pid_gamma, gamma_k0l, gamma_k0s;

//Int_t pip_bin = pion->GetXaxis()->FindBin(xx_piplus) ;
//Double_t gamma_sig_piplus = pi_track_gamma[pip_bin-1]*pi_pid_gamma[pip_bin-1] ;
//Int_t pim_bin = pion->GetXaxis()->FindBin(xx_piminus) ;
//Double_t gamma_sig_piminus = pi_track_gamma[pim_bin-1]*pi_pid_gamma[pim_bin-1] ;
//Double_t gamma_k0s = 1.+(out_ks[0]*TMath::Exp(-1.*out_ks[1]*xx_k0s) + out_ks[2]*xx_k0s + out_ks[3])/100. ; //1.+(out_ks[0] + out_ks[1]*xx_k0s + out_ks[2]/xx_k0s)/100. ;
/*
////////////////////////// K0Lpipi effi syst. ////////////////////////

double Gamma_Factor_k0lpipi(double xx_k0l, double xx_piplus, double xx_piminus) {

//K0L rec.
Double_t gamma_k0l = 1.+(TMath::Exp(out[0]+out[1]*xx_k0l) + out[2]*xx_k0l + out[3])/100. ;

//pi+/- track
piplus_trk_gamma = (TMath::Exp(-1.1890-23.4719*xx_piplus) + 0.0033*xx_piplus - 0.0014) + 1. ;
piminus_trk_gamma = (TMath::Exp(-1.1890-23.4719*xx_piminus) + 0.0033*xx_piminus - 0.0014) + 1. ;

//pi+/- PID
piplus_pid_gamma = (-.0031285 - .0201579 * xx_piplus - .0344877 * xx_piplus * xx_piplus) + 1. ;
piminus_pid_gamma = (-.0031285 - .0201579 * xx_piminus - .0344877 * xx_piminus * xx_piminus) + 1. ; 

return gamma_k0l*piplus_trk_gamma*piplus_pid_gamma*piminus_trk_gamma*piminus_pid_gamma ;
}

//////////////////// K0Spipi effi syst. ////////////////////

double Gamma_Factor_k0spipi(double xx_k0s, double xx_piplus, double xx_piminus) {

//K0S rec. 
if(xx_k0s<0.4) gamma_k0s = .0730142-.159545*xx_k0s + 1.;
else gamma_k0s = .00953762 + 1.;

//pi+/- track
piplus_trk_gamma = (TMath::Exp(-1.1890-23.4719*xx_piplus) + 0.0033*xx_piplus - 0.0014) + 1. ;
piminus_trk_gamma = (TMath::Exp(-1.1890-23.4719*xx_piminus) + 0.0033*xx_piminus - 0.0014) + 1. ;

//pi+/- PID
piplus_pid_gamma = (-.0031285 - .0201579 * xx_piplus - .0344877 * xx_piplus * xx_piplus) + 1. ;
piminus_pid_gamma = (-.0031285 - .0201579 * xx_piminus - .0344877 * xx_piminus * xx_piminus) + 1. ; 

return gamma_k0s*piplus_trk_gamma*piplus_pid_gamma*piminus_trk_gamma*piminus_pid_gamma ;
}
*/
// Gamma Fit parameters

//K0S
Double_t c0_k0s = 17.2306, c1_k0s = -11.3357, c2_k0s = 0.9741;
//Double_t c0_k0s = 7.3014; Double_t c1_k0s = -15.9545; Double_t c2_k0s = 0.9538 ;

//pi tracking
Double_t c0_pi_trk = 0.0182; Double_t c1_pi_trk = -0.1066; Double_t c2_pi_trk = 0.1888; Double_t c3_pi_trk = -0.1006;
//pi PID
Double_t c0_pi_pid = -0.003128; Double_t c1_pi_pid = 0.020158; Double_t c2_pi_pid = -0.034488;
//K0L
Double_t c0_k0l = 25.9620; Double_t c1_k0l = -10.8760; Double_t c2_k0l = 4.0073;

////////////// K0Lpipi effi syst. ////////////////

double Gamma_Factor_k0lpipi(double xx_k0l, double xx_piplus, double xx_piminus) {

//K0L rec.
Double_t gamma_k0l = 1. + (c0_k0l*TMath::Exp(c1_k0l*xx_k0l) + c2_k0l)/100. ;

//pi+/- track
piplus_trk_gamma  = 1. + (c0_pi_trk + c1_pi_trk*xx_piplus  +  c2_pi_trk*xx_piplus*xx_piplus   +  c3_pi_trk*xx_piplus*xx_piplus*xx_piplus) ;
piminus_trk_gamma = 1. + (c0_pi_trk + c1_pi_trk*xx_piminus +  c2_pi_trk*xx_piminus*xx_piminus +  c3_pi_trk*xx_piminus*xx_piminus*xx_piminus) ;

//pi+/- PID
piplus_pid_gamma  = 1. + (c0_pi_pid + c1_pi_pid * xx_piplus  + c2_pi_pid * xx_piplus * xx_piplus) ;
piminus_pid_gamma = 1. + (c0_pi_pid + c1_pi_pid * xx_piminus + c2_pi_pid * xx_piminus * xx_piminus) ;

return gamma_k0l*piplus_trk_gamma*piplus_pid_gamma*piminus_trk_gamma*piminus_pid_gamma ;
//return 1.;
}

///////////// K0Spipi effi syst. ///////////////////

double Gamma_Factor_k0spipi(double xx_k0s, double xx_piplus, double xx_piminus) {

//K0S rec.
gamma_k0s = 1. + (c0_k0s*TMath::Exp(c1_k0s*xx_k0s) + c2_k0s)/100. ;
//if(xx_k0s < 0.4) gamma_k0s = 1. + (c0_k0s + c1_k0s*xx_k0s)/100. ;
//else             gamma_k0s = 1. + (c2_k0s)/100. ;

//pi+/- track
piplus_trk_gamma  = 1. + (c0_pi_trk + c1_pi_trk*xx_piplus  +  c2_pi_trk*xx_piplus*xx_piplus   +  c3_pi_trk*xx_piplus*xx_piplus*xx_piplus) ;
piminus_trk_gamma = 1. + (c0_pi_trk + c1_pi_trk*xx_piminus +  c2_pi_trk*xx_piminus*xx_piminus +  c3_pi_trk*xx_piminus*xx_piminus*xx_piminus) ;

//pi+/- PID
piplus_pid_gamma  = 1. + (c0_pi_pid + c1_pi_pid * xx_piplus  + c2_pi_pid * xx_piplus * xx_piplus) ;
piminus_pid_gamma = 1. + (c0_pi_pid + c1_pi_pid * xx_piminus + c2_pi_pid * xx_piminus * xx_piminus) ;

return gamma_k0s*piplus_trk_gamma*piplus_pid_gamma*piminus_trk_gamma*piminus_pid_gamma ;
//return 1.;
}

////////////////////////////////////////////////////////////

int main()
{
  time_t my_time = time(NULL);

  double efficiency;
  int ii, jj;

  double mass, width, magni, phase, value, xxx, yyy, errxxx, erryyy, errmagni, errphase;
  string mode;

  Double_t mplus2, mminus2, pk0l, pk0s, ppip, ppim;
  Double_t Mplus2, Mminus2;

  std::map<string, Double_t> *mulfacbw_re = 0;
  std::map<string, Double_t> *mulfacbw_im = 0;
  std::map<string, Double_t> *Mulfacbw_re = 0;
  std::map<string, Double_t> *Mulfacbw_im = 0;
  std::map<string, Double_t> *lass_contribution_re = 0;
  std::map<string, Double_t> *lass_contribution_im = 0;
  std::map<string, Double_t> *LASS_contribution_re = 0;
  std::map<string, Double_t> *LASS_contribution_im = 0;

  std::vector<Double_t> *u1j_re = 0;
  std::vector<Double_t> *u1j_im = 0;
  std::vector<Double_t> *U1j_re = 0;
  std::vector<Double_t> *U1j_im = 0; 
  std::vector<Double_t> *BB = 0;
  std::vector<Double_t> *bb = 0;

  Int_t tag=0, norm_tag=0;
  Double_t mplus2dcs, mminus2dcs;
  Double_t Mplus2dcs, Mminus2dcs;

  std::map<string, Double_t> *mulfacbw_redcs = 0;
  std::map<string, Double_t> *mulfacbw_imdcs = 0;
  std::map<string, Double_t> *Mulfacbw_redcs = 0;
  std::map<string, Double_t> *Mulfacbw_imdcs = 0;
  std::map<string, Double_t> *lass_contribution_redcs = 0;
  std::map<string, Double_t> *lass_contribution_imdcs = 0;
  std::map<string, Double_t> *LASS_contribution_redcs = 0;
  std::map<string, Double_t> *LASS_contribution_imdcs = 0;

  std::vector<Double_t> *u1j_redcs = 0;
  std::vector<Double_t> *u1j_imdcs = 0;
  std::vector<Double_t> *U1j_redcs = 0;
  std::vector<Double_t> *U1j_imdcs = 0; 
  std::vector<Double_t> *BBdcs = 0;
  std::vector<Double_t> *bbdcs = 0;

  double Asum_re[5], Asum_im[5] ;

  TFile *filein_ks = TFile::Open("/home/anita/AmplitudeFitter/AmpFitter_01Sep2020/ForNormFactorRealData6C.root");
  TTree *tree1_ks;
  TTree *tree2_ks;

  filein_ks->GetObject("Norm", tree1_ks);
  tree1_ks->SetBranchAddress("norm_tag", &norm_tag);
  tree1_ks->SetBranchAddress("mplus2", &mplus2);				
  tree1_ks->SetBranchAddress("mminus2", &mminus2);	
  tree1_ks->SetBranchAddress("p_k0s", &pk0s);		
  tree1_ks->SetBranchAddress("p_pip", &ppip);		
  tree1_ks->SetBranchAddress("p_pim", &ppim);		
  tree1_ks->SetBranchAddress("mulfacbw_re", &mulfacbw_re);				//map
  tree1_ks->SetBranchAddress("mulfacbw_im", &mulfacbw_im);				//map
  tree1_ks->SetBranchAddress("u1j_re", &u1j_re);					//vector
  tree1_ks->SetBranchAddress("u1j_im", &u1j_im);					//vector
  tree1_ks->SetBranchAddress("bb", &bb);						//vector
  tree1_ks->SetBranchAddress("lass_contribution_re", &lass_contribution_re);	//map
  tree1_ks->SetBranchAddress("lass_contribution_im", &lass_contribution_im);	//map

  filein_ks->GetObject("datasample", tree2_ks);
  tree2_ks->SetBranchAddress("tag", &tag);
  tree2_ks->SetBranchAddress("Mplus2", &Mplus2);
  tree2_ks->SetBranchAddress("Mminus2", &Mminus2);
  tree2_ks->SetBranchAddress("Mulfacbw_re", &Mulfacbw_re);
  tree2_ks->SetBranchAddress("Mulfacbw_im", &Mulfacbw_im);
  tree2_ks->SetBranchAddress("U1j_re", &U1j_re);
  tree2_ks->SetBranchAddress("U1j_im", &U1j_im);
  tree2_ks->SetBranchAddress("BB", &BB);
  tree2_ks->SetBranchAddress("LASS_contribution_re", &LASS_contribution_re);
  tree2_ks->SetBranchAddress("LASS_contribution_im", &LASS_contribution_im);
  tree2_ks->SetBranchAddress("Asum_re", Asum_re);
  tree2_ks->SetBranchAddress("Asum_im", Asum_im);

  datanumber_ks = tree2_ks->GetEntries();
  cout<<"datanumber_ks: "<<datanumber_ks<<endl;

  TFile *filein_kl = TFile::Open("MulFactorRealDataEffiDCS.root");
  TTree *tree1_kl;
  TTree *tree2_kl;

  filein_kl->GetObject("Norm", tree1_kl);
  tree1_kl->SetBranchAddress("norm_tag", &norm_tag);
  tree1_kl->SetBranchAddress("mplus2", &mplus2);				
  tree1_kl->SetBranchAddress("mminus2", &mminus2);	
  tree1_kl->SetBranchAddress("p_k0l", &pk0l);	
  tree1_kl->SetBranchAddress("p_pip", &ppip);		
  tree1_kl->SetBranchAddress("p_pim", &ppim);		
  tree1_kl->SetBranchAddress("mulfacbw_re", &mulfacbw_re);				//map
  tree1_kl->SetBranchAddress("mulfacbw_im", &mulfacbw_im);				//map
  tree1_kl->SetBranchAddress("u1j_re", &u1j_re);					//vector
  tree1_kl->SetBranchAddress("u1j_im", &u1j_im);					//vector
  tree1_kl->SetBranchAddress("bb", &bb);						//vector
  tree1_kl->SetBranchAddress("lass_contribution_re", &lass_contribution_re);	//map
  tree1_kl->SetBranchAddress("lass_contribution_im", &lass_contribution_im);	//map

  filein_kl->GetObject("datasample", tree2_kl);
  tree2_kl->SetBranchAddress("tag", &tag);
  tree2_kl->SetBranchAddress("Mplus2", &Mplus2);
  tree2_kl->SetBranchAddress("Mminus2", &Mminus2);
  tree2_kl->SetBranchAddress("Mulfacbw_re", &Mulfacbw_re);
  tree2_kl->SetBranchAddress("Mulfacbw_im", &Mulfacbw_im);
  tree2_kl->SetBranchAddress("U1j_re", &U1j_re);
  tree2_kl->SetBranchAddress("U1j_im", &U1j_im);
  tree2_kl->SetBranchAddress("BB", &BB);
  tree2_kl->SetBranchAddress("LASS_contribution_re", &LASS_contribution_re);
  tree2_kl->SetBranchAddress("LASS_contribution_im", &LASS_contribution_im);
  tree2_kl->SetBranchAddress("Asum_re", Asum_re);
  tree2_kl->SetBranchAddress("Asum_im", Asum_im);

  datanumber_kl = tree2_kl->GetEntries();
  cout<<"datanumber_kl: "<<datanumber_kl<<endl;

  int flag_ar=0, flag_ph=0;

  string pp;

  ifstream fin4;
  fin4.open("/home/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/LASSFittedParameterInfoBelle2018.dat");

  while (1) {  
               if (!fin4.good()) break;
	       fin4 >> pp >> value ;
               lassparammap.insert(pair<string, double>(pp,value));
            }     
                   
  int i=0;        
  int j=0;        
  int k=0;        


  string cpresonance[10];

  k=0;


  int l=0;
  string modecf, modedcs, reso_cf[20], reso_dcs[20];
  ifstream fin5;
  fin5.open("/home/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/AllResonances.dat") ;
  while(1) { 
             fin5 >> modecf >> modedcs ;
	     if (!fin5.good()) break;
	     if(modecf[0] == '#') continue;
             reso_cf[l] = modecf ;
	     reso_dcs[l] = modedcs ;
             cout<<reso_cf[l]<<"\t"<<reso_dcs[l]<<endl;
             l++;
           }

  num_reso = l;
  cout<<"num_reso: "<<num_reso<<endl;
            
///////////////////////////////////////////////////////// Normalization ////////////////////////////////////////////////////////
                  
for(Int_t r=0; r<20; r++)
   { Z[r](0.0, 0.0); 
     Z_dcs[r](0., 0.);
for(int t=0; t<3; t++) {
     probZ_ks[r][t]=0.0; 
     probZ_dcs_ks[r][t]=0.0;
     probZ_kl[r][t]=0.0; 
     probZ_dcs_kl[r][t]=0.0;
     for(Int_t rr=0; rr<20; rr++)
        { reinterZ_ks[r][rr][t]=0.0;
          iminterZ_ks[r][rr][t]=0.0;

	  reinterZ_dcs_ks[r][rr][t]=0.0;
          iminterZ_dcs_ks[r][rr][t]=0.0;

	  reinterZ_cfdcs_ks[r][rr][t]=0.0;
          iminterZ_cfdcs_ks[r][rr][t]=0.0;

	  reinterZ_kl[r][rr][t]=0.0;
          iminterZ_kl[r][rr][t]=0.0;

	  reinterZ_dcs_kl[r][rr][t]=0.0;
          iminterZ_dcs_kl[r][rr][t]=0.0;

	  reinterZ_cfdcs_kl[r][rr][t]=0.0;
          iminterZ_cfdcs_kl[r][rr][t]=0.0;
        }
     for(Int_t l=0; l<5; l++)
        { reinterZval0_ks[r][l][t]=0.0;
          iminterZval0_ks[r][l][t]=0.0;
          reinterZval1_ks[r][l][t]=0.0;
          iminterZval1_ks[r][l][t]=0.0;

	  reinterZval0_dcs_ks[r][l][t]=0.0;
          iminterZval0_dcs_ks[r][l][t]=0.0;
          reinterZval1_dcs_ks[r][l][t]=0.0;
          iminterZval1_dcs_ks[r][l][t]=0.0;

	  reinterval0Z_ks[r][l][t]=0.0;
          iminterval0Z_ks[r][l][t]=0.0;
          reinterval1Z_ks[r][l][t]=0.0;
          iminterval1Z_ks[r][l][t]=0.0;

	  reinterZval0_kl[r][l][t]=0.0;
          iminterZval0_kl[r][l][t]=0.0;
          reinterZval1_kl[r][l][t]=0.0;
          iminterZval1_kl[r][l][t]=0.0;

	  reinterZval0_dcs_kl[r][l][t]=0.0;
          iminterZval0_dcs_kl[r][l][t]=0.0;
          reinterZval1_dcs_kl[r][l][t]=0.0;
          iminterZval1_dcs_kl[r][l][t]=0.0;

	  reinterval0Z_kl[r][l][t]=0.0;
          iminterval0Z_kl[r][l][t]=0.0;
          reinterval1Z_kl[r][l][t]=0.0;
          iminterval1Z_kl[r][l][t]=0.0;
        }
   }
   }

for(Int_t l=0; l<5; l++)
   { L[l](0.0, 0.0);  
     F[l](0.0, 0.0);  
for(int t=0; t<3; t++) { 
     mod2L_ks[l][t]=0.0;
     mod2F_ks[l][t]=0.0;
     mod2L_kl[l][t]=0.0;
     mod2F_kl[l][t]=0.0;
     for(Int_t k=0; k<5; k++)
        { reinterKM0_ks[l][k][t]=0.0;
          iminterKM0_ks[l][k][t]=0.0;
          reinterKM1_ks[l][k][t]=0.0;
          iminterKM1_ks[l][k][t]=0.0;
          reinterI_ks[l][k][t]  =0.0;
          iminterI_ks[l][k][t]  =0.0;

	  reinterKM0_kl[l][k][t]=0.0;
          iminterKM0_kl[l][k][t]=0.0;
          reinterKM1_kl[l][k][t]=0.0;
          iminterKM1_kl[l][k][t]=0.0;
          reinterI_kl[l][k][t]  =0.0;
          iminterI_kl[l][k][t]  =0.0;
        }
   }
   }

for(int t=0; t<3; t++) {
taggednormnumber_ks[t] = 0;
taggednormnumber_kl[t] = 0;
}

// pi+pi- channel
g[0][0]=0.22889;
g[1][0]=0.94128;
g[2][0]=0.36856;
g[3][0]=0.33650;
g[4][0]=0.18171;

// K+K- channel
g[0][1]=-0.55377;
g[1][1]=0.55095;
g[2][1]=0.23888;
g[3][1]=0.40907;
g[4][1]=-0.17558;

// 4pi channel
g[0][2]=0;
g[1][2]=0;
g[2][2]=0.55639;
g[3][2]=0.85679;
g[4][2]=-0.79658;

// eta eta channel
g[0][3]=-0.39899;
g[1][3]=0.39065;
g[2][3]=0.18340;
g[3][3]=0.19906;
g[4][3]=-0.00355;

//eta eta' channel
g[0][4]=-0.34639;
g[1][4]=0.31503;
g[2][4]=0.18681;
g[3][4]=-0.00984;
g[4][4]=0.22358;

std::cout<<"g[][] set..."<<std::endl;

ofstream fout;
fout.open("kskl_yields_check_111021.dat");

ofstream fout0 ;
ofstream fout1 ;
ofstream fout2 ;
ofstream fout3 ;
ofstream fout4 ;
fout0.open("r_t_check_111021.dat");
fout1.open("r_rr_t_check_111021.dat");
fout2.open("r_l_t_check_111021.dat");
fout3.open("l_t_check_111021.dat");
fout4.open("l_k_t_check_111021.dat");

pi_track_gamma[0]=1.0910;
pi_track_gamma[1]=1.0088;
pi_track_gamma[2]=0.9993;
pi_track_gamma[3]=0.9976;
pi_track_gamma[4]=0.9998;
pi_track_gamma[5]=1.0017;
pi_track_gamma[6]=1.0010;
pi_track_gamma[7]=1.0011;
pi_track_gamma[8]=1.0024;
pi_track_gamma[9]=1.0012; 

pi_pid_gamma[0]=0.9935;
pi_pid_gamma[1]=0.9990;
pi_pid_gamma[2]=0.9998;
pi_pid_gamma[3]=0.9992;
pi_pid_gamma[4]=0.9994;
pi_pid_gamma[5]=0.9985;
pi_pid_gamma[6]=0.9940;
pi_pid_gamma[7]=0.9914;
pi_pid_gamma[8]=0.9881;
pi_pid_gamma[9]=0.9888;

out_ks[0] = 2.89544 ;
out_ks[1] = -12.3569 ;
out_ks[2] = -0.283453 ; 
out_ks[3] = 1.22877 ; 

//Sum_i { |Sum_r A_r + A_KM|^2 * effi }				// i: events, r: resonances
//Sum_i { ( |Sum_r A_r|^2 + |A_KM|^2 + 2xRe[(Sum_r A_r)*A_KM] )xeffi }
//Part1 + Part2 + Part3

//------------------------------------------------------- Part1 --------------------------------------------------------------
//Sum_i{ ( |Sum_r A_r|^2}
//A_r^i = a_r e^{i phi_r} Z_r^i
// Sum_r{ |a_r|^2 ( Sum_i{|Z_r^i|^2} ) }   +   2 x Sum_i{Sum_b{Sum_a{Re(A*_b^i A_a^i) } } }

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// K0Spipi /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  for(Int_t j=0; j<tree1_ks->GetEntries(); j++)		
     { tree1_ks->GetEntry(j);	
       taggednormnumber_ks[norm_tag]++;
       normnumber_ks++;

       for(int i=0; i<(num_reso-2); i++) {	
       Z[i] = TComplex((*mulfacbw_re)[reso_cf[i]], (*mulfacbw_im)[reso_cf[i]]);
       Z_dcs[i] = TComplex((*mulfacbw_re)[reso_dcs[i]], (*mulfacbw_im)[reso_dcs[i]]);
       }
       
       Z[num_reso-2] = TComplex((*lass_contribution_re)["K0star1430minus"], (*lass_contribution_im)["K0star1430minus"]);	
       Z[num_reso-1] = TComplex((*lass_contribution_re)["K0star1430plus"], (*lass_contribution_im)["K0star1430plus"]);
       
       Z_dcs[num_reso-2] = TComplex((*lass_contribution_re)["K0star1430plus"], (*lass_contribution_im)["K0star1430plus"]);	
       Z_dcs[num_reso-1] = TComplex((*lass_contribution_re)["K0star1430minus"], (*lass_contribution_im)["K0star1430minus"]);

       for(Int_t r=0; r<num_reso; r++) {
       probZ_ks[r][norm_tag] = probZ_ks[r][norm_tag] + Z[r].Rho2() * Gamma_Factor_k0spipi(pk0s,ppip,ppim) ; 
       probZ_dcs_ks[r][norm_tag] = probZ_dcs_ks[r][norm_tag] + Z_dcs[r].Rho2() * Gamma_Factor_k0spipi(pk0s,ppip,ppim) ;	    
       }

       for(Int_t a=0; a<num_reso; a++) {
       for(Int_t b=0; b<num_reso; b++) {
       reinterZ_ks[a][b][norm_tag] = reinterZ_ks[a][b][norm_tag] + ((TComplex::Conjugate(Z[a]))*(Z[b])).Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);   
       iminterZ_ks[a][b][norm_tag] = iminterZ_ks[a][b][norm_tag] + ((TComplex::Conjugate(Z[a]))*(Z[b])).Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);   

       reinterZ_dcs_ks[a][b][norm_tag] = reinterZ_dcs_ks[a][b][norm_tag] + ((TComplex::Conjugate(Z_dcs[a]))*(Z_dcs[b])).Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);   
       iminterZ_dcs_ks[a][b][norm_tag] = iminterZ_dcs_ks[a][b][norm_tag] + ((TComplex::Conjugate(Z_dcs[a]))*(Z_dcs[b])).Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);   

       reinterZ_cfdcs_ks[a][b][norm_tag] = reinterZ_cfdcs_ks[a][b][norm_tag] + ((Z[a])*(TComplex::Conjugate(Z_dcs[b]))).Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 
       iminterZ_cfdcs_ks[a][b][norm_tag] = iminterZ_cfdcs_ks[a][b][norm_tag] + ((Z[a])*(TComplex::Conjugate(Z_dcs[b]))).Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 

       }
       }

//------------------------ Part2 -------------------------------

//Sum_i{|A_KM|^2}
//Sum_i{|value0|^2} + Sum_i{|value1|^2} + 2xSum_i{Re(value0* x value1)}

       u1j.clear();
       for(Int_t l=0; l<5; l++)
          { u1j.push_back(TComplex((*u1j_re)[l], (*u1j_im)[l]));
	    //U1jCal.push_back(TComplex((U1jreCal)[l], (U1jimCal)[l])) ;
          }

// Sum_i{|value0|^2}
   
       for(Int_t l=0; l<5; l++)
          { L[l] = TComplex(0.0, 0.0);
            F[l] = TComplex(0.0, 0.0);
            for(Int_t k=0; k<5; k++)
               { I[l][k] = TComplex(0.0, 0.0);
               }
          }
       
       for(Int_t l=0; l<5; l++)
          { L[0] = L[0] + (g[0][l]*u1j[l])/(*bb)[0]  ;			//beta[0]
            L[1] = L[1] + (g[1][l]*u1j[l])/(*bb)[1]  ;			//beta[1]
            L[2] = L[2] + (g[2][l]*u1j[l])/(*bb)[2]  ;			//beta[2]
            L[3] = L[3] + (g[3][l]*u1j[l])/(*bb)[3]  ;			//beta[3]
            L[4] = L[4] + (g[4][l]*u1j[l])/(*bb)[4]  ;			//beta[4]
          }

       for(Int_t k=0; k<5; k++)
          { mod2L_ks[k][norm_tag] = mod2L_ks[k][norm_tag] + L[k].Rho2() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);
          }

       for(Int_t a=0; a<4; a++)
          { for(Int_t b=a+1; b<5; b++)
               { reinterKM0_ks[a][b][norm_tag] = reinterKM0_ks[a][b][norm_tag] + ((TComplex::Conjugate(L[a]))*L[b]).Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);  
		 iminterKM0_ks[a][b][norm_tag] = iminterKM0_ks[a][b][norm_tag] + ((TComplex::Conjugate(L[a]))*L[b]).Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);  
	       }
          }

//Sum_i{|value1|^2}

       for(Int_t l=0; l<5; l++)
          { F[l] = u1j[l] * (*bb)[5];
            mod2F_ks[l][norm_tag] = mod2F_ks[l][norm_tag] + F[l].Rho2() * Gamma_Factor_k0spipi(pk0s,ppip,ppim) ;
          }       

       for(Int_t a=0; a<4; a++)
          { for(Int_t b=a+1; b<5; b++)
               { reinterKM1_ks[a][b][norm_tag] = reinterKM1_ks[a][b][norm_tag] + ((TComplex::Conjugate(F[a]))*F[b]).Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 	
		 iminterKM1_ks[a][b][norm_tag] = iminterKM1_ks[a][b][norm_tag] + ((TComplex::Conjugate(F[a]))*F[b]).Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 	
	       }
          }

//2xSum_i{Re(value0* x value1)}

       for(Int_t a=0; a<5; a++)
          { for(Int_t b=0; b<5; b++)
               { I[a][b] = u1j[b]*(*bb)[5]*(TComplex::Conjugate(L[a])) ;
                 reinterI_ks[a][b][norm_tag] = reinterI_ks[a][b][norm_tag] + I[a][b].Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);
                 iminterI_ks[a][b][norm_tag] = iminterI_ks[a][b][norm_tag] + I[a][b].Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);
	       }
          }       

//----------------------------------- Part3 ---------------------------------------------------
//Sum_i{2xRe[(Sum_r{A_r})*A_KM] }
//Sum_i{2xRe[(Sum_r{a_r e^{iphi_r} Z_r})* (value0 + value1)]}
       
       for(Int_t r=0; r<num_reso; r++) {
       for(Int_t l=0; l<5; l++) {
       reinterZval0_ks[r][l][norm_tag] = reinterZval0_ks[r][l][norm_tag] + ((TComplex::Conjugate(Z[r]))*L[l]).Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);   
       iminterZval0_ks[r][l][norm_tag] = iminterZval0_ks[r][l][norm_tag] + ((TComplex::Conjugate(Z[r]))*L[l]).Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);   
      
       reinterZval1_ks[r][l][norm_tag] = reinterZval1_ks[r][l][norm_tag] + ((TComplex::Conjugate(Z[r]))*F[l]).Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);   
       iminterZval1_ks[r][l][norm_tag] = iminterZval1_ks[r][l][norm_tag] + ((TComplex::Conjugate(Z[r]))*F[l]).Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim);  


       reinterZval0_dcs_ks[r][l][norm_tag] = reinterZval0_dcs_ks[r][l][norm_tag] + ((TComplex::Conjugate(Z_dcs[r]))*L[l]).Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 
       iminterZval0_dcs_ks[r][l][norm_tag] = iminterZval0_dcs_ks[r][l][norm_tag] + ((TComplex::Conjugate(Z_dcs[r]))*L[l]).Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 
      
       reinterZval1_dcs_ks[r][l][norm_tag] = reinterZval1_dcs_ks[r][l][norm_tag] + ((TComplex::Conjugate(Z_dcs[r]))*F[l]).Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 
       iminterZval1_dcs_ks[r][l][norm_tag] = iminterZval1_dcs_ks[r][l][norm_tag] + ((TComplex::Conjugate(Z_dcs[r]))*F[l]).Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

       reinterval0Z_ks[r][l][norm_tag] = reinterval0Z_ks[r][l][norm_tag] + (Z[r]*TComplex::Conjugate(L[l])).Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 
       iminterval0Z_ks[r][l][norm_tag] = iminterval0Z_ks[r][l][norm_tag] + (Z[r]*TComplex::Conjugate(L[l])).Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 
      
       reinterval1Z_ks[r][l][norm_tag] = reinterval1Z_ks[r][l][norm_tag] + (Z[r]*TComplex::Conjugate(F[l])).Re() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 
       iminterval1Z_ks[r][l][norm_tag] = iminterval1Z_ks[r][l][norm_tag] + (Z[r]*TComplex::Conjugate(F[l])).Im() * Gamma_Factor_k0spipi(pk0s,ppip,ppim); 

       }
       }


     }       

  std::cout<<"K0S done.."<<endl;

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// K0Lpipi ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

//Sum_i { |Sum_r A_r + A_KM|^2 * effi }				// i: events, r: resonances
//Sum_i { ( |Sum_r A_r|^2 + |A_KM|^2 + 2xRe[(Sum_r A_r)*A_KM] )xeffi }
//Part1 + Part2 + Part3

//------------------------------------------------------- Part1 --------------------------------------------------------------
//Sum_i{ ( |Sum_r A_r|^2}
//A_r^i = a_r e^{i phi_r} Z_r^i
// Sum_r{ |a_r|^2 ( Sum_i{|Z_r^i|^2} ) }   +   2 x Sum_i{Sum_b{Sum_a{Re(A*_b^i A_a^i) } } }


  for(Int_t j=0; j<tree1_kl->GetEntries(); j++)			
     { tree1_kl->GetEntry(j);	
       taggednormnumber_kl[norm_tag]++;
       normnumber_kl++;

       for(int i=0; i<(num_reso-2); i++)	//subtract however many LASS are included in the minimization 	/////// MODIFY /////////
          { Z[i] = TComplex((*mulfacbw_re)[reso_cf[i]], (*mulfacbw_im)[reso_cf[i]]);
            Z_dcs[i] = TComplex((*mulfacbw_re)[reso_dcs[i]], (*mulfacbw_im)[reso_dcs[i]]);
          }
       
       Z[num_reso-2] = TComplex((*lass_contribution_re)["K0star1430minus"], (*lass_contribution_im)["K0star1430minus"]);	
       Z[num_reso-1] = TComplex((*lass_contribution_re)["K0star1430plus"], (*lass_contribution_im)["K0star1430plus"]);
       
       Z_dcs[num_reso-2] = TComplex((*lass_contribution_re)["K0star1430plus"], (*lass_contribution_im)["K0star1430plus"]);	
       Z_dcs[num_reso-1] = TComplex((*lass_contribution_re)["K0star1430minus"], (*lass_contribution_im)["K0star1430minus"]);

       for(Int_t r=0; r<num_reso; r++) {
       probZ_kl[r][norm_tag] = probZ_kl[r][norm_tag] + Z[r].Rho2() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim) ; 
       probZ_dcs_kl[r][norm_tag] = probZ_dcs_kl[r][norm_tag] + Z_dcs[r].Rho2() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim) ;
       }

       for(Int_t a=0; a<num_reso; a++) {
       for(Int_t b=0; b<num_reso; b++) {
       reinterZ_kl[a][b][norm_tag] = reinterZ_kl[a][b][norm_tag] + ((TComplex::Conjugate(Z[a]))*(Z[b])).Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);   
       iminterZ_kl[a][b][norm_tag] = iminterZ_kl[a][b][norm_tag] + ((TComplex::Conjugate(Z[a]))*(Z[b])).Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);   

       reinterZ_dcs_kl[a][b][norm_tag] = reinterZ_dcs_kl[a][b][norm_tag] + ((TComplex::Conjugate(Z_dcs[a]))*(Z_dcs[b])).Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);   
       iminterZ_dcs_kl[a][b][norm_tag] = iminterZ_dcs_kl[a][b][norm_tag] + ((TComplex::Conjugate(Z_dcs[a]))*(Z_dcs[b])).Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);   

       reinterZ_cfdcs_kl[a][b][norm_tag] = reinterZ_cfdcs_kl[a][b][norm_tag] + ((Z[a])*(TComplex::Conjugate(Z_dcs[b]))).Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 
       iminterZ_cfdcs_kl[a][b][norm_tag] = iminterZ_cfdcs_kl[a][b][norm_tag] + ((Z[a])*(TComplex::Conjugate(Z_dcs[b]))).Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 

       }
       }

//------------------------ Part2 -------------------------------

//Sum_i{|A_KM|^2}
//Sum_i{|value0|^2} + Sum_i{|value1|^2} + 2xSum_i{Re(value0* x value1)}

       u1j.clear();
       for(Int_t l=0; l<5; l++)
          { u1j.push_back(TComplex((*u1j_re)[l], (*u1j_im)[l]));
	    //U1jCal.push_back(TComplex((U1jreCal)[l], (U1jimCal)[l])) ;
          }

// Sum_i{|value0|^2}
   
       for(Int_t l=0; l<5; l++)
          { L[l] = TComplex(0.0, 0.0);
            F[l] = TComplex(0.0, 0.0);
            for(Int_t k=0; k<5; k++)
               { I[l][k] = TComplex(0.0, 0.0);
               }
          }
       
       for(Int_t l=0; l<5; l++)
          { L[0] = L[0] + (g[0][l]*u1j[l])/(*bb)[0]  ;			//beta[0]
            L[1] = L[1] + (g[1][l]*u1j[l])/(*bb)[1]  ;			//beta[1]
            L[2] = L[2] + (g[2][l]*u1j[l])/(*bb)[2]  ;			//beta[2]
            L[3] = L[3] + (g[3][l]*u1j[l])/(*bb)[3]  ;			//beta[3]
            L[4] = L[4] + (g[4][l]*u1j[l])/(*bb)[4]  ;			//beta[4]
          }

       for(Int_t k=0; k<5; k++)
          { mod2L_kl[k][norm_tag] = mod2L_kl[k][norm_tag] + L[k].Rho2() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);			
	    //if(norm_tag==0 && (L[k].Rho2()<0. || L[k].Rho2()>100.)) std::cout<<"L["<<k<<"].Rho2(): "<<L[k].Rho2()<<endl;
          } 

       for(Int_t a=0; a<4; a++)
          { for(Int_t b=a+1; b<5; b++)
               { reinterKM0_kl[a][b][norm_tag] = reinterKM0_kl[a][b][norm_tag] + ((TComplex::Conjugate(L[a]))*L[b]).Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);  
		 iminterKM0_kl[a][b][norm_tag] = iminterKM0_kl[a][b][norm_tag] + ((TComplex::Conjugate(L[a]))*L[b]).Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 
	       }
          }

//Sum_i{|value1|^2}

       for(Int_t l=0; l<5; l++)
          { F[l] = u1j[l] * (*bb)[5];
            mod2F_kl[l][norm_tag] = mod2F_kl[l][norm_tag] + F[l].Rho2() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 		
          }       

       for(Int_t a=0; a<4; a++)
          { for(Int_t b=a+1; b<5; b++)
               { reinterKM1_kl[a][b][norm_tag] = reinterKM1_kl[a][b][norm_tag] + ((TComplex::Conjugate(F[a]))*F[b]).Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 	
		 iminterKM1_kl[a][b][norm_tag] = iminterKM1_kl[a][b][norm_tag] + ((TComplex::Conjugate(F[a]))*F[b]).Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 	
	       }
          }

//2xSum_i{Re(value0* x value1)}

       for(Int_t a=0; a<5; a++)
          { for(Int_t b=0; b<5; b++)
               { I[a][b] = u1j[b]*(*bb)[5]*(TComplex::Conjugate(L[a])) ;								
                 reinterI_kl[a][b][norm_tag] = reinterI_kl[a][b][norm_tag] + I[a][b].Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);
                 iminterI_kl[a][b][norm_tag] = iminterI_kl[a][b][norm_tag] + I[a][b].Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);
	       }
          }       

//----------------------------------- Part3 ---------------------------------------------------
//Sum_i{2xRe[(Sum_r{A_r})*A_KM] }
//Sum_i{2xRe[(Sum_r{a_r e^{iphi_r} Z_r})* (value0 + value1)]}
       
       for(Int_t r=0; r<num_reso; r++) {
       for(Int_t l=0; l<5; l++) {
       reinterZval0_kl[r][l][norm_tag] = reinterZval0_kl[r][l][norm_tag] + ((TComplex::Conjugate(Z[r]))*L[l]).Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);   
       iminterZval0_kl[r][l][norm_tag] = iminterZval0_kl[r][l][norm_tag] + ((TComplex::Conjugate(Z[r]))*L[l]).Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);   
      
       reinterZval1_kl[r][l][norm_tag] = reinterZval1_kl[r][l][norm_tag] + ((TComplex::Conjugate(Z[r]))*F[l]).Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);   
       iminterZval1_kl[r][l][norm_tag] = iminterZval1_kl[r][l][norm_tag] + ((TComplex::Conjugate(Z[r]))*F[l]).Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim);  


       reinterZval0_dcs_kl[r][l][norm_tag] = reinterZval0_dcs_kl[r][l][norm_tag] + ((TComplex::Conjugate(Z_dcs[r]))*L[l]).Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 
       iminterZval0_dcs_kl[r][l][norm_tag] = iminterZval0_dcs_kl[r][l][norm_tag] + ((TComplex::Conjugate(Z_dcs[r]))*L[l]).Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 
      
       reinterZval1_dcs_kl[r][l][norm_tag] = reinterZval1_dcs_kl[r][l][norm_tag] + ((TComplex::Conjugate(Z_dcs[r]))*F[l]).Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 
       iminterZval1_dcs_kl[r][l][norm_tag] = iminterZval1_dcs_kl[r][l][norm_tag] + ((TComplex::Conjugate(Z_dcs[r]))*F[l]).Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

       reinterval0Z_kl[r][l][norm_tag] = reinterval0Z_kl[r][l][norm_tag] + (Z[r]*TComplex::Conjugate(L[l])).Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 
       iminterval0Z_kl[r][l][norm_tag] = iminterval0Z_kl[r][l][norm_tag] + (Z[r]*TComplex::Conjugate(L[l])).Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 
      
       reinterval1Z_kl[r][l][norm_tag] = reinterval1Z_kl[r][l][norm_tag] + (Z[r]*TComplex::Conjugate(F[l])).Re() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 
       iminterval1Z_kl[r][l][norm_tag] = iminterval1Z_kl[r][l][norm_tag] + (Z[r]*TComplex::Conjugate(F[l])).Im() * Gamma_Factor_k0lpipi(pk0l,ppip,ppim); 

       }
       }


     }  

  std::cout<<"normnumber_kl: "<<normnumber_kl<<std::endl;

  std::cout<<"K0L done.."<<endl;

////////////////////////// Storing in ASCII files //////////////////////////////

fout << normnumber_ks << "\t" << taggednormnumber_ks[0] << "\t" << taggednormnumber_ks[1] << "\t" << taggednormnumber_ks[2] << endl;
fout << normnumber_kl << "\t" << taggednormnumber_kl[0] << "\t" << taggednormnumber_kl[1] << "\t" << taggednormnumber_kl[2] << endl;

for(int r=0; r<num_reso; r++) {
fout0 << probZ_ks[r][0] << "\t\t" << probZ_ks[r][1] << "\t\t" << probZ_ks[r][2] << endl;
}
for(int r=0; r<num_reso; r++) {
fout0 << probZ_dcs_ks[r][0] << "\t\t" << probZ_dcs_ks[r][1] << "\t\t" << probZ_dcs_ks[r][2] << endl;
}
for(int r=0; r<num_reso; r++) {
fout0 << probZ_kl[r][0] << "\t\t" << probZ_kl[r][1] << "\t\t" << probZ_kl[r][2] << endl;
}
for(int r=0; r<num_reso; r++) {
fout0 << probZ_dcs_kl[r][0] <<"\t\t" << probZ_dcs_kl[r][1] << "\t\t" << probZ_dcs_kl[r][2] << endl;
}

for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
fout1 << reinterZ_ks[r][rr][0] << "\t\t" << reinterZ_ks[r][rr][1] << "\t\t" << reinterZ_ks[r][rr][2] << "\t\t" << iminterZ_ks[r][rr][0] << "\t\t" << iminterZ_ks[r][rr][1] << "\t\t" << iminterZ_ks[r][rr][2] << endl;
}}
for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
fout1 << reinterZ_dcs_ks[r][rr][0] << "\t\t" << reinterZ_dcs_ks[r][rr][1] << "\t\t" << reinterZ_dcs_ks[r][rr][2] << "\t\t" << iminterZ_dcs_ks[r][rr][0] << "\t\t" << iminterZ_dcs_ks[r][rr][1] << "\t\t" << iminterZ_dcs_ks[r][rr][2] << endl;
}}
for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
fout1 << reinterZ_cfdcs_ks[r][rr][0] << "\t\t" << reinterZ_cfdcs_ks[r][rr][1] << "\t\t" << reinterZ_cfdcs_ks[r][rr][2] << "\t\t" << iminterZ_cfdcs_ks[r][rr][0] << "\t\t" << iminterZ_cfdcs_ks[r][rr][1] << "\t\t" << iminterZ_cfdcs_ks[r][rr][2] << endl;
}}
for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
fout1 << reinterZ_kl[r][rr][0] << "\t\t" << reinterZ_kl[r][rr][1] << "\t\t" << reinterZ_kl[r][rr][2] << "\t\t" << iminterZ_kl[r][rr][0] << "\t\t" << iminterZ_kl[r][rr][1] << "\t\t" << iminterZ_kl[r][rr][2] << endl;
}}
for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
fout1 << reinterZ_dcs_kl[r][rr][0] << "\t\t" << reinterZ_dcs_kl[r][rr][1] << "\t\t" << reinterZ_dcs_kl[r][rr][2] << "\t\t" << iminterZ_dcs_kl[r][rr][0] << "\t\t" << iminterZ_dcs_kl[r][rr][1] << "\t\t" << iminterZ_dcs_kl[r][rr][2] << endl;
}}
for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
fout1 << reinterZ_cfdcs_kl[r][rr][0] << "\t\t" << reinterZ_cfdcs_kl[r][rr][1] << "\t\t" << reinterZ_cfdcs_kl[r][rr][2] << "\t\t" << iminterZ_cfdcs_kl[r][rr][0] << "\t\t" << iminterZ_cfdcs_kl[r][rr][1] << "\t\t" << iminterZ_cfdcs_kl[r][rr][2] << endl;
}}

for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
fout2 << reinterZval0_ks[r][l][0] << "\t\t" << reinterZval0_ks[r][l][1] << "\t\t" << reinterZval0_ks[r][l][2] << "\t\t" << iminterZval0_ks[r][l][0] << "\t\t" << iminterZval0_ks[r][l][1] << "\t\t" << iminterZval0_ks[r][l][2] << "\t\t" << reinterZval1_ks[r][l][0] << "\t\t" << reinterZval1_ks[r][l][1] << "\t\t" << reinterZval1_ks[r][l][2] << "\t\t" << iminterZval1_ks[r][l][0] << "\t\t" << iminterZval1_ks[r][l][1] << "\t\t" << iminterZval1_ks[r][l][2] << endl;
}}
for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
fout2 << reinterZval0_dcs_ks[r][l][0] << "\t\t" << reinterZval0_dcs_ks[r][l][1] << "\t\t" << reinterZval0_dcs_ks[r][l][2] << "\t\t" << iminterZval0_dcs_ks[r][l][0] << "\t\t" << iminterZval0_dcs_ks[r][l][1] << "\t\t" << iminterZval0_dcs_ks[r][l][2] << "\t\t" << reinterZval1_dcs_ks[r][l][0] << "\t\t" << reinterZval1_dcs_ks[r][l][1] << "\t\t" << reinterZval1_dcs_ks[r][l][2] << "\t\t" << iminterZval1_dcs_ks[r][l][0] << "\t\t" << iminterZval1_dcs_ks[r][l][1] << "\t\t" << iminterZval1_dcs_ks[r][l][2] << endl;
}}
for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
fout2 << reinterval0Z_ks[r][l][0] << "\t\t" << reinterval0Z_ks[r][l][1] << "\t\t" << reinterval0Z_ks[r][l][2] << "\t\t" << iminterval0Z_ks[r][l][0] << "\t\t" << iminterval0Z_ks[r][l][1] << "\t\t" << iminterval0Z_ks[r][l][2] << "\t\t" << reinterval1Z_ks[r][l][0] << "\t\t" << reinterval1Z_ks[r][l][1] << "\t\t" << reinterval1Z_ks[r][l][2] << "\t\t" << iminterval1Z_ks[r][l][0] << "\t\t" << iminterval1Z_ks[r][l][1] << "\t\t" << iminterval1Z_ks[r][l][2] << endl;
}}
for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
fout2 << reinterZval0_kl[r][l][0] << "\t\t" << reinterZval0_kl[r][l][1] << "\t\t" << reinterZval0_kl[r][l][2] << "\t\t" << iminterZval0_kl[r][l][0] << "\t\t" << iminterZval0_kl[r][l][1] << "\t\t" << iminterZval0_kl[r][l][2] << "\t\t" << reinterZval1_kl[r][l][0] << "\t\t" << reinterZval1_kl[r][l][1] << "\t\t" << reinterZval1_kl[r][l][2] << "\t\t" << iminterZval1_kl[r][l][0] << "\t\t" << iminterZval1_kl[r][l][1] << "\t\t" << iminterZval1_kl[r][l][2] << endl;
}}
for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
fout2 << reinterZval0_dcs_kl[r][l][0] << "\t\t" << reinterZval0_dcs_kl[r][l][1] << "\t\t" << reinterZval0_dcs_kl[r][l][2] << "\t\t" << iminterZval0_dcs_kl[r][l][0] << "\t\t" << iminterZval0_dcs_kl[r][l][1] << "\t\t" << iminterZval0_dcs_kl[r][l][2] << "\t\t" << reinterZval1_dcs_kl[r][l][0] << "\t\t" << reinterZval1_dcs_kl[r][l][1] << "\t\t" << reinterZval1_dcs_kl[r][l][2] << "\t\t" << iminterZval1_dcs_kl[r][l][0] << "\t\t" << iminterZval1_dcs_kl[r][l][1] << "\t\t" << iminterZval1_dcs_kl[r][l][2] << endl;
}}
for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
fout2 << reinterval0Z_kl[r][l][0] << "\t\t" << reinterval0Z_kl[r][l][1] << "\t\t" << reinterval0Z_kl[r][l][2] << "\t\t" << iminterval0Z_kl[r][l][0] << "\t\t" << iminterval0Z_kl[r][l][1] << "\t\t" << iminterval0Z_kl[r][l][2] << "\t\t" << reinterval1Z_kl[r][l][0] << "\t\t" << reinterval1Z_kl[r][l][1] << "\t\t" << reinterval1Z_kl[r][l][2] << "\t\t" << iminterval1Z_kl[r][l][0] << "\t\t" << iminterval1Z_kl[r][l][1] << "\t\t" << iminterval1Z_kl[r][l][2] << endl;
}}

for(int l=0; l<5; l++) {
fout3 << mod2L_ks[l][0] << "\t\t" << mod2L_ks[l][1] << "\t\t" << mod2L_ks[l][2] << "\t\t" << mod2F_ks[l][0] << "\t\t" << mod2F_ks[l][1] << "\t\t" << mod2F_ks[l][2] << endl;
}
for(int l=0; l<5; l++) {
fout3 << mod2L_kl[l][0] << "\t\t" << mod2L_kl[l][1] << "\t\t" << mod2L_kl[l][2] << "\t\t" << mod2F_kl[l][0] << "\t\t" << mod2F_kl[l][1] << "\t\t" << mod2F_kl[l][2] << endl;
}

for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
fout4 << reinterKM0_ks[l][k][0] << "\t\t" << reinterKM0_ks[l][k][1] << "\t\t" << reinterKM0_ks[l][k][2] << "\t\t" << iminterKM0_ks[l][k][0] << "\t\t" << iminterKM0_ks[l][k][1] << "\t\t" << iminterKM0_ks[l][k][2] << endl;
}}
for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
fout4 << reinterKM1_ks[l][k][0] << "\t\t" << reinterKM1_ks[l][k][1] << "\t\t" << reinterKM1_ks[l][k][2] << "\t\t" << iminterKM1_ks[l][k][0] << "\t\t" << iminterKM1_ks[l][k][1] << "\t\t" << iminterKM1_ks[l][k][2] << endl;
}}
for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
fout4 << reinterI_ks[l][k][0] << "\t\t" << reinterI_ks[l][k][1] << "\t\t" << reinterI_ks[l][k][2] << "\t\t" << iminterI_ks[l][k][0] << "\t\t" << iminterI_ks[l][k][1] << "\t\t" << iminterI_ks[l][k][2] << endl;
}}
for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
fout4 << reinterKM0_kl[l][k][0] << "\t\t" << reinterKM0_kl[l][k][1] << "\t\t" << reinterKM0_kl[l][k][2] << "\t\t" << iminterKM0_kl[l][k][0] << "\t\t" << iminterKM0_kl[l][k][1] << "\t\t" << iminterKM0_kl[l][k][2] << endl;
}}
for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
fout4 << reinterKM1_kl[l][k][0] << "\t\t" << reinterKM1_kl[l][k][1] << "\t\t" << reinterKM1_kl[l][k][2] << "\t\t" << iminterKM1_kl[l][k][0] << "\t\t" << iminterKM1_kl[l][k][1] << "\t\t" << iminterKM1_kl[l][k][2] << endl;
}}
for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
fout4 << reinterI_kl[l][k][0] << "\t\t" << reinterI_kl[l][k][1] << "\t\t" << reinterI_kl[l][k][2] << "\t\t" << iminterI_kl[l][k][0] << "\t\t" << iminterI_kl[l][k][1] << "\t\t" << iminterI_kl[l][k][2] << endl;
}}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
  std::cout<<"Norm sample maps and vectors ready..."<<std::endl;

  time_t my_endtime = time(NULL);
  std::cout<<"Start time: "<<ctime(&my_time)<<endl;
  std::cout<<"End time: "<<ctime(&my_endtime)<<endl;

  return 0;

}


