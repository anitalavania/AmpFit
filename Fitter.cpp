#include <TMinuit.h>
#include <math.h>
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

map<string, double> massmap ;
map<string, double> widthmap ;
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

double k0sbkg_frac[3] ; //= 0.044493358 ; //Fixed K0S background fraction from Inclusive MC in the signal region of M_miss^2 distribution.
double NPbkg_frac[3] ; //= 0.064671507 ;
double chisq ;

double sb_pdf[110][110][3];
double SB_normfac[3];

//double rd[3], deltad[3], Rf[3] ;

double rd[3] = {0.05867, 0.0550, 0.0441} ;
double deltad[3] = {190.*(TMath::Pi()/180.), 161.*(TMath::Pi()/180.), 196.*(TMath::Pi()/180.)} ;
double Rf[3] = {1., .44, .79} ;


double r_kl = 1.0;
//double tan2thetaC=tan(13.04*TMath::Pi()/180.0)*tan(13.04*TMath::Pi()/180.0) ;
double tan2thetaC = (0.22650*0.22650)/(1.-(0.22650*0.22650)) ;	//sin(theta_C) = 0.22650 +/- 0.00048
double eitheta = 1.0;
double theta = 0. ;
//TComplex CP_mult;
//CP_mult = 1. - 2.*tan2thetaC*TComplex(r_kl*cos(theta), r_kl*sin(theta)) ;

double prob_dcs1=0., prob_dcs2=0., prob_dcs3=0., prob_dcs4=0.;

int fixflag_ar[50], fixflag_ph[50];
double effi[200][200][3];
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
int TAG_ks[50000], TAG_kl[50000], NPAR;

Double_t CMplus2_ks[200000], CMminus2_ks[200000];
//std::map<string, Double_t> CMulfacbw_re_ks[200000];
//std::map<string, Double_t> CMulfacbw_im_ks[200000];
//std::map<string, Double_t> CLASS_contribution_re_ks[200000];
//std::map<string, Double_t> CLASS_contribution_im_ks[200000];

std::vector<TComplex> CMulfacbw_ks[20000];
std::vector<TComplex> CMulfacbw_dcs_ks[20000];
std::vector<TComplex> CLASS_contribution_ks[20000];

std::vector <double> CU1j_re_ks ;
std::vector <double> CU1j_im_ks ; 
std::vector <double> CBB_ks[40000] ;

TComplex A0_ks[50000];
TComplex A1_ks[50000];
TComplex A2_ks[50000];
TComplex A3_ks[50000];
TComplex A4_ks[50000];

TComplex U1jCal_ks[50000][5];

Double_t CMplus2_kl[50000], CMminus2_kl[50000], Cpk0l[45000], Cpk0s[20000];
Double_t Cppip_kl[45000], Cppim_kl[45000], Cppip_ks[20000], Cppim_ks[20000] ;
//std::map<string, Double_t> CMulfacbw_re_kl[50000];
//std::map<string, Double_t> CMulfacbw_im_kl[50000];
//std::map<string, Double_t> CLASS_contribution_re_kl[50000];
//std::map<string, Double_t> CLASS_contribution_im_kl[50000];

std::vector<TComplex> CMulfacbw_kl[40000];
std::vector<TComplex> CMulfacbw_dcs_kl[40000];
std::vector<TComplex> CLASS_contribution_kl[40000];

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
TH2D *heffi = new TH2D("heffi", "heffi", 50, 0.3, 3.0, 50, 0.3, 3.0);
TH1F *histu1j = new TH1F("histu1j", "", 100, -2.0, 2.0);

double normprobcal_ks[3], normprobcal_kl[3], normprobcal_ksbkg[3];

TH1D *pion = new TH1D("pion gamma","",10,0.,1.) ;
Double_t pi_track_gamma[10], pi_pid_gamma[10] ;
double out[4] ;
Double_t piplus_trk_gamma, piminus_trk_gamma, piplus_pid_gamma, piminus_pid_gamma, gamma_k0l, gamma_k0s;

Double_t mc_yields[3][3] ;		//mc_yields[tag][event_type]
Double_t mc_bkgfrac[3][2], mc_bkgfrac_err[3][2];
double mc_np_frac[3], mc_np_frac_err[3] ;

// Gamma Fit parameters

//K0L
Double_t c0_k0l = 25.9620, c1_k0l = -10.8760, c2_k0l = 4.0073;
//pi tracking
Double_t c0_pi_trk = 0.0182, c1_pi_trk = -0.1066, c2_pi_trk = 0.1888, c3_pi_trk = -0.1006;
//pi PID
Double_t c0_pi_pid = -0.003128, c1_pi_pid = 0.020158, c2_pi_pid = -0.034488;

////////////// K0Lpipi effi syst. ////////////////

double Gamma_Factor_k0lpipi(double xx_k0l, double xx_piplus, double xx_piminus) {

//K0L rec.
Double_t gamma_k0l = 1. + (c0_k0l*TMath::Exp(c1_k0l*xx_k0l) + c2_k0l)/100. ;

//pi+/- track
piplus_trk_gamma = 1. + (c0_pi_trk + c1_pi_trk*xx_piplus +  c2_pi_trk*xx_piplus*xx_piplus +  c3_pi_trk*xx_piplus*xx_piplus*xx_piplus) ;
piminus_trk_gamma = 1. + (c0_pi_trk + c1_pi_trk*xx_piminus +  c2_pi_trk*xx_piminus*xx_piminus +  c3_pi_trk*xx_piminus*xx_piminus*xx_piminus) ;

//pi+/- PID
piplus_pid_gamma = 1. + (c0_pi_pid + c1_pi_pid * xx_piplus + c2_pi_pid * xx_piplus * xx_piplus) ;
piminus_pid_gamma = 1. + (c0_pi_pid + c1_pi_pid * xx_piminus + c2_pi_pid * xx_piminus * xx_piminus) ;

return gamma_k0l*piplus_trk_gamma*piplus_pid_gamma*piminus_trk_gamma*piminus_pid_gamma ;
}


///////////////////////////
///////// K0Spipi /////////
///////////////////////////

double prob_BWLASS_ks(TComplex *parameter, int tag)
{ 
  double prob=0.;
  for(int r=0; r<num_reso; r++)
     { 
       prob = prob + (parameter[r].Rho2())*probZ_ks[r][tag];
     }
  for(int a=0; a<(num_reso-1); a++)
     { for(int b=a+1; b<num_reso; b++)
          { prob = prob + 2.*(((TComplex::Conjugate(parameter[a])*parameter[b]).Re())*reinterZ_ks[a][b][tag] - ((TComplex::Conjugate(parameter[a])*parameter[b]).Im())*iminterZ_ks[a][b][tag]) ;
          }
     }

  return prob;
}

double prob_BWLASS_DCS_ks(TComplex *parameter, int tag)
{ 
  double prob=0.;
  for(int r=0; r<num_reso; r++)
     { prob = prob + (parameter[r].Rho2())*probZ_dcs_ks[r][tag];
     }
  for(int a=0; a<(num_reso-1); a++)
     { for(int b=a+1; b<num_reso; b++)
          { prob = prob + 2.*(((TComplex::Conjugate(parameter[a])*parameter[b]).Re())*reinterZ_dcs_ks[a][b][tag] - ((TComplex::Conjugate(parameter[a])*parameter[b]).Im())*iminterZ_dcs_ks[a][b][tag]) ;
          }
     }
  return prob;
}


double prob_KM_ks(TComplex *parameter, TComplex *beta, TComplex *fprod, int tag)
{
  double prob=0.;
//    |value0|^2
  for(Int_t k=0; k<5; k++) {
  prob += beta[k].Rho2() * mod2L_ks[k][tag];  
  }

  for(Int_t a=0; a<4; a++) {
  for(Int_t b=a+1; b<5; b++) {
  prob += 2.*(((TComplex::Conjugate(beta[a])*beta[b]).Re())*reinterKM0_ks[a][b][tag] - ((TComplex::Conjugate(beta[a])*beta[b]).Im())*iminterKM0_ks[a][b][tag]) ;
  }
  }
//    |value1|^2
  for(Int_t k=0; k<5; k++)
     { prob += fprod[k].Rho2() * mod2F_ks[k][tag];  
     }

  for(Int_t a=0; a<4; a++)
     { for(Int_t b=a+1; b<5; b++)
          { prob += 2.*(((TComplex::Conjugate(fprod[a])*fprod[b]).Re())*reinterKM1_ks[a][b][tag] - ((TComplex::Conjugate(fprod[a])*fprod[b]).Im())*iminterKM1_ks[a][b][tag]) ;
          }
     }
//    value0-value1 interference
  for(Int_t a=0; a<5; a++)
     { for(Int_t b=0; b<5; b++)
          { prob += 2.*(((TComplex::Conjugate(beta[a])*fprod[b]).Re())*reinterI_ks[a][b][tag] - ((TComplex::Conjugate(beta[a])*fprod[b]).Im())*iminterI_ks[a][b][tag]) ; 
          }
     }

  return prob;					
}

double prob_interference_ks(TComplex *parameter, TComplex *beta, TComplex *fprod, int tag)
{ 
  double prob=0.0;
  for(Int_t r=0; r<num_reso; r++) {
  for(Int_t l=0; l<5; l++) {
  prob += 2.*(((TComplex::Conjugate(parameter[r])*beta[l]).Re())*reinterZval0_ks[r][l][tag] - ((TComplex::Conjugate(parameter[r])*beta[l]).Im())*iminterZval0_ks[r][l][tag]) ; 
  prob += 2.*(((TComplex::Conjugate(parameter[r])*fprod[l]).Re())*reinterZval1_ks[r][l][tag] - ((TComplex::Conjugate(parameter[r])*fprod[l]).Im())*iminterZval1_ks[r][l][tag]) ;
  }
  }
  return prob;
}

double prob_interference_DCS_ks(TComplex *parameter, TComplex *beta, TComplex *fprod, int tag)
{ 
  double prob=0.0;
  for(Int_t r=0; r<num_reso; r++) {
  for(Int_t l=0; l<5; l++) {
  prob += 2.*(((TComplex::Conjugate(parameter[r])*beta[l]).Re())*reinterZval0_dcs_ks[r][l][tag] - ((TComplex::Conjugate(parameter[r])*beta[l]).Im())*iminterZval0_dcs_ks[r][l][tag]) ; 
  prob += 2.*(((TComplex::Conjugate(parameter[r])*fprod[l]).Re())*reinterZval1_dcs_ks[r][l][tag] - ((TComplex::Conjugate(parameter[r])*fprod[l]).Im())*iminterZval1_dcs_ks[r][l][tag]) ;
  }
  }
  return prob;
}

// Interference between CF and DCS (K0Spipi)

double CFDCS_BWBW_ks(TComplex *parameter, int tag, double *deltad)
{
double prob=0.;

for(int a=0; a<num_reso; a++) {
for(int b=0; b<num_reso; b++) {
prob += ( (parameter[a]*TComplex::Conjugate(parameter[b])*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Re()*reinterZ_cfdcs_ks[a][b][tag] - (parameter[a]*TComplex::Conjugate(parameter[b])*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Im()*iminterZ_cfdcs_ks[a][b][tag] ) ;
}
}

return prob;
}

double CFDCS_BWKM_ks(TComplex *parameter, TComplex *beta, TComplex *fprod, int tag, double *deltad)
{
double prob=0.;

for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
prob += (parameter[r]*TComplex::Conjugate(beta[l])*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Re()*reinterval0Z_ks[r][l][tag] - (parameter[r]*TComplex::Conjugate(beta[l])*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Im()*iminterval0Z_ks[r][l][tag] ;

prob += (parameter[r]*TComplex::Conjugate(fprod[l])*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Re()*reinterval1Z_ks[r][l][tag] - (parameter[r]*TComplex::Conjugate(fprod[l])*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Im()*iminterval1Z_ks[r][l][tag] ;
}
}
return prob;
}

double CFDCS_KMBW_ks(TComplex *parameter, TComplex *beta, TComplex *fprod, int tag, double *deltad)
{
double prob=0.;

for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
prob += (TComplex::Conjugate(parameter[r])*beta[l]*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Re()*reinterZval0_dcs_ks[r][l][tag] - (TComplex::Conjugate(parameter[r])*beta[l]*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Im()*iminterZval0_dcs_ks[r][l][tag] ;

prob += (TComplex::Conjugate(parameter[r])*fprod[l]*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Re()*reinterZval1_dcs_ks[r][l][tag] - (TComplex::Conjugate(parameter[r])*fprod[l]*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Im()*iminterZval1_dcs_ks[r][l][tag] ;
}
}
return prob;
}

////////////////////////////
///////// K0Lpipi /////////
//////////////////////////

double prob_BWLASS_kl(TComplex *parameter, int tag)
{ 
  double prob=0.;
  for(int r=0; r<num_reso; r++)
     { prob = prob + (parameter[r].Rho2())*probZ_kl[r][tag];
     }
  for(int a=0; a<(num_reso-1); a++)
     { for(int b=a+1; b<num_reso; b++)
          { prob = prob + 2.*(((TComplex::Conjugate(parameter[a])*parameter[b]).Re())*reinterZ_kl[a][b][tag] - ((TComplex::Conjugate(parameter[a])*parameter[b]).Im())*iminterZ_kl[a][b][tag]) ;
          }
     }
  return prob;
}

double prob_BWLASS_DCS_kl(TComplex *parameter, int tag)
{ 
  double prob=0.;
  for(int r=0; r<num_reso; r++)
     { prob = prob + (parameter[r].Rho2())*probZ_dcs_kl[r][tag];
     }
  for(int a=0; a<(num_reso-1); a++)
     { for(int b=a+1; b<num_reso; b++)
          { prob = prob + 2.*(((TComplex::Conjugate(parameter[a])*parameter[b]).Re())*reinterZ_dcs_kl[a][b][tag] - ((TComplex::Conjugate(parameter[a])*parameter[b]).Im())*iminterZ_dcs_kl[a][b][tag]) ;
          }
     }
  return prob;
}


double prob_KM_kl(TComplex *parameter, TComplex *beta, TComplex *fprod, int tag)
{
  double prob=0.;
//    |value0|^2
  for(Int_t k=0; k<5; k++) {
  prob += beta[k].Rho2() * mod2L_kl[k][tag];  
  }

  for(Int_t a=0; a<4; a++) {
  for(Int_t b=a+1; b<5; b++) {
  prob += 2.*(((TComplex::Conjugate(beta[a])*beta[b]).Re())*reinterKM0_kl[a][b][tag] - ((TComplex::Conjugate(beta[a])*beta[b]).Im())*iminterKM0_kl[a][b][tag]) ;
  }
  }
//    |value1|^2
  for(Int_t k=0; k<5; k++)
     { prob += fprod[k].Rho2() * mod2F_kl[k][tag];  
     }

  for(Int_t a=0; a<4; a++)
     { for(Int_t b=a+1; b<5; b++)
          { prob += 2.*(((TComplex::Conjugate(fprod[a])*fprod[b]).Re())*reinterKM1_kl[a][b][tag] - ((TComplex::Conjugate(fprod[a])*fprod[b]).Im())*iminterKM1_kl[a][b][tag]) ;
          }
     }
//    value0-value1 interference
  for(Int_t a=0; a<5; a++)
     { for(Int_t b=0; b<5; b++)
          { prob += 2.*(((TComplex::Conjugate(beta[a])*fprod[b]).Re())*reinterI_kl[a][b][tag] - ((TComplex::Conjugate(beta[a])*fprod[b]).Im())*iminterI_kl[a][b][tag]) ; 
          }
     }

  return prob;					
}

double prob_interference_kl(TComplex *parameter, TComplex *beta, TComplex *fprod, int tag, TComplex cpmult)
{ 
  double prob=0.0;
  for(Int_t r=0; r<num_reso; r++) {
  for(Int_t l=0; l<5; l++) {
  prob += 2.*(((TComplex::Conjugate(parameter[r])*cpmult*beta[l]).Re())*reinterZval0_kl[r][l][tag] - ((TComplex::Conjugate(parameter[r])*cpmult*beta[l]).Im())*iminterZval0_kl[r][l][tag]) ; 
  prob += 2.*(((TComplex::Conjugate(parameter[r])*cpmult*fprod[l]).Re())*reinterZval1_kl[r][l][tag] - ((TComplex::Conjugate(parameter[r])*cpmult*fprod[l]).Im())*iminterZval1_kl[r][l][tag]) ;
  }
  }
  return prob;
}

double prob_interference_DCS_kl(TComplex *parameter, TComplex *beta, TComplex *fprod, int tag, TComplex cpmult)
{ 
  double prob=0.0;
  for(Int_t r=0; r<num_reso; r++) {
  for(Int_t l=0; l<5; l++) {
  prob += 2.*(((TComplex::Conjugate(parameter[r])*cpmult*beta[l]).Re())*reinterZval0_dcs_kl[r][l][tag] - ((TComplex::Conjugate(parameter[r])*cpmult*beta[l]).Im())*iminterZval0_dcs_kl[r][l][tag]) ; 
  prob += 2.*(((TComplex::Conjugate(parameter[r])*cpmult*fprod[l]).Re())*reinterZval1_dcs_kl[r][l][tag] - ((TComplex::Conjugate(parameter[r])*cpmult*fprod[l]).Im())*iminterZval1_dcs_kl[r][l][tag]) ;
  }
  }
  return prob;
}

// Interference between CF and DCS (K0Lpipi)

double CFDCS_BWBW_kl(TComplex *parameter, int tag, double *deltad)
{
double prob=0.;

for(int a=0; a<num_reso; a++) {
for(int b=0; b<num_reso; b++) {
prob += ( (parameter[a]*TComplex::Conjugate(parameter[b])*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Re()*reinterZ_cfdcs_kl[a][b][tag] - (parameter[a]*TComplex::Conjugate(parameter[b])*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Im()*iminterZ_cfdcs_kl[a][b][tag] ) ;
}
}

return prob;
}

double CFDCS_BWKM_kl(TComplex *parameter, TComplex *beta, TComplex *fprod, int tag, double *deltad, TComplex cpmult)
{
double prob=0.;

for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
prob += (parameter[r]*TComplex::Conjugate(beta[l]*cpmult)*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Re()*reinterval0Z_kl[r][l][tag] - (parameter[r]*TComplex::Conjugate(beta[l]*cpmult)*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Im()*iminterval0Z_kl[r][l][tag] ;

prob += (parameter[r]*TComplex::Conjugate(fprod[l]*cpmult)*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Re()*reinterval1Z_kl[r][l][tag] - (parameter[r]*TComplex::Conjugate(fprod[l]*cpmult)*TComplex(cos(deltad[tag]), sin(deltad[tag]))).Im()*iminterval1Z_kl[r][l][tag] ;
}
}
return prob;
}

double CFDCS_KMBW_kl(TComplex *parameter, TComplex *beta, TComplex *fprod, int tag, double *deltad, TComplex cpmult)
{
double prob=0.;

for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
prob += (TComplex::Conjugate(parameter[r])*beta[l]*TComplex(cos(deltad[tag]), sin(deltad[tag]))*cpmult).Re()*reinterZval0_dcs_kl[r][l][tag] - (TComplex::Conjugate(parameter[r])*beta[l]*TComplex(cos(deltad[tag]), sin(deltad[tag]))*cpmult).Im()*iminterZval0_dcs_kl[r][l][tag] ;

prob += (TComplex::Conjugate(parameter[r])*fprod[l]*TComplex(cos(deltad[tag]), sin(deltad[tag]))*cpmult).Re()*reinterZval1_dcs_kl[r][l][tag] - (TComplex::Conjugate(parameter[r])*fprod[l]*TComplex(cos(deltad[tag]), sin(deltad[tag]))*cpmult).Im()*iminterZval1_dcs_kl[r][l][tag] ;
}
}
return prob;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

double likelihood(const double *par)			
{  

  double effiampl2 = 0.0, denomin=0.0, test_xx[2];
  double overallamp = 1.0;
  double overallphase = 0.0;
  TComplex param[20], param_dcs[20];
  TComplex param_kl[20], param_dcs_kl[20];
  int i=0, j=0;

  for(Int_t r=0; r<num_reso; r++) {
  if(fixflag_ar[r]==0 && fixflag_ph[r]==0) {
  param[r] = TComplex(par[i]*cos(par[i+1]), par[i]*sin(par[i+1])) ;
  i+=2;				
  }
  //else if(fixflag_ar[r]==0 && fixflag_ph[r]==2) {
  //param[r] = TComplex(par[i]*cos(fixedpar[j]), par[i]*sin(fixedpar[j])) ;  
  //i++; j++;						
  //}
  //else if(fixflag_ar[r]==2 && fixflag_ph[r]==0) {
  //param[r] = TComplex(fixedpar[j]*cos(par[i]), fixedpar[j]*sin(par[i])) ;  
  //i++; j++;
  //}
  else if(fixflag_ar[r]==2 && fixflag_ph[r]==2) {
  param[r] = TComplex(fixedpar[j]*cos(fixedpar[j+1]), fixedpar[j]*sin(fixedpar[j+1])) ;
  j+=2;
  }
  }

  //for(Int_t rr=0; rr<num_reso; rr++) cout<<"param["<<rr<<"]: "<<param[rr].Rho()<<"\t"<<param[rr].Theta()*(180./TMath::Pi())<<endl;
 
/////////// S-wave parameters ///////////

  TComplex _beta[5], _f[5], _beta_dcs[5], _f_dcs[5];

  int num_k = num_reso ;
  for(int pole=0; pole<4; pole++)
        {
                if(fixflag_ar[num_k]==0 && fixflag_ph[num_k]==0)
                {
                        _beta[pole] = TComplex(par[i]*cos(par[i+1]), par[i]*sin(par[i+1])) ;
                        i+=2;
                }
                else if(fixflag_ar[num_k]==0 && fixflag_ph[num_k]==2)
                {
                        _beta[pole] = TComplex(par[i]*cos(fixedpar[j+1]), par[i]*sin(fixedpar[j+1])) ;
                        i++; j++;
                }
                else if(fixflag_ar[num_k]==2 && fixflag_ph[num_k]==2)
                {
                        _beta[pole] = TComplex(fixedpar[j]*cos(fixedpar[j+1]), fixedpar[j]*sin(fixedpar[j+1])) ;
                        j+=2 ;
                }
                num_k++;
        }

  for(int pole=0; pole<4; pole++)
        {
                if(fixflag_ar[num_k]==0 && fixflag_ph[num_k]==0)
                {
                        _f[pole] = TComplex(par[i]*cos(par[i+1]), par[i]*sin(par[i+1])) ;
                        i+=2;
                }
                else if(fixflag_ar[num_k]==0 && fixflag_ph[num_k]==2)
                {
                        _f[pole] = TComplex(par[i]*cos(fixedpar[j+1]), par[i]*sin(fixedpar[j+1])) ;
                        i++; j++;
                }
                else if(fixflag_ar[num_k]==2 && fixflag_ph[num_k]==2)
                {
                        _f[pole] = TComplex(fixedpar[j]*cos(fixedpar[j+1]), fixedpar[j]*sin(fixedpar[j+1])) ;
                        j+=2 ;
                }
                num_k++;

        }


  _beta[4] = TComplex(0.0, 0.0);
  _f[4] = TComplex(0.0, 0.0);


  for(int t=0; t<3; t++) 
	{
  		normprobcal_ks[t] =  prob_BWLASS_ks(param, t) + prob_KM_ks(param, _beta, _f, t) + prob_interference_ks(param, _beta, _f, t) ;

  		normprobcal_ks[t] += rd[t]*rd[t]*(prob_BWLASS_DCS_ks(param, t) + prob_KM_ks(param, _beta, _f, t) + prob_interference_DCS_ks(param, _beta, _f, t)) ; 

  		normprobcal_ks[t] += -2.*rd[t]*Rf[t]*(CFDCS_BWBW_ks(param, t, deltad) + CFDCS_BWKM_ks(param, _beta, _f, t, deltad) + CFDCS_KMBW_ks(param, _beta, _f, t, deltad) + cos(deltad[t])*prob_KM_ks(param, _beta, _f, t)) ;

  		normprobcal_ks[t] = normprobcal_ks[t]/taggednormnumber_ks[t] ;		
  	}


  TComplex CP_mult[5];
  for(int r=0; r<5; r++) 
	{
  		//cout<<"CP mult par: "<<par[i]<<"\t"<<par[i+1]*(180./TMath::Pi())<<endl;
  		CP_mult[r] = (TComplex(1.,0.) - 2.*tan2thetaC*TComplex(par[i]*cos(par[i+1]), par[i]*sin(par[i+1]))) ;
  		i+=2;
  	}


  int w=0;
  for(Int_t r=0; r<num_reso; r++) {
  if(state[r]=="CF")       { param_kl[r] = 1.*param[r] ; } 
  else if(state[r]=="DCS") { param_kl[r] = -1.*param[r] ; } 
  else if(state[r]=="CP")  { 
  param_kl[r] = CP_mult[w]*param[r] ; 
  w++;
  } 
  }

  //NPbkg_frac = par[40] ;
  //k0sbkg_frac = par[41] ;

  NPbkg_frac[0] = par[NPAR-3] ;
  NPbkg_frac[1] = par[NPAR-2] ;
  NPbkg_frac[2] = par[NPAR-1] ;
  
  //k0sbkg_frac[0] = par[43] ;
  //k0sbkg_frac[1] = par[44] ;
  //k0sbkg_frac[2] = par[45] ;

  for(int t=0; t<3; t++) {
  normprobcal_kl[t] =  prob_BWLASS_kl(param_kl, t) + prob_KM_kl(param_kl, _beta, _f, t)*(CP_mult[4].Rho2()) + prob_interference_kl(param_kl, _beta, _f, t, CP_mult[4]) ;

  normprobcal_kl[t] += rd[t]*rd[t]*(prob_BWLASS_DCS_kl(param_kl, t) + prob_KM_kl(param_kl, _beta, _f, t)*(CP_mult[4].Rho2()) + prob_interference_DCS_kl(param_kl, _beta, _f, t, CP_mult[4])) ; 

  normprobcal_kl[t] += 2.*rd[t]*Rf[t]*(CFDCS_BWBW_kl(param_kl, t, deltad) + CFDCS_BWKM_kl(param_kl, _beta, _f, t, deltad, CP_mult[4]) + CFDCS_KMBW_kl(param_kl, _beta, _f, t, deltad, CP_mult[4]) + cos(deltad[t])*(CP_mult[4].Rho2())*prob_KM_kl(param_kl, _beta, _f, t)) ;

  normprobcal_kl[t] = normprobcal_kl[t]/taggednormnumber_kl[t] ;		
  }

  TComplex dummy(1.0, 0.0);

  for(int t=0; t<3; t++) {
  normprobcal_ksbkg[t] =  prob_BWLASS_kl(param, t) + prob_KM_kl(param, _beta, _f, t) + prob_interference_kl(param, _beta, _f, t, dummy) ;

  normprobcal_ksbkg[t] += rd[t]*rd[t]*(prob_BWLASS_DCS_kl(param, t) + prob_KM_kl(param, _beta, _f, t) + prob_interference_DCS_kl(param, _beta, _f, t, dummy)) ; 

  normprobcal_ksbkg[t] += -2.*rd[t]*Rf[t]*(CFDCS_BWBW_kl(param, t, deltad) + CFDCS_BWKM_kl(param, _beta, _f, t, deltad, dummy) + CFDCS_KMBW_kl(param, _beta, _f, t, deltad, dummy) + cos(deltad[t])*prob_KM_kl(param, _beta, _f, t)) ;

  normprobcal_ksbkg[t] = normprobcal_ksbkg[t]/taggednormnumber_kl[t] ;		
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Double_t ampl_sqr=0.0, DCS_contri=0., total_prob=0., ksbkg_prob=0. ;
  Double_t LL1=0., LogLL=0., LogLL_kl=0. ;
  Double_t pdf=0.;
  Double_t Interference=0.;

////////////////////////////////////////////////////////////////
////////////////////////// K0Spipi /////////////////////////////
////////////////////////////////////////////////////////////////

  for(Int_t i=0; i<datanumber_ks; i++)	
     { 

        TComplex acal(0.0, 0.0);
        TComplex acal_dcs(0.0, 0.0);
        total_prob=0.;

        value0 = TComplex(0.,0.);
        value1 = TComplex(0.,0.);

        value0 = _beta[0]*A0_ks[i]+_beta[1]*A1_ks[i]+_beta[2]*A2_ks[i]+_beta[3]*A3_ks[i]+_beta[4]*A4_ks[i]; 
        value1 = (_f[0]*U1jCal_ks[i][0]+_f[1]*U1jCal_ks[i][1]+_f[2]*U1jCal_ks[i][2]+_f[3]*U1jCal_ks[i][3]+_f[4]*U1jCal_ks[i][4]) * CBB_ks[i][5];
	
	for(int rr=0 ; rr<num_reso ; rr++)
	{
		acal += param[rr]*CMulfacbw_ks[i][rr] ;
		acal_dcs += param[rr]*CMulfacbw_dcs_ks[i][rr] ;
	}

	acal += (value0 + value1);
        acal_dcs += (value0 + value1);

        /////////////////////////////////////////////////////
        
        Interference = TComplex(1.0*cos(deltad[TAG_ks[i]]), -1.0*sin(deltad[TAG_ks[i]]))*acal_dcs*TComplex::Conjugate(acal) + TComplex(1.0*cos(deltad[TAG_ks[i]]), 1.0*sin(deltad[TAG_ks[i]]))*acal*TComplex::Conjugate(acal_dcs)   ;
        total_prob = acal.Rho2() + rd[TAG_ks[i]]*rd[TAG_ks[i]]*acal_dcs.Rho2() - rd[TAG_ks[i]]*Rf[TAG_ks[i]]*Interference ;

        ////////////////////////////////////////////////////

        if(total_prob != 0. && normprobcal_ks[TAG_ks[i]] != 0.) 
	{
		LL1 = total_prob/normprobcal_ks[TAG_ks[i]] ; 
        	LogLL = LogLL + log(LL1);             
        }
     }

/////////////////////////////////////////////////////////////////////
// (1-f_1-f_2)*P_{K0Lpipi} + f_1*P_{Peak K0Spipi bkg} + f_2*(P_{SB}) 
/////////////////////////////////////////////////////////////////////

  for(Int_t i=0; i<datanumber_kl; i++)		
     {  


        Int_t binx = histo->GetXaxis()->FindBin(CMminus2_kl[i]);
        Int_t biny = histo->GetYaxis()->FindBin(CMplus2_kl[i]);

        Int_t binm = heffi->GetXaxis()->FindBin(CMminus2_kl[i]);
        Int_t binp = heffi->GetYaxis()->FindBin(CMplus2_kl[i]);

        TComplex acal(0.0, 0.0);
        TComplex acal_dcs(0.0, 0.0);
        TComplex acal_ks(0.0, 0.0);
        TComplex acal_ks_dcs(0.0, 0.0);
        total_prob=0.;
        ksbkg_prob=0.;
        pdf=0.;

        value0 = TComplex(0.,0.);
        value1 = TComplex(0.,0.);

        value0 = _beta[0]*A0_kl[i]+_beta[1]*A1_kl[i]+_beta[2]*A2_kl[i]+_beta[3]*A3_kl[i]+_beta[4]*A4_kl[i]; 
        value1 = (_f[0]*U1jCal_kl[i][0]+_f[1]*U1jCal_kl[i][1]+_f[2]*U1jCal_kl[i][2]+_f[3]*U1jCal_kl[i][3]+_f[4]*U1jCal_kl[i][4]) * CBB_kl[i][5];
        
        ///////////////////////////////////////////////////////////////////////


        for(int rr=0 ; rr<num_reso ; rr++)
        {
                acal += param_kl[rr]*CMulfacbw_kl[i][rr] ;
                acal_dcs += param_kl[rr]*CMulfacbw_dcs_kl[i][rr] ;

                acal_ks += param[rr]*CMulfacbw_kl[i][rr] ;
                acal_ks_dcs += param[rr]*CMulfacbw_dcs_kl[i][rr] ;
        }

        acal += CP_mult[4]*(value0 + value1);
        acal_dcs += CP_mult[4]*(value0 + value1);
        acal_ks += (value0 + value1);
        acal_ks_dcs += (value0 + value1);


        Interference = TComplex(1.0*cos(deltad[TAG_kl[i]]), -1.0*sin(deltad[TAG_kl[i]]))*acal_dcs*TComplex::Conjugate(acal) + TComplex(1.0*cos(deltad[TAG_kl[i]]), 1.0*sin(deltad[TAG_kl[i]]))*acal*TComplex::Conjugate(acal_dcs)   ;
        total_prob = acal.Rho2() + rd[TAG_kl[i]]*rd[TAG_kl[i]]*acal_dcs.Rho2() + rd[TAG_kl[i]]*Rf[TAG_kl[i]]*Interference ;

        Interference = TComplex(1.0*cos(deltad[TAG_kl[i]]), -1.0*sin(deltad[TAG_kl[i]]))*acal_ks_dcs*TComplex::Conjugate(acal_ks) + TComplex(1.0*cos(deltad[TAG_kl[i]]), 1.0*sin(deltad[TAG_kl[i]]))*acal_ks*TComplex::Conjugate(acal_ks_dcs)   ;
        ksbkg_prob = acal_ks.Rho2() + rd[TAG_kl[i]]*rd[TAG_kl[i]]*acal_ks_dcs.Rho2() - rd[TAG_kl[i]]*Rf[TAG_kl[i]]*Interference ;

        ///////////////////////////////////////////////////////////////////////

        Double_t corr_effi = effi[binm][binp][TAG_kl[i]] * Gamma_Factor_k0lpipi(Cpk0l[i],Cppip_kl[i],Cppim_kl[i]) ;

        if(normprobcal_kl[TAG_kl[i]]!=0.)    pdf = (1.-k0sbkg_frac[TAG_kl[i]]-NPbkg_frac[TAG_kl[i]])*(total_prob/normprobcal_kl[TAG_kl[i]]) ;

        if(normprobcal_ksbkg[TAG_kl[i]]!=0.) pdf += k0sbkg_frac[TAG_kl[i]]*(ksbkg_prob/normprobcal_ksbkg[TAG_kl[i]]) ; 

        if(corr_effi!=0.)                    pdf += NPbkg_frac[TAG_kl[i]]*(sb_pdf[binx][biny][TAG_kl[i]]/(corr_effi*SB_normfac[TAG_kl[i]])) ;	

        if(pdf>0.) LogLL_kl += log(pdf);

     }

  chisq =  pow((NPbkg_frac[0] - mc_np_frac[0])/mc_np_frac_err[0],2) ;
  chisq += pow((NPbkg_frac[1] - mc_np_frac[1])/mc_np_frac_err[1],2) ;
  chisq += pow((NPbkg_frac[2] - mc_np_frac[2])/mc_np_frac_err[2],2) ;

  //cout<<"LogLL: "<<LogLL<<endl;
  //cout<<"LogLL_kl: "<<LogLL_kl<<endl;
  //cout<<"chisq: "<<chisq<<endl;
  //cout<<"LogLL_ksbkg: "<<LogLL_ksbkg<<endl;
  //cout<<"LogLL_NPbkg: "<<LogLL_NPbkg<<endl;

  return -2.*(LogLL + LogLL_kl) + chisq ;

}


int main()
{
  const char * minName = "Minuit2" ;
  const char *algoName = "" ;
  int randomSeed = -1 ;

  time_t my_time = time(NULL);

  int ii, jj;

  double mass, width, magni, phase, value, xxx, yyy, errxxx, erryyy, errmagni, errphase;
  string mode;

  Double_t pdf1,pdf2,pdf3,pdf4,pdf5,pdf6;

  ifstream fin2_t0;
  fin2_t0.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/Data_SB_PDF_KPi_mBC.txt");
  while(1) {
  if(!fin2_t0.good()) break;
  fin2_t0 >> ii >> jj >> pdf1 >> pdf2  ;
  sb_pdf[ii][jj][0] = pdf2 ;
  }

  ifstream fin2_t1;
  fin2_t1.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/Data_SB_PDF_K3Pi_mBC.txt");
  while(1) {
  if(!fin2_t1.good()) break;
  fin2_t1 >> ii >> jj >> pdf1 >> pdf2 ;
  sb_pdf[ii][jj][1] = pdf2 ;
  }

  ifstream fin2_t2;
  fin2_t2.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/Data_SB_PDF_KPiPi0_mBC.txt");
  while(1) {
  if(!fin2_t2.good()) break;
  fin2_t2 >> ii >> jj >> pdf1 >> pdf2 ;
  sb_pdf[ii][jj][2] = pdf2 ;
  }

  cout<<"sb_pdf[41][57][2]: "<<sb_pdf[41][57][2]<<endl;

  std::cout<<"SB PDF done."<<endl;

//////////// sideband normalization factor /////////////

  SB_normfac[0] = 0.000547359 ;  
  SB_normfac[1] = 0.00137712 ;  
  SB_normfac[2] = 0.00112104 ;  

  double efficiency[3];

  TFile* fileIn2  = new TFile("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/effiprofile_PHSP.root");
  TTree* t2;
  t2 = (TTree*)fileIn2->Get("Effi");
  t2->SetBranchAddress("efficiency",efficiency);
  t2->SetBranchAddress("ii",&ii);				//ii and jj run from 1 to 50
  t2->SetBranchAddress("jj",&jj);

  for(int m=0; m<t2->GetEntries(); m++)
     { t2->GetEntry(m);
       for(int t=0; t<3; t++) {
       if(efficiency[t]>=0.0 && efficiency[t]<=1.0) effi[ii][jj][t] = efficiency[t];
       else effi[ii][jj][t] = 0.0;
       //effi[ii][jj] = 1.0;
       }
     }


  mc_np_frac[0] = 0.057792787 ;
  mc_np_frac[1] = 0.093381042 ;
  mc_np_frac[2] = 0.088815787 ;
  //mc_np_frac[0] = 0.0596 ;
  //mc_np_frac[1] = 0.0584 ;
  //mc_np_frac[2] = 0.0648 ;
  mc_np_frac_err[0] = 0.0016 ;
  mc_np_frac_err[1] = 0.0034 ;
  mc_np_frac_err[2] = 0.0023 ;

  Double_t mplus2, mminus2;
  Double_t Mplus2, Mminus2, P_K0L, P_K0S;
  Double_t P_PIP, P_PIM;

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

  TFile *filein_ks = TFile::Open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/ForNormFactorRealData6C.root");
  TTree *tree1_ks;
  TTree *tree2_ks;

  filein_ks->GetObject("Norm", tree1_ks);
  tree1_ks->SetBranchAddress("norm_tag", &norm_tag);
  tree1_ks->SetBranchAddress("mplus2", &mplus2);				
  tree1_ks->SetBranchAddress("mminus2", &mminus2);		
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
  tree2_ks->SetBranchAddress("P_K0S", &P_K0S);
  tree2_ks->SetBranchAddress("P_PIP", &P_PIP);
  tree2_ks->SetBranchAddress("P_PIM", &P_PIM);

  datanumber_ks = tree2_ks->GetEntries();
  cout<<"datanumber_ks: "<<datanumber_ks<<endl;

  //TFile *filein_kl = TFile::Open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/MulFactorRealDataEffiDCS.root");
  TFile *filein_kl = TFile::Open("/Users/anita/github/check.root");
  TTree *tree1_kl;
  TTree *tree2_kl;

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
  tree2_kl->SetBranchAddress("P_K0L", &P_K0L);
  tree2_kl->SetBranchAddress("P_PIP", &P_PIP);
  tree2_kl->SetBranchAddress("P_PIM", &P_PIM);

  datanumber_kl = tree2_kl->GetEntries();
  cout<<"datanumber_kl: "<<datanumber_kl<<endl;

  ifstream parin0;
  parin0.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/InitialParameters.dat");
  int freeparcounter=0;
  string proxystate ;
  int flag_ar=0, flag_ph=0;

  while(1) {
                parin0 >> mode >> proxystate >> flag_ar >> xxx >> errxxx >> flag_ph >> yyy >> erryyy ;
                if (!parin0.good()) break;
                if(mode[0] == '#') continue;

                if(flag_ar==0 && flag_ph==0)           freeparcounter+=2 ;
                else if(flag_ar==0 && flag_ph==2)      freeparcounter+=1 ;
                else if(flag_ar==2 && flag_ph==0)      freeparcounter+=1 ;

           }

  std::cout<<"freeparcounter: "<<freeparcounter<<std::endl;

  const int npar = freeparcounter+10+3 ;
  std::cout<<"npar: "<<npar<<std::endl;

  double par[npar];
  string parName[npar];
  double stepSize[npar];

  NPAR = npar;

  string pp, OmittedComponent;

  ifstream fin4;
  fin4.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/LASSFittedParameterInfoBelle2018.dat");

  while (1) {  
               if (!fin4.good()) break;
	       fin4 >> pp >> value ;
               lassparammap.insert(pair<string, double>(pp,value));
            }     
                   
  ifstream parin1;  
  parin1.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/InitialParameters.dat");
  int i=0;        
  int j=0;        
  int k=0;        

  while(1) {	  parin1 >> mode >> state[k] >> flag_ar >> xxx >> errxxx >> flag_ph >> yyy >> erryyy ;
		  if (!parin1.good()) break;
                  if(mode[0] == '#') continue; 
                  fixflag_ar[k] = flag_ar;
                  fixflag_ph[k] = flag_ph;
		  resonance[k] = mode;
                  
		  amplmap[mode] = xxx; phasemap[mode] = yyy*(TMath::Pi()/180.);
		  erramplmap[mode] = errxxx; errphasemap[mode] = erryyy*(TMath::Pi()/180.) ;  

	          if(flag_ar==0) {
		  par[i] = xxx; stepSize[i] = errxxx ; parName[i]="ar_"+mode;
                  i++;
                  if(flag_ph==0) {
                  par[i] = yyy*(TMath::Pi()/180.); stepSize[i] = erryyy*(TMath::Pi()/180.) ; parName[i]="phir_"+mode;
                  i++; 
                  }
                  else if(flag_ph==2) {
                  fixedpar[j] = yyy*(TMath::Pi()/180.);
                  j++;
                  }
                  }
		  else if(flag_ar==2) {
		  fixedpar[j] = xxx;
                  j++;
		  if(flag_ph==0) {
                  par[i] = yyy*(TMath::Pi()/180.);
		  stepSize[i] = erryyy*(TMath::Pi()/180.) ; parName[i]="phir_"+mode;
		  i++;
                  }
                  else if(flag_ph==2) {
                  fixedpar[j] = yyy*(TMath::Pi()/180.);
                  j++;
                  if(j==4) OmittedComponent = mode ;
                  }
		  }
                  cout<<resonance[k]<<"\t"<<xxx<<"\t"<<yyy<<endl;
		  k++;
           }

  num_reso = k - 8 ;
  cout<<"num_reso: "<<num_reso<<endl;
            
  std::cout<<"OmittedComponent: "<<OmittedComponent<<std::endl;

  string cpresonance[10];

  k=0;
  Int_t fixflag_r[5], fixflag_dlt[5];

  ifstream parin2;  
  parin2.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/InitialParameters_rdelta.dat");
  while(1) {      parin2 >> mode >> flag_ar >> xxx >> errxxx >> flag_ph >> yyy >> erryyy ;
	          if (!parin2.good()) break;
                  if(mode[0] == '#') continue;
                  fixflag_r[k] = flag_ar;
                  fixflag_dlt[k] = flag_ph;

	          if(flag_ar==0) {
		  par[i] = xxx; stepSize[i] = errxxx ; parName[i]="r_"+mode;
                  i++;
                  if(flag_ph==0) {
                  par[i] = yyy*(TMath::Pi()/180.); stepSize[i] = erryyy*(TMath::Pi()/180.) ; parName[i]="delta_"+mode;
                  i++; 
                  }
                  else if(flag_ph==2) {
                  fixedpar[j] = yyy*(TMath::Pi()/180.);
                  j++;
                  }
                  }
		  else if(flag_ar==2) {
		  fixedpar[j] = xxx;
                  j++;
		  if(flag_ph==0) {
                  par[i] = yyy*(TMath::Pi()/180.);
		  stepSize[i] = erryyy*(TMath::Pi()/180.) ; parName[i]="delta_"+mode;
		  i++;
                  }
                  else if(flag_ph==2) {
                  fixedpar[j] = yyy*(TMath::Pi()/180.);
                  j++;
                  }
		  }      
		  cpresonance[k] = "rdelta_"+mode;
	          cout<<mode<<"\t"<<xxx<<"\t"<<yyy<<endl;
		  k++;
            }

  std::cout<<"npar-3: "<<npar-3<<endl;

  par[npar-3]      = 0.0578 ;
  par[npar-2]      = 0.0933 ;
  par[npar-1]      = 0.0888 ;

  parName[npar-3]  = "NPbkg_frac_kpi ";
  parName[npar-2]  = "NPbkg_frac_k3pi";
  parName[npar-1]  = "NPbkg_frac_kpipi0";

  stepSize[npar-3] = 0.0005 ;    
  stepSize[npar-2] = 0.0005 ;    
  stepSize[npar-1] = 0.0005 ;    

  k0sbkg_frac[0] = 0.047909717 ; //0.048004682 ;
  k0sbkg_frac[1] = 0.040662811 ;
  k0sbkg_frac[2] = 0.048558820 ;


  int l=0;
  string modecf, modedcs, reso_cf[20], reso_dcs[20];
  ifstream fin5;
  fin5.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/AllResonances.dat") ;
  while(1) { 
             fin5 >> modecf >> modedcs ;
	     if (!fin5.good()) break;
	     if(modecf[0] == '#') continue;
             reso_cf[l] = modecf ;
	     reso_dcs[l] = modedcs ;
             cout<<reso_cf[l]<<"\t"<<reso_dcs[l]<<endl;
             l++;
           }


  allnum_reso = l; 				

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

ifstream yield_fin;
yield_fin.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/kskl_yields.dat");

ifstream norm_fin0;
ifstream norm_fin1;
ifstream norm_fin2;
ifstream norm_fin3;
ifstream norm_fin4;

norm_fin0.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/r_t.dat");
norm_fin1.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/r_rr_t.dat");
norm_fin2.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/r_l_t.dat");
norm_fin3.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/l_t.dat");
norm_fin4.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/l_k_t.dat");

yield_fin >> normnumber_ks >> taggednormnumber_ks[0] >> taggednormnumber_ks[1] >> taggednormnumber_ks[2] ;
yield_fin >> normnumber_kl >> taggednormnumber_kl[0] >> taggednormnumber_kl[1] >> taggednormnumber_kl[2] ;

cout<<"normnumber_ks: "<<normnumber_ks<<endl;
cout<<"normnumber_kl: "<<normnumber_kl<<endl;

for(int r=0; r<num_reso; r++) {
norm_fin0 >> probZ_ks[r][0] >> probZ_ks[r][1] >> probZ_ks[r][2] ;
}
for(int r=0; r<num_reso; r++) {
norm_fin0 >> probZ_dcs_ks[r][0] >> probZ_dcs_ks[r][1] >> probZ_dcs_ks[r][2] ;
}
for(int r=0; r<num_reso; r++) {
norm_fin0 >> probZ_kl[r][0] >> probZ_kl[r][1] >> probZ_kl[r][2] ;
}
for(int r=0; r<num_reso; r++) {
norm_fin0 >> probZ_dcs_kl[r][0] >> probZ_dcs_kl[r][1] >> probZ_dcs_kl[r][2] ;
}


for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
norm_fin1 >> reinterZ_ks[r][rr][0] >> reinterZ_ks[r][rr][1] >> reinterZ_ks[r][rr][2] >> iminterZ_ks[r][rr][0] >> iminterZ_ks[r][rr][1] >> iminterZ_ks[r][rr][2] ;
}}
for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
norm_fin1 >> reinterZ_dcs_ks[r][rr][0] >> reinterZ_dcs_ks[r][rr][1] >> reinterZ_dcs_ks[r][rr][2] >> iminterZ_dcs_ks[r][rr][0] >> iminterZ_dcs_ks[r][rr][1] >> iminterZ_dcs_ks[r][rr][2] ;
}}
for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
norm_fin1 >> reinterZ_cfdcs_ks[r][rr][0] >> reinterZ_cfdcs_ks[r][rr][1] >> reinterZ_cfdcs_ks[r][rr][2] >> iminterZ_cfdcs_ks[r][rr][0] >> iminterZ_cfdcs_ks[r][rr][1] >> iminterZ_cfdcs_ks[r][rr][2] ;
}}
for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
norm_fin1 >> reinterZ_kl[r][rr][0] >> reinterZ_kl[r][rr][1] >> reinterZ_kl[r][rr][2] >> iminterZ_kl[r][rr][0] >> iminterZ_kl[r][rr][1] >> iminterZ_kl[r][rr][2] ;
}}
for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
norm_fin1 >> reinterZ_dcs_kl[r][rr][0] >> reinterZ_dcs_kl[r][rr][1] >> reinterZ_dcs_kl[r][rr][2] >> iminterZ_dcs_kl[r][rr][0] >> iminterZ_dcs_kl[r][rr][1] >> iminterZ_dcs_kl[r][rr][2] ;
}}
for(int r=0; r<num_reso; r++) {
for(int rr=0; rr<num_reso; rr++) {
norm_fin1 >> reinterZ_cfdcs_kl[r][rr][0] >> reinterZ_cfdcs_kl[r][rr][1] >> reinterZ_cfdcs_kl[r][rr][2] >> iminterZ_cfdcs_kl[r][rr][0] >> iminterZ_cfdcs_kl[r][rr][1] >> iminterZ_cfdcs_kl[r][rr][2] ;
}}

for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
norm_fin2 >> reinterZval0_ks[r][l][0] >> reinterZval0_ks[r][l][1] >> reinterZval0_ks[r][l][2] >> iminterZval0_ks[r][l][0] >> iminterZval0_ks[r][l][1] >> iminterZval0_ks[r][l][2] >> reinterZval1_ks[r][l][0] >> reinterZval1_ks[r][l][1] >> reinterZval1_ks[r][l][2] >> iminterZval1_ks[r][l][0] >> iminterZval1_ks[r][l][1] >> iminterZval1_ks[r][l][2] ;
}}
for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
norm_fin2 >> reinterZval0_dcs_ks[r][l][0] >> reinterZval0_dcs_ks[r][l][1] >> reinterZval0_dcs_ks[r][l][2] >> iminterZval0_dcs_ks[r][l][0] >> iminterZval0_dcs_ks[r][l][1] >> iminterZval0_dcs_ks[r][l][2] >> reinterZval1_dcs_ks[r][l][0] >> reinterZval1_dcs_ks[r][l][1] >> reinterZval1_dcs_ks[r][l][2] >> iminterZval1_dcs_ks[r][l][0] >> iminterZval1_dcs_ks[r][l][1] >> iminterZval1_dcs_ks[r][l][2] ;
}}
for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
norm_fin2 >> reinterval0Z_ks[r][l][0] >> reinterval0Z_ks[r][l][1] >> reinterval0Z_ks[r][l][2] >> iminterval0Z_ks[r][l][0] >> iminterval0Z_ks[r][l][1] >> iminterval0Z_ks[r][l][2] >> reinterval1Z_ks[r][l][0] >> reinterval1Z_ks[r][l][1] >> reinterval1Z_ks[r][l][2] >> iminterval1Z_ks[r][l][0] >> iminterval1Z_ks[r][l][1] >> iminterval1Z_ks[r][l][2] ;
}}
for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
norm_fin2 >> reinterZval0_kl[r][l][0] >> reinterZval0_kl[r][l][1] >> reinterZval0_kl[r][l][2] >> iminterZval0_kl[r][l][0] >> iminterZval0_kl[r][l][1] >> iminterZval0_kl[r][l][2] >> reinterZval1_kl[r][l][0] >> reinterZval1_kl[r][l][1] >> reinterZval1_kl[r][l][2] >> iminterZval1_kl[r][l][0] >> iminterZval1_kl[r][l][1] >> iminterZval1_kl[r][l][2] ;
}}
for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
norm_fin2 >> reinterZval0_dcs_kl[r][l][0] >> reinterZval0_dcs_kl[r][l][1] >> reinterZval0_dcs_kl[r][l][2] >> iminterZval0_dcs_kl[r][l][0] >> iminterZval0_dcs_kl[r][l][1] >> iminterZval0_dcs_kl[r][l][2] >> reinterZval1_dcs_kl[r][l][0] >> reinterZval1_dcs_kl[r][l][1] >> reinterZval1_dcs_kl[r][l][2] >> iminterZval1_dcs_kl[r][l][0] >> iminterZval1_dcs_kl[r][l][1] >> iminterZval1_dcs_kl[r][l][2] ;
}}
for(int r=0; r<num_reso; r++) {
for(int l=0; l<5; l++) {
norm_fin2 >> reinterval0Z_kl[r][l][0] >> reinterval0Z_kl[r][l][1] >> reinterval0Z_kl[r][l][2] >> iminterval0Z_kl[r][l][0] >> iminterval0Z_kl[r][l][1] >> iminterval0Z_kl[r][l][2] >> reinterval1Z_kl[r][l][0] >> reinterval1Z_kl[r][l][1] >> reinterval1Z_kl[r][l][2] >> iminterval1Z_kl[r][l][0] >> iminterval1Z_kl[r][l][1] >> iminterval1Z_kl[r][l][2] ;
}}

for(int l=0; l<5; l++) {
norm_fin3 >> mod2L_ks[l][0] >> mod2L_ks[l][1] >> mod2L_ks[l][2] >> mod2F_ks[l][0] >> mod2F_ks[l][1] >> mod2F_ks[l][2] ;
}
for(int l=0; l<5; l++) {
norm_fin3 >> mod2L_kl[l][0] >> mod2L_kl[l][1] >> mod2L_kl[l][2] >> mod2F_kl[l][0] >> mod2F_kl[l][1] >> mod2F_kl[l][2] ;
}

for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
norm_fin4 >> reinterKM0_ks[l][k][0] >> reinterKM0_ks[l][k][1] >> reinterKM0_ks[l][k][2] >> iminterKM0_ks[l][k][0] >> iminterKM0_ks[l][k][1] >> iminterKM0_ks[l][k][2] ;
}}
for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
norm_fin4 >> reinterKM1_ks[l][k][0] >> reinterKM1_ks[l][k][1] >> reinterKM1_ks[l][k][2] >> iminterKM1_ks[l][k][0] >> iminterKM1_ks[l][k][1] >> iminterKM1_ks[l][k][2] ;
}}
for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
norm_fin4 >> reinterI_ks[l][k][0] >> reinterI_ks[l][k][1] >> reinterI_ks[l][k][2] >> iminterI_ks[l][k][0] >> iminterI_ks[l][k][1] >> iminterI_ks[l][k][2] ;
}}
for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
norm_fin4 >> reinterKM0_kl[l][k][0] >> reinterKM0_kl[l][k][1] >> reinterKM0_kl[l][k][2] >> iminterKM0_kl[l][k][0] >> iminterKM0_kl[l][k][1] >> iminterKM0_kl[l][k][2] ;
}}
for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
norm_fin4 >> reinterKM1_kl[l][k][0] >> reinterKM1_kl[l][k][1] >> reinterKM1_kl[l][k][2] >> iminterKM1_kl[l][k][0] >> iminterKM1_kl[l][k][1] >> iminterKM1_kl[l][k][2] ;
}}
for(int l=0; l<5; l++) {
for(int k=0; k<5; k++) {
norm_fin4 >> reinterI_kl[l][k][0] >> reinterI_kl[l][k][1] >> reinterI_kl[l][k][2] >> iminterI_kl[l][k][0] >> iminterI_kl[l][k][1] >> iminterI_kl[l][k][2] ;
}}
    
std::cout<<"Norm sample maps and vectors ready..."<<std::endl;

//////////////////////////////////////////////// DataSample - K0Spipi ///////////////////////////////////////////////////


  for(Int_t j=0; j<datanumber_ks; j++)	
     { tree2_ks->GetEntry(j);
       TAG_ks[j] = tag;
       CMplus2_ks[j] = Mplus2;   
       CMminus2_ks[j] = Mminus2;
       Cpk0s[j] = P_K0S;
       Cppip_ks[j] = P_PIP;
       Cppim_ks[j] = P_PIM;
       //CMulfacbw_re_ks[j] = *Mulfacbw_re;				//Array of map
       //CMulfacbw_im_ks[j] = *Mulfacbw_im;				//Array of map
       //CLASS_contribution_re_ks[j] = *LASS_contribution_re;	//Array of map
       //CLASS_contribution_im_ks[j] = *LASS_contribution_im;	//Array of map

       for(int rr=0 ; rr<num_reso-2 ; rr++)
	{
		CMulfacbw_ks[j].push_back(TComplex( (*Mulfacbw_re)[reso_cf[rr]] , (*Mulfacbw_im)[reso_cf[rr]] )) ;
		CMulfacbw_dcs_ks[j].push_back(TComplex( (*Mulfacbw_re)[reso_dcs[rr]] , (*Mulfacbw_im)[reso_dcs[rr]] )) ;

	}	

       CMulfacbw_ks[j].push_back(TComplex( (*LASS_contribution_re)[reso_cf[num_reso-2]] , (*LASS_contribution_im)[reso_cf[num_reso-2]] ));
       CMulfacbw_dcs_ks[j].push_back(TComplex( (*LASS_contribution_re)[reso_dcs[num_reso-2]], (*LASS_contribution_im)[reso_dcs[num_reso-2]] ));
       CMulfacbw_ks[j].push_back(TComplex( (*LASS_contribution_re)[reso_cf[num_reso-1]], (*LASS_contribution_im)[reso_cf[num_reso-1]] ));
       CMulfacbw_dcs_ks[j].push_back(TComplex( (*LASS_contribution_re)[reso_dcs[num_reso-1]], (*LASS_contribution_im)[reso_dcs[num_reso-1]] ));

       CU1j_re_ks = *U1j_re;					//array of vector
       CU1j_im_ks = *U1j_im;					//array of vector
       CBB_ks[j] = *BB;						//array of vector

       A0_ks[j] = TComplex(Asum_re[0], Asum_im[0]) ;
       A1_ks[j] = TComplex(Asum_re[1], Asum_im[1]) ;
       A2_ks[j] = TComplex(Asum_re[2], Asum_im[2]) ;
       A3_ks[j] = TComplex(Asum_re[3], Asum_im[3]) ;
       A4_ks[j] = TComplex(Asum_re[4], Asum_im[4]) ;

       U1jCal_ks[j][0] = TComplex(CU1j_re_ks[0], CU1j_im_ks[0]);
       U1jCal_ks[j][1] = TComplex(CU1j_re_ks[1], CU1j_im_ks[1]);
       U1jCal_ks[j][2] = TComplex(CU1j_re_ks[2], CU1j_im_ks[2]);
       U1jCal_ks[j][3] = TComplex(CU1j_re_ks[3], CU1j_im_ks[3]);
       U1jCal_ks[j][4] = TComplex(CU1j_re_ks[4], CU1j_im_ks[4]);

     }

//////////////////////////////////////////// DataSample - K0Lpipi //////////////////////////////////////////////////


  for(Int_t j=0; j<datanumber_kl; j++)	
     { tree2_kl->GetEntry(j);
       TAG_kl[j] = tag;
       CMplus2_kl[j] = Mplus2;   
       CMminus2_kl[j] = Mminus2;
       Cpk0l[j] = P_K0L;
       Cppip_kl[j] = P_PIP;
       Cppim_kl[j] = P_PIM;
       //CMulfacbw_re_kl[j] = *Mulfacbw_re;				//Array of map
       //CMulfacbw_im_kl[j] = *Mulfacbw_im;				//Array of map
       //CLASS_contribution_re_kl[j] = *LASS_contribution_re;	//Array of map
       //CLASS_contribution_im_kl[j] = *LASS_contribution_im;	//Array of map

       for(int rr=0 ; rr<num_reso-2 ; rr++)
        {
                CMulfacbw_kl[j].push_back(TComplex( (*Mulfacbw_re)[reso_cf[rr]] , (*Mulfacbw_im)[reso_cf[rr]] )) ;
                CMulfacbw_dcs_kl[j].push_back(TComplex( (*Mulfacbw_re)[reso_dcs[rr]] , (*Mulfacbw_im)[reso_dcs[rr]] )) ;
        }

       CMulfacbw_kl[j].push_back(TComplex( (*LASS_contribution_re)[reso_cf[num_reso-2]], (*LASS_contribution_im)[reso_cf[num_reso-2]] ));
       CMulfacbw_dcs_kl[j].push_back(TComplex( (*LASS_contribution_re)[reso_dcs[num_reso-2]], (*LASS_contribution_im)[reso_dcs[num_reso-2]] ));
       CMulfacbw_kl[j].push_back(TComplex( (*LASS_contribution_re)[reso_cf[num_reso-1]], (*LASS_contribution_im)[reso_cf[num_reso-1]] ));
       CMulfacbw_dcs_kl[j].push_back(TComplex( (*LASS_contribution_re)[reso_dcs[num_reso-1]], (*LASS_contribution_im)[reso_dcs[num_reso-1]] ));


       CU1j_re_kl = *U1j_re;					//array of vector
       CU1j_im_kl = *U1j_im;					//array of vector
       CBB_kl[j] = *BB;						//array of vector

       A0_kl[j] = TComplex(Asum_re[0], Asum_im[0]) ;
       A1_kl[j] = TComplex(Asum_re[1], Asum_im[1]) ;
       A2_kl[j] = TComplex(Asum_re[2], Asum_im[2]) ;
       A3_kl[j] = TComplex(Asum_re[3], Asum_im[3]) ;
       A4_kl[j] = TComplex(Asum_re[4], Asum_im[4]) ;

       U1jCal_kl[j][0] = TComplex(CU1j_re_kl[0], CU1j_im_kl[0]);
       U1jCal_kl[j][1] = TComplex(CU1j_re_kl[1], CU1j_im_kl[1]);
       U1jCal_kl[j][2] = TComplex(CU1j_re_kl[2], CU1j_im_kl[2]);
       U1jCal_kl[j][3] = TComplex(CU1j_re_kl[3], CU1j_im_kl[3]);
       U1jCal_kl[j][4] = TComplex(CU1j_re_kl[4], CU1j_im_kl[4]);

     }

  std::cout<<"Data sample maps and vectors ready..."<<std::endl;

/////////////////////////////////////////////////////// Minimization ////////////////////////////////////////////////////

   ROOT::Math::Minimizer* min =  ROOT::Math::Factory::CreateMinimizer(minName, algoName);

   min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2 
   min->SetMaxIterations(10000000);  // for GSL 
   min->SetTolerance(0.001);
   min->SetPrintLevel(2);

   ROOT::Math::Functor f(&likelihood,npar); 
    
   min->SetFunction(f);
 
   for(int j=0; j<npar; j++) { 
   min->SetVariable(j, parName[j].c_str(), par[j], stepSize[j]);
   }

   min->SetVariableLimits(npar-3,0.,1.); 
   min->SetVariableLimits(npar-2,0.,1.); 
   min->SetVariableLimits(npar-1,0.,1.); 

   // minimization
   min->Minimize(); 

  const double *outpar_mig = min->X();
  const double *err_mig = min->Errors();


  string head_c = "CovMat_" ; 
  string head_a = "ArDeltar_" ; 
  string head_u = "UspinBkgFrac_" ; 
  string tail = "050622.dat" ;

  string filename_c = head_c + tail ;
  string filename_a = head_a + tail ;
  string filename_u = head_u + tail ;

   double covmat_mig[npar][npar];
   min->GetCovMatrix(*covmat_mig) ;
   ofstream fop_mig;
   fop_mig.open(filename_c) ; //, std::ios::app);
   for(int i=0; i<npar; i++) {
   for(int j=0; j<npar; j++) {
   //cout<<covmat[i][j]<<"\t" ;
   fop_mig << covmat_mig[i][j] << "\t" ;
   }
   cout<<endl;
   fop_mig << endl;
   }
   fop_mig << endl;

   i=0; j=0;
   ofstream fout_mig;
   fout_mig.open(filename_a) ; //,std::ios::app); 

   for( int r = 0 ; r < (num_reso+8) ; r++ )
        {
                if(fixflag_ar[r]==0 && fixflag_ph[r]==0)
                {
                        fout_mig.setf(ios::left);
                        fout_mig << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << outpar_mig[i] << setw(18) << setprecision(4) << err_mig[i] << setw(18) << setprecision(4) << outpar_mig[i+1]*(180./TMath::Pi()) << setw(18) << setprecision(4) <<err_mig[i+1]*(180./TMath::Pi())<<std::endl;
                        cout << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << outpar_mig[i] << setw(18) << setprecision(4) << err_mig[i] << setw(18) << setprecision(4) << outpar_mig[i+1]*(180./TMath::Pi()) << setw(18) << setprecision(4) <<err_mig[i+1]*(180./TMath::Pi())<<std::endl;
                        i+=2;
                }
                else if(fixflag_ar[r]==0 && fixflag_ph[r]==2)
                {
                        fout_mig.setf(ios::left);
                        fout_mig << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << outpar_mig[i] << setw(18) << setprecision(4) << err_mig[i] << setw(18) << setprecision(4) << fixedpar[j]*(180./TMath::Pi()) << setw(18) << setprecision(4) << 0.0*(180./TMath::Pi()) << std::endl;
                        cout << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << outpar_mig[i] << setw(18) << setprecision(4) << err_mig[i] << setw(18) << setprecision(4) << fixedpar[j]*(180./TMath::Pi()) << setw(18) << setprecision(4) << 0.0*(180./TMath::Pi()) << std::endl;
                        i++; j++;
                }
                else if(fixflag_ar[r]==2 && fixflag_ph[r]==0)
                {
                        fout_mig.setf(ios::left);
                        fout_mig << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << fixedpar[j] << setw(18) << setprecision(4) << 0.0 << setw(18) << setprecision(4) << outpar_mig[i]*(180./TMath::Pi()) << setw(18) << setprecision(4) << err_mig[i]*(180./TMath::Pi()) << std::endl;
                        cout << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << fixedpar[j] << setw(18) << setprecision(4) << 0.0 << setw(18) << setprecision(4) << outpar_mig[i]*(180./TMath::Pi()) << setw(18) << setprecision(4) << err_mig[i]*(180./TMath::Pi()) << std::endl;
                        i++; j++;
                }
                else if(fixflag_ar[r]==2 && fixflag_ph[r]==2)
                {
                        fout_mig.setf(ios::left);
                        fout_mig << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << fixedpar[j] << setw(18) << setprecision(4) << 0.0 << setw(18) << setprecision(4) << fixedpar[j+1]*(180./TMath::Pi())<< setw(18) << setprecision(4) <<0.0*(180./TMath::Pi()) << std::endl;
                        cout << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << fixedpar[j] << setw(18) << setprecision(4) << 0.0 << setw(18) << setprecision(4) << fixedpar[j+1]*(180./TMath::Pi())<< setw(18) << setprecision(4) <<0.0*(180./TMath::Pi()) << std::endl;
                        j+=2;
                }
        }
   cout << endl;

  ofstream fout2_mig;
  fout2_mig.open(filename_u) ; //,std::ios::app) ;
  for(int r=0; r<5; r++) 
	{
		if(fixflag_r[r]==0 && fixflag_dlt[r]==0) 
		{
			fout2_mig.setf(ios::left);
			fout2_mig << setw(18) << cpresonance[r] << setw(18) << setprecision(4) << outpar_mig[i] << setw(18) << setprecision(4) << err_mig[i] << setw(18) << setprecision(4) << outpar_mig[i+1]*(180./TMath::Pi()) << setw(18) << setprecision(4) << err_mig[i+1]*(180./TMath::Pi()) << std::endl;
			cout << setw(18) << cpresonance[r] << setw(18) << setprecision(4) << outpar_mig[i] << setw(18) << setprecision(4) << err_mig[i] << setw(18) << setprecision(4) << outpar_mig[i+1]*(180./TMath::Pi()) << setw(18) << setprecision(4) << err_mig[i+1]*(180./TMath::Pi()) << std::endl;
  		i+=2;
  		}

		if(fixflag_r[r]==2 && fixflag_dlt[r]==2) 
		{
			fout2_mig.setf(ios::left);
			fout2_mig << setw(18) << cpresonance[r] << setw(18) << setprecision(4) << fixedpar[j] << setw(18) << setprecision(4) << 0.0 << setw(18) << setprecision(4) << fixedpar[j+1]*(180./TMath::Pi()) << setw(18) << setprecision(4) << 0.0*(180./TMath::Pi()) << std::endl;
			cout << setw(18) << cpresonance[r] << setw(18) << setprecision(4) << fixedpar[j] << setw(18) << setprecision(4) << 0.0 << setw(18) << setprecision(4) << fixedpar[j+1]*(180./TMath::Pi()) << setw(18) << setprecision(4) << 0.0*(180./TMath::Pi()) << std::endl;
			j+=2;
		}

  	}

  std::cout<<"i in printing: "<<i<<std::endl;

  fout2_mig.setf(ios::left);
  fout2_mig << "NPbkg_frac_kpi"   << setw(24) << setprecision(4) << outpar_mig[npar-3] << setw(18) << setprecision(4) << err_mig[npar-3] << "\t\t" << 0. << "\t\t" << 0. << std::endl;
  fout2_mig.setf(ios::left);
  fout2_mig << "NPbkg_frac_k3pi"  << setw(24) << setprecision(4) << outpar_mig[npar-2] << setw(18) << setprecision(4) << err_mig[npar-2] << "\t\t" << 0. << "\t\t" << 0. << std::endl;
  fout2_mig.setf(ios::left);
  fout2_mig << "NPbkg_frac_kpipi0"<< setw(24) << setprecision(4) << outpar_mig[npar-1] << setw(18) << setprecision(4) << err_mig[npar-1] << "\t\t" << 0. << "\t\t" << 0. << std::endl;

  fout2_mig.setf(ios::left);
  fout2_mig << "k0sbkg_frac_kpi    " << setw(24) << setprecision(4) << k0sbkg_frac[0] << "\t\t" << 0. << "\t\t" << 0. << "\t\t" << 0. << std::endl;
  fout2_mig.setf(ios::left);
  fout2_mig << "k0sbkg_frac_k3pi   " << setw(24) << setprecision(4) << k0sbkg_frac[1] << "\t\t" << 0. << "\t\t" << 0. << "\t\t" << 0. << std::endl;
  fout2_mig.setf(ios::left);
  fout2_mig << "k0sbkg_frac_kpipi0 " << setw(24) << setprecision(4) << k0sbkg_frac[2] << "\t\t" << 0. << "\t\t" << 0. << "\t\t" << 0. << std::endl;

  fout2_mig << endl;


  //HESSE minimization
  
  min->Hesse();

  int hesse_cov_status = min->CovMatrixStatus();
  cout<<"Covariance matrix status (anita): "<<hesse_cov_status<<endl;

   double covmat[npar][npar];
   ofstream fout2;
   ofstream fout;
   string resoName[50];
   const double *outpar ;
   const double *err ;
   ofstream fop;
   double elow,eup = 0;

   if(hesse_cov_status != 3) goto label ;

   min->GetCovMatrix(*covmat) ;

   fop.open(filename_c);
   cout<<"Covariance matrix..."<<endl;
   for(int i=0; i<npar; i++) {
   for(int j=0; j<npar; j++) {
        cout<<covmat[i][j]<<"\t" ;
        fop << covmat[i][j] << "\t" ;
   }
   cout<<endl;
   fop << endl;
   }
   fop << endl;
   outpar = min->X();
   err = min->Errors();

   i=0; j=0;
   fout.open(filename_a) ; //,std::ios::app); 

   for( int r = 0 ; r < (num_reso+8) ; r++ )
        {
                if(fixflag_ar[r]==0 && fixflag_ph[r]==0)
                {
                        fout.setf(ios::left);
                        fout << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << outpar[i] << setw(18) << setprecision(4) << err[i] << setw(18) << setprecision(4) << outpar[i+1]*(180./TMath::Pi()) << setw(18) << setprecision(4) <<err[i+1]*(180./TMath::Pi())<<std::endl;
                        i+=2;
                }
                else if(fixflag_ar[r]==0 && fixflag_ph[r]==2)
                {
                        fout.setf(ios::left);
                        fout << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << outpar[i] << setw(18) << setprecision(4) << err[i] << setw(18) << setprecision(4) << fixedpar[j]*(180./TMath::Pi()) << setw(18) << setprecision(4) << 0.0*(180./TMath::Pi()) << std::endl;
                        i++; j++;
                }
                else if(fixflag_ar[r]==2 && fixflag_ph[r]==0)
                {
                        fout.setf(ios::left);
                        fout << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << fixedpar[j] << setw(18) << setprecision(4) << 0.0 << setw(18) << setprecision(4) << outpar[i]*(180./TMath::Pi()) << setw(18) << setprecision(4) << err[i]*(180./TMath::Pi()) << std::endl;
                        i++; j++;
                }
                else if(fixflag_ar[r]==2 && fixflag_ph[r]==2)
                {
                        fout.setf(ios::left);
                        fout << setw(18) << resonance[r] << setw(18) << state[r] << setw(18) << setprecision(4) << fixedpar[j] << setw(18) << setprecision(4) << 0.0 << setw(18) << setprecision(4) << fixedpar[j+1]*(180./TMath::Pi())<< setw(18) << setprecision(4) <<0.0*(180./TMath::Pi()) << std::endl;
                        j+=2;
                }
        }
   fout << endl;


  fout2.open(filename_u) ; //,std::ios::app) ;
  for(int r=0; r<5; r++)
        {
                if(fixflag_r[r]==0 && fixflag_dlt[r]==0)
                {
                        fout2.setf(ios::left);
                        fout2 << setw(18) << cpresonance[r] << setw(18) << setprecision(4) << outpar[i] << setw(18) << setprecision(4) << err[i] << setw(18) << setprecision(4) << outpar[i+1]*(180./TMath::Pi()) << setw(18) << setprecision(4) << err[i+1]*(180./TMath::Pi()) << std::endl;
                i+=2;
                }

                if(fixflag_r[r]==2 && fixflag_dlt[r]==2)
                {
                        fout2.setf(ios::left);
                        fout2 << setw(18) << cpresonance[r] << setw(18) << setprecision(4) << fixedpar[j] << setw(18) << setprecision(4) << 0.0 << setw(18) << setprecision(4) << fixedpar[j+1]*(180./TMath::Pi()) << setw(18) << setprecision(4) << 0.0*(180./TMath::Pi()) << std::endl;
                        j+=2;
                }

        }


  fout2.setf(ios::left);
  fout2 << "NPbkg_frac_kpi"   << setw(24) << setprecision(4) << outpar[npar-3] << setw(18) << setprecision(4) << err[npar-3] << "\t\t" << 0. << "\t\t" << 0. << std::endl;
  fout2.setf(ios::left);
  fout2 << "NPbkg_frac_k3pi"  << setw(24) << setprecision(4) << outpar[npar-2] << setw(18) << setprecision(4) << err[npar-2] << "\t\t" << 0. << "\t\t" << 0. << std::endl;
  fout2.setf(ios::left);
  fout2 << "NPbkg_frac_kpipi0"<< setw(24) << setprecision(4) << outpar[npar-1] << setw(18) << setprecision(4) << err[npar-1] << "\t\t" << 0. << "\t\t" << 0. << std::endl;

  fout2.setf(ios::left);
  fout2 << "k0sbkg_frac_kpi    " << setw(24) << setprecision(4) << k0sbkg_frac[0] << "\t\t" << 0. << "\t\t" << 0. << "\t\t" << 0. << std::endl;
  fout2.setf(ios::left);
  fout2 << "k0sbkg_frac_k3pi   " << setw(24) << setprecision(4) << k0sbkg_frac[1] << "\t\t" << 0. << "\t\t" << 0. << "\t\t" << 0. << std::endl;
  fout2.setf(ios::left);
  fout2 << "k0sbkg_frac_kpipi0 " << setw(24) << setprecision(4) << k0sbkg_frac[2] << "\t\t" << 0. << "\t\t" << 0. << "\t\t" << 0. << std::endl;

  fout2 << endl;


  label:





  time_t my_endtime = time(NULL);
  std::cout<<"Start time: "<<ctime(&my_time)<<endl;
  std::cout<<"End time: "<<ctime(&my_endtime)<<endl;

  return 0;

}


