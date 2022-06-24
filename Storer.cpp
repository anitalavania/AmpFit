#include "TMinuit.h"
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
#include "TH2F.h"
#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TVirtualFitter.h"
#include "TRandom3.h"
#include "amp_func.h"
//#include "omp.h"

using namespace std; 
  
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
bool Gen;
bool Fit;

double g[5][5];

double mass, width, magni, phase, value, xxx, yyy, errxxx, erryyy, errmagni, errphase;
string mode, param;
double effi[60][60];
Double_t sk0spip[5000000], sk0spim[5000000], k0lmom[5000000], sK0SPIP[50000], sK0SPIM[50000], K0Lmom[50000];
Double_t PIPmom[50000], PIMmom[50000], pipmom[5000000], pimmom[5000000];
Int_t TAG[50000], norm_TAG[5000000];

int datanumber, normnumber;

class unicorn	
{
	public:

	string horse;

	void animal()
	{
		std::cout << "Animal name is "<<horse<<endl;
	}
};


TComplex amplitude_BW(double *phsp, double mR, double gammaR, int spin, string reso)
{ 

  double pi180inv = TMath::Pi() / 180.;

  TComplex matrixEl;

  double mm2 = phsp[0];
  double mp2 = phsp[1];

  const double mD0 = 1.86483;
  const double mKl = 0.49761;
  const double mPi = 0.13957;

  double R_r = 1.5; // "resonance radius" for Blatt-Weisskopf barrier factors.
  double R_D = 5.0; // "D meson radius" for Blatt-Weisskopf barrier factors.

  double fR,fD;

  double mpippim2 = mD0 * mD0 + mKl * mKl + (2 * mPi * mPi) - mp2 - mm2;

  double mA, mB, mC, mAB, mBC, mAC;

  if(reso=="k0spim") 
    { mA = mKl; mB = mC = mPi;
      mAB = sqrt(mm2);      mAC = sqrt(mp2); mBC = sqrt(mpippim2); 
    }
  else if(reso=="pippim") 
    { mA = mB = mPi; mC = mKl;
      mAB = sqrt(mpippim2); mAC = sqrt(mp2); mBC = sqrt(mm2);       
    }
  else if(reso=="k0spip") 
    { mA = mKl; mB = mC = mPi;
      mAB = sqrt(mp2);      mAC = sqrt(mm2); mBC = sqrt(mpippim2); 
    }

  double mD = mD0;
  Int_t power;
  double ZL=0.0;
  double a1=0., a2=0., a3=0.;

//One of the resonance's daughter particles' momentum in resonance rest frame
  double q = sqrt((((mAB*mAB-mA*mA-mB*mB)*(mAB*mAB-mA*mA-mB*mB)/4.0) - mA*mA*mB*mB)/(mAB*mAB));			        
  double q0 = sqrt((((mR*mR-mA*mA-mB*mB)*(mR*mR-mA*mA-mB*mB)/4.0) - mA*mA*mB*mB)/(mR*mR));				

//Spectator particle's momentum in resonance rest frame
  double p = sqrt( (((mD*mD-mAB*mAB-mC*mC)*(mD*mD-mAB*mAB-mC*mC)/4.0) - mAB*mAB*mC*mC)/(mD*mD));			
  double p0 = (((mD*mD-mR*mR-mC*mC)*(mD*mD-mR*mR-mC*mC)/4.0) - mR*mR*mC*mC)/(mD*mD);				      

  //double p = sqrt( (((mD*mD-mAB*mAB+mC*mC)*(mD*mD-mAB*mAB+mC*mC)/4.0) - mD*mD*mC*mC))/(mAB);			
  //double p0 = (((mD*mD-mR*mR+mC*mC)*(mD*mD-mR*mR+mC*mC)/4.0) - mD*mD*mC*mC)/(mR*mR);				      
  
  if ( p0>0 ) { p0=sqrt(p0); } else {p0=0.;}

//Spectator particle's momentum in resonance rest frame (According to P. Del Amo Sanchez et al.)
//  double p = sqrt( (mD*mD+mC*mC-mAB*mAB)*(mD*mD+mC*mC-mAB*mAB) - 4.0*mD*mD*mC*mC ) / (2.0*mAB) ;
//  double p0 = sqrt( (mD*mD+mC*mC-mR*mR)*(mD*mD+mC*mC-mR*mR) - 4.0*mD*mD*mC*mC ) / (2.0*mR) ;
//  if ( p0>0 ) { p0=sqrt(p0); } else {p0=0;}
  switch (spin) {
    case 0:
		fR=1.0;
		fD=1.0;
		power=1;
                ZL = 1.;
		break;
    case 1:
		fR=sqrt(1.0+R_r*R_r*q0*q0)/sqrt(1.0+R_r*R_r*q*q);
		fD=sqrt(1.0+R_D*R_D*p0*p0)/sqrt(1.0+R_D*R_D*p*p);
		power=3;
		ZL = (mAC*mAC-mBC*mBC+((mD*mD-mC*mC)*(mB*mB-mA*mA)/(mAB*mAB))) ;
		//ZL = (mBC*mBC-mAC*mAC+((mD*mD-mC*mC)*(mA*mA-mB*mB)/(mAB*mAB))) ;
		break;
    case 2:
		fR=sqrt(9.0 + 3.0*R_r*R_r*q0*q0 + R_r*R_r*q0*q0*R_r*R_r*q0*q0)/sqrt(9.0 + 3.0*R_r*R_r*q*q + R_r*R_r*q*q*R_r*R_r*q*q);
		fD=sqrt(9.0 + 3.0*R_D*R_D*p0*p0 + R_D*R_D*p0*p0*R_D*R_D*p0*p0)/sqrt(9.0 + 3.0*R_D*R_D*p*p + R_D*R_D*p*p*R_D*R_D*p*p);
		power=5;
                ZL = pow(mBC*mBC-mAC*mAC+(mD*mD-mC*mC)*(mA*mA-mB*mB)/(mAB*mAB),2)-
		1./3.*(mAB*mAB-2.*(mD*mD+mC*mC)+pow(mD*mD-mC*mC,2)/(mAB*mAB))*
		(mAB*mAB-2.*(mA*mA+mB*mB)+pow(mA*mA-mB*mB,2)/(mAB*mAB));
		break;

    default:
		std::cout << "Incorrect spin in Resonance.cc" << std::endl;
	}

/*
  switch (spin) {
    case 0:
		fR=1.0;
		fD=1.0;
		power=1;
                ZL = 1.;
		break;
    case 1:
		fR=sqrt(1.0+R_r*R_r*q0*q0)/sqrt(1.0+R_r*R_r*q*q);
		fD=sqrt(1.0+R_D*R_D*p0*p0)/sqrt(1.0+R_D*R_D*p*p);
		//fD=sqrt(1.0+(R_D+p0)*(R_D+p0))/sqrt(1.0+(R_D+p)*(R_D+p)) ;						//According to P. Del Amo Sanchez et al.
		power=3;
		ZL = (mAC*mAC-mBC*mBC+((mD*mD-mC*mC)*(mB*mB-mA*mA)/(mAB*mAB))) ;
		//ZL = (mBC*mBC-mAC*mAC-((mD*mD-mC*mC)*(mB*mB-mA*mA)/(mAB*mAB))) ;					//According to P. Del Amo Sanchez et al.
		break;
    case 2:
		fR=sqrt(9.0+3.0*R_r*R_r*q0*q0+R_r*R_r*q0*q0*R_r*R_r*q0*q0)/sqrt(9.0+3.0*R_r*R_r*q*q+R_r*R_r*q*q*R_r*R_r*q*q);
		fD=sqrt(9.0+3.0*R_D*R_D*p0*p0+R_D*R_D*p0*p0*R_D*R_D*p0*p0)/sqrt(9.0+3.0*R_D*R_D*p*p+R_D*R_D*p*p*R_D*R_D*p*p);
                //fD=sqrt(9.0+3.0*(R_D+p0)*(R_D+p0)+pow(R_D+p,4))/sqrt(9.0+3.0*(R_D+p)*(R_D+p)+pow(R_D+p,4)) ;		//According to P. Del Amo Sanchez et al.
		power=5;
		//a1 = mBC*mBC-mAC*mAC+((mD*mD-mC*mC)*(mA*mA-mB*mB)/(mAB*mAB)) ;						//According to P. Del Amo Sanchez et al.
                //a2 = mAB*mAB-2.*mD*mD-2.*mC*mC+((mD*mD-mC*mC)*(mD*mD-mC*mC))/mAB*mAB ;					//According to P. Del Amo Sanchez et al.
                //a3 = mAB*mAB-2.*mA*mA-2.*mB*mB+((mA*mA-mB*mB)*(mA*mA-mB*mB))/mAB*mAB ;					//According to P. Del Amo Sanchez et al.

		//ZL = a1*a1 - (a2*a3)/3. ;										//According to P. Del Amo Sanchez et al.

                ZL = pow(mBC*mBC-mAC*mAC+(mD*mD-mC*mC)*(mA*mA-mB*mB)/(mAB*mAB),2) - 1./3.*(mAB*mAB-2.*(mD*mD+mC*mC)+pow(mD*mD-mC*mC,2)/(mAB*mAB))*(mAB*mAB-2.*(mA*mA+mB*mB)+pow(mA*mA-mB*mB,2)/(mAB*mAB));
		break;

    default:
		std::cout << "Incorrect spin in Resonance.cc" << std::endl;
	}
*/  
// Compute the running width.
  double gammaAB = gammaR*pow(q/q0,power)*(mR/mAB)*fR*fR;

  matrixEl = (fR*fD*ZL/(mR*mR-mAB*mAB-TComplex(0.0,mR*gammaAB)));

  return matrixEl;

}
/*
vector<double> amplitude_kstar(double *phsp, double mR, string reso)
{

  double pi180inv = TMath::Pi() / 180.;

  TComplex matrixEl;

  double mm2 = phsp[0];
  double mp2 = phsp[1];

  const double mD0 = 1.86483;
  const double mKl = 0.49761;
  const double mPi = 0.13957;

  double R_r = 1.5; // "resonance radius" for Blatt-Weisskopf barrier factors.
  double R_D = 5.0; // "D meson radius" for Blatt-Weisskopf barrier factors.

  double fR,fD;

  double mpippim2 = mD0 * mD0 + mKl * mKl + (2 * mPi * mPi) - mp2 - mm2;

  double mA, mB, mC, mAB, mBC, mAC;

  if(reso=="k0spim")
    { mA = mKl; mB = mC = mPi;
      mAB = sqrt(mm2);      mAC = sqrt(mp2); mBC = sqrt(mpippim2);
    }
  else if(reso=="k0spip")
    { mA = mKl; mB = mC = mPi;
      mAB = sqrt(mp2);      mAC = sqrt(mm2); mBC = sqrt(mpippim2);
    }

  double mD = mD0;
  Int_t power;
  double ZL=0.0;
  double a1=0., a2=0., a3=0.;

//One of the resonance's daughter particles' momentum in resonance rest frame
  double q = sqrt((((mAB*mAB-mA*mA-mB*mB)*(mAB*mAB-mA*mA-mB*mB)/4.0) - mA*mA*mB*mB)/(mAB*mAB));
  double q0 = sqrt((((mR*mR-mA*mA-mB*mB)*(mR*mR-mA*mA-mB*mB)/4.0) - mA*mA*mB*mB)/(mR*mR));

//Spectator particle's momentum in resonance rest frame
  double p = sqrt( (((mD*mD-mAB*mAB-mC*mC)*(mD*mD-mAB*mAB-mC*mC)/4.0) - mAB*mAB*mC*mC)/(mD*mD));
  double p0 = (((mD*mD-mR*mR-mC*mC)*(mD*mD-mR*mR-mC*mC)/4.0) - mR*mR*mC*mC)/(mD*mD);

  if ( p0>0 ) { p0=sqrt(p0); } else {p0=0.;}

  fR=sqrt(1.0+R_r*R_r*q0*q0)/sqrt(1.0+R_r*R_r*q*q);
  fD=sqrt(1.0+R_D*R_D*p0*p0)/sqrt(1.0+R_D*R_D*p*p);
  power=3;
  ZL = (mAC*mAC-mBC*mBC+((mD*mD-mC*mC)*(mB*mB-mA*mA)/(mAB*mAB))) ;

// Compute the running width.
  double gammaAB = pow(q/q0,power)*(mR/mAB)*fR*fR;

  std::vector<double> kstr_temp ;

  kstr_temp.push_back(fR*fD*ZL) ;
  kstr_temp.push_back(mR*gammaAB) ;

  return kstr_temp;

}
*/

vector<double> amplitude_Kmatrix(double *phsp, int jj)
{ 

        const double mD0 = 1.86483;
        const double mKl = 0.49761;
        const double mPi = 0.13957;

        double mm2 = phsp[0]; double mp2 = phsp[1];

        //cout<<"mm2: "<<mm2<<",\t mp2: "<<mp2<<endl;

        double mAB2 = mD0 * mD0 + mKl * mKl + 2 * mPi * mPi - mp2 - mm2;			//m_pippim^2
        double mAB = abs(sqrt(mAB2));

	Double_t s = mAB*mAB;

	// Define the complex coupling constants
	Double_t g[5][5]; // g[Physical pole]Decay channel]

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

	// Define masses of the physical poles (in GeV)
	Double_t ma[5];

	ma[0]=0.651;
	ma[1]=1.20360;
	ma[2]=1.55817;
	ma[3]=1.21000;
	ma[4]=1.82206;

	// Define variables
	TComplex n11,n12,n13,n14,n15,n21,n22,n23,n24,n25,n31,n32,n33,n34,n35,n41,n42,n43,n44,n45,n51,n52,n53,n54,n55;
	Double_t rho1sq,rho2sq,rho4sq,rho5sq;
	TComplex rho1,rho2,rho3,rho4,rho5;
	TComplex rho[5];
	TComplex pole,SVT,Adler;
	TComplex det;
	TComplex i[5][5];
	Double_t f[5][5];

	// pi+, K+, eta, and eta' PDG masses
	Double_t mpi=0.13957;
	Double_t mK=0.493677;
	Double_t meta=0.54775;
	Double_t metap=0.95778;

	// Init matrices and vectors with zeros
	TComplex K[5][5];
	for(Int_t k=0;k<5;k++) {
		for(Int_t l=0;l<5;l++) {
			i[k][l]=TComplex(0.,0.);
			K[k][l]=TComplex(0.,0.);
			f[k][l]=0.;
		}
		rho[k]=0.;
	}

	// Fill scattering data values
	Double_t s_scatt=-3.92637;
	Double_t sa=1.0;
	Double_t sa_0=-0.15;

	// f_scattering
	f[0][0]=0.23399;
	f[0][1]=0.15044;
	f[0][2]=-0.20545;
	f[0][3]=0.32825;
	f[0][4]=0.35412;

	f[1][0]=f[0][1];
	f[2][0]=f[0][2];
	f[3][0]=f[0][3];
	f[4][0]=f[0][4];

	// Compute phase space factors
	rho1sq=(1.0-(pow((mpi+mpi),2)/s));
	if(rho1sq >=0.) {
		rho1=TComplex(sqrt(rho1sq),0.);
	}
	else{
		rho1=TComplex(0.,sqrt(-rho1sq));
	}
	rho[0]=rho1;

	rho2sq=(1.0-(pow((mK+mK),2)/s));
	if(rho2sq >=0.) {
		rho2=TComplex(sqrt(rho2sq),0.);
	}
	else{
		rho2=TComplex(0.,sqrt(-rho2sq));
	}

	rho[1]=rho2;

	rho3=TComplex(0.,0.);

	if(s<=1) {
		Double_t real = 1.2274+0.00370909/(s*s) - (0.111203)/(s) - 6.39017*s +16.8358*s*s - 21.8845*s*s*s + 11.3153*s*s*s*s;
		Double_t cont32=sqrt(1.0-(16.0*mpi*mpi));
		rho3=TComplex(cont32*real,0.);
	}
	else{
		rho3=TComplex(sqrt(1.0-(16.0*mpi*mpi/s)),0.);
	}
	rho[2]=rho3;

	rho4sq=(1.0-(pow((meta+meta),2)/s));
	if(rho4sq>=0.) {
		rho4=TComplex(sqrt(rho4sq),0.);
	}
	else{
		rho4=TComplex(0.,sqrt(-rho4sq));
	}
	rho[3]=rho4;

	rho5sq=(1.0-(pow((meta+metap),2)/s));
	if(rho5sq >=0.) {
		rho5=TComplex(sqrt(rho5sq),0.);
	}
	else{
		rho5=TComplex(0.,sqrt(-rho5sq));
	}
	rho[4]=rho5;

	// Sum over the poles
	for(Int_t k=0;k<5;k++) {
		for(Int_t l=0;l<5;l++) {
			for (Int_t pole_index=0;pole_index<5;pole_index++) {
				Double_t A=g[pole_index][k]*g[pole_index][l];
				Double_t B=ma[pole_index]*ma[pole_index]-s;
				K[k][l]=K[k][l]+TComplex(A/B,0.);
			}
		}
	}

	for(Int_t k=0;k<5;k++) {
		for(Int_t l=0;l<5;l++) {
			Double_t C=f[k][l]*(1.0-s_scatt);
			Double_t D=(s-s_scatt);
			K[k][l]=K[k][l]+TComplex(C/D,0.);
		}
	}

	for(Int_t k=0;k<5;k++) {
		for(Int_t l=0;l<5;l++) {
			Double_t E=(s-(sa*mpi*mpi*0.5))*(1.0-sa_0);
			Double_t F=(s-sa_0);
			K[k][l]=K[k][l]*TComplex(E/F,0.);
		}
	}

	n11=TComplex(1.,0.)-TComplex(0.,1.)*K[0][0]*rho[0];
	n12=TComplex(0.,0.)-TComplex(0.,1.)*K[0][1]*rho[1];
	n13=TComplex(0.,0.)-TComplex(0.,1.)*K[0][2]*rho[2];
	n14=TComplex(0.,0.)-TComplex(0.,1.)*K[0][3]*rho[3];
	n15=TComplex(0.,0.)-TComplex(0.,1.)*K[0][4]*rho[4];

	n21=TComplex(0.,0.)-TComplex(0.,1.)*K[1][0]*rho[0];
	n22=TComplex(1.,0.)-TComplex(0.,1.)*K[1][1]*rho[1];
	n23=TComplex(0.,0.)-TComplex(0.,1.)*K[1][2]*rho[2];
	n24=TComplex(0.,0.)-TComplex(0.,1.)*K[1][3]*rho[3];
	n25=TComplex(0.,0.)-TComplex(0.,1.)*K[1][4]*rho[4];

	n31=TComplex(0.,0.)-TComplex(0.,1.)*K[2][0]*rho[0];
	n32=TComplex(0.,0.)-TComplex(0.,1.)*K[2][1]*rho[1];
	n33=TComplex(1.,0.)-TComplex(0.,1.)*K[2][2]*rho[2];
	n34=TComplex(0.,0.)-TComplex(0.,1.)*K[2][3]*rho[3];
	n35=TComplex(0.,0.)-TComplex(0.,1.)*K[2][4]*rho[4];

	n41=TComplex(0.,0.)-TComplex(0.,1.)*K[3][0]*rho[0];
	n42=TComplex(0.,0.)-TComplex(0.,1.)*K[3][1]*rho[1];
	n43=TComplex(0.,0.)-TComplex(0.,1.)*K[3][2]*rho[2];
	n44=TComplex(1.,0.)-TComplex(0.,1.)*K[3][3]*rho[3];
	n45=TComplex(0.,0.)-TComplex(0.,1.)*K[3][4]*rho[4];

	n51=TComplex(0.,0.)-TComplex(0.,1.)*K[4][0]*rho[0];
	n52=TComplex(0.,0.)-TComplex(0.,1.)*K[4][1]*rho[1];
	n53=TComplex(0.,0.)-TComplex(0.,1.)*K[4][2]*rho[2];
	n54=TComplex(0.,0.)-TComplex(0.,1.)*K[4][3]*rho[3];
	n55=TComplex(1.,0.)-TComplex(0.,1.)*K[4][4]*rho[4];

	// Compute the determinant
	det = (n15*n24*n33*n42*n51 - n14*n25*n33*n42*n51 - n15*n23*n34*n42*n51 +
			n13*n25*n34*n42*n51 + n14*n23*n35*n42*n51 - n13*n24*n35*n42*n51 -
			n15*n24*n32*n43*n51 + n14*n25*n32*n43*n51 + n15*n22*n34*n43*n51 -
			n12*n25*n34*n43*n51 - n14*n22*n35*n43*n51 + n12*n24*n35*n43*n51 +
			n15*n23*n32*n44*n51 - n13*n25*n32*n44*n51 - n15*n22*n33*n44*n51 +
			n12*n25*n33*n44*n51 + n13*n22*n35*n44*n51 - n12*n23*n35*n44*n51 -
			n14*n23*n32*n45*n51 + n13*n24*n32*n45*n51 + n14*n22*n33*n45*n51 -
			n12*n24*n33*n45*n51 - n13*n22*n34*n45*n51 + n12*n23*n34*n45*n51 -
			n15*n24*n33*n41*n52 + n14*n25*n33*n41*n52 + n15*n23*n34*n41*n52 -
			n13*n25*n34*n41*n52 - n14*n23*n35*n41*n52 + n13*n24*n35*n41*n52 +
			n15*n24*n31*n43*n52 - n14*n25*n31*n43*n52 - n15*n21*n34*n43*n52 +
			n11*n25*n34*n43*n52 + n14*n21*n35*n43*n52 - n11*n24*n35*n43*n52 -
			n15*n23*n31*n44*n52 + n13*n25*n31*n44*n52 + n15*n21*n33*n44*n52 -
			n11*n25*n33*n44*n52 - n13*n21*n35*n44*n52 + n11*n23*n35*n44*n52 +
			n14*n23*n31*n45*n52 - n13*n24*n31*n45*n52 - n14*n21*n33*n45*n52 +
			n11*n24*n33*n45*n52 + n13*n21*n34*n45*n52 - n11*n23*n34*n45*n52 +
			n15*n24*n32*n41*n53 - n14*n25*n32*n41*n53 - n15*n22*n34*n41*n53 +
			n12*n25*n34*n41*n53 + n14*n22*n35*n41*n53 - n12*n24*n35*n41*n53 -
			n15*n24*n31*n42*n53 + n14*n25*n31*n42*n53 + n15*n21*n34*n42*n53 -
			n11*n25*n34*n42*n53 - n14*n21*n35*n42*n53 + n11*n24*n35*n42*n53 +
			n15*n22*n31*n44*n53 - n12*n25*n31*n44*n53 - n15*n21*n32*n44*n53 +
			n11*n25*n32*n44*n53 + n12*n21*n35*n44*n53 - n11*n22*n35*n44*n53 -
			n14*n22*n31*n45*n53 + n12*n24*n31*n45*n53 + n14*n21*n32*n45*n53 -
			n11*n24*n32*n45*n53 - n12*n21*n34*n45*n53 + n11*n22*n34*n45*n53 -
			n15*n23*n32*n41*n54 + n13*n25*n32*n41*n54 + n15*n22*n33*n41*n54 -
			n12*n25*n33*n41*n54 - n13*n22*n35*n41*n54 + n12*n23*n35*n41*n54 +
			n15*n23*n31*n42*n54 - n13*n25*n31*n42*n54 - n15*n21*n33*n42*n54 +
			n11*n25*n33*n42*n54 + n13*n21*n35*n42*n54 - n11*n23*n35*n42*n54 -
			n15*n22*n31*n43*n54 + n12*n25*n31*n43*n54 + n15*n21*n32*n43*n54 -
			n11*n25*n32*n43*n54 - n12*n21*n35*n43*n54 + n11*n22*n35*n43*n54 +
			n13*n22*n31*n45*n54 - n12*n23*n31*n45*n54 - n13*n21*n32*n45*n54 +
			n11*n23*n32*n45*n54 + n12*n21*n33*n45*n54 - n11*n22*n33*n45*n54 +
			n14*n23*n32*n41*n55 - n13*n24*n32*n41*n55 - n14*n22*n33*n41*n55 +
			n12*n24*n33*n41*n55 + n13*n22*n34*n41*n55 - n12*n23*n34*n41*n55 -
			n14*n23*n31*n42*n55 + n13*n24*n31*n42*n55 + n14*n21*n33*n42*n55 -
			n11*n24*n33*n42*n55 - n13*n21*n34*n42*n55 + n11*n23*n34*n42*n55 +
			n14*n22*n31*n43*n55 - n12*n24*n31*n43*n55 - n14*n21*n32*n43*n55 +
			n11*n24*n32*n43*n55 + n12*n21*n34*n43*n55 - n11*n22*n34*n43*n55 -
			n13*n22*n31*n44*n55 + n12*n23*n31*n44*n55 + n13*n21*n32*n44*n55 -
			n11*n23*n32*n44*n55 - n12*n21*n33*n44*n55 + n11*n22*n33*n44*n55);


	// The 1st row of the inverse matrix {(I-iKp)^-1}_0j
	i[0][0] = (n25*n34*n43*n52 -
			n24*n35*n43*n52 - n25*n33*n44*n52 + n23*n35*n44*n52 +
			n24*n33*n45*n52 - n23*n34*n45*n52 - n25*n34*n42*n53 +
			n24*n35*n42*n53 + n25*n32*n44*n53 - n22*n35*n44*n53 -
			n24*n32*n45*n53 + n22*n34*n45*n53 + n25*n33*n42*n54 -
			n23*n35*n42*n54 - n25*n32*n43*n54 + n22*n35*n43*n54 +
			n23*n32*n45*n54 - n22*n33*n45*n54 - n24*n33*n42*n55 +
			n23*n34*n42*n55 + n24*n32*n43*n55 - n22*n34*n43*n55 -
			n23*n32*n44*n55 + n22*n33*n44*n55)/det;

	i[0][1] = (-n15*n34*n43*n52 +
			n14*n35*n43*n52 + n15*n33*n44*n52 - n13*n35*n44*n52 -
			n14*n33*n45*n52 + n13*n34*n45*n52 + n15*n34*n42*n53 -
			n14*n35*n42*n53 - n15*n32*n44*n53 + n12*n35*n44*n53 +
			n14*n32*n45*n53 - n12*n34*n45*n53 - n15*n33*n42*n54 +
			n13*n35*n42*n54 + n15*n32*n43*n54 - n12*n35*n43*n54 -
			n13*n32*n45*n54 + n12*n33*n45*n54 + n14*n33*n42*n55 -
			n13*n34*n42*n55 - n14*n32*n43*n55 + n12*n34*n43*n55 +
			n13*n32*n44*n55 - n12*n33*n44*n55)/det;

	i[0][2] = (n15*n24*n43*n52 -
			n14*n25*n43*n52 - n15*n23*n44*n52 + n13*n25*n44*n52 +
			n14*n23*n45*n52 - n13*n24*n45*n52 - n15*n24*n42*n53 +
			n14*n25*n42*n53 + n15*n22*n44*n53 - n12*n25*n44*n53 -
			n14*n22*n45*n53 + n12*n24*n45*n53 + n15*n23*n42*n54 -
			n13*n25*n42*n54 - n15*n22*n43*n54 + n12*n25*n43*n54 +
			n13*n22*n45*n54 - n12*n23*n45*n54 - n14*n23*n42*n55 +
			n13*n24*n42*n55 + n14*n22*n43*n55 - n12*n24*n43*n55 -
			n13*n22*n44*n55 + n12*n23*n44*n55)/det;

	i[0][3] = (-n15*n24*n33*n52 +
			n14*n25*n33*n52 + n15*n23*n34*n52 - n13*n25*n34*n52 -
			n14*n23*n35*n52 + n13*n24*n35*n52 + n15*n24*n32*n53 -
			n14*n25*n32*n53 - n15*n22*n34*n53 + n12*n25*n34*n53 +
			n14*n22*n35*n53 - n12*n24*n35*n53 - n15*n23*n32*n54 +
			n13*n25*n32*n54 + n15*n22*n33*n54 - n12*n25*n33*n54 -
			n13*n22*n35*n54 + n12*n23*n35*n54 + n14*n23*n32*n55 -
			n13*n24*n32*n55 - n14*n22*n33*n55 + n12*n24*n33*n55 +
			n13*n22*n34*n55 - n12*n23*n34*n55)/det;

	i[0][4] = (n15*n24*n33*n42 -
			n14*n25*n33*n42 - n15*n23*n34*n42 + n13*n25*n34*n42 +
			n14*n23*n35*n42 - n13*n24*n35*n42 - n15*n24*n32*n43 +
			n14*n25*n32*n43 + n15*n22*n34*n43 - n12*n25*n34*n43 -
			n14*n22*n35*n43 + n12*n24*n35*n43 + n15*n23*n32*n44 -
			n13*n25*n32*n44 - n15*n22*n33*n44 + n12*n25*n33*n44 +
			n13*n22*n35*n44 - n12*n23*n35*n44 - n14*n23*n32*n45 +
			n13*n24*n32*n45 + n14*n22*n33*n45 - n12*n24*n33*n45 -
			n13*n22*n34*n45 + n12*n23*n34*n45)/det;

/*
        det = n11*(n22*(n33*n44-n34*n43) -n23*(n32*n44-n34*n42) +n24*(n32*n43-n33*n42)) - n12*(n21*(n33*n44-n34*n43) - n23*(n31*n44-n34*n41) + n24*(n31*n43-n33*n41)) + n13*(n21*(n32*n44-n34*n42) - n22*(n31*n44-n34*n41) + n24*(n31*n42-n32*n41)) - n14*(n21*(n32*n43-n33*n42) - n22*(n31*n43-n33*n41) + n23*(n31*n42-n32*n41)) ;

	i[0][0] = (n22*(n33*n44-n34*n43) - n23*(n32*n44-n34*n42) + n24*(n32*n43-n33*n42))/det ;

	i[0][1] = -1.0*( n12*(n33*n44-n34*n43) - n13*(n32*n44-n34*n42) + n14*(n32*n43-n33*n42) )/det ;

	i[0][2] = (n12*(n23*n44-n24*n43) - n13*(n22*n44-n24*n42) + n14*(n22*n43-n23*n42))/det ;

	i[0][3] = -1.0*( n12*(n23*n34-n24*n33) - n13*(n22*n34-n24*n32) + n14*(n22*n33-n23*n32) )/det ;
*/
//--------------------------------------------------------------------------------------------------------------------------------
	double _s0prod = -0.07;
        std::vector<double> U1j;

        if(det.Re()==0. && det.Im()==0.) 
	{
		for(Int_t j=0;j<5;j++) 
		{ 
			U1j.push_back(9999.); 
			U1j.push_back(9999.);
		}
        }

        else 
	{
		for(Int_t j=0;j<5;j++) 
		{ 
			U1j.push_back(i[0][j].Re()); 
			U1j.push_back(i[0][j].Im());
		}
        }

        for(Int_t pole_index=0;pole_index<5;pole_index++) 
	{
		U1j.push_back(ma[pole_index]*ma[pole_index]-s);
        }
	
        U1j.push_back((1-_s0prod)/(s-_s0prod));

	return U1j;

}

TComplex amplitude_LASS(double *phsp, double mR, double gammaR, map<string, double> scalarparam, string reso)
{
        double mab2=0.0;
        double mm2 = phsp[0];
        double mp2 = phsp[1];
        if(reso == "k0spim") mab2 = mm2;
        else if(reso == "k0spip") mab2 = mp2;

	double pi180inv = TMath::Pi() / 180.;

	double s = mab2;

	//cout<<"m_KPi^{2}(RS/WS): "<<s<<endl;

        const double mD0 = 1.86483;
        const double mKl = 0.49761;
        const double mPi = 0.13957;

	double _a = scalarparam["a"];
	double _r = scalarparam["r"];
	double _R = scalarparam["R"]; // Amplitude magnitude of the resonant term
	double _F = scalarparam["F"]; // Amplitude magnitude of the non-resonant term
	double _phiR = scalarparam["phiR"]; // Phase of the resonant term
	double _phiF = scalarparam["phiF"]; // Phase of the non-resonant term

	double fR=1.0; // K*0(1430) has spin zero
	int power=1; // Power is 1 for spin zero

	double mAB=abs(sqrt(mab2)); // (_p4_d1+_p4_d2).mass();

	double mA=mKl; // _p4_d1.mass();
	double mB=mPi; // _p4_d2.mass();
        double mC=mPi;
        double mD = mD0;

	double pAB=sqrt( (((mAB*mAB-mA*mA-mB*mB)*(mAB*mAB-mA*mA-mB*mB)/4.0) - mA*mA*mB*mB)/(mAB*mAB));
	double q=pAB;

	double pR=sqrt( (((mR*mR-mA*mA-mB*mB)*(mR*mR-mA*mA-mB*mB)/4.0) - mA*mA*mB*mB)/(mR*mR));
        double q0=pR;

//         double q = sqrt( (((mD*mD-mAB*mAB-mC*mC)*(mD*mD-mAB*mAB-mC*mC)/4.0) - mAB*mAB*mC*mC)/(mD*mD));			
//         double q0 = (((mD*mD-mR*mR-mC*mC)*(mD*mD-mR*mR-mC*mC)/4.0) - mR*mR*mC*mC)/(mD*mD);				        
//        if ( q0>0 ) { q0=sqrt(q0); } else {q0=0;}


	// Running width.
	double g = gammaR*pow(q/q0,power)*(mR/mAB)*fR*fR;

	TComplex propagator_relativistic_BreitWigner = 1./(mR*mR - mAB*mAB - TComplex(0.,mR*g));

	// Non-resonant phase shift
	Double_t cot_deltaF = 1.0/(_a*q) + 0.5*_r*q;
	Double_t qcot_deltaF = 1.0/_a + 0.5*_r*q*q;

	// Compute resonant part
	TComplex expi2deltaF = TComplex(qcot_deltaF, q)/ TComplex(qcot_deltaF, -q);

	TComplex resonant_term_T = _R * TComplex(cos(_phiR + 2 * _phiF), sin(_phiR + 2 * _phiF)) * propagator_relativistic_BreitWigner * mR * gammaR * mR / q0 * expi2deltaF;

//        TComplex resonant_term_T = _R * TComplex(cos(_phiR + 2 * _phiF), sin(_phiR + 2 * _phiF)) * propagator_relativistic_BreitWigner * expi2deltaF * ((mR*g*sin(_phiR)) + (mR*mR - s)*cos(_phiR)) ; 		//From 2018 paper

	// Compute non-resonant part
	TComplex non_resonant_term_F = _F * TComplex(cos(_phiF), sin(_phiF)) * (cos(_phiF) + cot_deltaF * sin(_phiF)) * sqrt(s) / TComplex(qcot_deltaF, -q);

//	TComplex non_resonant_term_F = _F * TComplex(cos(_phiF), sin(_phiF)) * (cos(_phiF) + cot_deltaF * sin(_phiF)) * TComplex(cot_deltaF, 1.0) / (1.0 + pow(cot_deltaF, 2));			//From 2018 paper

	// Add non-resonant and resonant terms
	TComplex LASS_contribution = non_resonant_term_F + resonant_term_T;

	return LASS_contribution;

}



bool inDalitz_01_02(Double_t x, Double_t y) {
	static const Double_t PI_mass = 0.13957;
	static const Double_t K0_mass = 0.49761;
	static const Double_t mD0 = 1.86483;

	double msquared01 = x; 
	double msquared02 = y; 
	double msquared12 = mD0*mD0 + K0_mass*K0_mass + 2*PI_mass*PI_mass - msquared01 - msquared02;

	Double_t local_ma = K0_mass;
	Double_t local_mb = PI_mass;
	Double_t local_mc = PI_mass;

	Double_t local_xmin = pow(local_ma + local_mb,2);
	Double_t local_xmax = pow(mD0 - local_mc,2);

	Double_t ebab = (x - local_ma*local_ma + local_mb*local_mb)/(2.0*sqrt(x));
	Double_t ecab = (mD0*mD0 - x - local_mc*local_mc)/(2.0*sqrt(x));

	Double_t yhi = pow(ebab+ecab,2) - pow( sqrt(ebab*ebab-local_mb*local_mb)-sqrt(ecab*ecab-local_mc*local_mc) ,2);
	Double_t ylo = pow(ebab+ecab,2) - pow( sqrt(ebab*ebab-local_mb*local_mb)+sqrt(ecab*ecab-local_mc*local_mc) ,2);

	bool inDal = false;

	if ((local_xmin <= x) && (x <= local_xmax) && (ylo <= msquared12) && (msquared12 <= yhi)) { inDal = true; }

	return inDal;
}


int main()
{

  bool GenFit = false;					//"true" if fitting a generated sample; "false" for data
  bool acceptance = true;				//"true" if applying efficiency with SMC rec. events using MC integration

  double efficiency;
  int ii, jj;
  int normnumberReal=0, datanumberReal=0;
  int normnumberGen=0, datanumberGen=0;

  ifstream fin1;
  fin1.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/MassWidth_PDG.dat");

  while (1) {  
               if (!fin1.good()) break;
	       fin1 >> mode >> mass >> width ;
               massmap.insert(pair<string, double>(mode,mass));
               widthmap.insert(pair<string, double>(mode,width));
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

for(int k=0; k<80000; k++)
   { TAG[k]=500;
   }

if(GenFit == false)
{
  double Klpx, Klpy, Klpz, KlE, pippx, pippy, pippz, pipE, pimpx, pimpy, pimpz, pimE;
  int tag=0;

  TFile* fileIn1  = new TFile("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/klpipivsTagsData_LS.root");
  TTree* t1;
  t1 = (TTree*)fileIn1->Get("Mom");
  t1->SetBranchAddress("tag",&tag);
  t1->SetBranchAddress("Kl_px",&Klpx);
  t1->SetBranchAddress("Kl_py",&Klpy);
  t1->SetBranchAddress("Kl_pz",&Klpz);
  t1->SetBranchAddress("Kl_E",&KlE);
  t1->SetBranchAddress("Pip_px",&pippx);
  t1->SetBranchAddress("Pip_py",&pippy);
  t1->SetBranchAddress("Pip_pz",&pippz);
  t1->SetBranchAddress("Pip_E",&pipE);
  t1->SetBranchAddress("Pim_px",&pimpx);
  t1->SetBranchAddress("Pim_py",&pimpy);
  t1->SetBranchAddress("Pim_pz",&pimpz);
  t1->SetBranchAddress("Pim_E",&pimE);

  int cc=0;
  for(int i=0; i<t1->GetEntries(); i++)
     {      t1->GetEntry(i);

       	    TLorentzVector K0L(Klpx, Klpy, Klpz, KlE);
            TLorentzVector PIP(pippx, pippy, pippz, pipE);
            TLorentzVector PIM(pimpx, pimpy, pimpz, pimE);
	    if(inDalitz_01_02((K0L+PIM).M2(), (K0L+PIP).M2()) == true)
              { sK0SPIP[cc] = (K0L+PIP).M2(); //cout<<"sK0SPIP["<<i<<"]: "<<sK0SPIP[i]<<endl;
                sK0SPIM[cc] = (K0L+PIM).M2();
                TAG[cc] = tag;
                K0Lmom[cc] = K0L.Rho();
                PIPmom[cc] = PIP.Rho();
		PIMmom[cc] = PIM.Rho();
                cc++;
              }
     }
       
  datanumberReal = cc;
  cout<<"datanumber: "<<datanumberReal<<endl;

}

  ifstream fin4;
  fin4.open("/Users/anita/AmplitudeFitter/D0toK0Lpipi/AmplitudeFit/LASSFittedParameterInfoBelle2018.dat");

  while (1) {  
               if (!fin4.good()) break;
	       fin4 >> param >> value ;
               lassparammap.insert(pair<string, double>(param,value));
            }

  Int_t norm_tag;
  Double_t mplus2, mminus2, pk0l, ppip, ppim;
  std::map<string, Double_t> mulfacbw_re;
  std::map<string, Double_t> mulfacbw_im;
  std::vector<Double_t> u1j_re;
  std::vector<Double_t> u1j_im;
  std::vector<Double_t> bb;
  std::map<string, Double_t> lass_contribution_re;
  std::map<string, Double_t> lass_contribution_im;
  
  Int_t tag;
  Double_t Mplus2, Mminus2, Pk0l, Ppip, Ppim;
  std::map<string, Double_t> Mulfacbw_re;
  std::map<string, Double_t> Mulfacbw_im;
  std::vector<Double_t> U1j_re;
  std::vector<Double_t> U1j_im;
  std::vector<Double_t> BB;
  std::map<string, Double_t> LASS_contribution_re;
  std::map<string, Double_t> LASS_contribution_im;
  TComplex Asum(0.0, 0.0) ;
  double Asum_re[5], Asum_im[5] ;

  Double_t Pmplus2, Pmminus2;
  std::map<string, Double_t> Pmulfacbw_re;
  std::map<string, Double_t> Pmulfacbw_im;
  std::vector<Double_t> Pu1j_re;
  std::vector<Double_t> Pu1j_im;
  std::vector<Double_t> Pbb;
  std::map<string, Double_t> Plass_contribution_re;
  std::map<string, Double_t> Plass_contribution_im;

  Double_t P2mpippim2;

  Double_t P2mplus2, P2mminus2;
  std::map<string, Double_t> P2mulfacbw_re;
  std::map<string, Double_t> P2mulfacbw_im;
  std::vector<Double_t> P2u1j_re;
  std::vector<Double_t> P2u1j_im;
  std::vector<Double_t> P2bb;
  std::map<string, Double_t> P2lass_contribution_re;
  std::map<string, Double_t> P2lass_contribution_im;

  TFile* fileout = new TFile("check.root", "RECREATE");

  TTree* tr = new TTree("datasample","");
  tr->Branch("tag",&tag,"tag/I");
  tr->Branch("Mplus2",&Mplus2,"Mplus2/D");
  tr->Branch("Mminus2",&Mminus2,"Mminus2/D");
  tr->Branch("Mulfacbw_re",&Mulfacbw_re);
  tr->Branch("Mulfacbw_im",&Mulfacbw_im);
  tr->Branch("U1j_re", &U1j_re);
  tr->Branch("U1j_im", &U1j_im);
  tr->Branch("BB", &BB);
  tr->Branch("LASS_contribution_re", &LASS_contribution_re);
  tr->Branch("LASS_contribution_im", &LASS_contribution_im);
  tr->Branch("Asum_re",Asum_re,"Asum_re[5]/D");
  tr->Branch("Asum_im",Asum_im,"Asum_im[5]/D");
  tr->Branch("P_K0L",&Pk0l,"Pk0l/D");
  tr->Branch("P_PIP",&Ppip,"Ppip/D");
  tr->Branch("P_PIM",&Ppim,"Ppim/D");
  

  double test_xx[2];

  TComplex pass(0.0, 0.0);
  std::vector<double> passKM;
  //std::vector<double> pass_kstar;

  double u1j_re_limit[5][2], u1j_im_limit[5][2] ;
  u1j_re_limit[0][0] =  0.   ; u1j_re_limit[0][1] = 1.    ;
  u1j_re_limit[1][0] = -0.29 ; u1j_re_limit[1][1] = 0.12  ;
  u1j_re_limit[2][0] = -0.17 ; u1j_re_limit[2][1] = 0.065 ;
  u1j_re_limit[3][0] = -0.66 ; u1j_re_limit[3][1] = 0.1   ;
  u1j_re_limit[4][0] = -1.36 ; u1j_re_limit[4][1] = 0.18  ;

  u1j_im_limit[0][0] = -0.58  ; u1j_im_limit[0][1] = 0.58 ;
  u1j_im_limit[1][0] =  0.00  ; u1j_im_limit[1][1] = 0.28 ;
  u1j_im_limit[2][0] = -0.135 ; u1j_im_limit[2][1] = 0.10 ;
  u1j_im_limit[3][0] = -0.13  ; u1j_im_limit[3][1] = 0.40 ;
  u1j_im_limit[4][0] = -0.36  ; u1j_im_limit[4][1] = 0.80 ;


  amp_func bw_instance;
  //unicorn obj;
  //obj.horse = "donkey";
  //obj.animal();


  //for(Int_t j=0; j<datanumberReal; j++)
  for(Int_t j=0 ; j<100 ; j++)
     { 
	Mplus2 = sK0SPIP[j];
        Mminus2 = sK0SPIM[j];
        Pk0l = K0Lmom[j];
        Ppip = PIPmom[j];
        Ppim = PIMmom[j];
        tag = TAG[j];
        Mulfacbw_re.clear();
        Mulfacbw_im.clear();
        LASS_contribution_re.clear();
        LASS_contribution_im.clear();
        U1j_re.clear();
        U1j_im.clear();
        BB.clear();
        passKM.clear();

        test_xx[0] = sK0SPIM[j]; test_xx[1] = sK0SPIP[j];
        if(inDalitz_01_02(test_xx[0], test_xx[1]) == false) continue;

        //pass = amplitude_BW(test_xx, massmap["Rho770"], widthmap["Rho770"], 1, "pippim");
	pass = bw_instance.amplitude_BW(test_xx, massmap["Rho770"], widthmap["Rho770"], 1, "pippim");
        Mulfacbw_re.insert(pair<string, double>("Rho770", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("Rho770", pass.Im()));

        pass = amplitude_BW(test_xx, massmap["Kstar892minus"], widthmap["Kstar892minus"], 1, "k0spim");
        Mulfacbw_re.insert(pair<string, double>("Kstar892minus", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("Kstar892minus", pass.Im()));

        pass = amplitude_BW(test_xx, massmap["K2star1430minus"], widthmap["K2star1430minus"], 2, "k0spim");
        Mulfacbw_re.insert(pair<string, double>("K2star1430minus", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("K2star1430minus", pass.Im()));

        pass = amplitude_BW(test_xx, massmap["ftwo1270"], widthmap["ftwo1270"], 2, "pippim");
        Mulfacbw_re.insert(pair<string, double>("ftwo1270", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("ftwo1270", pass.Im()));

        pass = amplitude_BW(test_xx, massmap["Rho1450"], widthmap["Rho1450"], 1, "pippim");
        Mulfacbw_re.insert(pair<string, double>("Rho1450", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("Rho1450", pass.Im()));

        pass = amplitude_BW(test_xx, massmap["Kstar892plus"], widthmap["Kstar892plus"], 1, "k0spip");
        Mulfacbw_re.insert(pair<string, double>("Kstar892plus", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("Kstar892plus", pass.Im()));

        pass = amplitude_BW(test_xx, massmap["Kstar1680minus"], widthmap["Kstar1680minus"], 1, "k0spim");
        Mulfacbw_re.insert(pair<string, double>("Kstar1680minus", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("Kstar1680minus", pass.Im()));

        pass = amplitude_BW(test_xx, massmap["Kstar1680plus"], widthmap["Kstar1680plus"], 1, "k0spip");
        Mulfacbw_re.insert(pair<string, double>("Kstar1680plus", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("Kstar1680plus", pass.Im()));

        pass = amplitude_BW(test_xx, massmap["Omega782"], widthmap["Omega782"], 1, "pippim");
        Mulfacbw_re.insert(pair<string, double>("Omega782", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("Omega782", pass.Im()));

        pass = amplitude_BW(test_xx, massmap["Kstar1410minus"], widthmap["Kstar1410minus"], 1, "k0spim");
        Mulfacbw_re.insert(pair<string, double>("Kstar1410minus", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("Kstar1410minus", pass.Im()));

        pass = amplitude_BW(test_xx, massmap["K2star1430plus"], widthmap["K2star1430plus"], 2, "k0spip");
        Mulfacbw_re.insert(pair<string, double>("K2star1430plus", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("K2star1430plus", pass.Im()));

        pass = amplitude_BW(test_xx, massmap["Kstar1410plus"], widthmap["Kstar1410plus"], 1, "k0spip");
        Mulfacbw_re.insert(pair<string, double>("Kstar1410plus", pass.Re()));
        Mulfacbw_im.insert(pair<string, double>("Kstar1410plus", pass.Im()));

        pass = amplitude_LASS(test_xx, massmap["K0star1430minus"], widthmap["K0star1430minus"], lassparammap, "k0spim");
        LASS_contribution_re.insert(pair<string, double>("K0star1430minus", pass.Re()));
        LASS_contribution_im.insert(pair<string, double>("K0star1430minus", pass.Im()));

        pass = amplitude_LASS(test_xx, massmap["K0star1430plus"], widthmap["K0star1430plus"], lassparammap, "k0spip");
        LASS_contribution_re.insert(pair<string, double>("K0star1430plus", pass.Re()));
        LASS_contribution_im.insert(pair<string, double>("K0star1430plus", pass.Im()));

        passKM = amplitude_Kmatrix(test_xx, j) ;

        bool reject=false;
        for(Int_t k=0; k<5; k++)
           { U1j_re.push_back(passKM[2*k]);
             U1j_im.push_back(passKM[2*k+1]);
             if(U1j_re[k]==9999. || U1j_im[k]==9999. || U1j_re[k]==0. || U1j_im[k]==0.) reject=true;
           }

        for(int k=10; k<=15; k++)
           { BB.push_back(passKM[k]);
           }

        for(int pole_index=0;pole_index<5;pole_index++) 
         {
        		Asum_re[pole_index] = 0. ;
        		Asum_im[pole_index] = 0. ;
        		for(Int_t l=0;l<5;l++) 
         	{
        			Asum_re[pole_index] += (g[pole_index][l]*U1j_re[l])/BB[pole_index] ;
        			Asum_im[pole_index] += (g[pole_index][l]*U1j_im[l])/BB[pole_index] ;
        		}
        }
        if(reject==true) continue;
        tr->Fill();

     }


 fileout->Write();
 //fileout->Close();

 return 0;


}


