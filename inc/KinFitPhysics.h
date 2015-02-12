#ifndef __KinFitPhysics_h__
#define __KinFitPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "GTreePluto.h"
#include "PPhysics.h"
#include "KinFit.h"
#include "TRandom3.h"

class KinFitPhysics : public PPhysics
{
private:
    typedef std::vector<std::vector<TLorentzVector>> vvP4;
    typedef vvP4::iterator outer_it;
    typedef std::vector<TLorentzVector>::iterator inner_it;

	GH1* time;
	GH1* time_cut;
	GH1* time_2g;
	GH1* time_2g_cut;

	GH1* IM;
	GH1* IM_2g;

	GH1* MM;
	GH1* MM_2g;

	TH1* TaggerAccScal;

	// fit related stuff
	TH1* n_particles;
	TH1* invM_2g;
	TH1* missM;
	TH1* expectProt;

	TH1* photon1PullE;
	TH1* photon1PullTheta;
	TH1* photon1PullPhi;
	TH1* photon2PullE;
	TH1* photon2PullTheta;
	TH1* photon2PullPhi;
	TH1* protonPullE;
	TH1* protonPullTheta;
	TH1* protonPullPhi;

	TH1* chisq_hist;
	TH1* prob_hist;
	TH1* pi0_fitted;
	TH1* nIter;
	TH1* fitStatus;

protected:
	// some general constants
	const unsigned int N_FINAL_STATE = 3;
	const double MASS_ETAPRIME = 957.78;
	//const double MASS_PROTON = 938.272;
	//const double MASS_ETA = 547.862;
	const double MASS_PIZERO = 134.9766;
	const double R2D = 180./3.14159265359;
	const double MAX_DOUBLE = 0x1.FFFFFFFFFFFFFp1023;

	// Class used for kinematic fitting
	KinFit kinFit;

	virtual Bool_t Start();
	virtual void ProcessEvent();
	virtual void ProcessScalerRead();
	virtual Bool_t Write();

public:
	KinFitPhysics();
	virtual ~KinFitPhysics();
	virtual Bool_t Init(const char* configfile);

};

#endif  // __KinFitPhysics_h__
