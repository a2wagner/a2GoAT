#ifndef __KinFitPhysics_h__
#define __KinFitPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <exception>

#include "GTreeManager.h"
#include "GTreePluto.h"
#include "PPhysics.h"
#include "KinFit.h"
#include "TRandom3.h"
#include "particle_t.h"

class KinFitPhysics : public PPhysics
{
private:
	typedef std::vector<TLorentzVector> p4vector;
	typedef std::vector<p4vector> vvP4;
	typedef vvP4::iterator outer_it;
	typedef p4vector::iterator inner_it;

	typedef std::vector<particle_t> particle_vector;
	typedef particle_vector::iterator particle_it;

	class no_final_state_exception : public std::exception
	{
		virtual const char* what() const throw()
		{
			return "No final state particle! Something went wrong...";
		}
	} noFS;

	GH1* time;
	GH1* time_cut;
	GH1* time_2g;
	GH1* time_2g_cut;

	GH1* IM;
	GH1* IM_2g;

	GH1* MM;
	GH1* MM_2g;

	TH1* TaggerAccScal;

	// used to store if the current file contains MC data
	bool MC;
	// manage MC trees
	bool pluto_tree, geant_tree;
	GTreeA2Geant* geant;
	GTreePluto* pluto;

	// store infos of collected particles
	unsigned int nParticles, nParticlesCB, nParticlesTAPS;

	// fit related stuff
	TH1* n_particles;
	TH1* invM_2g;
	TH1* missM;
	TH1* expectProt;

	TH2* pi0_fit_steps;
	TH2* coplanarity;

	TH2* missM_pi0_vs_pTheta;
	TH2* invM_gg_vs_beamE;
	TH2* p_kinE_vs_pTheta_true;
	TH2* invM_gg_smeared_vs_invM_gg_fit;

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
	TH1* pi0_smeared;
	TH1* pi0_fitted;
	TH1* copl_smeared;
	TH1* copl_fitted;
	TH1* nIter_converged;
	TH1* nIter_all;
	TH1* fitStatus;

	TH1* photon1P_Dsmeared;
	TH1* photon2P_Dsmeared;
	TH1* protonP_Dsmeared;
	TH1* photon1P_Dfitted;
	TH1* photon2P_Dfitted;
	TH1* protonP_Dfitted;

protected:
	// some general constants
	const unsigned int N_FINAL_STATE = 3;
	const double MASS_ETAPRIME = 957.78;
	//const double MASS_PROTON = 938.272;
	//const double MASS_ETA = 547.862;
	const double MASS_PIZERO = 134.9766;
	const double R2D = 180./3.14159265359;
	const double MAX_DOUBLE = 0x1.FFFFFFFFFFFFFp1023;

	// class used for kinematic fitting
	KinFit kinFit;

	// methods to get relative errors
	double sigma_E(const TLorentzVector* const p);
	double sigma_theta(const TLorentzVector* const p);
	double sigma_phi(const TLorentzVector* const p);

	// collect particles
	void GetParticles();
	void GetTrueParticles(particle_vector* particles, particle_t* beam);

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
