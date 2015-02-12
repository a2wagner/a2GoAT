#ifndef __KinFitPhysics_h__
#define __KinFitPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <algorithm>

#include "GTreeManager.h"
#include "GTreePluto.h"
#include "PPhysics.h"
#include "KinFit.h"
#include "TRandom3.h"

/* prepare a special enum class where one can check if a value was defined in the enum */
template <typename T>
struct enum_traits {
	static constexpr bool is_enum = std::is_enum<T>::value;
	static constexpr bool sorted = false;
};

template<typename T, size_t N>
T *endof(T (&ra)[N]) {
	return ra + N;
}

template<typename T, typename ValType>
T check(ValType v) {
	typedef enum_traits<T> traits;
	const T *first = traits::enumerators;
	const T *last = endof(traits::enumerators);
	if (traits::sorted) {  // probably premature optimization
		if (std::binary_search(first, last, v))
			return T(v);
	} else if (std::find(first, last, v) != last) {
		return T(v);
	}
	throw false;
}

// "enhanced" definition of enum
enum particle_id {
	rootino = 0,
	photon = 1,
	positron,
	electron,
	antimuon = 5,
	muon,
	piplus = 8,
	piminus,
	neutron = 13,
	proton = 14
};

/*template<>
struct enum_traits<particle_id> {
	static constexpr bool is_enum = std::is_enum<particle_id>::value;
	static constexpr bool sorted = true;
	static const particle_id enumerators[];
};

constexpr particle_id enum_traits<particle_id>::enumerators[] = {
	rootino,
	photon,
	positron,
	electron,
	antimuon,
	muon,
	piplus,
	piminus,
	proton
};*/


template<>
struct enum_traits<particle_id> {
	static constexpr bool is_enum = std::is_enum<particle_id>::value;
	static constexpr bool sorted = true;
	static constexpr particle_id enumerators[] = {
		rootino,
		photon,
		positron,
		electron,
		antimuon,
		muon,
		piplus,
		piminus,
		neutron,
		proton
	};
};

//typedef std::pair<TLorentzVector, particle_id> particle_t;
struct particle_t {
	TLorentzVector p4;
	particle_id id;
	particle_t() : p4(), id(rootino) {}
	particle_t(TLorentzVector p4, particle_id id) : p4(p4), id(id) {}
	bool is_final_state()
	{
		return check<particle_id>(id);
	}
};

class KinFitPhysics : public PPhysics
{
private:

	typedef std::vector<TLorentzVector> p4vector;
	typedef std::vector<p4vector> vvP4;
	typedef vvP4::iterator outer_it;
	typedef p4vector::iterator inner_it;

	typedef std::vector<particle_t> particle_vector;
	typedef particle_vector::iterator particle_it;

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
