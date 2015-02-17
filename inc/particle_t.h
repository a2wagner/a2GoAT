#ifndef __PARTICLE_T__
#define __PARTICLE_T__

#include <algorithm>
#include <math.h>

#include "TLorentzVector.h"

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
	return T(0);
}

template<typename T, typename ValType>
T charged(ValType v) {
	typedef enum_traits<T> traits;
	const T *first = traits::charged;
	const T *last = endof(traits::charged);
	if (traits::sorted) {  // probably premature optimization
		if (std::binary_search(first, last, v))
			return T(v);
	} else if (std::find(first, last, v) != last) {
		return T(v);
	}
	return T(0);
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
	static constexpr particle_id charged[] = {
		positron,
		electron,
		antimuon,
		muon,
		piplus,
		piminus,
		proton
	};
};

//typedef std::pair<TLorentzVector, particle_id> particle_t;
struct particle_t {
	TLorentzVector p4;
	particle_id id;
	particle_t() : p4(), id(rootino) {}
	particle_t(TLorentzVector p4) : p4(p4), id(rootino) {}
	particle_t(particle_id id) : p4(), id(id) {}
	particle_t(TLorentzVector p4, particle_id id) : p4(p4), id(id) {}
	bool is_final_state()
	{
		return check<particle_id>(id);
	}
	bool is_charged()
	{
		return charged<particle_id>(id);
	}
	void transform(double factor)
	{
		double pt = p4.Pt()*factor, phi = p4.Phi(), e = p4.E()*factor;
		p4.SetXYZT(pt*cos(phi), pt*sin(phi), pt*sinh(p4.Eta()), e);
	}
	void GeV2MeV(){ transform(1000.); }
	void MeV2GeV(){ transform(.0001); }

	/* Some frequently used methods from TLorentzVector */
	Double_t E()     { return p4.E(); }
	Double_t P()     { return p4.P(); }
	Double_t Px()    { return p4.Px(); }
	Double_t Py()    { return p4.Py(); }
	Double_t Pz()    { return p4.Pz(); }
	Double_t Pt()    { return p4.Pt(); }
	Double_t M()     { return p4.M(); }
	Double_t M2()    { return p4.M2(); }
	Double_t Theta() { return p4.Theta(); }
	Double_t Phi()   { return p4.Phi(); }
	TVector3 Vect()  { return p4.Vect(); }
	void Print() { p4.Print(); }

	/* Add a few basic operators */
	particle_t operator + (const particle_t& p) const {
		return particle_t(p4+p.p4);
	}
	particle_t operator + (const TLorentzVector& p) const {
		return particle_t(p4+p);
	}
	particle_t operator += (const particle_t& p) const {
		return particle_t(p4+p.p4);
	}
	particle_t operator += (const TLorentzVector& p) const {
		return particle_t(p4+p);
	}
	particle_t operator - (const particle_t& p) const {
		return particle_t(p4-p.p4);
	}
	particle_t operator - (const TLorentzVector& p) const {
		return particle_t(p4-p);
	}
	particle_t operator -= (const particle_t& p) const {
		return particle_t(p4-p.p4);
	}
	particle_t operator -= (const TLorentzVector& p) const {
		return particle_t(p4-p);
	}
	/* Implicit cast to TLorentzVector for easier usage (like TLorentzVector = particle_t + particle_t */
	operator TLorentzVector() { return p4; }
};

#endif  // __PARTICLE_T__
