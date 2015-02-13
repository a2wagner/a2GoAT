#ifndef __PARTICLE_T__
#define __PARTICLE_T__

#include <algorithm>

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

// better to avoid linker errors if this is defined in a header and included in multiple files
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

#endif  // __PARTICLE_T__
