#ifndef __KinFit_h__
#define __KinFit_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "KinFitter/TKinFitter.h"
#include "KinFitter/TFitParticlePThetaPhi.h"
#include "KinFitter/TFitConstraintEp.h"
#include "KinFitter/TFitConstraintM.h"

class KinFit : public TKinFitter
{
private:
	double chi2, ndf, probability;

protected:

public:
	int fillMatrixDiagonal(TMatrixD* m, Double_t* e, const int rows, const int cols);
	int fillSquareMatrixDiagonal(TMatrixD* m, Double_t* e, const int rows);
	double getChi2(){ return chi2; };
	double getNdf(){ return ndf; };
	double getProbability(){ return probability; };

	KinFit();
	virtual ~KinFit();

	int fit();
	void reset();

	template<typename T>
	inline void fill_array(T* t, std::initializer_list<T> il)
	{
		std::copy(il.begin(), il.end(), t);
	}

};

#endif  // __KinFit_h__
