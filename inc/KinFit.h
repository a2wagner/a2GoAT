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

protected:

public:
	int fillMatrixDiagonal(TMatrixD* m, Double_t* e, int rows, int cols);
	int fillSquareMatrixDiagonal(TMatrixD* m, Double_t* e, int rows);
	KinFit();
	virtual ~KinFit();

};

#endif  // __KinFit_h__
