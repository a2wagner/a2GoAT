/**
 * KinFit
 *
 * Wrapper class serving several methods to easily use the KinFitter
 */

#include "KinFit.h"

KinFit::KinFit() : chi2(0.), ndf(0.), probability(0.)
{
}

KinFit::~KinFit()
{
}

int KinFit::fillMatrixDiagonal(TMatrixD* m, Double_t* e, const int rows, const int cols)
{
	if (!m->IsValid()) {
		fprintf(stderr, "Error: Matrix is not valid!\n");
		return 1;
	}
	// set matrix to zero and resize it
	m->Zero();
	m->ResizeTo(rows, cols);
	// get pointer to matrix elements, iterate over them and add array entry if diagonal element
	Double_t *ep = m->GetMatrixArray();
	int idx = 0;
	//memset(ep, 0, m->GetNoElements*sizeof(Double_t));
	memset(ep, 0, rows*cols*sizeof(Double_t));
	for (Int_t i = 0; i < rows; i++)
		for (Int_t j = 0; j < cols; j++)
			*ep++ = (i == j ? e[idx++] : 0.);

	return 0;
}

int KinFit::fillSquareMatrixDiagonal(TMatrixD* m, Double_t* e, const int rows)
{
	return fillMatrixDiagonal(m, e, rows, rows);
}

int KinFit::fit()
{
	int ret = TKinFitter::fit();

	chi2 = getS();
	ndf = getNDF();
	probability = TMath::Prob(chi2, ndf);

	return ret;
}

void KinFit::reset()
{
	TKinFitter::reset();
	chi2 = 0.;
	ndf = 0.;
	probability = 0.;
}

