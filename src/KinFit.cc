/**
 * KinFit
 *
 * Wrapper class serving several methods to easily use the KinFitter
 */

#include "KinFit.h"

KinFit::KinFit()
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

double KinFit::getProbability()
{
	return TMath::Prob(getS(), getNDF());
}

