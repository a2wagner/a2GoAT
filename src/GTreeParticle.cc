#include "GTreeParticle.h"
#include "GTreeManager.h"


using namespace std;


GTreeParticle::GTreeParticle(GTreeManager *Manager, const TString& _Name)    :
    GTree(Manager,_Name),
    nParticles(0),
    particles(new TClonesArray("TLorentzVector", 32000))
{
}

GTreeParticle::~GTreeParticle()
{
}

void    GTreeParticle::SetBranchAdresses()
{
    tree_in->SetBranchAddress("nParticles",&nParticles);
    tree_in->SetBranchAddress("particles.", &particles);
    tree_in->SetBranchAddress("rawIndex",daughters);
}

void    GTreeParticle::SetBranches()
{
    tree_out->Branch("nParticles",&nParticles,"nParticles/i");
    tree_out->Branch("particles.", &particles, 32000, 0);
    tree_out->Branch("rawIndex",daughters,"rawIndex[nParticles]/I");
}
