#include "GTreePluto.h"
#include "PParticle.h"

GTreePluto::GTreePluto(GTreeManager* Manager):
    GTree(Manager, TString("data")),
    PlutoMCTrue(NULL),
    plutoID(-1),
    plutoRandomID(-1)
{}


GTreePluto::~GTreePluto()
{}

GTreePluto::ParticleList GTreePluto::GetFinalState() const
{
    ParticleList list;

    for( Long64_t i=0;i<PlutoMCTrue->GetEntries(); ++i) {

        const PParticle* const particle = dynamic_cast<const PParticle*>((*PlutoMCTrue)[i]);

        if(particle && particle->GetDaughterIndex() == -1)
            list.push_back(particle);
    }

    return list;
}

GTreePluto::ParticleList GTreePluto::GetAllParticles() const
{
    ParticleList list;

    for( Long64_t i=0;i<PlutoMCTrue->GetEntries(); ++i) {

        const PParticle* const particle = dynamic_cast<const PParticle*>((*PlutoMCTrue)[i]);

        if(particle)
            list.push_back(particle);
    }

    return list;
}

PParticle* GTreePluto::GetMCTrue(const int idx) const
{
    if (!PlutoMCTrue) {
        //return NULL;
        fprintf(stderr, "ERROR: No Pluto Tree in current file!\n");
        exit(1);
    }

    return dynamic_cast<PParticle*>(PlutoMCTrue->At(idx));
}

TLorentzVector GTreePluto::GetTrueP4(const int idx) const
{
    if (!PlutoMCTrue) {
        fprintf(stderr, "ERROR: No Pluto Tree in current file!\n");
        exit(1);
    }

    return dynamic_cast<PParticle*>(PlutoMCTrue->At(idx))->Vect4();
}

void    GTreePluto::SetBranchAdresses()
{
    inputTree->SetBranchAddress("Particles", 	&PlutoMCTrue);
    inputTree->SetBranchAddress("plutoID", 	&plutoID);
    inputTree->SetBranchAddress("plutoRandomID", 	&plutoRandomID);
}

void    GTreePluto::SetBranches()
{
    outputTree->Branch("Particles", PlutoMCTrue);
    outputTree->Branch("plutoID", &plutoID);
    outputTree->Branch("plutoRandomID", &plutoRandomID);
}
