#ifndef __SASCHAPHYSICS_H__
#define __SASCHAPHYSICS_H__

#include "AntPhysics.h"

#include <APLCON.hpp>
#include "plot/Histogram.h"


#include <vector>
#include <map>
#include <random>

#include "TH1D.h"
#include "TH2D.h"

template<typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> v = v1;
    v.insert(v.end(),v2.begin(),v2.end());
    return v;
}

namespace ant {
namespace analysis {

class SaschaPhysics: public Physics {

protected:

    // lightweight structure for linking to fitter
    struct FitParticle {
        void SetFromVector(const TLorentzVector& p_)
        {
            Ek = p_.E()-p_.M();
            Theta = p_.Theta();
            Phi = p_.Phi();
            Ek_Sigma = 0.02*Ek*pow(Ek,-0.36);
            Theta_Sigma = 2.5*TMath::DegToRad();
            if (Theta > 20*TMath::DegToRad() && Theta < 160*TMath::DegToRad())
                Phi_Sigma = Theta_Sigma/sin(Theta);
            else
                Phi_Sigma = 1*TMath::DegToRad();
        }

        static TLorentzVector Make(const std::vector<double>& EkThetaPhi, const Double_t m);
        static TLorentzVector Make(const FitParticle& p, const Double_t m)
        {
            return Make(std::vector<double>{p.Ek, p.Theta, p.Phi}, m);
        }

        std::vector<double*> Link()
        {
            return {std::addressof(Ek),
                    std::addressof(Theta),
                    std::addressof(Phi)};
        }
        std::vector<double*> LinkSigma()
        {
            return {std::addressof(Ek_Sigma),
                    std::addressof(Theta_Sigma),
                    std::addressof(Phi_Sigma)};
        }

        void Smear();

        double Ek;
        double Ek_Sigma;
        double Theta;
        double Theta_Sigma;
        double Phi;
        double Phi_Sigma;
    };

    template<typename T>
    bool contains(const std::vector<T>& vec, const T& val)
    {
        return std::find(vec.begin(), vec.end(), val) != vec.end();
    }

    template<typename T>
    bool contains(const std::vector<const T*>& vec, const T& val)
    {
        return std::find(vec.begin(), vec.end(), &val) != vec.end();
    }

    typedef std::vector<ant::Particle> particle_vector;

    const map<short, string> component = {{0, "Energy"}, {1, "Theta"}, {2, "Phi"}};

    // choose here what you want to do
    // please also provide GoAT trees with matching MC true information...
    static constexpr bool includeCoplanarityConstraint = true;
    static constexpr bool includeIMconstraint = false;
    static constexpr bool includeVertexFit = true;
    static constexpr size_t nFinalState = 4;
    const double IM = ParticleTypeDatabase::EtaPrime.Mass();


    size_t nParticles, nParticlesCB, nParticlesTAPS;

    TH2D* banana;
    TH1D* particles;
    TH1D* tagger;
    TH1D* ntagged;
    TH1D* cbesum;
    // different checks
    TH2D* lepton_energies;
    TH2D* lepton_energies_true;
    TH2D* photon_energy_vs_opening_angle;
    TH2D* photon_energy_vs_opening_angle_true;
    TH2D* theta_vs_clusters;
    TH2D* opening_angle_vs_q2;
    TH2D* opening_angle_vs_E_high;
    TH2D* opening_angle_vs_E_low;
    TH2D* dEvE;
    TH2D* crystals_vs_ecl_charged;
    TH2D* crystals_vs_ecl_uncharged;
    TH2D* crystals_vs_ecl_charged_candidates;
    TH2D* energy_vs_momentum_z_balance;

    TH1D* opening_angle_leptons;
    TH1D* opening_angle_leptons_true;
    TH1D* energy_lepton1;
    TH1D* energy_lepton1_true;
    TH1D* energy_lepton2;
    TH1D* energy_lepton2_true;
    TH1D* energy_photon;
    TH1D* energy_photon_true;
    // proton checks
    TH1D* proton_energy;
    TH1D* proton_energy_true;
    TH1D* proton_energy_fit;
    TH1D* proton_energy_delta;
    TH1D* proton_angle_TAPS_expected;

    TH1D* coplanarity;
    TH1D* missing_mass;

    std::map<const ParticleTypeDatabase::Type*, TH1D*> numParticleType;

    TH1D* chisquare;
    TH1D* probability;
    TH1D* iterations;
    std::map<std::string, TH1D*> pulls;


    TH1D* im_true;
    TH1D* im_smeared;
    TH1D* im_fit;

    TH1D* vertex_z_after;
    TH1D* vertex_z_before;

    TH1D* q2_dist_before;
    TH1D* q2_dist_after;

    TH1D* coplanarity_fit;
    TH1D* missing_mass_fit;
    TH2D* energy_vs_momentum_z_balance_fit;


    void FillIM(TH1D* h, const std::vector<FitParticle>& final_state);

    // collect particles
    void GetParticles(const ant::Event& event, particle_vector& particles);
    void GetTrueParticles(const ant::Event& event, particle_vector& particles);

    APLCON fitter;
    FitParticle beam;
    std::vector<FitParticle> final_state;
    FitParticle proton;


public:
    SaschaPhysics(const mev_t energy_scale=1000.0);
    virtual ~SaschaPhysics() {}
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}
#endif  // __SASCHAPHYSICS_H__
