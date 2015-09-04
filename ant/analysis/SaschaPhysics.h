#ifndef __SASCHAPHYSICS_H__
#define __SASCHAPHYSICS_H__

#include "AntPhysics.h"

#include <APLCON.hpp>
#include "plot/Histogram.h"
#include "base/interval.h"


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

    // Class to group histograms
    class HistList {
    protected:
        static constexpr unsigned int max_hist_per_canvas = 20;

    public:
        std::string pref;  // prefix to label whole group of histograms
        mutable std::map<std::string, TH1*> h;  // container for histograms by name (without prefix)
        std::map<std::string, std::string> h_title;  // container for histogram titles by name (without prefix)

        // Add 1D histogram
        void AddHistogram(const std::string& name,      // short specifier for histogram
                          const std::string& title,     // descriptive title for histogram
                          const std::string& x_label,   // x axis label
                          const std::string& y_label,   // y axis label
                          const int x_bins_n,           // number of bins in x
                          const double x_bins_low,      // lower bound of x axis
                          const double x_bins_up        // upper bound of x axis
                          );
        void AddHistogram(const std::string &name,
                          const std::string &title,
                          const std::string &x_label,
                          const std::string &y_label,
                          const BinSettings &bins);

        // Add 2D histogram
        void AddHistogram(const std::string& name,      // short specifier for histogram
                          const std::string& title,     // descriptive title for histogram
                          const std::string& x_label,   // x axis label
                          const std::string& y_label,   // y axis label
                          const int x_bins_n,           // number of bins in x
                          const double x_bins_low,      // lower bound of x axis
                          const double x_bins_up,       // upper bound of y axis
                          const int y_bins_n,           // number of bins in y
                          const double y_bins_low,      // lower bound of y axis
                          const double y_bins_up        // upper bound of y axis
                          );
        void AddHistogram(const string &name,
                          const string &title,
                          const string &x_label,
                          const string &y_label,
                          const BinSettings &x_bins,
                          const BinSettings &y_bins);

        HistList(const std::string& prefix, const mev_t energy_scale = 1000.0);

        void Draw();
        void AddScaled(const HistList& h2, const Double_t f = 1.);

        HistList& operator*= (const Double_t factor);
        HistList operator= (const HistList& other);

        TH1* operator[] (const std::string& key)
        {
            return h[key];
        }

        const TH1* operator[] (const std::string& key) const
        {
            return h[key];
        }
    };

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

    // histograms after applying cuts
    TH1D* im_cut;
    TH1D* im_fit_cut;
    TH1D* q2_dist_cut;
    TH1D* q2_dist_fit_cut;
    TH1D* coplanarity_cut;
    TH1D* coplanarity_fit_cut;
    TH1D* proton_angle_TAPS_expected_cut;
    TH1D* missing_mass_cut;
    TH1D* missing_mass_fit_cut;
    TH2D* energy_vs_momentum_z_balance_cut;
    TH2D* energy_vs_momentum_z_balance_fit_cut;

    TH2D* dEvE_cut;
    TH2D* crystals_vs_ecl_cut;
    TH2D* crystals_vs_ecl_charged_cut;
    TH2D* crystals_vs_ecl_uncharged_cut;

    // invM spectra for different q^2 ranges
    static constexpr double im_q2_mev_steps = 50.;
    std::vector<TH1D*> im_q2;
    TH1D* im_q2_0_50;
    TH1D* im_q2_50_100;
    TH1D* im_q2_100_150;
    TH1D* im_q2_150_200;
    TH1D* im_q2_200_250;
    TH1D* im_q2_250_300;
    TH1D* im_q2_300_350;
    TH1D* im_q2_350_400;
    TH1D* im_q2_400_450;
    TH1D* im_q2_450_500;
    TH1D* im_q2_500_550;
    TH1D* im_q2_550_600;
    TH1D* im_q2_600_650;
    TH1D* im_q2_650_700;
    TH1D* im_q2_700_750;
    TH1D* im_q2_750_800;
    TH1D* im_q2_800_850;
    TH1D* im_q2_850_900;


    void FillIM(TH1D* h, const std::vector<FitParticle>& final_state);

    // collect particles
    void GetParticles(const ant::Event& event, particle_vector& particles);
    void GetTrueParticles(const ant::Event& event, particle_vector& particles);

    APLCON fitter;
    FitParticle beam;
    std::vector<FitParticle> final_state;
    FitParticle proton;


public:
    SaschaPhysics(const mev_t energy_scale = 1000.0);
    virtual ~SaschaPhysics() {}
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}
#endif  // __SASCHAPHYSICS_H__
