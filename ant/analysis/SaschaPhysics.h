#ifndef __SASCHAPHYSICS_H__
#define __SASCHAPHYSICS_H__

#include "AntPhysics.h"

#include <APLCON.hpp>
#include "plot/Histogram.h"
#include "base/interval.h"


#include <vector>
#include <map>
#include <random>
#include <fstream>
#include <type_traits>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

template<typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> v = v1;
    v.insert(v.end(),v2.begin(),v2.end());
    return v;
}

template<typename T>
const size_t sum_vector(const std::vector<T>& v, typename std::enable_if<std::is_arithmetic<T>::value>::type* = 0)
{
    size_t sum = 0;
    for_each(v.cbegin(), v.cend(), [&sum](T n){ sum += n; });
    return sum;
}

template<typename T>
const size_t non_zero_entries(const std::vector<T>& v, typename std::enable_if<std::is_arithmetic<T>::value>::type* = 0)
{
    return count_if(v.cbegin(), v.cend(), [](T i){ return i; });
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
        void SetFromParticle(const Particle& p, const vector<double>& sigmas)
        {
            Ek = p.Ek();
            Theta = p.Theta();
            Phi = p.Phi();

            Ek_Sigma = sigmas.at(0);
            Theta_Sigma = sigmas.at(1)*TMath::DegToRad();
            Phi_Sigma = sigmas.at(2)*TMath::DegToRad();
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
        void AddHistogram(const std::string &name,
                          const std::string &title,
                          const std::string &x_label,
                          const std::string &y_label,
                          const BinSettings &x_bins,
                          const BinSettings &y_bins);

        // Add 3D histogram
        void AddHistogram(const std::string& name,      // short specifier for histogram
                          const std::string& title,     // descriptive title for histogram
                          const std::string& x_label,   // x axis label
                          const std::string& y_label,   // y axis label
                          const std::string& z_label,   // z axis label
                          const int x_bins_n,           // number of bins in x
                          const double x_bins_low,      // lower bound of x axis
                          const double x_bins_up,       // upper bound of y axis
                          const int y_bins_n,           // number of bins in y
                          const double y_bins_low,      // lower bound of y axis
                          const double y_bins_up,       // upper bound of y axis
                          const int z_bins_n,           // number of bins in z
                          const double z_bins_low,      // lower bound of z axis
                          const double z_bins_up);      // upper bound of z axis
        void AddHistogram(const std::string &name,
                          const std::string &title,
                          const std::string &x_label,
                          const std::string &y_label,
                          const std::string &z_label,
                          const BinSettings &x_bins,
                          const BinSettings &y_bins,
                          const BinSettings &z_bins);

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

    HistList prompt;
    HistList random;
    HistList diff;

    interval<double> prompt_window;
    interval<double> random_window1;
    interval<double> random_window2;

    typedef std::vector<ant::Particle> particle_vector;

    const map<short, string> component = {{0, "Energy"}, {1, "Theta"}, {2, "Phi"}};

    // choose here what you want to do
    // please also provide GoAT trees with matching MC true information...
    static constexpr bool includeCoplanarityConstraint = true;
    static constexpr bool includeIMconstraint = false;
    static constexpr bool includeVertexFit = true;
    static constexpr size_t nFinalState = 4;
    //const double IM = ParticleTypeDatabase::EtaPrime.Mass();
    static constexpr double IM = 957.78;
    // threshold for cluster energies
    static constexpr double CLUSTER_TRESH = 25.;
    // threshold for CB energy sum
    static constexpr double CB_ESUM = 500.;

    size_t nParticles, nParticlesCB, nParticlesTAPS;

    std::map<const ParticleTypeDatabase::Type*, TH1D*> numParticleType;

    // store true particle information in case of MC (if used)
    particle_vector true_particles;

    // histogram to keep track of efficencies
    TH1D* accepted_events;
    // histograms before prompt / random handling
    TH1D* cb_esum;
    TH2D* pid;
    TH1D* tagger_spectrum;
    TH1D* tagger_time;
    TH1D* particle_types;
    TH1D* n_part;
    TH1D* n_cluster_cb;
    TH1D* n_cluster_taps;
    TH1D* etap_chi2;
    TH2D* etap_chi2_vs_q2;
    // different MC true information
    TH2D* eplus_true_theta_vs_energy;
    TH2D* eminus_true_theta_vs_energy;
    TH2D* photon_true_theta_vs_energy;
    TH2D* proton_true_theta_vs_energy;
    TH2D* true_open_angle_vs_q2;
    TH1D* opening_angle_leptons_true;
    TH2D* opening_angle_vs_photon_energy_true;
    TH2D* lepton_energies_true;
    TH1D* energy_lepton1_true;
    TH1D* energy_lepton2_true;
    TH2D* n_cluster_cb_vs_q2;
    TH2D* n_cluster_taps_vs_q2;
    TH2D* n_cluster_cb_vs_open_angle;
    TH2D* n_cluster_taps_vs_open_angle;
    TH2D* expected_proton_diff_vs_q2;
    TH2D* expected_proton_diff_vs_q2_rebin;
    TH1D* protons_found;

    std::map<std::string, TH1*> pulls_prompt;
    std::map<std::string, TH1*> pulls_random;
    std::map<std::string, TH1*> pulls_diff;

    // invM spectra for different q^2 ranges
    static constexpr double im_q2_mev_steps = 50.;
    static constexpr double im_q2_upper_bound = 900.;
    std::vector<TH1*> im_q2_prompt;
    std::vector<TH1*> im_q2_random;
    std::vector<TH1*> im_q2_diff;

    // corrections, uncertainties
    std::vector<double> cb_gain;
    std::vector<double> taps_gain;
    std::vector<double> cb_time_correction;
    TH2D* cb_energy_correction;
    TH2D* taps_energy_correction;
    TH1D* cb_theta_correction;
    TH1D* taps_theta_correction;
    TH2D* photon_energy_uncertainties;
    TH2D* photon_theta_uncertainties;
    TH2D* photon_phi_uncertainties;
    TH2D* proton_energy_uncertainties;
    TH2D* proton_theta_uncertainties;
    TH2D* proton_phi_uncertainties;

    // read histogram from file and return pointer to the object
    TH1* read_hist(const char*, const char*) const;
    void read_file(const char*, std::vector<double>&, const int = 0);

    void set_fit_particle(const Particle&, FitParticle&);
    // uncertainties
    void get_uncertainties(const Particle&, vector<double>&);
    double get_sigma_energy(const Particle&) const;
    double get_sigma_theta(const Particle&) const;
    double get_sigma_phi(const Particle&) const;
    double get_sigma(const Particle&, TH2* const) const;
    // corrections
    void apply_time_correction(const TrackList&, const TrackList&);
    void apply_energy_correction(const TrackList&, const TrackList&);
    void apply_theta_correction(const TrackList&, const TrackList&);

    void FillIM(TH1* h, const std::vector<FitParticle>& final_state);

    // collect particles
    void GetTracks(const Event& event, TrackList& tracksCB, TrackList& tracksTAPS);
    bool IdentifyTracks(const TrackList& tracksCB, const TrackList& tracksTAPS, particle_vector& particles);
    void GetParticles(const ant::Event& event, particle_vector& particles);
    void GetTrueParticles(const ant::Event& event, particle_vector& particles);

    void sort_particles(particle_vector&);
    void fill_MC_true(const ParticleList&, const size_t, const size_t);
    const particle_vector get_MC_true_particles(const ParticleList&);
    const Particle get_true_particle(const ParticleList&, const size_t);
    const Particle get_true_positron(const ParticleList&);
    const Particle get_true_electron(const ParticleList&);
    const Particle get_true_photon(const ParticleList&);
    const Particle get_true_proton(const ParticleList&);
    void proton_tests(const TrackList&, const TaggerHitPtr);

    APLCON fitter;
    FitParticle beam;
    std::vector<FitParticle> final_state;
    FitParticle proton;
    // eta' fitter
    APLCON etap_fit;
    std::vector<FitParticle> etap_fs;


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
