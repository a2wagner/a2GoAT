#pragma once

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

template<typename T>
void shift_right(std::vector<T>& v)
{
    std::rotate(v.begin(), v.end() -1, v.end());
}

namespace ant {
namespace analysis {

class etaDalitz_CalibPID: public Physics {

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
    static constexpr bool includeCoplanarityConstraint = false;
    static constexpr bool includeVertexFit = false;
    static constexpr size_t nFinalState = 4;
    //const double IM = ParticleTypeDatabase::EtaPrime.Mass();
    static constexpr double IM = 547.85;
    // threshold for cluster energies
    static constexpr double CLUSTER_TRESH = 25.;
    // threshold for CB energy sum
    static constexpr double CB_ESUM = 500.;
    // threshold for PID energy to count particle as charged
    static constexpr double PID_THRESH = .5;
    // distance to TAPS [cm]
    static constexpr double TAPS_DISTANCE = 145;
    // radius CB [cm]
    static constexpr double RADIUS_CB = 25.4;
    // limits for target if z vertex is used in KinFit
    static constexpr double TARGET_MIN = -10.;
    static constexpr double TARGET_MAX = 10.;
    // threshold to check if double value should be treated as zero
    static constexpr double EPSILON = 2*std::numeric_limits<double>::epsilon();

    size_t nParticles, nParticlesCB, nParticlesTAPS;

    // store true particle information in case of MC (if used)
    particle_vector true_particles;

    // histogram to keep track of efficencies
    TH1D* steps;
    // histograms before prompt / random handling
    TH1D* cb_esum;
    TH2D* pid;
    TH1D* tagger_spectrum;
    TH1D* tagger_time;
    TH1D* n_part;
    TH1D* cluster_time;
    TH1D* cb_time;
    TH1D* taps_time;
    TH1D* cb_time_avg;
    TH1D* cb_avg_tagger_time_diff;
    TH1D* cb_avg_taps_time_diff;
    TH2D* energy_vs_cb_avg_taps_time_diff;
    // eta Dalitz PID calibration related
    TH2* eegPID = nullptr;
    TH1D* etaIM = nullptr;
    TH1D* etaIM_fit = nullptr;
    TH1D* etaIM_cand = nullptr;
    TH1D* etaIM_final = nullptr;
    TH1D* MM = nullptr;
    TH1D* hCopl = nullptr;
    TH1D* hCopl_final = nullptr;
    TH1D* hChi2 = nullptr;
    TH1D* hProb = nullptr;
    TH1D* hIter = nullptr;
    TH1D* pTheta = nullptr;
    TH1D* protonVeto = nullptr;

    std::map<std::string, TH1*> pulls_prompt;
    std::map<std::string, TH1*> pulls_random;
    std::map<std::string, TH1*> pulls_diff;

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
    void GetParticles(const ant::Event& event, particle_vector& particles);
    void GetTrueParticles(const ant::Event& event, particle_vector& particles);

    double calculate_energy_weighted_cb_time_average(const TrackList&) const;
    void sort_particles(particle_vector&);
    const particle_vector get_MC_true_particles(const ParticleList&);
    const Particle get_true_particle(const ParticleList&, const size_t);

    size_t get_histogram_neigbours(TH1* const, const int, vector<double>&, const int depth = 1) const;
    double local_bin_average(TH1* const, const int) const;

    APLCON fitter;
    FitParticle beam;
    std::vector<FitParticle> final_state;
    FitParticle proton;


public:
    etaDalitz_CalibPID(const mev_t energy_scale = 1000.0);
    virtual ~etaDalitz_CalibPID() {}
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}
