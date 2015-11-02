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
    static constexpr double cluster_thresh = 25.;
    // threshold for CB energy sum
    static constexpr double cb_esum = 500.;

    size_t nParticles, nParticlesCB, nParticlesTAPS;

    std::map<const ParticleTypeDatabase::Type*, TH1D*> numParticleType;

    TH1D* accepted_events;

    std::map<std::string, TH1*> pulls_prompt;
    std::map<std::string, TH1*> pulls_random;
    std::map<std::string, TH1*> pulls_diff;

    // invM spectra for different q^2 ranges
    static constexpr double im_q2_mev_steps = 50.;
    static constexpr double im_q2_upper_bound = 900.;
    std::vector<TH1*> im_q2_prompt;
    std::vector<TH1*> im_q2_random;
    std::vector<TH1*> im_q2_diff;

    void FillIM(TH1* h, const std::vector<FitParticle>& final_state);

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
