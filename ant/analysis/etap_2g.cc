#include "etap_2g.h"
#include "Particle.h"
#include "Track.h"
#include "plot/root_draw.h"
#include <string>
#include "utils/combinatorics.h"
#include <vector>
#include <numeric>
#include <functional>
#include <APLCON.hpp>
#include <iomanip>

using namespace std;
using namespace ant;


TLorentzVector analysis::etap_2g::FitParticle::Make(const std::vector<double>& EkThetaPhi, const Double_t m)
{
    const double E = EkThetaPhi[0] + m;
    const Double_t p = sqrt(E*E - m*m);
    TVector3 pv(1,0,0);
    pv.SetMagThetaPhi(p, EkThetaPhi[1], EkThetaPhi[2]);
    TLorentzVector l(pv, E);
    return l;
}

void analysis::etap_2g::FitParticle::Smear()
{
}

void analysis::etap_2g::FillIM(TH1 *h, const std::vector<analysis::etap_2g::FitParticle>& final_state)
{
    TLorentzVector sum(0,0,0,0);
    for (auto p = final_state.cbegin(); p != final_state.cend()-1; ++p)
        sum += FitParticle::Make(*p, ParticleTypeDatabase::Photon.Mass());
    h->Fill(sum.M());
}

void analysis::etap_2g::HistList::AddHistogram(const string &name,
                                               const string &title,
                                               const string &x_label,
                                               const string &y_label,
                                               const int x_bins_n,
                                               const double x_bins_low,
                                               const double x_bins_up)
{
    // setup one dimensional histogram TH1D
    h_title[name] = title;
    h[name] = HistogramFactory::Default().Make1D(
                title,
                x_label,
                y_label,
                BinSettings(x_bins_n, x_bins_low, x_bins_up),
                pref+name);
}

void analysis::etap_2g::HistList::AddHistogram(const string &name,
                                               const string &title,
                                               const string &x_label,
                                               const string &y_label,
                                               const BinSettings &bins)
{
    // setup one dimensional histogram TH1D
    h_title[name] = title;
    h[name] = HistogramFactory::Default().Make1D(
                title,
                x_label,
                y_label,
                bins,
                pref+name);
}

void analysis::etap_2g::HistList::AddHistogram(const string &name,
                                               const string &title,
                                               const string &x_label,
                                               const string &y_label,
                                               const int x_bins_n,
                                               const double x_bins_low,
                                               const double x_bins_up,
                                               const int y_bins_n,
                                               const double y_bins_low,
                                               const double y_bins_up)
{
    // setup two dimensional histogram TH2D
    h_title[name] = title;
    h[name] = HistogramFactory::Default().Make2D(
                title,
                x_label,
                y_label,
                BinSettings(x_bins_n, x_bins_low, x_bins_up),
                BinSettings(y_bins_n, y_bins_low, y_bins_up),
                pref+name);
}

void analysis::etap_2g::HistList::AddHistogram(const string &name,
                                               const string &title,
                                               const string &x_label,
                                               const string &y_label,
                                               const BinSettings &x_bins,
                                               const BinSettings &y_bins)
{
    // setup two dimensional histogram TH2D
    h_title[name] = title;
    h[name] = HistogramFactory::Default().Make2D(
                title,
                x_label,
                y_label,
                x_bins,
                y_bins,
                pref+name);
}

void analysis::etap_2g::HistList::AddHistogram(const string &name,
                                               const string &title,
                                               const string &x_label,
                                               const string &y_label,
                                               const string &z_label,
                                               const int x_bins_n,
                                               const double x_bins_low,
                                               const double x_bins_up,
                                               const int y_bins_n,
                                               const double y_bins_low,
                                               const double y_bins_up,
                                               const int z_bins_n,
                                               const double z_bins_low,
                                               const double z_bins_up)
{
    // setup three dimensional histogram TH3D
    h_title[name] = title;
    h[name] = HistogramFactory::Default().Make3D(
                title,
                x_label,
                y_label,
                z_label,
                BinSettings(x_bins_n, x_bins_low, x_bins_up),
                BinSettings(y_bins_n, y_bins_low, y_bins_up),
                BinSettings(z_bins_n, z_bins_low, z_bins_up),
                pref+name);
}

void analysis::etap_2g::HistList::AddHistogram(const string &name,
                                               const string &title,
                                               const string &x_label,
                                               const string &y_label,
                                               const string &z_label,
                                               const BinSettings &x_bins,
                                               const BinSettings &y_bins,
                                               const BinSettings &z_bins)
{
    // setup three dimensional histogram TH3D
    h_title[name] = title;
    h[name] = HistogramFactory::Default().Make3D(
                title,
                x_label,
                y_label,
                z_label,
                x_bins,
                y_bins,
                z_bins,
                pref+name);
}

void analysis::etap_2g::HistList::Draw()
{
    if (!h.size())
        return;

    std::vector<canvas*> cv;
    size_t n = h.size() / max_hist_per_canvas;
    if (h.size() % max_hist_per_canvas)
        n++;
    cv.resize(n);

    size_t i = 0, j;
    auto it = h.begin();
    for (auto& c : cv) {
        c = new canvas("etap_2g: Overview " + pref + std::to_string(++i));
        j = 0;
        while (j++ < max_hist_per_canvas && it != h.end()) {
            TH2D* h2 = dynamic_cast<TH2D*>(it->second);
            if (h2 != nullptr)
                *c << drawoption("colz");
            *c << it->second;
            ++it;
        }
        *c << endc;
    }
}

analysis::etap_2g::HistList &analysis::etap_2g::HistList::operator*=(const Double_t factor)
{
    for (auto i : h)
        i.second->Scale(factor);

    return *this;
}

analysis::etap_2g::HistList analysis::etap_2g::HistList::operator=(const analysis::etap_2g::HistList &other)
{
    for (auto i : h) {
        TH1* h = i.second;
        h->Reset();
        h->Add(other[i.first]);
    }

    return *this;
}

void analysis::etap_2g::HistList::AddScaled(const analysis::etap_2g::HistList &h2, const Double_t f)
{
    for (auto i : h)
        i.second->Add(h2.h.at(i.first), f);
}

analysis::etap_2g::HistList::HistList(const string &prefix, const mev_t energy_scale)
{
    pref = prefix + '_';

    const BinSettings energy_bins(1600, 0, energy_scale*1.6);
    const BinSettings tagger_bins(200, 1400, 1600);
    const BinSettings veto_bins(1000, 0, 10);
    const BinSettings chisquare_bins(100, 0, 30);
    const BinSettings probability_bins(100, 0, 1);
    const BinSettings iterations_bins(15, 0, 15);
    const BinSettings im_bins(200, IM-100, IM+100);
    const BinSettings vertex_bins(200, -10, 10);
    const BinSettings phi_bins(720, 0, 360);
    const BinSettings energy_range_bins(1600, -1600, 1600);
    const BinSettings q2_bins(2000, 0, 1000);
    const BinSettings count_bins(50, 0, 50);

    AddHistogram("tagger_energy", "Tagger Energy", "Photon Beam Energy [MeV]", "#", tagger_bins);
    AddHistogram("n_tagged", "Tagger Hits", "Tagger Hits / event", "#", count_bins);

    AddHistogram("q2_dist_before", "q^{2} Distribution before KinFit", "q^{2} [MeV]", "#", q2_bins);
    AddHistogram("q2_dist_after", "q^{2} Distribution after KinFit", "q^{2} [MeV]", "#", q2_bins);

    AddHistogram("coplanarity", "Coplanarity #eta' proton", "Coplanarity [#circ]", "#", phi_bins);
    AddHistogram("missing_mass", "Missing Mass Proton", "m_{miss} [MeV]", "#", energy_bins);
    AddHistogram("energy_vs_momentum_z_balance", "Energy vs. p_{z} balance", "p_{z} [GeV]", "Energy [GeV]",
                 energy_range_bins, energy_range_bins);

    // make fitter histograms
    AddHistogram("chisquare", "ChiSquare", "#chi^{2}", "#", chisquare_bins);
    AddHistogram("probability", "Probability", "Probability", "#", probability_bins);
    AddHistogram("iterations", "Number of iterations", "Iterations", "#", iterations_bins);

    stringstream fs;
    fs << "gg";
    AddHistogram("im_true", "IM "+fs.str()+" true", "IM", "#", im_bins);
    AddHistogram("im_smeared", "IM "+fs.str()+" smeared", "IM", "#", im_bins);
    AddHistogram("im_fit", "IM "+fs.str()+" fit", "IM", "#", im_bins);

    AddHistogram("vertex_z_before", "Vertex Z Before", "v_z [cm]", "#", vertex_bins);
    AddHistogram("vertex_z_after", "Vertex Z After", "v_z [cm]", "#", vertex_bins);

    AddHistogram("coplanarity_fit", "Coplanarity #eta' proton fitted", "Coplanarity [#circ]", "#", phi_bins);
    AddHistogram("missing_mass_fit", "Missing Mass Proton fitted", "m_{miss} [MeV]", "#", energy_bins);
    AddHistogram("energy_vs_momentum_z_balance_fit", "Energy vs. p_{z} balance fitted", "p_{z} [GeV]", "Energy [GeV]",
                 energy_range_bins, energy_range_bins);

    // histograms after cuts
    AddHistogram("im_cut", "IM "+fs.str()+" after cuts", "IM", "#", im_bins);
    AddHistogram("im_fit_cut", "IM "+fs.str()+" fit after cuts", "IM", "#", im_bins);
    AddHistogram("q2_dist_cut", "q^{2} Distribution after Cuts", "q^{2} [GeV]", "#", q2_bins);
    AddHistogram("q2_dist_fit_cut", "q^{2} Distribution fit after Cuts", "q^{2} [GeV]", "#", q2_bins);
    AddHistogram("coplanarity_cut", "Coplanarity #eta' proton after Cuts", "Coplanarity [#circ]", "#", phi_bins);
    AddHistogram("coplanarity_fit_cut", "Coplanarity #eta' proton fit after Cuts", "Coplanarity [#circ]", "#", phi_bins);
    AddHistogram("missing_mass_cut", "Missing Mass Proton after Cuts", "m_{miss} [MeV]", "#", energy_bins);
    AddHistogram("missing_mass_fit_cut", "Missing Mass Proton fit after Cuts", "m_{miss} [MeV]", "#", energy_bins);
    AddHistogram("energy_vs_momentum_z_balance_cut", "Energy vs. p_{z} balance after Cuts", "p_{z} [GeV]", "Energy [GeV]",
                 energy_range_bins, energy_range_bins);
    AddHistogram("energy_vs_momentum_z_balance_fit_cut", "Energy vs. p_{z} balance fit after Ctus", "p_{z} [GeV]",
                 "Energy [GeV]", energy_range_bins, energy_range_bins);

    AddHistogram("dEvE_cut", "dE vs. E Cut", "E_{Crystals} [MeV]", "dE_{Veto} [MeV]", energy_bins, veto_bins);

    // invM spectra for different q^2 ranges
    int start_range = 0;
    while (start_range < int(im_q2_upper_bound)) {
        char title[40];
        char name[20];
        sprintf(name, "im_q2_%d_%d", int(start_range), int(start_range + im_q2_mev_steps));
        sprintf(title, "IM %d MeV < q^{2} < %d MeV", int(start_range), int(start_range + im_q2_mev_steps));
        AddHistogram(name, title, "IM [MeV]", "#", im_bins);
        start_range += im_q2_mev_steps;
    }
}

ant::analysis::etap_2g::etap_2g(const mev_t energy_scale) :
    Physics("etap_2g"),
    prompt("prompt", energy_scale),
    random("random", energy_scale),
    diff("diff", energy_scale),
    prompt_window(-6, 6),
    random_window1(-40, -15),
    random_window2(15, 40),
    fitter("etap_2g"),
    final_state(nFinalState),
    detectors(nFinalState + 1)  // beam photon needs to be considered
{
    cout << "Eta' 2gamma Physics:\n";
    cout << "Prompt window: " << prompt_window << " ns\n";
    cout << "Random window 1: " << random_window1 << " ns\n";
    cout << "Random window 2: " << random_window2 << " ns\n";

    // save current directory for later usage (ROOT changes it while reading stuff)
    const char* dir = gDirectory->GetPath();

    // read in histograms with uncertainties and corrections as well as correction factors from files
    read_file("files/data.2014.07/CB_lingain.txt", cb_gain);
    read_file("files/data.2014.07/TAPS_lingain.txt", taps_gain);
    read_file("files/data.2014.07/time_corr_CB.txt", cb_time_correction);
    cb_energy_correction = static_cast<TH2D*>(read_hist("files/corr/CB_e_corr.root", "g_peak_E_CB"));
    taps_energy_correction = static_cast<TH2D*>(read_hist("files/corr/TAPS_e_corr.root", "g_peak_E_TAPS"));
    cb_theta_correction = static_cast<TH1D*>(read_hist("files/corr/CB_th_corr.root", "photon_dtheta_v_theta_CB_pfx"));
    taps_theta_correction = static_cast<TH1D*>(read_hist("files/corr/TAPS_th_corr.root", "photon_dtheta_v_theta_TAPS_pfx"));
    photon_energy_uncertainties = static_cast<TH2D*>(read_hist("files/APLCONunc/photon_uncertainties_vz.root", "E"));
    photon_theta_uncertainties = static_cast<TH2D*>(read_hist("files/APLCONunc/photon_uncertainties_vz.root", "theta"));
    photon_phi_uncertainties = static_cast<TH2D*>(read_hist("files/APLCONunc/photon_uncertainties_vz.root", "phi"));
    //proton_energy_uncertainties = static_cast<TH2D*>(read_hist("files/APLCONunc/proton_uncertainties_vz.root", "E"));
    proton_theta_uncertainties = static_cast<TH2D*>(read_hist("files/APLCONunc/proton_uncertainties_vz.root", "p_theta"));
    proton_phi_uncertainties = static_cast<TH2D*>(read_hist("files/APLCONunc/proton_uncertainties_vz.root", "p_phi"));

    // histogram to count number of different particle types
    for (auto& t : ParticleTypeDatabase::DetectableTypes())
        numParticleType[t]= HistFac.makeTH1D("Number of " + t->PrintName(),
                                      "number of " + t->PrintName() + "/ event",
                                      "", BinSettings(16,0,16), "n_" + t->PrintName());

    // histogram to count how many events passed the different selection criteria
    accepted_events = HistFac.makeTH1D("Accepted events", "condition", "#events", BinSettings(10), "accepted_events");

    // some histograms filled before tagger time window information is used
    const BinSettings energy_bins(1600, 0, energy_scale*1.6);
    const BinSettings veto_bins(1000, 0, 10);
    const BinSettings energy_tagger(200, 1400, 1600);
    const BinSettings time_tagger(1300, -650, 650);
    const BinSettings particle_bins(10);
    const BinSettings time_bins(1000, -50, 50);
    cb_esum = HistFac.makeTH1D("CB Energy Sum", "E [MeV]", "#", energy_bins, "cb_esum");
    pid = HistFac.makeTH2D("PID Bananas", "Cluster Energy [MeV]", "Veto Energy [MeV]", energy_bins, veto_bins, "pid");
    tagger_spectrum = HistFac.makeTH1D("Tagger Spectrum", "Photon Beam Energy [MeV]", "#", energy_tagger, "tagger_spectrum");
    tagger_time = HistFac.makeTH1D("Tagger Time", "Tagger Time [ns]", "#", time_tagger, "tagger_time");
    particle_types = HistFac.makeTH1D("Identified particles", "Particle Type", "#", particle_bins, "particle_types");
    n_part = HistFac.makeTH1D("Number of particles", "#particles", "#", particle_bins, "n_part");
    n_cluster_cb = HistFac.makeTH1D("Number of Clusters in CB", "#particles", "#", particle_bins, "n_cluster_cb");
    n_cluster_taps = HistFac.makeTH1D("Number of Clusters in TAPS", "#particles", "#", particle_bins, "n_cluster_taps");
    cluster_time = HistFac.makeTH1D("Cluster Time", "time [ns]", "#", time_tagger, "cluster_time");
    cb_time = HistFac.makeTH1D("CB Time", "time [ns]", "#", time_tagger, "cb_time");
    taps_time = HistFac.makeTH1D("TAPS Time", "time [ns]", "#", time_tagger, "taps_time");
    cb_time_avg = HistFac.makeTH1D("Energy weighted CB time average", "avg time [ns]", "#", time_bins, "cb_time_avg");
    cb_avg_tagger_time_diff = HistFac.makeTH1D("Time difference CB average and Tagger", "time [ns]", "#", time_tagger,
                                               "cb_avg_tagger_time_diff;");
    cb_avg_taps_time_diff = HistFac.makeTH1D("Time difference CB average and TAPS elemets", "time [ns]", "#", time_bins,
                                             "cb_avg_taps_time_diff");
    energy_vs_cb_avg_taps_time_diff = HistFac.makeTH2D("Energy vs. Time difference CB average and TAPS elemets", "time [ns]",
                                                       "E [MeV]", time_bins, energy_bins, "energy_vs_cb_avg_taps_time_diff");

    // prepare invM histograms for different q2 ranges
    int start_range = 0;
    while (start_range < int(im_q2_upper_bound)) {
        char name[20];
        sprintf(name, "im_q2_%d_%d", int(start_range), int(start_range + im_q2_mev_steps));
        im_q2_prompt.push_back(prompt[name]);
        im_q2_random.push_back(random[name]);
        im_q2_diff.push_back(diff[name]);
        start_range += im_q2_mev_steps;
    }

    // setup fitter for eta' 2g decay
    fitter.LinkVariable("Beam", beam.Link(), beam.LinkSigma());
    fitter.LinkVariable("Photon1", final_state[0].Link(), final_state[0].LinkSigma());
    fitter.LinkVariable("Photon2", final_state[1].Link(), final_state[1].LinkSigma());
    fitter.LinkVariable("Proton", final_state[2].Link(), final_state[2].LinkSigma());

    vector<string> all_names = {"Beam", "Photon1", "Photon2", "Proton"};
    if (includeVertexFit) {
        all_names.push_back("v_z");
        //fitter.AddUnmeasuredVariable("v_z"); // default value 0
        auto v_z_settings = APLCON::Variable_Settings_t::Default;
        v_z_settings.Limit.High = TARGET_MAX;
        v_z_settings.Limit.Low = TARGET_MIN;
        fitter.AddMeasuredVariable("v_z", 0., 2.3, v_z_settings);  // default value 0
    }

    auto recalculate_cluster = [this] (const vector<vector<double>>& particles, size_t id) -> vector<double>
    {
        constexpr double ln2 = 0.693;
        // nuclear properties from http://pdg.lbl.gov/2015/AtomicNuclearProperties/
        // general calorimetry infos: http://www.kip.uni-heidelberg.de/~coulon/Lectures/Detectors/Free_PDFs/Lecture9.pdf
        // radiation lengths
        constexpr double X0_NaI = 2.588;
        constexpr double X0_BaF2 = 2.026;
        constexpr double X0_PbWO4 = 0.89;
        // critical energies
        constexpr double Ec_NaI = 13.3;
        constexpr double Ec_BaF2 = 13.7;
        constexpr double Ec_PbWO4 = 9.6;
        // last element in particles is v_z (scalar has dimension 1)
        const double v_z = particles.back()[0];

        // get correct constants for the material
        double X0, Ec;
        if (detectors.at(id) & ant::detector_t::NaI) {
            X0 = X0_NaI;
            Ec = Ec_NaI;
        } else if (detectors.at(id) & ant::detector_t::BaF2) {
            X0 = X0_BaF2;
            Ec = Ec_BaF2;
        } else if (detectors.at(id) & ant::detector_t::PbWO4) {
            X0 = X0_PbWO4;
            Ec = Ec_PbWO4;
        } else  // no known material
            return particles[id];

        // shower depth ~ ln(E/E_c)
        const double L = log(particles[id][0]/Ec);
        // E = E_0 exp(-x/X_0)
        const double r = X0*L/ln2;  // shower depth maximum
        const double theta = particles[id][1];
        const double R = detectors.at(id) & ant::detector_t::anyCB ? RADIUS_CB + r : TAPS_DISTANCE/cos(theta) + r;
        const double theta_new = atan2(R*sin(theta), R*cos(theta) - v_z);
        // return new particle values
        vector<double> part = particles[id];
        part[1] = theta_new;
        return part;
    };

    // Constraint: Incoming 4-vector = Outgoing 4-vector
    auto EnergyMomentumBalance = [recalculate_cluster] (vector<vector<double>>& particles) -> vector<double>
    {
        const TLorentzVector target(0,0,0, ParticleTypeDatabase::Proton.Mass());
        TLorentzVector diff(0,0,0,0);
        if (includeVertexFit) {
            // assume first particle is beam photon
            diff = target + FitParticle::Make(particles[0], ParticleTypeDatabase::Photon.Mass());
            /* from here all clusters are recalculated by taking a shifted z vertex position into account */
            // assume second and third particle outgoing photons
            diff -= FitParticle::Make(recalculate_cluster(particles, 1), ParticleTypeDatabase::Photon.Mass());
            diff -= FitParticle::Make(recalculate_cluster(particles, 2), ParticleTypeDatabase::Photon.Mass());
            // assume last particle outgoing proton
            diff -= FitParticle::Make(recalculate_cluster(particles, 3), ParticleTypeDatabase::Proton.Mass());
        } else {
            // assume first particle is beam photon
            TLorentzVector diff = target + FitParticle::Make(particles[0], ParticleTypeDatabase::Photon.Mass());
            // assume second and third particle outgoing photons
            diff -= FitParticle::Make(particles[1], ParticleTypeDatabase::Photon.Mass());
            diff -= FitParticle::Make(particles[2], ParticleTypeDatabase::Photon.Mass());
            // assume last particle outgoing proton
            diff -= FitParticle::Make(particles[3], ParticleTypeDatabase::Proton.Mass());
        }

        return {diff.X(), diff.Y(), diff.Z(), diff.T()};
    };
    fitter.AddConstraint("EnergyMomentumBalance", all_names, EnergyMomentumBalance);

    // Constraint: Coplanarity between eta' and recoil proton
    auto CoplanarityConstraint = [recalculate_cluster] (const vector<vector<double>>& particles) -> double
    {
        TLorentzVector etap(0,0,0,0);
        TLorentzVector proton;
        if (includeVertexFit) {
            etap = FitParticle::Make(recalculate_cluster(particles, 1), ParticleTypeDatabase::Photon.Mass());
            etap += FitParticle::Make(recalculate_cluster(particles, 2), ParticleTypeDatabase::Photon.Mass());
            proton = FitParticle::Make(recalculate_cluster(particles, 3), ParticleTypeDatabase::Proton.Mass());
        } else {
            etap = FitParticle::Make(particles[1], ParticleTypeDatabase::Photon.Mass());
            etap += FitParticle::Make(particles[2], ParticleTypeDatabase::Photon.Mass());
            proton = FitParticle::Make(particles[3], ParticleTypeDatabase::Proton.Mass());
        }
        return abs(etap.Phi() - proton.Phi())*TMath::RadToDeg() - 180.;
    };
    if (includeCoplanarityConstraint)
        fitter.AddConstraint("CoplanarityConstraint", all_names, CoplanarityConstraint);

    // Constraint: Invariant mass of nPhotons equals constant IM,
    // make lambda catch also this with [&] specification
    auto RequireIM = [recalculate_cluster] (const vector<vector<double>>& particles) -> double
    {
        TLorentzVector sum(0,0,0,0);
        if (includeVertexFit) {
            sum += FitParticle::Make(recalculate_cluster(particles, 1), ParticleTypeDatabase::Photon.Mass());
            sum += FitParticle::Make(recalculate_cluster(particles, 2), ParticleTypeDatabase::Photon.Mass());
        } else {
            sum += FitParticle::Make(particles[1], ParticleTypeDatabase::Photon.Mass());
            sum += FitParticle::Make(particles[2], ParticleTypeDatabase::Photon.Mass());
        }

        return sum.M() - IM;
    };
    if (includeIMconstraint)
        fitter.AddConstraint("RequireIM", all_names, RequireIM);

    // Constraint: Vertex position in z direction: v_z (positive if upstream)
    // if the particle originated from (0,0,v_z) instead of origin,
    // the corrected angle theta' is given by
    // tan(theta') = (R sin(theta))/(R cos(theta) - v_z)
    // R is the CB radius, 10in aka 25.4cm

/*    auto VertexConstraint = [&] (vector<vector<double>>& particles) -> double
    {
        constexpr double R = 25.4;
        // last element in particles is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = particles.back()[0];
        particles.resize(particles.size()-1); // get rid of last element
        // correct each photon's theta angle,
        // then calculate invariant mass of all photons
        TLorentzVector sum(0,0,0,0);
        for (auto i = particles.begin()+1; i != particles.end()-1; ++i) {
            const double theta = (*i)[1]; // second element is theta
            const double theta_p = std::atan2( R*sin(theta), R*cos(theta) - v_z);
            (*i)[1] = theta_p;
            sum += FitParticle::Make(*i, ParticleTypeDatabase::Photon.Mass());
        }
        return sum.M() - IM;
    };

    if (includeVertexFit) {
        fitter.AddUnmeasuredVariable("v_z"); // default value 0
        fitter.AddConstraint("VertexConstraint", all_names + std::vector<string>{"v_z"}, VertexConstraint);
    }

    static_assert(!(includeIMconstraint && includeVertexFit), "Do not enable Vertex and IM Fit at the same time");
*/
    // create pull histograms
    const BinSettings pull_bins(50,-3,3);
    // change back to former directory to store pulls
    gDirectory->cd(dir);
    gDirectory->mkdir("Pulls");
    gDirectory->cd("Pulls");
    for (const auto& varname : fitter.VariableNames()) {
        string title(varname);
        size_t pos = title.find("[");
        if (pos != string::npos) {
            const string prop = " " + component.at(atoi(&varname.at(pos+1)));
            title.replace(pos, 3, prop);
        }
        prompt.AddHistogram("pull_" + varname, "Pull " + title, "Pull", "#", pull_bins);
        pulls_prompt[varname] = prompt["pull_" + varname];
        random.AddHistogram("pull_" + varname, "Pull " + title, "Pull", "#", pull_bins);
        pulls_random[varname] = random["pull_" + varname];
        diff.AddHistogram("pull_" + varname, "Pull " + title, "Pull", "#", pull_bins);
        pulls_diff[varname] = diff["pull_" + varname];
    }

    APLCON::Fit_Settings_t settings = fitter.GetSettings();
    settings.MaxIterations = 10;
    fitter.SetSettings(settings);

    cout.precision(3);
    APLCON::PrintFormatting::Width = 11;
}

void ant::analysis::etap_2g::ProcessEvent(const ant::Event &event)
{
    const TLorentzVector target(0., 0., 0., ParticleTypeDatabase::Proton.Mass());

    const bool MC = !(event.MCTrue().Particles().GetAll().empty());
    if (MC)
        true_particles.clear();

    accepted_events->Fill("all events", 1);

    cb_esum->Fill(event.Reconstructed().TriggerInfos().CBEenergySum());

    if (MC && event.Reconstructed().TriggerInfos().CBEenergySum() < CB_ESUM)
        return;

    if (MC)
        accepted_events->Fill("CB Esum", 1);

    for (const auto& track : event.Reconstructed().Tracks()) {
        pid->Fill(track->ClusterEnergy(), track->VetoEnergy());
        // fill time information
        cluster_time->Fill(track->Time());
        if (track->Detector() & detector_t::anyCB)
            cb_time->Fill(track->Time());
        else
            taps_time->Fill(track->Time());
    }

    // determine #clusters per event per detector
    /*size_t nCB = 0, nTAPS = 0;
    std::for_each(event.Reconstructed().Tracks().begin(), event.Reconstructed().Tracks().end(),
                  [&nCB, &nTAPS](const auto& track){ track->Detector() & detector_t::anyCB ? nCB++ : nTAPS++; });
    prompt["n_cluster_cb"]->Fill(nCB);
    prompt["n_cluster_taps"]->Fill(nTAPS);*/

    // get the tracks in CB and TAPS
    TrackList tracksCB;
    TrackList tracksTAPS;
    GetTracks(event, tracksCB, tracksTAPS);
    //printf("Particles found: %zu CB, %zu TAPS\n", tracksCB.size(), tracksTAPS.size());
    n_cluster_cb->Fill(tracksCB.size());
    n_cluster_taps->Fill(tracksTAPS.size());

    const double cb_time_avg_energy_weighted = calculate_energy_weighted_cb_time_average(tracksCB);
    cb_time_avg->Fill(cb_time_avg_energy_weighted);
    for (const auto& track: tracksTAPS) {
        cb_avg_taps_time_diff->Fill(cb_time_avg_energy_weighted - track->Time());
        energy_vs_cb_avg_taps_time_diff->Fill(cb_time_avg_energy_weighted - track->Time(), track->ClusterEnergy());
    }

    n_part->Fill(event.Reconstructed().Tracks().size());
    for (auto& particle : event.Reconstructed().Particles().GetAll()) {
        particle_types->Fill(particle->Type().PrintName().c_str(), 1);
    }

    for (auto& t : ParticleTypeDatabase::DetectableTypes()) {
        try {
            numParticleType.at(t)->Fill(event.Reconstructed().Particles().Get(*t).size());
        } catch (...) {}
    }

    particle_vector particles;
    if (!IdentifyTracks(tracksCB, tracksTAPS, particles))
        return;
    accepted_events->Fill("identification", 1);

    TaggerHistList tagger_hits;
    size_t prompt_n_tagger = 0, random_n_tagger = 0;
    //static size_t count = 0;
    constexpr bool trueMC = true;  // use MC true
    if (trueMC && MC)
        tagger_hits = event.MCTrue().TaggerHits();
    else
        tagger_hits = event.Reconstructed().TaggerHits();

    for (const auto& taggerhit : tagger_hits) {
        tagger_spectrum->Fill(taggerhit->PhotonEnergy());
        tagger_time->Fill(taggerhit->Time());

        const double cb_avg_tagger_time = taggerhit->Time() - cb_time_avg_energy_weighted;
        cb_avg_tagger_time_diff->Fill(cb_avg_tagger_time);

        // determine if the event is in the prompt or random window
        bool is_prompt = false;
        if (prompt_window.Contains(taggerhit->Time()))
            is_prompt = true;
        else if (random_window1.Contains(taggerhit->Time()) || random_window1.Contains(taggerhit->Time()))
            is_prompt = false;
        else
            continue;

        accepted_events->Fill("Tagger time", 1);

        is_prompt ? prompt_n_tagger++ : random_n_tagger++;

        // make sure the correct histogram will be filled
        HistList& h = is_prompt ? prompt : random;

        h["tagger_energy"]->Fill(taggerhit->PhotonEnergy());

//        // find the photons and one proton
//        size_t foundPhotons = 0;
//        for(const auto& p : event.MCTrue().Particles().GetAll()) {
//            if(p->Type() == ParticleTypeDatabase::Proton) {
//                proton.SetFromVector(*p);
//            }
//            else if(foundPhotons<nPhotons && p->Type() == ParticleTypeDatabase::Photon) {
//                photons[foundPhotons].SetFromVector(*p);
//                foundPhotons++;
//           }
//        }
//        if(foundPhotons != nPhotons)
//            continue;

        nParticles = nParticlesCB = nParticlesTAPS = 0;
        //particle_vector particles;
        /*GetParticles(event, particles);
        if (nParticles != nFinalState)
            continue;

        accepted_events->Fill("n_FS", 1);*/

        // check if we have the right amount of leptons, photons and protons
        unsigned short ng = 0, np = 0;
        for (auto it = particles.begin(); it != particles.end(); ++it) {
            if (it->Type() == ParticleTypeDatabase::Photon)
                ng++;
            else if (it->Type() == ParticleTypeDatabase::Proton)
                np++;
        }
        if (ng != 2 || np != 1)
            continue;

        accepted_events->Fill("#part", 1);

        // sort order of particles
        sort_particles(particles);

        TLorentzVector proton = particles.back();
        TLorentzVector etap(0., 0., 0., 0.);
        for (auto it = particles.cbegin(); it != particles.cend()-1; ++it)
            etap += *it;

/*        // perform cut on PID elements to suppress conversion lepton pairs
        element_index_t pid1 = particles.at(0).Tracks().front()->CentralVeto(),
                        pid2 = particles.at(1).Tracks().front()->CentralVeto();
        if (pid1 == pid2)
            continue;
        accepted_events->Fill("PID Cut", 1);*/

        TLorentzVector balanceP4 = taggerhit->PhotonBeam() + target;
        //for (const auto& p : particles)
        //    balanceP4 -= p;
        balanceP4 -= etap + proton;
        h["energy_vs_momentum_z_balance"]->Fill(balanceP4.Pz(), balanceP4.E());

        const double copl = abs(etap.Phi() - particles.back().Phi())*TMath::RadToDeg();
        h["coplanarity"]->Fill(copl);
        TLorentzVector missingProton = taggerhit->PhotonBeam() + target - etap;
        double missM = missingProton.M();
        h["missing_mass"]->Fill(missM);

        //std::copy(particles.begin(), particles.end(), final_state);
//        std::transform(particles.begin(), particles.end(), final_state.begin(),
//                       [](Particle& p) -> double { return FitParticle().SetFromVector(p); });
        auto fs_it = final_state.begin();
        auto det_it = detectors.begin();
        det_it++;  // skip first element which is the beam photon --> not considered when recalculating clusters due to v_z
        for (const auto& p : particles) {
            //(it++)->SetFromVector(p);
            set_fit_particle(p, *fs_it++);
            *det_it++ = p.Tracks().front()->Detector();
        }

        double q2_before = (particles[0] + particles[1]).M();
        h["q2_dist_before"]->Fill(q2_before);

        // set proton energy sigma to zero to indicate it's unmeasured
        final_state[2].Ek_Sigma = 0;

//        std::transform(myv1.begin(), myv1.end(), myv1.begin(), [](double d) -> double { return d * 3; });
//        for_each(begin(myv1), end(myv1), [](double& a) { a *= 3; });

        beam.SetFromVector(taggerhit->PhotonBeam());
        // set beam sigmas, energy 2/sqrt(3)
        beam.Ek_Sigma = 1.1547;
        beam.Theta_Sigma = .0001;
        beam.Phi_Sigma = .0001;

        FillIM(h["im_true"], final_state);

//        // smear the MC true data
//        proton.Smear();
//        for(auto& photon : photons)
//            photon.Smear();
//        beam.Smear();

//        FillIM(h["im_smeared"], photons);

/*        if (missM < 880. && missM > 1130.)
            continue;
        accepted_events->Fill("miss mass", 1);*/

/*        // Cut on missing mass of the proton
        if (missM < 900. && missM > 990.)
            continue;
        // Cut on Veto energy
        if (highVeto)
            continue;
        // Cut on opening angle between TAPS cluster and expected proton
        // (mainly for suppressing bad events which will be unnecassarily fitted, affects only a few background events)
        if (openAngle_p_TAPS_expected > 5.)
            continue;
*/
        // let APLCON do the work
        const APLCON::Result_t& result = fitter.DoFit();

        //cout << result << endl;

        if (result.Status != APLCON::Result_Status_t::Success) {
            //cout << result << endl;
            continue;
        }

        accepted_events->Fill("KinFit", 1);

        for (const auto& it_map : result.Variables) {
            const string& varname = it_map.first;
            const APLCON::Result_Variable_t& var = it_map.second;
            if (is_prompt)
                pulls_prompt.at(varname)->Fill(var.Pull);
            else
                pulls_random.at(varname)->Fill(var.Pull);
        }
        h["chisquare"]->Fill(result.ChiSquare);
        h["probability"]->Fill(result.Probability);
        h["iterations"]->Fill(result.NIterations);

        if (includeVertexFit) {
            h["vertex_z_after"]->Fill(result.Variables.at("v_z").Value.After);
            h["vertex_z_before"]->Fill(result.Variables.at("v_z").Value.Before);
        }

        FillIM(h["im_fit"], final_state);

        double q2_after = (FitParticle::Make(final_state[0], ParticleTypeDatabase::eMinus.Mass())
                + FitParticle::Make(final_state[1], ParticleTypeDatabase::eMinus.Mass())).M();
        h["q2_dist_after"]->Fill(q2_after);

        TLorentzVector balanceP4_fit = target + FitParticle::Make(beam, ParticleTypeDatabase::Photon.Mass());
        // create the fitted final state particles
        TLorentzVector proton_fit = FitParticle::Make(final_state[2], ParticleTypeDatabase::Proton.Mass());
        TLorentzVector etap_fit(0., 0., 0., 0.);
        etap_fit += FitParticle::Make(final_state[0], ParticleTypeDatabase::Photon.Mass());
        etap_fit += FitParticle::Make(final_state[1], ParticleTypeDatabase::Photon.Mass());
        balanceP4_fit -= etap_fit;
        balanceP4_fit -= proton_fit;
        h["energy_vs_momentum_z_balance_fit"]->Fill(balanceP4_fit.Pz(), balanceP4_fit.E());

        const double copl_fit = abs(etap_fit.Phi() - particles.back().Phi())*TMath::RadToDeg();
        h["coplanarity_fit"]->Fill(copl_fit);
        TLorentzVector missingProtonFit = target + FitParticle::Make(beam, ParticleTypeDatabase::Photon.Mass()) - etap_fit;
        h["missing_mass_fit"]->Fill(missingProtonFit.M());

        // Cut on chi^2
        if (result.ChiSquare > 10.)
            return;

        accepted_events->Fill("#chi^{2} cut", 1);

        // fill the invM histograms for different q2 ranges
        try {
            if (q2_after < im_q2_upper_bound) {
                if (is_prompt)
                    im_q2_prompt.at(static_cast<int>(q2_after/im_q2_mev_steps))->Fill(etap_fit.M());
                else
                    im_q2_random.at(static_cast<int>(q2_after/im_q2_mev_steps))->Fill(etap_fit.M());
            }
        } catch (const std::out_of_range& oor) {
            std::cerr << "Out of Range error: " << oor.what() << '\n';
            cerr << "q2_after: " << q2_after << ", im_q2_upper_bound: " << im_q2_upper_bound << endl;
        }

        h["im_cut"]->Fill(etap.M());
        FillIM(h["im_fit_cut"], final_state);
        h["q2_dist_cut"]->Fill(q2_before);
        h["q2_dist_fit_cut"]->Fill(q2_after);
        h["coplanarity_cut"]->Fill(copl);
        h["missing_mass_cut"]->Fill(missM);
        h["missing_mass_fit_cut"]->Fill(missingProtonFit.M());
        h["coplanarity_fit_cut"]->Fill(copl_fit);
        h["energy_vs_momentum_z_balance_cut"]->Fill(balanceP4.Pz(), balanceP4.E());
        h["energy_vs_momentum_z_balance_fit_cut"]->Fill(balanceP4_fit.Pz(), balanceP4_fit.E());

        for (const auto& p : particles) {
            TrackPtr tr = p.Tracks().front();
            h["dEvE_cut"]->Fill(tr->ClusterEnergy(), tr->VetoEnergy());
        }
    }
    prompt["n_tagged"]->Fill(prompt_n_tagger);
    random["n_tagged"]->Fill(random_n_tagger);
}

void ant::analysis::etap_2g::Finish()
{
    double factor = -prompt_window.Length() / (random_window1.Length() + random_window2.Length());
    diff = prompt;
    diff.AddScaled(random, factor);
}

void ant::analysis::etap_2g::ShowResult()
{
    canvas c("etap_2g: Overview");
    c << drawoption("colz") << pid
      << padoption::set(padoption_t::Legend)
      << particle_types
      << padoption::unset(padoption_t::Legend)
      << tagger_spectrum << diff["n_tagged"] << cb_esum << endc;

    canvas c_pulls("etap_2g: Pulls");
    c_pulls << padoption::set(padoption_t::LogY);
    for (auto& p : pulls_diff)
        c_pulls << p.second;
    c_pulls << endc;

    canvas c_fitter("etap_2g: Fitter");
    c_fitter <<  diff["chisquare"] <<  diff["probability"] << diff["iterations"]
             <<  diff["im_true"] <<  diff["im_smeared"] <<  diff["im_fit"]
             <<  diff["vertex_z_before"] <<  diff["vertex_z_after"] << endc;

    // draw all random subtracted histograms
    //diff.Draw();
}

void ant::analysis::etap_2g::sort_particles(particle_vector& particles)
{
//    for (const auto& p : particles)
//        cout << p.Type() << endl;

    // swap the proton and the last vector entry if the proton is not the last particle
    for (auto it = particles.begin(); it != particles.end()-1; ++it)
        if (it->Type() == ParticleTypeDatabase::Proton) {
            std::iter_swap(it, particles.end()-1);
            break;
        }
}

const vector<Particle> ant::analysis::etap_2g::get_MC_true_particles(const ParticleList& particles)
{
    if (!true_particles.empty())
        return true_particles;

    Particle proton_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
    for (const auto& p : particles) {
        if (p->Type() == ParticleTypeDatabase::Photon)
            true_particles.push_back(Particle(*p));
        else if (p->Type() == ParticleTypeDatabase::Proton)
            proton_true = Particle(*p);
    }
    if (proton_true.Type() == ParticleTypeDatabase::Neutron || true_particles.size() != nFinalState-1)
        accepted_events->Fill("wrong MC", 1);

    true_particles.push_back(proton_true);
    return true_particles;
}

const Particle ant::analysis::etap_2g::get_true_particle(const ParticleList& particles, const size_t pos)
{
    if (true_particles.empty())
        return get_MC_true_particles(particles).at(pos);
    else
        return true_particles.at(pos);
}

double ant::analysis::etap_2g::calculate_energy_weighted_cb_time_average(const TrackList& tracks) const
{
    double time = 0, e_sum = 0;
    for (const auto& track : tracks) {
        time += track->Time();
        e_sum += track->ClusterEnergy();
    }
    return time/e_sum;
}

void ant::analysis::etap_2g::GetTracks(const Event& event, TrackList& tracksCB, TrackList& tracksTAPS)
{
    for (const auto& track : event.Reconstructed().Tracks()) {
        // ignore particles below set cluster energy threshold
        if (track->ClusterEnergy() < CLUSTER_TRESH)
            continue;
        if (track->Detector() & detector_t::anyCB)
            tracksCB.push_back(track);
        else if (track->Detector() & detector_t::anyTAPS)
            tracksTAPS.push_back(track);
    }
}

bool ant::analysis::etap_2g::IdentifyTracks(const TrackList& tracksCB, const TrackList& tracksTAPS,
                                            particle_vector& particles)
{
    bool success = false;
    size_t nCB = tracksCB.size(), nCharged = 0;
    if (nCB + tracksTAPS.size() != nFinalState)
        return success;
    if (nCB == 2) {
        for (const auto& track : tracksCB) {
            particles.push_back(Particle(ParticleTypeDatabase::Photon, track));
            if (track->Detector() & detector_t::PID)
                nCharged++;
        }
        particles.push_back(Particle(ParticleTypeDatabase::Proton, tracksTAPS.front()));
        if (nCharged)
            return success;
        success = true;
    } else if (nCB == 1) {
        particles.push_back(Particle(ParticleTypeDatabase::Photon, tracksCB.front()));
        if (tracksCB.front()->Detector() & detector_t::PID)
            nCharged++;
        if (nCharged)
            return success;
        //TODO: Correct for path lengths!
        double time_track0 = abs(tracksTAPS.at(0)->Time()), time_track1 = abs(tracksTAPS.at(1)->Time());
        unsigned short proton_id = 0;
        if (time_track1 > time_track0)
            proton_id = 1;
        particles.push_back(Particle(ParticleTypeDatabase::Photon, tracksTAPS.at(1-proton_id)));
        particles.push_back(Particle(ParticleTypeDatabase::Proton, tracksTAPS.at(proton_id)));
        success = true;
    }

    return success;
}

void ant::analysis::etap_2g::GetParticles(const ant::Event& event, particle_vector& particles)
{
    // example of how to collect particles in a user defined way
    for (const auto& track : event.Reconstructed().Tracks()) {
        // ignore particles below set cluster energy threshold
        if (track->ClusterEnergy() < CLUSTER_TRESH)
            continue;
        if (track->Detector() & ant::detector_t::NaI) {
            nParticlesCB++;
            //if (track->VetoEnergy() > 0.)  // PID entry? --> charged
            if (track->Detector() & detector_t::PID)
                particles.push_back(Particle(ParticleTypeDatabase::eMinus, track));
            else
                particles.push_back(Particle(ParticleTypeDatabase::Photon, track));
        } else if (track->Detector() & ant::detector_t::BaF2 || track->Detector() & ant::detector_t::PbWO4) {
            nParticlesTAPS++;
            //if (track->VetoEnergy() > 0.)  // Veto entry? --> charged
            if (track->Detector() & detector_t::Veto)
                particles.push_back(Particle(ParticleTypeDatabase::Proton, track));
            else
                particles.push_back(Particle(ParticleTypeDatabase::Photon, track));
        }
    //std::cout << track->Detector() << std::endl;
    }
    nParticles = nParticlesCB + nParticlesTAPS;
}

void ant::analysis::etap_2g::GetTrueParticles(const ant::Event& event, particle_vector& particles)
{
    for (const auto& p : event.MCTrue().Particles().GetAll())
        if (contains(ParticleTypeDatabase::MCFinalStateTypes(), p->Type()))
            particles.push_back(Particle(*p));
    nParticles = particles.size();
}

TH1* ant::analysis::etap_2g::read_hist(const char* file, const char* hist_name) const
{
    TFile f(file);
    if (!f.IsOpen()) {
        cerr << "Couldn't open file " << file << endl;
        exit(1);
    }
    TH2* h = static_cast<TH2*>(f.Get(hist_name));
    if (!h) {
        cerr << "Couldn't find histogram '" << hist_name << "' in file " << file << endl;
        exit(1);
    }
    h->SetDirectory(0);  // dissociate histogram from the file
    f.Close();
    return h;
}

void ant::analysis::etap_2g::read_file(const char* file, std::vector<double>& values, const int fill_size)
{
    values.clear();
    ifstream ifs(file);

    if (ifs) {
        string line;
        string buffer;
        stringstream ss;

        getline(ifs, line);
        ss << line;

        while (getline(ss, buffer, '\t'))
            values.push_back(stod(buffer));

        ifs.close();
    } else {
        cerr << "File " << file << " couldn't be found" << endl;
        if (!fill_size)
            exit(1);
        else {
            cout << "Use " << fill_size << " zeros instead" << endl;
            values.resize(fill_size);
        }
    }
}

void ant::analysis::etap_2g::set_fit_particle(const Particle& p, FitParticle& part)
{
    vector<double> sigmas;
    get_uncertainties(p, sigmas);
    part.SetFromParticle(p, sigmas);
}

void ant::analysis::etap_2g::get_uncertainties(const Particle& p, vector<double>& sigmas)
{
    sigmas.clear();
    sigmas.push_back(get_sigma_energy(p));
    sigmas.push_back(get_sigma_theta(p));
    sigmas.push_back(get_sigma_phi(p));
}

double ant::analysis::etap_2g::get_sigma_energy(const Particle& p) const
{
    if (p.Type() == ParticleTypeDatabase::Proton)
        return 0;//get_sigma(p, proton_energy_uncertainties);
    else
        return get_sigma(p, photon_energy_uncertainties);
}

double ant::analysis::etap_2g::get_sigma_theta(const Particle& p) const
{
    if (p.Type() == ParticleTypeDatabase::Proton)
        return get_sigma(p, proton_theta_uncertainties);
    else
        return get_sigma(p, photon_theta_uncertainties);
}

double ant::analysis::etap_2g::get_sigma_phi(const Particle& p) const
{
    if (p.Type() == ParticleTypeDatabase::Proton)
        return get_sigma(p, proton_phi_uncertainties);
    else
        return get_sigma(p, photon_phi_uncertainties);
}

double ant::analysis::etap_2g::get_sigma(const Particle& p, TH2* const h) const
{
    int bin;
    if (p.Type() == ParticleTypeDatabase::Proton)
        bin = h->FindBin(p.Tracks().front()->CentralCrystal());
    else
        bin = h->FindBin(p.Ek(), p.Theta()*TMath::RadToDeg());
    double sigma = static_cast<double>(h->GetBinContent(bin));
    if (sigma < EPSILON)
        sigma = local_bin_average(h, bin);

    return sigma;
}

// select ring of neighbouring bins, depth defines how many rings
// store valid bin contents in passed vector (not under-/overflow, above epsilon (threshold for empty bins))
// returns size of possible neighbours
size_t ant::analysis::etap_2g::get_histogram_neigbours(TH1* const h,
                                                             const int bin,
                                                             vector<double>& values,
                                                             const int depth) const
{
    if (depth < 1) {
        cerr << "Depth of neighbouring crystals has to be greater than 0!" << endl;
        return 0;
    }

    values.clear();
    const int dim = h->GetDimension();
    const int overflow_bin_x = h->GetNbinsX() + 1;
    size_t possible_neighbours = 0;
    if (dim == 1) {
        for (int x = bin - depth; x <= bin + depth; x++) {
            possible_neighbours++;
            if (x < 1)
                continue;
            if (x >= overflow_bin_x)
                continue;

            double content = static_cast<double>(h->GetBinContent(x));
            if (content > EPSILON)
                values.push_back(content);
        }
    } else if (dim == 2) {
        const int bins_x = h->GetNbinsX() + 2;
        const int bins = bins_x*(h->GetNbinsY()+2);
        for (int x = -depth; x <= depth; x++)
            for (int y = -depth; y <= depth; y++) {
                possible_neighbours++;
                int i = (bin + x) + y*(bins_x);
                if (i < 1 || i >= bins)
                    continue;
                if (h->IsBinUnderflow(i) || h->IsBinOverflow(i))
                    continue;

                double content = static_cast<double>(h->GetBinContent(i));
                if (content > EPSILON)
                    values.push_back(content);
            }
    } else
        cerr << "histogram dimension " << dim << " not supported" << endl;

    return possible_neighbours;
}

// get a local bin average based on neighbouring bins
// amount of neighbours taken into account is based on how many valid bins surround the bin of interest
double ant::analysis::etap_2g::local_bin_average(TH1* const h, const int bin) const
{
    if (!h) {
        cerr << "No valid histogram given!" << endl;
        return 0;
    }

    vector<double> values;
    size_t neighbours = get_histogram_neigbours(h, bin, values);
    if (values.size()/neighbours < .3)
        get_histogram_neigbours(h, bin, values, 2);
    if (!values.size()) {
        cerr << "Couldn't find any neighbouring bins to calculate average" << endl;
        return 0;
    }

    return sum_vector(values)/values.size();
}

void ant::analysis::etap_2g::apply_time_correction(const TrackList& tracksCB, const TrackList& tracksTAPS)
{
    for (const auto track : tracksCB)
        track->SetTime(track->Time() - cb_time_correction.at(track->CentralCrystal()));
}

void ant::analysis::etap_2g::apply_energy_correction(const TrackList& tracksCB, const TrackList& tracksTAPS)
{
    double e, delta;

    for (const auto track : tracksCB) {
        e = track->ClusterEnergy()*cb_gain.at(track->CentralCrystal());
        delta = static_cast<double>(cb_energy_correction->GetBinContent(cb_energy_correction->FindBin(e, track->Theta())));
        track->SetClusterEnergy(e*(1 - delta));
    }

    for (const auto track : tracksTAPS) {
        e = track->ClusterEnergy()*taps_gain.at(track->CentralCrystal());
        delta = static_cast<double>(taps_energy_correction->GetBinContent(taps_energy_correction->FindBin(e, track->Theta())));
        track->SetClusterEnergy(e*(1 - delta));
    }
}

void ant::analysis::etap_2g::apply_theta_correction(const TrackList& tracksCB, const TrackList& tracksTAPS)
{
    double delta;

    for (const auto track : tracksCB) {
        delta = static_cast<double>(cb_theta_correction->GetBinContent(cb_theta_correction->FindBin(track->Theta())));
        track->SetTheta(track->Theta() - delta);
    }

    for (const auto track : tracksTAPS) {
        delta = static_cast<double>(taps_theta_correction->GetBinContent(taps_theta_correction->FindBin(track->Theta())));
        track->SetTheta(track->Theta() - delta);
    }
}
