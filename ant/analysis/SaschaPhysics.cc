#include "SaschaPhysics.h"
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

TLorentzVector analysis::SaschaPhysics::FitParticle::Make(const std::vector<double>& EkThetaPhi, const Double_t m)
{
    const double E = EkThetaPhi[0] + m;
    const Double_t p = sqrt( E*E - m*m );
    TVector3 pv(1,0,0);
    pv.SetMagThetaPhi(p, EkThetaPhi[1], EkThetaPhi[2]);
    TLorentzVector l(pv, E);
    return l;
}

void analysis::SaschaPhysics::FitParticle::Smear()
{
}

void analysis::SaschaPhysics::FillIM(TH1 *h, const std::vector<analysis::SaschaPhysics::FitParticle>& final_state)
{
    TLorentzVector sum(0,0,0,0);
//    for (auto p = final_state.cbegin(); p != final_state.cend()-1; ++p) {
//        sum += FitParticle::Make(*p, ParticleTypeDatabase::Photon.Mass());
//    }
    sum += FitParticle::Make(final_state[0], ParticleTypeDatabase::eMinus.Mass());
    sum += FitParticle::Make(final_state[1], ParticleTypeDatabase::eMinus.Mass());
    sum += FitParticle::Make(final_state[2], ParticleTypeDatabase::Photon.Mass());
    h->Fill(sum.M());
}

void analysis::SaschaPhysics::HistList::AddHistogram(const string &name,
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

void analysis::SaschaPhysics::HistList::AddHistogram(const string &name,
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

void analysis::SaschaPhysics::HistList::AddHistogram(const string &name,
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

void analysis::SaschaPhysics::HistList::AddHistogram(const string &name,
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

void analysis::SaschaPhysics::HistList::AddHistogram(const string &name,
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

void analysis::SaschaPhysics::HistList::AddHistogram(const string &name,
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

void analysis::SaschaPhysics::HistList::Draw()
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
        c = new canvas("SaschaPhysics: Overview " + pref + std::to_string(++i));
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

analysis::SaschaPhysics::HistList &analysis::SaschaPhysics::HistList::operator*=(const Double_t factor)
{
    for (auto i : h)
        i.second->Scale(factor);

    return *this;
}

analysis::SaschaPhysics::HistList analysis::SaschaPhysics::HistList::operator=(const analysis::SaschaPhysics::HistList &other)
{
    for (auto i : h) {
        TH1* h = i.second;
        h->Reset();
        h->Add(other[i.first]);
    }

    return *this;
}

void analysis::SaschaPhysics::HistList::AddScaled(const analysis::SaschaPhysics::HistList &h2, const Double_t f)
{
    for (auto i : h)
        i.second->Add(h2.h.at(i.first), f);
}

analysis::SaschaPhysics::HistList::HistList(const string &prefix, const mev_t energy_scale)
{
    pref = prefix + '_';

    const BinSettings energy_bins(1600,0,energy_scale*1.6);
    const BinSettings tagger_bins(1500,300,1800);
    const BinSettings ntaggerhits_bins(15);
    const BinSettings veto_bins(1000,0,10);
    const BinSettings particle_bins(10,0,10);
    const BinSettings particlecount_bins(16,0,16);
    const BinSettings chisquare_bins(100,0,30);
    const BinSettings probability_bins(100,0,1);
    const BinSettings iterations_bins(15,0,15);
    const BinSettings im_bins(200,IM-100,IM+100);
    const BinSettings vertex_bins(200,-10,10);
    const BinSettings theta_bins(720, 0, 180);
    const BinSettings phi_bins(720, 0, 360);
    const BinSettings angle_bins(500, 0, 50);
    const BinSettings energy_range_bins(1600, -1600, 1600);
    const BinSettings q2_bins(2000, 0, 1000);
    const BinSettings count_bins(50, 0, 50);
    const BinSettings tof_bins(1000, -50, 50);

    AddHistogram("pid", "PID Bananas", "Cluster Energy [MeV]", "Veto Energy [MeV]", energy_bins, veto_bins);
    AddHistogram("particle_types", "Identified particles", "Particle Type", "#", particle_bins);
    AddHistogram("n_part", "Number of particles", "#particles", "#", particle_bins);
    AddHistogram("n_cluster_cb", "Number of Clusters in CB", "#particles", "#", particle_bins);
    AddHistogram("n_cluster_taps", "Number of Clusters in TAPS", "#particles", "#", particle_bins);
    AddHistogram("tagger_spectrum", "Tagger Spectrum", "Photon Beam Energy [MeV]", "#", tagger_bins);
    AddHistogram("tagger_time", "Tagger Time", "Tagger Time [ns]", "#", energy_range_bins);
    AddHistogram("tagger_energy", "Tagger Energy", "Photon Beam Energy [MeV]", "#", tagger_bins);
    AddHistogram("n_tagged", "Tagger Hits", "Tagger Hits / event", "#", count_bins);
    AddHistogram("cb_esum", "CB Energy Sum", "E [MeV]", "#", energy_bins);

    AddHistogram("proton_id", "TOF vs. Cluster Energy vs. Cluster Size", "Cluster size", "Cluster Energy [GeV]",
                 "TOF [ns]", iterations_bins, BinSettings(500, 0, 1000), BinSettings(400, -100, 100));
    AddHistogram("tof_proton", "TOF proton", "#", "TOF [ns]", tof_bins);
    AddHistogram("E_vs_tof", "Energy vs. TOF", "TOF [ns]", "Energy [MeV]", tof_bins, BinSettings(1000, 0, 1000));
    AddHistogram("E_vs_tof_veto", "Energy vs. TOF w/ Veto", "TOF [ns]", "Energy [MeV]", tof_bins, BinSettings(1000, 0, 1000));
    AddHistogram("E/cl_vs_tof", "Energy/Cluster Size vs. TOF Veto", "TOF [ns]", "Energy/Cluster Size [MeV]", tof_bins,
                 BinSettings(1000, 0, 1000));
    AddHistogram("E/cl_vs_tof_veto", "Energy/Cluster Size vs. TOF w/ Veto", "TOF [ns]", "Energy/Cluster Size [MeV]",
                 tof_bins, BinSettings(1000, 0, 1000));
    AddHistogram("Ecenter/cl_vs_tof", "Central Cluster Energy/Cluster Size vs. TOF Veto", "TOF [ns]",
                 "Energy/Cluster Size [MeV]", tof_bins, BinSettings(1000, 0, 1000));
    AddHistogram("Ecenter/cl_vs_tof_veto", "Central Cluster Energy/Cluster Size vs. TOF w/ Veto", "TOF [ns]",
                 "Energy/Cluster Size [MeV]", tof_bins, BinSettings(1000, 0, 1000));
    AddHistogram("proton_id_BaF", "Energy/Cluster Size vs. TOF", "TOF [ns]", "Energy/Cluster Size [MeV]", tof_bins,
                 BinSettings(800, 0, 800));
    AddHistogram("proton_id_PbWO", "Energy/Cluster Size vs. TOF", "TOF [ns]", "Energy/Cluster Size [MeV]", tof_bins,
                 BinSettings(800, 0, 800));
    AddHistogram("proton_energies", "Short Energy vs. Cluster Energy", "Cluster [MeV]", "Short [MeV]",
                 energy_bins, theta_bins);
    AddHistogram("proton_psa", "PSA Angle vs. PSA Radius", "PSA Radius [MeV]", "PSA Angle [#circ]",
                 energy_bins, theta_bins);

    AddHistogram("q2_dist_before", "q^{2} Distribution before KinFit", "q^{2} [MeV]", "#", q2_bins);
    AddHistogram("q2_dist_after", "q^{2} Distribution after KinFit", "q^{2} [MeV]", "#", q2_bins);

    // im vs q2 histograms before and after several cuts
    AddHistogram("q2_im_base_selection", "q^{2} vs. IM - Base Selection", "IM [MeV]", "q^{2} [MeV]", im_bins, q2_bins);
    AddHistogram("q2_im_pid_cut", "q^{2} vs. IM - PID Cut", "IM [MeV]", "q^{2} [MeV]", im_bins, q2_bins);
    AddHistogram("q2_im_miss_mass", "q^{2} vs. IM - Missing Mass", "IM [MeV]", "q^{2} [MeV]", im_bins, q2_bins);
    AddHistogram("q2_im_before_fit", "q^{2} vs. IM - Before KinFit", "IM [MeV]", "q^{2} [MeV]", im_bins, q2_bins);
    AddHistogram("q2_im_after_fit", "q^{2} vs. IM - After KinFit", "IM [MeV]", "q^{2} [MeV]", im_bins, q2_bins);
    AddHistogram("q2_im_chi2_cut", "q^{2} vs. IM - #chi^{2} Cut", "IM [MeV]", "q^{2} [MeV]", im_bins, q2_bins);

    // different checks
    AddHistogram("lepton_energies", "Lepton Energies", "E(lepton1) [MeV]", "E(lepton2) [MeV]", energy_bins, energy_bins);
    AddHistogram("lepton_energies_true", "True Lepton Energies", "E(lepton1)_{true} [MeV]", "E(lepton2)_{true} [MeV]",
                 energy_bins, energy_bins);
    AddHistogram("photon_energy_vs_opening_angle", "Photon energy vs. opening angle", "Opening angle [#circ]",
                 "E_{#gamma} [MeV]", theta_bins, energy_bins);
    AddHistogram("photon_energy_vs_opening_angle_true", "True photon energy vs. opening angle", "Opening angle [#circ]",
                 "E_{#gamma} [MeV]", theta_bins, energy_bins);
    AddHistogram("theta_vs_clusters", "Theta vs. #Clusters", "#Clusters", "#vartheta [#circ]", particlecount_bins, theta_bins);
    AddHistogram("opening_angle_vs_q2", "Opening Angle Leptons vs. q^{2}", "q^{2} [GeV^{2}]", "Opening angle [#circ]",
                 q2_bins, theta_bins);
    AddHistogram("opening_angle_vs_E_high", "Opening Angle Leptons vs. E_{high}(Lepton)", "E [MeV]", "Opening angle [#circ]",
                 energy_bins, theta_bins);
    AddHistogram("opening_angle_vs_E_low", "Opening Angle Leptons vs. E_{low}(Lepton)", "E [MeV]", "Opening angle [#circ]",
                 energy_bins, theta_bins);
    AddHistogram("dEvE", "dE vs. E", "E_{Crystals} [MeV]", "dE_{Veto} [MeV]", energy_bins, veto_bins);
    AddHistogram("crystals_vs_ecl_charged", "#Crystals vs. Cluster Energy", "Cluster Energy [GeV]", "#Crystals",
                 q2_bins, count_bins);
    AddHistogram("crystals_vs_ecl_uncharged", "#Crystals vs. Cluster Energy", "Cluster Energy [GeV]", "#Crystals",
                 q2_bins, count_bins);
    AddHistogram("crystals_vs_ecl_charged_candidates", "#Crystals vs. Cluster Energy (Candidates)", "Cluster Energy [GeV]",
                 "#Crystals", q2_bins, count_bins);
    AddHistogram("energy_vs_momentum_z_balance", "Energy vs. p_{z} balance", "p_{z} [GeV]", "Energy [GeV]",
                 energy_range_bins, energy_range_bins);

    AddHistogram("opening_angle_leptons", "Lepton opening angle", "Opening angle [#circ]", "#", theta_bins);
    AddHistogram("opening_angle_leptons_true", "Lepton opening angle true", "Opening angle [#circ]", "#", theta_bins);
    AddHistogram("energy_lepton1", "Energy 1st lepton", "E [MeV]", "#", energy_bins);
    AddHistogram("energy_lepton1_true", "True energy 1st lepton", "E_{true} [MeV]", "#", energy_bins);
    AddHistogram("energy_lepton2", "Energy 2nd lepton", "E [MeV]", "#", energy_bins);
    AddHistogram("energy_lepton2_true", "True energy 2nd lepton", "E_{true} [MeV]", "#", energy_bins);
    AddHistogram("energy_photon", "Energy photon", "E [MeV]", "#", energy_bins);
    AddHistogram("energy_photon_true", "True energy photon", "E_{true} [MeV]", "#", energy_bins);
    // proton checks
    AddHistogram("proton_energy", "Energy proton", "E [MeV]", "#", energy_bins);
    AddHistogram("proton_energy_true", "Energy proton true", "E [MeV]", "#", energy_bins);
    AddHistogram("proton_energy_fit", "Energy proton fitted", "E [MeV]", "#", energy_bins);
    AddHistogram("proton_energy_delta", "#DeltaE_{proton} reconstructed - fitted", "E [MeV]", "#", energy_range_bins);
    AddHistogram("proton_angle_TAPS_expected", "Opening Angle reconstr. Cluster_{TAPS} - expected proton",
                 "opening angle [#circ]", "#", angle_bins);

    AddHistogram("coplanarity", "Coplanarity #eta' proton", "Coplanarity [#circ]", "#", phi_bins);
    AddHistogram("missing_mass", "Missing Mass Proton", "m_{miss} [MeV]", "#", energy_bins);

    // make fitter histograms
    AddHistogram("chisquare", "ChiSquare", "#chi^{2}", "#", chisquare_bins);
    AddHistogram("probability", "Probability", "Probability", "#", probability_bins);
    AddHistogram("iterations", "Number of iterations", "Iterations", "#", iterations_bins);

    stringstream fs;
    fs << "e+e-g";
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
    AddHistogram("proton_angle_TAPS_expected_cut", "Opening Angle reconstr. Cluster_{TAPS} - expected proton after Cuts",
                 "opening angle [#circ]", "#", angle_bins);
    AddHistogram("energy_vs_momentum_z_balance_cut", "Energy vs. p_{z} balance after Cuts", "p_{z} [GeV]", "Energy [GeV]",
                 energy_range_bins, energy_range_bins);
    AddHistogram("energy_vs_momentum_z_balance_fit_cut", "Energy vs. p_{z} balance fit after Ctus", "p_{z} [GeV]",
                 "Energy [GeV]", energy_range_bins, energy_range_bins);

    AddHistogram("dEvE_cut", "dE vs. E Cut", "E_{Crystals} [MeV]", "dE_{Veto} [MeV]", energy_bins, veto_bins);
    AddHistogram("crystals_vs_ecl_cut", "#Crystals vs. Cluster Energy Cut", "Cluster Energy [GeV]", "#Crystals",
                 q2_bins, count_bins);
    AddHistogram("crystals_vs_ecl_charged_cut", "#Crystals vs. Cluster Energy Cut", "Cluster Energy [GeV]",
                 "#Crystals", q2_bins, count_bins);
    AddHistogram("crystals_vs_ecl_uncharged_cut", "#Crystals vs. Cluster Energy Cut", "Cluster Energy [GeV]",
                 "#Crystals", q2_bins, count_bins);

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

ant::analysis::SaschaPhysics::SaschaPhysics(const mev_t energy_scale) :
    Physics("SaschaPhysics"),
    prompt("prompt", energy_scale),
    random("random", energy_scale),
    diff("diff", energy_scale),
    prompt_window(-6, 6),
    random_window1(-40, -15),
    random_window2(15, 40),
    fitter("SaschaPhysics"),
    final_state(nFinalState)
{
    cout << "Eta' Dalitz Physics:\n";
    cout << "Prompt window: " << prompt_window << " ns\n";
    cout << "Random window 1: " << random_window1 << " ns\n";
    cout << "Random window 2: " << random_window2 << " ns\n";

    // histogram to count number of different particle types
    for (auto& t : ParticleTypeDatabase::DetectableTypes())
        numParticleType[t]= HistFac.makeTH1D("Number of " + t->PrintName(),
                                      "number of " + t->PrintName() + "/ event",
                                      "", BinSettings(16,0,16), "n_" + t->PrintName());

    // histogram to count how many events passed the different selection criteria
    accepted_events = HistFac.makeTH1D("Accepted events", "condition", "#events", BinSettings(10), "accepted_events");

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

    // setup fitter for Dalitz decay
    fitter.LinkVariable("Beam", beam.Link(), beam.LinkSigma());
    fitter.LinkVariable("Lepton1", final_state[0].Link(), final_state[0].LinkSigma());
    fitter.LinkVariable("Lepton2", final_state[1].Link(), final_state[1].LinkSigma());
    fitter.LinkVariable("Photon", final_state[2].Link(), final_state[2].LinkSigma());
    fitter.LinkVariable("Proton", final_state[3].Link(), final_state[3].LinkSigma());

    vector<string> all_names = {"Beam", "Lepton1", "Lepton2", "Photon", "Proton"};

    // Constraint: Incoming 4-vector = Outgoing 4-vector
    auto EnergyMomentumBalance = [] (const vector<vector<double>>& particles) -> vector<double>
    {
        const TLorentzVector target(0,0,0, ParticleTypeDatabase::Proton.Mass());
        // assume first particle is beam photon
        TLorentzVector diff = target + FitParticle::Make(particles[0], ParticleTypeDatabase::Photon.Mass());
        // assume second and third particle outgoing leptons
        diff -= FitParticle::Make(particles[1], ParticleTypeDatabase::eMinus.Mass());
        diff -= FitParticle::Make(particles[2], ParticleTypeDatabase::eMinus.Mass());
        // assume fourth particle outgoing photon
        diff -= FitParticle::Make(particles[3], ParticleTypeDatabase::Photon.Mass());
        // assume last particle outgoing proton
        diff -= FitParticle::Make(particles[4], ParticleTypeDatabase::Proton.Mass());

        return {diff.X(), diff.Y(), diff.Z(), diff.T()};
    };
    fitter.AddConstraint("EnergyMomentumBalance", all_names, EnergyMomentumBalance);

    // Constraint: Coplanarity between eta' and recoil proton
    auto CoplanarityConstraint = [] (const vector<vector<double>>& particles) -> double
    {
        TLorentzVector etap = FitParticle::Make(particles[1], ParticleTypeDatabase::eMinus.Mass());
        etap += FitParticle::Make(particles[2], ParticleTypeDatabase::eMinus.Mass());
        etap += FitParticle::Make(particles[3], ParticleTypeDatabase::Photon.Mass());
        TLorentzVector proton = FitParticle::Make(particles[4], ParticleTypeDatabase::Proton.Mass());
        return abs(etap.Phi() - proton.Phi())*TMath::RadToDeg() - 180.;
    };
    if (includeCoplanarityConstraint)
        fitter.AddConstraint("CoplanarityConstraint", all_names, CoplanarityConstraint);

    // Constraint: Invariant mass of nPhotons equals constant IM,
    // make lambda catch also this with [&] specification
    auto RequireIM = [&] (const vector<vector<double>>& particles) -> double
    {
        TLorentzVector sum(0,0,0,0);
        sum += FitParticle::Make(particles[1], ParticleTypeDatabase::eMinus.Mass());
        sum += FitParticle::Make(particles[2], ParticleTypeDatabase::eMinus.Mass());
        sum += FitParticle::Make(particles[3], ParticleTypeDatabase::Photon.Mass());

        return sum.M() - IM;
    };
    if (includeIMconstraint)
        fitter.AddConstraint("RequireIM", all_names, RequireIM);

    // Constraint: Vertex position in z direction: v_z (positive if upstream)
    // if the particle originated from (0,0,v_z) instead of origin,
    // the corrected angle theta' is given by
    // tan(theta') = (R sin(theta))/(R cos(theta) - v_z)
    // R is the CB radius, 10in aka 25.4cm

    auto VertexConstraint = [&] (vector<vector<double>>& particles) -> double
    {
        constexpr double R = 25.4;
        // last element in particles is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = particles.back()[0];
        particles.resize(particles.size()-1); // get rid of last element
        // correct each particle's theta angle,
        // then calculate invariant mass of all particles
        TLorentzVector sum(0,0,0,0);
        size_t n = 0;
        for (auto i = particles.begin()+1; i != particles.end()-1; ++i) {
            const double theta = (*i)[1]; // second element is theta
            const double theta_p = std::atan2( R*sin(theta), R*cos(theta) - v_z);
            (*i)[1] = theta_p;
            if (n < 2)  // first two particles are the electron and positron
                sum += FitParticle::Make(*i, ParticleTypeDatabase::eMinus.Mass());
            else  // last particle is the photon
                sum += FitParticle::Make(*i, ParticleTypeDatabase::Photon.Mass());
            n++;
        }
        return sum.M() - IM;
    };

    if (includeVertexFit) {
        fitter.AddUnmeasuredVariable("v_z"); // default value 0
        fitter.AddConstraint("VertexConstraint", all_names + std::vector<string>{"v_z"}, VertexConstraint);
    }

    static_assert(!(includeIMconstraint && includeVertexFit), "Do not enable Vertex and IM Fit at the same time");

    // create pull histograms
    const BinSettings pull_bins(50,-3,3);
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
    settings.MaxIterations = 8;
    fitter.SetSettings(settings);

    cout.precision(3);
    APLCON::PrintFormatting::Width = 11;
}

void ant::analysis::SaschaPhysics::ProcessEvent(const ant::Event &event)
{
    accepted_events->Fill("all events", 1);

    prompt["cb_esum"]->Fill(event.Reconstructed().TriggerInfos().CBEenergySum());

    if (event.Reconstructed().TriggerInfos().CBEenergySum() < cb_esum)
        return;

    accepted_events->Fill("CB Esum", 1);

    for (const auto& track : event.Reconstructed().Tracks()) {
        prompt["pid"]->Fill(track->ClusterEnergy(), track->VetoEnergy());
    }

    // determine #clusters per event per detector
    size_t nCB = 0, nTAPS = 0;
    std::for_each(event.Reconstructed().Tracks().begin(), event.Reconstructed().Tracks().end(),
                  [&nCB, &nTAPS](const auto& track){ track->Detector() & detector_t::anyCB ? nCB++ : nTAPS++; });
    prompt["n_cluster_cb"]->Fill(nCB);
    prompt["n_cluster_taps"]->Fill(nTAPS);

    for (auto& particle : event.Reconstructed().Particles().GetAll()) {
        prompt["particle_types"]->Fill(particle->Type().PrintName().c_str(), 1);
    }

    for (auto& t : ParticleTypeDatabase::DetectableTypes()) {
        try {
            numParticleType.at(t)->Fill(event.Reconstructed().Particles().Get(*t).size());
        } catch (...) {}
    }

    TaggerHistList tagger_hits;
    size_t prompt_n_tagger = 0, random_n_tagger = 0;
    //static size_t count = 0;
    const bool MC = false;
    if (MC)
        tagger_hits = event.MCTrue().TaggerHits();
    else
        tagger_hits = event.Reconstructed().TaggerHits();

    for (const auto& taggerhit : tagger_hits) {
        prompt["tagger_spectrum"]->Fill(taggerhit->PhotonEnergy());
        prompt["tagger_time"]->Fill(taggerhit->Time());

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

        // proton tests (TOF, PSA)
        for (auto& track : event.Reconstructed().Tracks()) {
            double tof = taggerhit->Time() - track->Time();
            double shortE = track->ShortEnergy(), clusterE = track->ClusterEnergy();
            const double r2d = TMath::RadToDeg();
            if (is_prompt) {
                if (track->Detector() & ant::detector_t::anyTAPS) {
                    dynamic_cast<TH3D*>(prompt["proton_id"])->Fill(track->ClusterSize(), clusterE, tof);
                    prompt["proton_energies"]->Fill(clusterE, shortE);
                    prompt["E_vs_tof"]->Fill(tof, clusterE);
                    prompt["E/cl_vs_tof"]->Fill(tof, clusterE/track->ClusterSize());
                    if (track->Detector() & detector_t::Veto) {
                        prompt["E_vs_tof_veto"]->Fill(tof, clusterE);
                        prompt["E/cl_vs_tof_veto"]->Fill(tof, clusterE/track->ClusterSize());
                    }
                    if (track->Theta()*TMath::RadToDeg() > 5.)
                        prompt["proton_id_BaF"]->Fill(tof, clusterE/track->ClusterSize());
                    else
                        prompt["proton_id_PbWO"]->Fill(tof, clusterE/track->ClusterSize());
                    //prompt["Ecenter/cl_vs_tof"]->Fill(tof, track->CentralCrystal()/track->ClusterSize());
                    //prompt["Ecenter/cl_vs_tof_veto",]->Fill(tof, track->CentralCrystal()/track->ClusterSize());
                    prompt["proton_psa"]->Fill(sqrt(shortE*shortE + clusterE*clusterE), atan2(clusterE, shortE)*r2d);
                } else if (track->Detector() & detector_t::anyCB) {
                    dynamic_cast<TH3D*>(random["proton_id"])->Fill(track->ClusterSize(), clusterE, tof);
                    random["proton_energies"]->Fill(clusterE, shortE);
                    random["E_vs_tof"]->Fill(tof, clusterE);
                    random["E/cl_vs_tof"]->Fill(tof, clusterE/track->ClusterSize());
                    if (track->Detector() & detector_t::PID) {
                        random["E_vs_tof_veto"]->Fill(tof, clusterE);
                        random["E/cl_vs_tof_veto"]->Fill(tof, clusterE/track->ClusterSize());
                    }
                    //random["Ecenter/cl_vs_tof"]->Fill(tof, track->CentralCrystal()/track->ClusterSize());
                    //random["Ecenter/cl_vs_tof_veto",]->Fill(tof, track->CentralCrystal()/track->ClusterSize());
                    random["proton_psa"]->Fill(sqrt(shortE*shortE + clusterE*clusterE), atan2(clusterE, shortE)*r2d);
                }
            }
        }

        // make sure the correct histogram will be filled
        HistList& h = is_prompt ? prompt : random;

        h["n_part"]->Fill(event.Reconstructed().Particles().GetAll().size());
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
        particle_vector particles;
        GetParticles(event, particles);
        if (nParticles != nFinalState)
            continue;

        accepted_events->Fill("n_FS", 1);

        // check if we have the right amount of leptons, photons and protons
        unsigned short ne = 0, ng = 0, np = 0;
        for (auto it = particles.begin(); it != particles.end(); ++it) {
            if (it->Type() == ParticleTypeDatabase::Photon)
                ng++;
            else if (it->Type() == ParticleTypeDatabase::eCharged)
                ne++;
            else if (it->Type() == ParticleTypeDatabase::Proton)
                np++;
        }
        if (ne != 2 || ng != 1 || np != 1)
            continue;

        accepted_events->Fill("#part", 1);

        // check theta distribution in dependence of the number of clusters
        TrackPtr tr;
        bool highVeto = false;
        for (const auto& p : particles) {
            tr = p.Tracks().front();
            if (tr->VetoEnergy() > 3.)
                highVeto = true;
            h["theta_vs_clusters"]->Fill(nParticles, p.Theta()*TMath::RadToDeg());
            h["crystals_vs_ecl_charged_candidates"]->Fill(tr->ClusterEnergy(), tr->ClusterSize());
            h["dEvE"]->Fill(tr->ClusterEnergy(), tr->VetoEnergy());
            if (p.Type().Charged())
                h["crystals_vs_ecl_charged"]->Fill(tr->ClusterEnergy(), tr->ClusterSize());
            else
                h["crystals_vs_ecl_uncharged"]->Fill(tr->ClusterEnergy(), tr->ClusterSize());
        }

//        cout << "We have the right 4 particles\nNow sort them, order before:\n";
//        for (const auto& p : particles)
//            cout << p.Type() << endl;

        // swap the proton and the last vector entry if the proton is not the last particle
        for (auto it = particles.begin(); it != particles.end()-1; ++it)
            if (it->Type() == ParticleTypeDatabase::Proton) {
                std::iter_swap(it, particles.end()-1);
                break;
            }
        // swap the photon at the third position
        for (auto it = particles.begin(); it != particles.end()-2; ++it)
            if (it->Type() == ParticleTypeDatabase::Photon) {
                std::iter_swap(it, particles.end()-2);
                break;
            }

        TLorentzVector proton = particles.back();
        TLorentzVector etap(0., 0., 0., 0.);
        for (auto it = particles.cbegin(); it != particles.cend()-1; ++it)
            etap += *it;
        double q2 = (particles.at(0) + particles.at(1)).M();
        h["q2_im_base_selection"]->Fill(etap.M(), q2);

        // perform cut on PID elements to suppress conversion lepton pairs
        element_index_t pid1 = particles.at(0).Tracks().front()->CentralVeto(),
                        pid2 = particles.at(1).Tracks().front()->CentralVeto();
        if (pid1 == pid2)
            continue;
        h["q2_im_pid_cut"]->Fill(etap.M(), q2);
        accepted_events->Fill("PID Cut", 1);

        h["tof_proton"]->Fill(taggerhit->Time() - particles.back().Tracks().front()->Time());
        //prompt["E_vs_tof"]->Fill(taggerhit->Time() - particles.back().Tracks().front()->Time(),
        //                         particles.back().Tracks().front()->ClusterEnergy());

        const TLorentzVector target(0., 0., 0., ParticleTypeDatabase::Proton.Mass());
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

        h["proton_energy"]->Fill(proton.T());
        double openAngle_p_TAPS_expected = proton.Angle(missingProton.Vect())*TMath::RadToDeg();
        h["proton_angle_TAPS_expected"]->Fill(openAngle_p_TAPS_expected);

        //std::copy(particles.begin(), particles.end(), final_state);
//        std::transform(particles.begin(), particles.end(), final_state.begin(),
//                       [](Particle& p) -> double { return FitParticle().SetFromVector(p); });
        auto it = final_state.begin();
        for (const auto& p : particles)
            (it++)->SetFromVector(p);

        if (!(event.MCTrue().Particles().GetAll().empty())) {
            // use the first particle in the particles vector as a placeholder
            Particle ePlus_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
            Particle eMinus_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
            Particle photon_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
            Particle proton_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
            for (const auto& p : event.MCTrue().Particles().GetAll()) {
                if (p->Type() == ParticleTypeDatabase::Photon)
                    photon_true = Particle(*p);
                else if (p->Type() == ParticleTypeDatabase::eMinus)
                    eMinus_true = Particle(*p);
                else if (p->Type() == ParticleTypeDatabase::ePlus)
                    ePlus_true = Particle(*p);
                else if (p->Type() == ParticleTypeDatabase::Proton)
                    proton_true = Particle(*p);
            }

            double lepton_open_angle_true = ePlus_true.Angle(eMinus_true.Vect())*TMath::RadToDeg();
            double en_lep1_true, en_lep2_true;
            if (eMinus_true.Ek() > ePlus_true.Ek()) {
                en_lep1_true = eMinus_true.Ek();
                en_lep2_true = ePlus_true.Ek();
            } else {
                en_lep1_true = ePlus_true.Ek();
                en_lep2_true = eMinus_true.Ek();
            }
            if ((ePlus_true.Type() == ParticleTypeDatabase::eCharged)
                    && (eMinus_true.Type() == ParticleTypeDatabase::eCharged)) {
                h["opening_angle_leptons_true"]->Fill(lepton_open_angle_true);
                h["photon_energy_vs_opening_angle_true"]->Fill(lepton_open_angle_true, photon_true.Ek());
                h["lepton_energies_true"]->Fill(en_lep1_true, en_lep2_true);
                h["energy_lepton1_true"]->Fill(en_lep1_true);
                h["energy_lepton2_true"]->Fill(en_lep2_true);
            }
            h["energy_photon_true"]->Fill(photon_true.Ek());
            h["proton_energy_true"]->Fill(proton_true.Ek());
        }

        double q2_before = (particles[0] + particles[1]).M();
        double lepton_open_angle = particles[0].Angle(particles[1].Vect())*TMath::RadToDeg();
        double en_lep1 = particles[0].T();
        double en_lep2 = particles[1].T();
        if (en_lep1 < en_lep2) {
            en_lep1 = en_lep2;
            en_lep2 = particles[0].T();
        }
        h["opening_angle_leptons"]->Fill(lepton_open_angle);
        h["lepton_energies"]->Fill(en_lep1, en_lep2);
        h["energy_lepton1"]->Fill(en_lep1);
        h["energy_lepton2"]->Fill(en_lep2);
        h["energy_photon"]->Fill(particles[2].E());
        h["photon_energy_vs_opening_angle"]->Fill(lepton_open_angle, particles[2].E());
        h["opening_angle_vs_q2"]->Fill(q2_before, lepton_open_angle);
        h["opening_angle_vs_E_high"]->Fill(en_lep1 > en_lep2 ? en_lep1 : en_lep2, lepton_open_angle);
        h["opening_angle_vs_E_low"]->Fill(en_lep1 < en_lep2 ? en_lep1 : en_lep2, lepton_open_angle);
        h["q2_dist_before"]->Fill(q2_before);

        // set proton energy sigma to zero to indicate it's unmeasured
        final_state[3].Ek_Sigma = 0;

//        std::transform(myv1.begin(), myv1.end(), myv1.begin(), [](double d) -> double { return d * 3; });
//        for_each(begin(myv1), end(myv1), [](double& a) { a *= 3; });

        beam.SetFromVector(taggerhit->PhotonBeam());
        // set beam energy sigma to 2 MeV
        beam.Ek_Sigma = 2.;

        FillIM(h["im_true"], final_state);

//        // smear the MC true data
//        proton.Smear();
//        for(auto& photon : photons)
//            photon.Smear();
//        beam.Smear();

//        FillIM(h["im_smeared"], photons);

        if (missM < 880. && missM > 1130.)
            continue;
        accepted_events->Fill("miss mass", 1);

        h["q2_im_miss_mass"]->Fill(etap.M(), q2);
        h["q2_im_before_fit"]->Fill(etap.M(), q2);
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
        TLorentzVector proton_fit = FitParticle::Make(final_state[3], ParticleTypeDatabase::Proton.Mass());
        TLorentzVector etap_fit(0., 0., 0., 0.);
        etap_fit += FitParticle::Make(final_state[0], ParticleTypeDatabase::eMinus.Mass());
        etap_fit += FitParticle::Make(final_state[1], ParticleTypeDatabase::eMinus.Mass());
        etap_fit += FitParticle::Make(final_state[2], ParticleTypeDatabase::Photon.Mass());
        balanceP4_fit -= etap_fit;
        balanceP4_fit -= proton_fit;
        h["energy_vs_momentum_z_balance_fit"]->Fill(balanceP4_fit.Pz(), balanceP4_fit.E());
        h["q2_im_after_fit"]->Fill(etap_fit.M(), q2_after);

        const double copl_fit = abs(etap_fit.Phi() - particles.back().Phi())*TMath::RadToDeg();
        h["coplanarity_fit"]->Fill(copl_fit);
        TLorentzVector missingProtonFit = target + FitParticle::Make(beam, ParticleTypeDatabase::Photon.Mass()) - etap_fit;
        h["missing_mass_fit"]->Fill(missingProtonFit.M());

        h["proton_energy_fit"]->Fill(proton_fit.T());
        h["proton_energy_delta"]->Fill(proton.T() - proton_fit.T());

        // Cut on chi^2
        if (result.ChiSquare > 10.)
            return;
        h["q2_im_chi2_cut"]->Fill(etap_fit.M(), q2_after);

        accepted_events->Fill("#chi^{2} cut", 1);

        // fill the invM histograms for different q2 ranges
        if (q2_after < im_q2_upper_bound) {
            if (is_prompt)
                im_q2_prompt.at(static_cast<int>(q2_after/im_q2_mev_steps))->Fill(etap_fit.M());
            else
                im_q2_random.at(static_cast<int>(q2_after/im_q2_mev_steps))->Fill(etap_fit.M());
        }

        h["im_cut"]->Fill(etap.M());
        FillIM(h["im_fit_cut"], final_state);
        h["q2_dist_cut"]->Fill(q2_before);
        h["q2_dist_fit_cut"]->Fill(q2_after);
        h["coplanarity_cut"]->Fill(copl);
        h["missing_mass_cut"]->Fill(missM);
        h["missing_mass_fit_cut"]->Fill(missingProtonFit.M());
        h["coplanarity_fit_cut"]->Fill(copl_fit);
        h["proton_angle_TAPS_expected_cut"]->Fill(openAngle_p_TAPS_expected);
        h["energy_vs_momentum_z_balance_cut"]->Fill(balanceP4.Pz(), balanceP4.E());
        h["energy_vs_momentum_z_balance_fit_cut"]->Fill(balanceP4_fit.Pz(), balanceP4_fit.E());

        for (const auto& p : particles) {
            tr = p.Tracks().front();
            h["crystals_vs_ecl_cut"]->Fill(tr->ClusterEnergy(), tr->ClusterSize());
            h["dEvE_cut"]->Fill(tr->ClusterEnergy(), tr->VetoEnergy());
            if (p.Type().Charged())
                h["crystals_vs_ecl_charged_cut"]->Fill(tr->ClusterEnergy(), tr->ClusterSize());
            else
                h["crystals_vs_ecl_uncharged_cut"]->Fill(tr->ClusterEnergy(), tr->ClusterSize());
        }
    }
    prompt["n_tagged"]->Fill(prompt_n_tagger);
    random["n_tagged"]->Fill(random_n_tagger);
}

void ant::analysis::SaschaPhysics::Finish()
{
    double factor = -prompt_window.Length() / (random_window1.Length() + random_window2.Length());
    diff = prompt;
    diff.AddScaled(random, factor);
}

void ant::analysis::SaschaPhysics::ShowResult()
{
    canvas c("SaschaPhysics: Overview");
    c << drawoption("colz") << diff["pid"]
      << padoption::set(padoption_t::Legend)
      << diff["particle_types"]
      << padoption::unset(padoption_t::Legend)
      << diff["tagger_spectrum"] << diff["n_tagged"] << diff["cb_esum"] << endc;

    canvas c_pulls("SaschaPhysics: Pulls");
    c_pulls << padoption::set(padoption_t::LogY);
    for (auto& p : pulls_diff)
        c_pulls << p.second;
    c_pulls << endc;

    canvas c_fitter("SaschaPhysics: Fitter");
    c_fitter <<  diff["chisquare"] <<  diff["probability"] << diff["iterations"]
             <<  diff["im_true"] <<  diff["im_smeared"] <<  diff["im_fit"]
             <<  diff["vertex_z_before"] <<  diff["vertex_z_after"] << endc;

    // draw all random subtracted histograms
    //diff.Draw();
}

void ant::analysis::SaschaPhysics::GetParticles(const ant::Event& event, particle_vector& particles)
{
    // example of how to collect particles in a user defined way
    for (const auto& track : event.Reconstructed().Tracks()) {
        // ignore particles below set cluster energy threshold
        if (track->ClusterEnergy() < cluster_thresh)
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

void ant::analysis::SaschaPhysics::GetTrueParticles(const ant::Event& event, particle_vector& particles)
{
    for (const auto& p : event.MCTrue().Particles().GetAll())
        if (contains(ParticleTypeDatabase::MCFinalStateTypes(), p->Type()))
            particles.push_back(Particle(*p));
    nParticles = particles.size();
}
