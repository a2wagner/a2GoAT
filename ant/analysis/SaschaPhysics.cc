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
    const Double_t p = sqrt(E*E - m*m);
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

    const BinSettings energy_bins(1600, 0, energy_scale*1.6);
    const BinSettings tagger_bins(200, 1400, 1600);
    const BinSettings veto_bins(1000, 0, 10);
    const BinSettings chisquare_bins(100, 0, 30);
    const BinSettings probability_bins(100, 0, 1);
    const BinSettings iterations_bins(15, 0, 15);
    const BinSettings im_bins(200, IM-100, IM+100);
    const BinSettings vertex_bins(200, -10, 10);
    const BinSettings theta_bins(720, 0, 180);
    const BinSettings phi_bins(720, 0, 360);
    const BinSettings angle_bins(500, 0, 50);
    const BinSettings energy_range_bins(1600, -1600, 1600);
    const BinSettings q2_bins(2000, 0, 1000);
    const BinSettings count_bins(50, 0, 50);
    const BinSettings tof_bins(1000, -50, 50);

    AddHistogram("tagger_energy", "Tagger Energy", "Photon Beam Energy [MeV]", "#", tagger_bins);
    AddHistogram("n_tagged", "Tagger Hits", "Tagger Hits / event", "#", count_bins);

    // proton tests
    //AddHistogram("proton_id", "TOF vs. Cluster Energy vs. Cluster Size", "Cluster size", "Cluster Energy [GeV]",
    //             "TOF [ns]", iterations_bins, BinSettings(500, 0, 1000), BinSettings(400, -100, 100));
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
    AddHistogram("photon_energy_vs_opening_angle", "Photon energy vs. opening angle", "Opening angle [#circ]",
                 "E_{#gamma} [MeV]", theta_bins, energy_bins);
    AddHistogram("theta_vs_clusters", "Theta vs. #Clusters", "#Clusters", "#vartheta [#circ]", BinSettings(10), theta_bins);
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
    AddHistogram("energy_photon_cm", "Energy Photon CM", "E [MeV]", "#", energy_bins);

    AddHistogram("opening_angle_leptons", "Lepton opening angle", "Opening angle [#circ]", "#", theta_bins);
    AddHistogram("energy_lepton1", "Energy 1st lepton", "E [MeV]", "#", energy_bins);
    AddHistogram("energy_lepton2", "Energy 2nd lepton", "E [MeV]", "#", energy_bins);
    AddHistogram("energy_photon", "Energy photon", "E [MeV]", "#", energy_bins);
    // proton checks
    AddHistogram("proton_energy", "Energy proton", "E [MeV]", "#", energy_bins);
    AddHistogram("proton_energy_fit", "Energy proton fitted", "E [MeV]", "#", energy_bins);
    AddHistogram("proton_energy_delta", "#DeltaE_{proton} reconstructed - fitted", "E [MeV]", "#", energy_range_bins);
    AddHistogram("proton_angle_TAPS_expected", "Opening Angle reconstr. Cluster_{TAPS} - expected proton",
                 "opening angle [#circ]", "#", angle_bins);
    AddHistogram("protons_found", "#Protons in predicted cone", "#protons", "#", BinSettings(5));
    AddHistogram("nTAPS_vs_protons_found", "#clusters in TAPS vs. #protons in predicted cone", "#protons", "#clusters TAPS",
                 BinSettings(5), BinSettings(5));
    AddHistogram("protons_diff_clusters", "#different TAPS clusters identified as proton", "#clusters", "#",BinSettings(5));

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
    cb_time_window(-20, 20),
    taps_time_window(-20, 20),//(-10, 10), //TODO: for now bigger (TOF), check what is best
    fitter("SaschaPhysics"),
    final_state(nFinalState),
    detectors(nFinalState + 1),  // beam photon needs to be considered
    etap_fit("eta' fitter"),
    etap_fs(3)
{
    cout << "Eta' Dalitz Physics:\n";
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
    const BinSettings chi2_bins(150, 0, 50);
    const BinSettings q2_bins(2000, 0, 1000);
    const BinSettings time_bins(1000, -50, 50);
    const BinSettings im_bins(200, IM-100, IM+100);
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
    etap_chi2 = HistFac.makeTH1D("#chi^{2} Fit #eta' (particle selection)", "#chi^{2}", "#", chi2_bins, "etap_chi2");
    etap_chi2_vs_q2 = HistFac.makeTH2D("#chi^{2} Fit #eta' vs. q^{2} (particle selection)", "q^{2} [MeV]", "#chi^{2}",
                                       q2_bins, chi2_bins, "etap_chi2_vs_q2");
    q2_vs_im_candidates = HistFac.makeTH2D("q^{2} vs. IM - Candidates", "IM [MeV]", "q^{2} [MeV]", im_bins, q2_bins,
                                           "q2_vs_im_candidates");
    // pi0 background tests
    sub_pi0_cand_im = HistFac.makeTH1D("#pi^{0} Candidate Invariant Mass", "IM [MeV]", "#", q2_bins, "sub_pi0_cand_im");
    sub_pi0_cand_im_cut = HistFac.makeTH1D("#pi^{0} Candidate Invariant Mass Cut", "IM [MeV]", "#", q2_bins,
                                           "sub_pi0_cand_im_cut");
    q2_vs_im_pi0_cut = HistFac.makeTH2D("q^{2} vs. IM - Candidates", "IM [MeV]", "q^{2} [MeV]", im_bins, q2_bins,
                                        "q2_vs_im_pi0_cut");

    // true MC checks
    const BinSettings true_energy_bins(1200, 0, energy_scale*1.2);
    const BinSettings theta_bins(720, 0, 180);
    eplus_true_theta_vs_energy = HistFac.makeTH2D("True e^{+} #vartheta vs. Energy", "Energy [MeV]", "#vartheta [#circ]",
                                                  true_energy_bins, theta_bins, "eplus_true_theta_vs_energy");
    eminus_true_theta_vs_energy = HistFac.makeTH2D("True e^{-} #vartheta vs. Energy", "Energy [MeV]", "#vartheta [#circ]",
                                                   true_energy_bins, theta_bins, "eminus_true_theta_vs_energy");
    photon_true_theta_vs_energy = HistFac.makeTH2D("True #gamma #vartheta vs. Energy", "Energy [MeV]", "#vartheta [#circ]",
                                                   true_energy_bins, theta_bins, "photon_true_theta_vs_energy");
    proton_true_theta_vs_energy = HistFac.makeTH2D("True p #vartheta vs. Energy", "Energy [MeV]", "#vartheta [#circ]",
                                                   BinSettings(600, 0, 600), BinSettings(100, 0, 25),
                                                   "proton_true_theta_vs_energy");
    true_open_angle_vs_q2 = HistFac.makeTH2D("True Lepton Opening Angle vs. q^{2}", "q^{2} [MeV]", "Angle [#circ]",
                                             q2_bins, theta_bins, "true_open_angle_vs_q2");
    opening_angle_leptons_true = HistFac.makeTH1D("True Lepton Opening Angle", "Opening angle [#circ]", "#", theta_bins,
                                                  "opening_angle_leptons_true");
    opening_angle_vs_photon_energy_true = HistFac.makeTH2D("True Lepton Opening Angle vs. #gamma Energy",
                                                           "E_{#gamma} [MeV]", "Angle [#circ]", true_energy_bins, theta_bins,
                                                           "opening_angle_vs_photon_energy_true");
    lepton_energies_true = HistFac.makeTH2D("True Lepton Energies", "E(lepton1)_{true} [MeV]", "E(lepton2)_{true} [MeV]",
                                            true_energy_bins, true_energy_bins, "lepton_energies_true");
    energy_lepton1_true = HistFac.makeTH1D("True energy higher energetic lepton", "E_{true} [MeV]", "#", true_energy_bins,
                                           "energy_lepton1_true");
    energy_lepton2_true = HistFac.makeTH1D("True energy lower energetic lepton", "E_{true} [MeV]", "#", true_energy_bins,
                                           "energy_lepton2_true");
    energy_photon_cm_true = HistFac.makeTH1D("True Photon CM Energy", "E_{true} [MeV]", "#", true_energy_bins,
                                             "energy_photon_cm_true");
    n_cluster_cb_vs_q2 = HistFac.makeTH2D("Number of Clusters in CB vs. q^{2}", "q^{2} [MeV]", "#clusters",
                                          q2_bins, BinSettings(6), "n_cluster_cb_vs_q2");
    n_cluster_taps_vs_q2 = HistFac.makeTH2D("Number of Clusters in TAPS vs. q^{2}", "q^{2} [MeV]", "#clusters",
                                            q2_bins, BinSettings(5), "n_cluster_taps_vs_q2");
    n_cluster_cb_vs_open_angle = HistFac.makeTH2D("True p #vartheta vs. Energy", "Angle [#circ]", "#clusters",
                                                  theta_bins, BinSettings(6), "n_cluster_cb_vs_open_angle");
    n_cluster_taps_vs_open_angle = HistFac.makeTH2D("True p #vartheta vs. Energy", "Angle [#circ]", "#clusters",
                                                    theta_bins, BinSettings(5), "n_cluster_taps_vs_open_angle");
    // for Achims proton id test
    expected_proton_diff_vs_q2 = HistFac.makeTH2D("Expected Proton Difference vs. q^{2}", "q^{2} [MeV]", "Angle [#circ]",
                                                  q2_bins, BinSettings(100, 0, 25), "expected_proton_diff_vs_q2");
    expected_proton_diff_vs_q2_rebin = HistFac.makeTH2D("Expected Proton Difference vs. q^{2}", "q^{2} [MeV]", "Angle [#circ]",
                                                        BinSettings(20, 0, 1000), BinSettings(100, 0, 25),
                                                        "expected_proton_diff_vs_q2_rebin");
    protons_found = HistFac.makeTH1D("#Protons in predicted cone", "protons", "#", BinSettings(5), "protons_found");


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
    if (includeVertexFit) {
        all_names.push_back("v_z");
        //fitter.AddUnmeasuredVariable("v_z"); // default value 0
        auto v_z_settings = APLCON::Variable_Settings_t::Default;
        v_z_settings.Limit.High = TARGET_MAX;
        v_z_settings.Limit.Low = TARGET_MIN;
        fitter.AddMeasuredVariable("v_z", 4., 2.3, v_z_settings);  // default value 0
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
    auto EnergyMomentumBalance = [recalculate_cluster] (const vector<vector<double>>& particles) -> vector<double>
    {
        const TLorentzVector target(0,0,0, ParticleTypeDatabase::Proton.Mass());
        TLorentzVector diff(0,0,0,0);
        if (includeVertexFit) {
            // assume first particle is beam photon
            diff = target + FitParticle::Make(particles[0], ParticleTypeDatabase::Photon.Mass());
            /* from here all clusters are recalculated by taking a shifted z vertex position into account */
            // assume second and third particle outgoing leptons
            diff -= FitParticle::Make(recalculate_cluster(particles, 1), ParticleTypeDatabase::eMinus.Mass());
            diff -= FitParticle::Make(recalculate_cluster(particles, 2), ParticleTypeDatabase::eMinus.Mass());
            // assume fourth particle outgoing photon
            diff -= FitParticle::Make(recalculate_cluster(particles, 3), ParticleTypeDatabase::Photon.Mass());
            // assume last particle outgoing proton
            diff -= FitParticle::Make(recalculate_cluster(particles, 4), ParticleTypeDatabase::Proton.Mass());
        } else {
            // assume first particle is beam photon
            diff = target + FitParticle::Make(particles[0], ParticleTypeDatabase::Photon.Mass());
            // assume second and third particle outgoing leptons
            diff -= FitParticle::Make(particles[1], ParticleTypeDatabase::eMinus.Mass());
            diff -= FitParticle::Make(particles[2], ParticleTypeDatabase::eMinus.Mass());
            // assume fourth particle outgoing photon
            diff -= FitParticle::Make(particles[3], ParticleTypeDatabase::Photon.Mass());
            // assume last particle outgoing proton
            diff -= FitParticle::Make(particles[4], ParticleTypeDatabase::Proton.Mass());
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
            etap = FitParticle::Make(recalculate_cluster(particles, 1), ParticleTypeDatabase::eMinus.Mass());
            etap += FitParticle::Make(recalculate_cluster(particles, 2), ParticleTypeDatabase::eMinus.Mass());
            etap += FitParticle::Make(recalculate_cluster(particles, 3), ParticleTypeDatabase::Photon.Mass());
            proton = FitParticle::Make(recalculate_cluster(particles, 4), ParticleTypeDatabase::Proton.Mass());
        } else {
            etap = FitParticle::Make(particles[1], ParticleTypeDatabase::eMinus.Mass());
            etap += FitParticle::Make(particles[2], ParticleTypeDatabase::eMinus.Mass());
            etap += FitParticle::Make(particles[3], ParticleTypeDatabase::Photon.Mass());
            proton = FitParticle::Make(particles[4], ParticleTypeDatabase::Proton.Mass());
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
            sum += FitParticle::Make(recalculate_cluster(particles, 1), ParticleTypeDatabase::eMinus.Mass());
            sum += FitParticle::Make(recalculate_cluster(particles, 2), ParticleTypeDatabase::eMinus.Mass());
            sum += FitParticle::Make(recalculate_cluster(particles, 3), ParticleTypeDatabase::Photon.Mass());
        } else {
            sum += FitParticle::Make(particles[1], ParticleTypeDatabase::eMinus.Mass());
            sum += FitParticle::Make(particles[2], ParticleTypeDatabase::eMinus.Mass());
            sum += FitParticle::Make(particles[3], ParticleTypeDatabase::Photon.Mass());
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
    settings.MaxIterations = 8;
    fitter.SetSettings(settings);

    cout.precision(3);
    APLCON::PrintFormatting::Width = 11;

    // eta' fitter
    etap_fit.LinkVariable("Lepton1", etap_fs[0].Link(), etap_fs[0].LinkSigma());
    etap_fit.LinkVariable("Lepton2", etap_fs[1].Link(), etap_fs[1].LinkSigma());
    etap_fit.LinkVariable("Photon", etap_fs[2].Link(), etap_fs[2].LinkSigma());
    vector<string> etap_names = {"Lepton1", "Lepton2", "Photon"};
    // Constraint: Invariant mass of nPhotons equals constant IM,
    // make lambda catch also this with [&] specification
    auto etapIM = [&] (const vector<vector<double>>& particles) -> double
    {
        TLorentzVector sum(0,0,0,0);
        sum += FitParticle::Make(particles[0], ParticleTypeDatabase::eMinus.Mass());
        sum += FitParticle::Make(particles[1], ParticleTypeDatabase::eMinus.Mass());
        sum += FitParticle::Make(particles[2], ParticleTypeDatabase::Photon.Mass());

        return sum.M() - IM;
    };
    etap_fit.AddConstraint("IM", etap_names, etapIM);
    // iteration steps
    APLCON::Fit_Settings_t etap_settings = etap_fit.GetSettings();
    etap_settings.MaxIterations = 10;
    etap_fit.SetSettings(etap_settings);
}

void ant::analysis::SaschaPhysics::ProcessEvent(const ant::Event &event)
{
    const TLorentzVector target(0., 0., 0., ParticleTypeDatabase::Proton.Mass());
    TLorentzVector etap;
    double q2;

    const bool MC = !(event.MCTrue().Particles().GetAll().empty());
    if (MC)
        true_particles.clear();

    // do a q2 preselection
    double q2_true = -1.;
    if (MC) {
        Particle ePlus_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
        Particle eMinus_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
        for (const auto& p : event.MCTrue().Particles().GetAll()) {
            if (p->Type() == ParticleTypeDatabase::eMinus)
                eMinus_true = Particle(*p);
            else if (p->Type() == ParticleTypeDatabase::ePlus)
                ePlus_true = Particle(*p);
        }
        //double lepton_open_angle_true = ePlus_true.Angle(eMinus_true.Vect())*TMath::RadToDeg();
        q2_true = (ePlus_true + eMinus_true).M();
        if (q2_true < 50.)
            return;
    }

    accepted_events->Fill("all events", 1);

    cb_esum->Fill(event.Reconstructed().TriggerInfos().CBEenergySum());

    if (MC && event.Reconstructed().TriggerInfos().CBEenergySum() < CB_ESUM)
        return;

    if (MC)
        accepted_events->Fill("MC CB Esum", 1);

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

    // get the tracks in CB and TAPS, time cuts on clusters are applied
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

    // fill true MC information
    if (MC)
        fill_MC_true(event.MCTrue().Particles().GetAll(), tracksCB.size(), tracksTAPS.size());

    if (tracksCB.size() + tracksTAPS.size() != nFinalState)
        return;

    accepted_events->Fill("#part FS", 1);

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

    sort_particles(particles);
    etap = TLorentzVector(0., 0., 0., 0.);
    for (auto it = particles.cbegin(); it != particles.cend()-1; ++it)
        etap += *it;
    q2 = (particles.at(0) + particles.at(1)).M();
    q2_vs_im_candidates->Fill(etap.M(), q2);

    /* Achim's proposed test */
/*    if (tracksCB.size() != 3)
        return;
    accepted_events->Fill("#tracks CB", 1);
    size_t nCharged = 0;
    for (const auto& track : tracksCB) {
        if (track->Detector() & detector_t::PID) {
            particles.push_back(Particle(ParticleTypeDatabase::eMinus, track));
            nCharged++;
        } else
            particles.push_back(Particle(ParticleTypeDatabase::Photon, track));
    }
    if (nCharged != 2)
        return;
    accepted_events->Fill("2 leptons", 1);
    // swap the photon at the third (last) position
    for (auto it = particles.begin(); it != particles.end()-1; ++it)
        if (it->Type() == ParticleTypeDatabase::Photon) {
            std::iter_swap(it, particles.end()-1);
            break;
        }
    TLorentzVector tmp_etap(0., 0., 0., 0.);
    for (auto it = particles.cbegin(); it != particles.cend(); ++it)
        tmp_etap += *it;

    auto etap_it = etap_fs.begin();
    for (const auto& p : particles)
        (etap_it++)->SetFromVector(p);

    const APLCON::Result_t& etap_res = etap_fit.DoFit();
    if (etap_res.Status != APLCON::Result_Status_t::Success)
        return;
    accepted_events->Fill("fit #eta'", 1);
    // Cut on chi^2
    etap_chi2->Fill(etap_res.ChiSquare);
    etap_chi2_vs_q2->Fill(q2_true, etap_res.ChiSquare);
    if (etap_res.ChiSquare > 15.)
        return;
    accepted_events->Fill("#chi^{2} #eta'", 1);
    TLorentzVector fit_etap(0., 0., 0., 0.);
    fit_etap += FitParticle::Make(etap_fs[0], ParticleTypeDatabase::eMinus.Mass());
    fit_etap += FitParticle::Make(etap_fs[1], ParticleTypeDatabase::eMinus.Mass());
    fit_etap += FitParticle::Make(etap_fs[2], ParticleTypeDatabase::Photon.Mass());
*/    /* finished proton test, next test tagger hit combinations */

    /* test for pion background */
    //interval<double> pion_cut(113., 158.5);  // 2sigma
    interval<double> pion_cut(102., 170.);  // 3sigma
    TLorentzVector pi0;
    std::vector<FitParticle> pi0_cand(3);
    auto part_it = pi0_cand.begin();
    for (const auto& p : particles)
        if (p.Type() == ParticleTypeDatabase::eCharged || p.Type() == ParticleTypeDatabase::Photon)
            (part_it++)->SetFromVector(p);
    //const std::vector<std::array<size_t, 2>> pi0_comb = {{0, 1}, {0, 2}, {1, 2}};
    const std::vector<std::array<size_t, 2>> pi0_comb = {{0, 2}, {1, 2}};  // third position is the photon
    for (const auto comb : pi0_comb) {
        pi0 = TLorentzVector(0., 0., 0., 0.);
        for (const auto idx : comb)
            pi0 += FitParticle::Make(pi0_cand.at(idx), ParticleTypeDatabase::Photon.Mass());
        sub_pi0_cand_im->Fill(pi0.M());
    }
    for (const auto comb : pi0_comb) {
        pi0 = TLorentzVector(0., 0., 0., 0.);
        for (const auto idx : comb)
            pi0 += FitParticle::Make(pi0_cand.at(idx), ParticleTypeDatabase::Photon.Mass());
        // apply an anti pion cut
        if (pion_cut.Contains(pi0.M()))
            return;
    }
    for (const auto comb : pi0_comb) {
        pi0 = TLorentzVector(0., 0., 0., 0.);
        for (const auto idx : comb)
            pi0 += FitParticle::Make(pi0_cand.at(idx), ParticleTypeDatabase::Photon.Mass());
        sub_pi0_cand_im_cut->Fill(pi0.M());
    }
    q2_vs_im_pi0_cut->Fill(etap.M(), q2);
    accepted_events->Fill("anti #pi^{0} cut", 1);

    /* process tagger hits */

    TaggerHistList tagger_hits;
    size_t prompt_n_tagger = 0, random_n_tagger = 0;
    //static size_t count = 0;
    constexpr bool trueMC = true;  // use MC true
    if (trueMC && MC)
        tagger_hits = event.MCTrue().TaggerHits();
    else
        tagger_hits = event.Reconstructed().TaggerHits();

    size_t proton_count = 0;
/*    vector<size_t> taps_cluster_as_proton_prompt;
    vector<size_t> taps_cluster_as_proton_random;
    taps_cluster_as_proton_prompt.resize(tracksTAPS.size());
    taps_cluster_as_proton_random.resize(tracksTAPS.size());
    constexpr double proton_cone = 4.;*/
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

        // plot proton tests (TOF, PSA) in case of prompt hit
        if (is_prompt)
            proton_tests(event.Reconstructed().Tracks(), taggerhit);

        /* Achim's proposed test, continued */
        //double best_chi2 = 1000.;
/*        Particle proton_candidate(ParticleTypeDatabase::Proton, 0., 0., 0.);
        TLorentzVector expected_proton = taggerhit->PhotonBeam() + target - fit_etap;
        if (MC) {
            Particle true_proton = get_true_proton(event.MCTrue().Particles().GetAll());
            expected_proton_diff_vs_q2->Fill(q2_true, expected_proton.Angle(true_proton.Vect())*TMath::RadToDeg());
            expected_proton_diff_vs_q2_rebin->Fill(q2_true, expected_proton.Angle(true_proton.Vect())*TMath::RadToDeg());
        }
        bool proton_found = false;
        size_t i = 0;
        for (const auto& track : tracksTAPS) {
            proton_candidate = Particle(ParticleTypeDatabase::Proton, track);
            if (expected_proton.Angle(proton_candidate.Vect())*TMath::RadToDeg() < proton_cone) {
                if (proton_found)
                    accepted_events->Fill("2nd proton", 1);
                proton_found = true;
                //break;
                proton_count++;
                if (is_prompt)
                    taps_cluster_as_proton_prompt[i]++;
                else
                    taps_cluster_as_proton_random[i]++;
            }
            i++;
        }
        if (!proton_found)
            continue;
        accepted_events->Fill("proton found", 1);
        particles.push_back(proton_candidate);
*/        /* finished proton test */

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
        //TODO: check if the following is needed anymore
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

        // sort order of particles
        sort_particles(particles);

        TLorentzVector proton = particles.back();
        etap = TLorentzVector(0., 0., 0., 0.);
        for (auto it = particles.cbegin(); it != particles.cend()-1; ++it)
            etap += *it;
        q2 = (particles.at(0) + particles.at(1)).M();
        h["q2_im_base_selection"]->Fill(etap.M(), q2);
        TLorentzVector* photon_cm = static_cast<TLorentzVector*>(particles.at(2).Clone());
        photon_cm->Boost(-etap.BoostVector());
        h["energy_photon_cm"]->Fill(photon_cm->E());

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
        auto fs_it = final_state.begin();
        auto det_it = detectors.begin();
        det_it++;  // skip first element which is the beam photon --> not considered when recalculating clusters due to v_z
        for (const auto& p : particles) {
            //(it++)->SetFromVector(p);
            set_fit_particle(p, *fs_it++);
            *det_it++ = p.Tracks().front()->Detector();
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
    protons_found->Fill(proton_count);
/*    const size_t protons_found_prompt = sum_vector(taps_cluster_as_proton_prompt);
    const size_t protons_found_random = sum_vector(taps_cluster_as_proton_random);
    prompt["protons_found"]->Fill(protons_found_prompt);
    random["protons_found"]->Fill(protons_found_random);
    prompt["nTAPS_vs_protons_found"]->Fill(protons_found_prompt, tracksTAPS.size());
    random["nTAPS_vs_protons_found"]->Fill(protons_found_random, tracksTAPS.size());
    const size_t diff_clusters_prompt = non_zero_entries(taps_cluster_as_proton_prompt);
    const size_t diff_clusters_random = non_zero_entries(taps_cluster_as_proton_random);
    prompt["protons_diff_clusters"]->Fill(diff_clusters_prompt);
    random["protons_diff_clusters"]->Fill(diff_clusters_random);
    prompt["proton_diff_clusters_vs_q2"]->Fill(q2_true, diff_clusters_prompt);
    random["proton_diff_clusters_vs_q2"]->Fill(q2_true, diff_clusters_random);*/
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
    c << drawoption("colz") << pid
      << padoption::set(padoption_t::Legend)
      << particle_types
      << padoption::unset(padoption_t::Legend)
      << tagger_spectrum << diff["n_tagged"] << cb_esum << endc;

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

void ant::analysis::SaschaPhysics::sort_particles(particle_vector& particles)
{
//    cout << "We have the right 4 particles\nNow sort them, order before:\n";
//    for (const auto& p : particles)
//        cout << p.Type() << endl;

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
}

const vector<Particle> ant::analysis::SaschaPhysics::get_MC_true_particles(const ParticleList& particles)
{
    if (!true_particles.empty())
        return true_particles;

    Particle ePlus_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
    Particle eMinus_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
    Particle photon_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
    Particle proton_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
    for (const auto& p : particles) {
        if (p->Type() == ParticleTypeDatabase::Photon)
            photon_true = Particle(*p);
        else if (p->Type() == ParticleTypeDatabase::eMinus)
            eMinus_true = Particle(*p);
        else if (p->Type() == ParticleTypeDatabase::ePlus)
            ePlus_true = Particle(*p);
        else if (p->Type() == ParticleTypeDatabase::Proton)
            proton_true = Particle(*p);
    }
    if (ePlus_true.Type() == ParticleTypeDatabase::Neutron
            || eMinus_true.Type() == ParticleTypeDatabase::Neutron
            || photon_true.Type() == ParticleTypeDatabase::Neutron
            || proton_true.Type() == ParticleTypeDatabase::Neutron)
        accepted_events->Fill("wrong MC", 1);

    true_particles = {ePlus_true, eMinus_true, photon_true, proton_true};
    return true_particles;
}

const Particle ant::analysis::SaschaPhysics::get_true_particle(const ParticleList& particles, const size_t pos)
{
    if (true_particles.empty())
        return get_MC_true_particles(particles).at(pos);
    else
        return true_particles.at(pos);
}

const Particle ant::analysis::SaschaPhysics::get_true_positron(const ParticleList& particles)
{
    return get_true_particle(particles, 0);
}

const Particle ant::analysis::SaschaPhysics::get_true_electron(const ParticleList& particles)
{
    return get_true_particle(particles, 1);
}

const Particle ant::analysis::SaschaPhysics::get_true_photon(const ParticleList& particles)
{
    return get_true_particle(particles, 2);
}

const Particle ant::analysis::SaschaPhysics::get_true_proton(const ParticleList& particles)
{
    return get_true_particle(particles, 3);
}

void ant::analysis::SaschaPhysics::fill_MC_true(const ParticleList& particles, const size_t tracksCB, const size_t tracksTAPS)
{
    const particle_vector true_part = get_MC_true_particles(particles);
    Particle ePlus_true = true_part.at(0);
    Particle eMinus_true = true_part.at(1);
    Particle photon_true = true_part.at(2);
    Particle proton_true = true_part.at(3);
    TLorentzVector etap_true = ePlus_true + eMinus_true + photon_true;

    double lepton_open_angle_true = ePlus_true.Angle(eMinus_true.Vect())*TMath::RadToDeg();
    double q2 = (ePlus_true + eMinus_true).M();
    eplus_true_theta_vs_energy->Fill(ePlus_true.Ek(), ePlus_true.Theta()*TMath::RadToDeg());
    eminus_true_theta_vs_energy->Fill(eMinus_true.Ek(), eMinus_true.Theta()*TMath::RadToDeg());
    photon_true_theta_vs_energy->Fill(photon_true.Ek(), photon_true.Theta()*TMath::RadToDeg());
    proton_true_theta_vs_energy->Fill(proton_true.Ek(), proton_true.Theta()*TMath::RadToDeg());
    true_open_angle_vs_q2->Fill(q2, lepton_open_angle_true);
    n_cluster_cb_vs_q2->Fill(q2, tracksCB);
    n_cluster_taps_vs_q2->Fill(q2, tracksTAPS);
    n_cluster_cb_vs_open_angle->Fill(lepton_open_angle_true, tracksCB);
    n_cluster_taps_vs_open_angle->Fill(lepton_open_angle_true, tracksTAPS);
    opening_angle_leptons_true->Fill(lepton_open_angle_true);
    opening_angle_vs_photon_energy_true->Fill(photon_true.Ek(), lepton_open_angle_true);

    double en_lep1_true, en_lep2_true;
    if (eMinus_true.Ek() > ePlus_true.Ek()) {
        en_lep1_true = eMinus_true.Ek();
        en_lep2_true = ePlus_true.Ek();
    } else {
        en_lep1_true = ePlus_true.Ek();
        en_lep2_true = eMinus_true.Ek();
    }
    lepton_energies_true->Fill(en_lep1_true, en_lep2_true);
    energy_lepton1_true->Fill(en_lep1_true);
    energy_lepton2_true->Fill(en_lep2_true);

    TLorentzVector* photon_cm = static_cast<TLorentzVector*>(photon_true.Clone());
    photon_cm->Boost(-etap_true.BoostVector());
    energy_photon_cm_true->Fill(photon_cm->E());
}

void ant::analysis::SaschaPhysics::proton_tests(const TrackList& tracks, const TaggerHitPtr taggerhit)
{
    for (const auto& track : tracks) {
        double tof = taggerhit->Time() - track->Time();
        double shortE = track->ShortEnergy(), clusterE = track->ClusterEnergy();
        const double r2d = TMath::RadToDeg();

        if (track->Detector() & ant::detector_t::anyTAPS) {
            //dynamic_cast<TH3D*>(prompt["proton_id"])->Fill(track->ClusterSize(), clusterE, tof);
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
            //dynamic_cast<TH3D*>(random["proton_id"])->Fill(track->ClusterSize(), clusterE, tof);
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

double ant::analysis::SaschaPhysics::calculate_energy_weighted_cb_time_average(const TrackList& tracks) const
{
    double time = 0, e_sum = 0;
    for (const auto& track : tracks) {
        time += track->Time();
        e_sum += track->ClusterEnergy();
    }
    return time/e_sum;
}

void ant::analysis::SaschaPhysics::GetTracks(const Event& event, TrackList& tracksCB, TrackList& tracksTAPS)
{
    for (const auto& track : event.Reconstructed().Tracks()) {
        // ignore particles below set cluster energy threshold
        if (track->ClusterEnergy() < CLUSTER_TRESH)
            continue;
        if (track->Detector() & detector_t::anyCB && cb_time_window.Contains(track->Time()))
            tracksCB.push_back(track);
        else if (track->Detector() & detector_t::anyTAPS && taps_time_window.Contains(track->Time()))
            tracksTAPS.push_back(track);
    }
}

bool ant::analysis::SaschaPhysics::IdentifyTracks(const TrackList& tracksCB, const TrackList& tracksTAPS,
                                                  particle_vector& particles)
{
    bool success = false;
    size_t nCB = tracksCB.size(), nCharged = 0;
    if (nCB + tracksTAPS.size() != nFinalState)
        return success;
    if (nCB == 3) {
        for (const auto& track : tracksCB) {
            if (track->Detector() & detector_t::PID)
                particles.push_back(Particle(ParticleTypeDatabase::eMinus, track));
            else
                particles.push_back(Particle(ParticleTypeDatabase::Photon, track));
        }
        particles.push_back(Particle(ParticleTypeDatabase::Proton, tracksTAPS.front()));
        success = true;
    } else if (nCB == 2) {
        for (const auto& track : tracksCB) {
            if (track->Detector() & detector_t::PID) {
                particles.push_back(Particle(ParticleTypeDatabase::eMinus, track));
                nCharged++;
            } else
                particles.push_back(Particle(ParticleTypeDatabase::Photon, track));
        }
        if (nCharged == 0)
            return success;
        //TODO: Correct for path lengths!
        double time_track0 = abs(tracksTAPS.at(0)->Time()), time_track1 = abs(tracksTAPS.at(1)->Time());
        unsigned short proton_id = 0;
        if (time_track1 > time_track0)
            proton_id = 1;
        if (nCharged < 2)
            particles.push_back(Particle(ParticleTypeDatabase::eMinus, tracksTAPS.at(1-proton_id)));
        else
            particles.push_back(Particle(ParticleTypeDatabase::Photon, tracksTAPS.at(1-proton_id)));
        particles.push_back(Particle(ParticleTypeDatabase::Proton, tracksTAPS.at(proton_id)));
        success = true;
    }

    return success;
}

void ant::analysis::SaschaPhysics::GetParticles(const ant::Event& event, particle_vector& particles)
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

void ant::analysis::SaschaPhysics::GetTrueParticles(const ant::Event& event, particle_vector& particles)
{
    for (const auto& p : event.MCTrue().Particles().GetAll())
        if (contains(ParticleTypeDatabase::MCFinalStateTypes(), p->Type()))
            particles.push_back(Particle(*p));
    nParticles = particles.size();
}

TH1* ant::analysis::SaschaPhysics::read_hist(const char* file, const char* hist_name) const
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

void ant::analysis::SaschaPhysics::read_file(const char* file, std::vector<double>& values, const int fill_size)
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

void ant::analysis::SaschaPhysics::set_fit_particle(const Particle& p, FitParticle& part)
{
    vector<double> sigmas;
    get_uncertainties(p, sigmas);
    part.SetFromParticle(p, sigmas);
}

void ant::analysis::SaschaPhysics::get_uncertainties(const Particle& p, vector<double>& sigmas)
{
    sigmas.clear();
    sigmas.push_back(get_sigma_energy(p));
    sigmas.push_back(get_sigma_theta(p));
    sigmas.push_back(get_sigma_phi(p));
}

double ant::analysis::SaschaPhysics::get_sigma_energy(const Particle& p) const
{
    if (p.Type() == ParticleTypeDatabase::Proton)
        return 0;//get_sigma(p, proton_energy_uncertainties);
    else
        return get_sigma(p, photon_energy_uncertainties);
}

double ant::analysis::SaschaPhysics::get_sigma_theta(const Particle& p) const
{
    if (p.Type() == ParticleTypeDatabase::Proton)
        return get_sigma(p, proton_theta_uncertainties);
    else
        return get_sigma(p, photon_theta_uncertainties);
}

double ant::analysis::SaschaPhysics::get_sigma_phi(const Particle& p) const
{
    if (p.Type() == ParticleTypeDatabase::Proton)
        return get_sigma(p, proton_phi_uncertainties);
    else
        return get_sigma(p, photon_phi_uncertainties);
}

double ant::analysis::SaschaPhysics::get_sigma(const Particle& p, TH2* const h) const
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
size_t ant::analysis::SaschaPhysics::get_histogram_neigbours(TH1* const h,
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
double ant::analysis::SaschaPhysics::local_bin_average(TH1* const h, const int bin) const
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

void ant::analysis::SaschaPhysics::apply_time_correction(const TrackList& tracksCB, const TrackList& tracksTAPS)
{
    for (const auto track : tracksCB)
        track->SetTime(track->Time() - cb_time_correction.at(track->CentralCrystal()));
}

void ant::analysis::SaschaPhysics::apply_energy_correction(const TrackList& tracksCB, const TrackList& tracksTAPS)
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

void ant::analysis::SaschaPhysics::apply_theta_correction(const TrackList& tracksCB, const TrackList& tracksTAPS)
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
