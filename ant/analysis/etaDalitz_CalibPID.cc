#include "etaDalitz_CalibPID.h"
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


TLorentzVector analysis::etaDalitz_CalibPID::FitParticle::Make(const std::vector<double>& EkThetaPhi, const Double_t m)
{
    const double E = EkThetaPhi[0] + m;
    const Double_t p = sqrt(E*E - m*m);
    TVector3 pv(1,0,0);
    pv.SetMagThetaPhi(p, EkThetaPhi[1], EkThetaPhi[2]);
    TLorentzVector l(pv, E);
    return l;
}

void analysis::etaDalitz_CalibPID::FitParticle::Smear()
{
}

void analysis::etaDalitz_CalibPID::FillIM(TH1 *h, const std::vector<analysis::etaDalitz_CalibPID::FitParticle>& final_state)
{
    TLorentzVector sum(0,0,0,0);
    for (auto p = final_state.cbegin(); p != final_state.cend()-1; ++p)
        sum += FitParticle::Make(*p, ParticleTypeDatabase::Photon.Mass());
    h->Fill(sum.M());
}

void analysis::etaDalitz_CalibPID::HistList::AddHistogram(const string &name,
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

void analysis::etaDalitz_CalibPID::HistList::AddHistogram(const string &name,
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

void analysis::etaDalitz_CalibPID::HistList::AddHistogram(const string &name,
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

void analysis::etaDalitz_CalibPID::HistList::AddHistogram(const string &name,
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

void analysis::etaDalitz_CalibPID::HistList::AddHistogram(const string &name,
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

void analysis::etaDalitz_CalibPID::HistList::AddHistogram(const string &name,
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

void analysis::etaDalitz_CalibPID::HistList::Draw()
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
        c = new canvas("etaDalitz_CalibPID: Overview " + pref + std::to_string(++i));
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

analysis::etaDalitz_CalibPID::HistList &analysis::etaDalitz_CalibPID::HistList::operator*=(const Double_t factor)
{
    for (auto i : h)
        i.second->Scale(factor);

    return *this;
}

analysis::etaDalitz_CalibPID::HistList analysis::etaDalitz_CalibPID::HistList::operator=(const analysis::etaDalitz_CalibPID::HistList &other)
{
    for (auto i : h) {
        TH1* h = i.second;
        h->Reset();
        h->Add(other[i.first]);
    }

    return *this;
}

void analysis::etaDalitz_CalibPID::HistList::AddScaled(const analysis::etaDalitz_CalibPID::HistList &h2, const Double_t f)
{
    for (auto i : h)
        i.second->Add(h2.h.at(i.first), f);
}

analysis::etaDalitz_CalibPID::HistList::HistList(const string &prefix, const mev_t energy_scale)
{
    pref = prefix + '_';

    const BinSettings energy_bins(1600, 0, energy_scale*1.6);
    const BinSettings tagger_bins(200, 1400, 1600);
    const BinSettings veto_bins(1000, 0, 10);
    const BinSettings chisquare_bins(100, 0, 30);
    const BinSettings probability_bins(100, 0, 1);
    const BinSettings iterations_bins(15, 0, 15);
    const BinSettings im_bins(200, IM-100, IM+100);
    const BinSettings phi_bins(720, 0, 360);
    const BinSettings energy_range_bins(1600, -1600, 1600);
    const BinSettings q2_bins(2000, 0, 1000);
    const BinSettings count_bins(50, 0, 50);

    AddHistogram("tagger_energy", "Tagger Energy", "Photon Beam Energy [MeV]", "#", tagger_bins);
    AddHistogram("n_tagged", "Tagger Hits", "Tagger Hits / event", "#", count_bins);

    AddHistogram("coplanarity", "Coplanarity #eta' proton", "Coplanarity [#circ]", "#", phi_bins);
    AddHistogram("missing_mass", "Missing Mass Proton", "m_{miss} [MeV]", "#", energy_bins);

    stringstream fs;
    fs << "eeg";
    AddHistogram("im_true", "IM "+fs.str()+" true", "IM", "#", im_bins);
    AddHistogram("im_smeared", "IM "+fs.str()+" smeared", "IM", "#", im_bins);
    AddHistogram("im_fit", "IM "+fs.str()+" fit", "IM", "#", im_bins);

    AddHistogram("coplanarity_fit", "Coplanarity #eta' proton fitted", "Coplanarity [#circ]", "#", phi_bins);
    AddHistogram("missing_mass_fit", "Missing Mass Proton fitted", "m_{miss} [MeV]", "#", energy_bins);
}

ant::analysis::etaDalitz_CalibPID::etaDalitz_CalibPID(const mev_t energy_scale) :
    Physics("etaDalitz_CalibPID"),
    prompt("prompt", energy_scale),
    random("random", energy_scale),
    diff("diff", energy_scale),
    prompt_window(-6, 6),
    random_window1(-40, -15),
    random_window2(15, 40),
    fitter("etaDalitz_CalibPID"),
    final_state(nFinalState)
{
    cout << "Eta Dalitz CalibPID Physics:\n";
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

    // histogram to count how many events passed the different selection criteria
    steps = HistFac.makeTH1D("Accepted events", "step", "#", BinSettings(10), "steps");

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
    n_part = HistFac.makeTH1D("Number of particles", "#particles", "#", particle_bins, "n_part");
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

    const BinSettings pid_channels(24);
    const BinSettings veto_energy(1000, 0, 10);
    const BinSettings energy(1200);
    eegPID = HistFac.makeTH2D("PID 2 charged 1 neutral", "PID Energy [MeV]", "#", veto_energy, pid_channels, "eegPID");
    etaIM = HistFac.makeTH1D("IM #eta all comb", "IM [MeV]", "#", energy, "etaIM");
    etaIM_fit = HistFac.makeTH1D("IM #eta fitted", "IM [MeV]", "#", energy, "etaIM_fit");
    etaIM_cand = HistFac.makeTH1D("IM #eta candidates", "IM [MeV]", "#", energy, "etaIM_cand");
    etaIM_final = HistFac.makeTH1D("IM #eta final", "IM [MeV]", "#", energy, "etaIM_final");
    MM = HistFac.makeTH1D("Missing Mass proton", "MM [MeV]", "#", BinSettings(1600), "MM");
    hCopl = HistFac.makeTH1D("Coplanarity #eta - proton all comb", "coplanarity [#circ]", "#", BinSettings(720, -180, 180), "hCopl");
    hCopl_final = HistFac.makeTH1D("Coplanarity #eta - proton final", "coplanarity [#circ]", "#", BinSettings(720, -180, 180), "hCopl_final");
    hChi2 = HistFac.makeTH1D("#chi^{2}", "#chi^{2}", "#", BinSettings(500, 0, 100), "hChi2");
    hProb = HistFac.makeTH1D("Probability", "probability", "#", BinSettings(500, 0, 1), "hProb");
    hIter = HistFac.makeTH1D("# Iterations", "#iterations", "#", BinSettings(20), "hIter");
    pTheta = HistFac.makeTH1D("#vartheta proton candidate", "#vartheta_{p} [#circ]", "#", BinSettings(720, 0, 180), "pTheta");
    protonVeto = HistFac.makeTH1D("Veto energy identified proton", "Veto [MeV]", "#", veto_energy, "protonVeto");

    // setup fitter for eta Dalitz decay
    fitter.LinkVariable("Beam", beam.Link(), beam.LinkSigma());
    fitter.LinkVariable("Photon1", final_state[0].Link(), final_state[0].LinkSigma());
    fitter.LinkVariable("Photon2", final_state[1].Link(), final_state[1].LinkSigma());
    fitter.LinkVariable("Photon3", final_state[2].Link(), final_state[2].LinkSigma());
    fitter.LinkVariable("Proton", final_state[3].Link(), final_state[3].LinkSigma());

    vector<string> all_names = {"Beam", "Photon1", "Photon2", "Photon3", "Proton"};
    if (includeVertexFit) {
        all_names.push_back("v_z");
        //fitter.AddUnmeasuredVariable("v_z"); // default value 0
        auto v_z_settings = APLCON::Variable_Settings_t::Default;
        v_z_settings.Limit.High = TARGET_MAX;
        v_z_settings.Limit.Low = TARGET_MIN;
        fitter.AddMeasuredVariable("v_z", 0., 2.3, v_z_settings);  // default value 0
    }

    // Constraint: Incoming 4-vector = Outgoing 4-vector
    auto EnergyMomentumBalance = [] (vector<vector<double>>& particles) -> vector<double>
    {
        const TLorentzVector target(0,0,0, ParticleTypeDatabase::Proton.Mass());
        // assume first particle is beam photon
        TLorentzVector diff = target + FitParticle::Make(particles[0], ParticleTypeDatabase::Photon.Mass());
        // assume second to fourth particle outgoing photons (eeg)
        diff -= FitParticle::Make(particles[1], ParticleTypeDatabase::Photon.Mass());
        diff -= FitParticle::Make(particles[2], ParticleTypeDatabase::Photon.Mass());
        diff -= FitParticle::Make(particles[3], ParticleTypeDatabase::Photon.Mass());
        // assume last particle outgoing proton
        diff -= FitParticle::Make(particles[4], ParticleTypeDatabase::Proton.Mass());

        return {diff.X(), diff.Y(), diff.Z(), diff.T()};
    };
    fitter.AddConstraint("EnergyMomentumBalance", all_names, EnergyMomentumBalance);

    // Constraint: Coplanarity between eta' and recoil proton
    auto CoplanarityConstraint = [] (const vector<vector<double>>& particles) -> double
    {
        TLorentzVector etap(0,0,0,0);
        TLorentzVector proton;
        etap = FitParticle::Make(particles[1], ParticleTypeDatabase::Photon.Mass());
        etap += FitParticle::Make(particles[2], ParticleTypeDatabase::Photon.Mass());
        etap += FitParticle::Make(particles[3], ParticleTypeDatabase::Photon.Mass());
        proton = FitParticle::Make(particles[4], ParticleTypeDatabase::Proton.Mass());

        return abs(etap.Phi() - proton.Phi())*TMath::RadToDeg() - 180.;
    };
    if (includeCoplanarityConstraint)
        fitter.AddConstraint("CoplanarityConstraint", all_names, CoplanarityConstraint);

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

void ant::analysis::etaDalitz_CalibPID::ProcessEvent(const ant::Event &event)
{
    const TLorentzVector target(0., 0., 0., ParticleTypeDatabase::Proton.Mass());

    const bool MC = !(event.MCTrue().Particles().GetAll().empty());
    if (MC)
        true_particles.clear();

    steps->Fill("all events", 1);

    cb_esum->Fill(event.Reconstructed().TriggerInfos().CBEenergySum());
/*
    if (MC && event.Reconstructed().TriggerInfos().CBEenergySum() < CB_ESUM)
        return;

    if (MC)
        steps->Fill("CB Esum", 1);
*/
    for (const auto& track : event.Reconstructed().Tracks()) {
        pid->Fill(track->ClusterEnergy(), track->VetoEnergy());
        // fill time information
        cluster_time->Fill(track->Time());
        if (track->Detector() & detector_t::anyCB)
            cb_time->Fill(track->Time());
        else
            taps_time->Fill(track->Time());
    }

    const double cb_time_avg_energy_weighted = calculate_energy_weighted_cb_time_average(event.Reconstructed().Tracks());
    cb_time_avg->Fill(cb_time_avg_energy_weighted);

    n_part->Fill(event.Reconstructed().Tracks().size());



    TrackList tracks = event.Reconstructed().Tracks();

    if (tracks.size() != 4)
        return;
    steps->Fill("#tracks", 1);

    TLorentzVector eta;
    TLorentzVector proton;
    //static constexpr double ETA_IM = 547.83;
    //static constexpr double ETA_SIGMA = 50.;
    //const interval<double> eta_im({ETA_IM-ETA_SIGMA, ETA_IM+ETA_SIGMA});
    const interval<double> coplanarity({-25, 25});
    particle_vector particles;
    TrackList comb;
    for (auto t : tracks)
        comb.emplace_back(t);

    // require at least 2 candidates with PID/Veto entries
    if (std::count_if(comb.begin(), comb.end(), [](TrackPtr t){ return t->VetoEnergy(); }) < 2)
        return;
    steps->Fill("#Veto", 1);

    TLorentzVector missing;
    const interval<double> mm({ParticleTypeDatabase::Proton.Mass()-150., ParticleTypeDatabase::Proton.Mass()+150.});
    double min_chi2 = std::numeric_limits<double>::infinity();
    size_t best_comb = tracks.size();

    for (const auto& taggerhit : event.Reconstructed().TaggerHits()) {
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

        steps->Fill("Tagger time", 1);

        // make sure the correct histogram will be filled
        HistList& h = is_prompt ? prompt : random;

        h["tagger_energy"]->Fill(taggerhit->PhotonEnergy());

        for (size_t i = 0; i < tracks.size(); i++) {  // loop to test all different combinations
            // ensure the possible proton candidate is kinematically allowed
            if (comb.back()->Theta() * TMath::RadToDeg() > 60.) {
                shift_right(comb);
                continue;
            }
            steps->Fill("proton #vartheta", 1);

            // require 2 PID entries for the eta candidate
            if (std::count_if(comb.begin(), comb.end()-1, [](TrackPtr t){ return t->VetoEnergy(); }) < 2) {
                shift_right(comb);
                continue;
            }
            steps->Fill("2 PIDs", 1);

            particles.clear();
            eta.SetXYZT(0,0,0,0);
            for (size_t j = 0; j < comb.size()-1; j++) {
                particles.emplace_back(Particle(ParticleTypeDatabase::Photon, comb.at(j)));
                eta += particles.back();
            }
            particles.emplace_back(Particle(ParticleTypeDatabase::Proton, comb.back()));
            proton = particles.back();
            etaIM->Fill(eta.M());

            const double copl = abs(eta.Phi() - proton.Phi())*TMath::RadToDeg() - 180.;
            hCopl->Fill(copl);
            if (!coplanarity.Contains(copl)) {
                shift_right(comb);
                continue;
            }
            steps->Fill("coplanarity", 1);

            missing = taggerhit->PhotonBeam() + target - eta;
            MM->Fill(missing.M());
            if (!mm.Contains(missing.M())) {
                shift_right(comb);
                continue;
            }
            steps->Fill("missing mass", 1);

            // at this point, the event and combination looks good
            // check this combination via kinematic fit

            auto fs_it = final_state.begin();
            for (const auto& p : particles)
                set_fit_particle(p, *fs_it++);

            // set proton energy sigma to zero to indicate it's unmeasured
            final_state[3].Ek_Sigma = 0;

            beam.SetFromVector(taggerhit->PhotonBeam());
            // set beam sigmas, energy 2/sqrt(3)
            beam.Ek_Sigma = 1.1547;
            beam.Theta_Sigma = .0001;
            beam.Phi_Sigma = .0001;

            // let APLCON do the work
            const APLCON::Result_t& result = fitter.DoFit();

            //cout << result << endl;

            if (result.Status != APLCON::Result_Status_t::Success) {
                shift_right(comb);
                continue;
            }
            steps->Fill("KinFit", 1);

            const double chi2 = result.ChiSquare;
            const double prob = result.Probability;
            const int iterations = result.NIterations;

            for (const auto& it_map : result.Variables) {
                const string& varname = it_map.first;
                const APLCON::Result_Variable_t& var = it_map.second;
                if (is_prompt)
                    pulls_prompt.at(varname)->Fill(var.Pull);
                else
                    pulls_random.at(varname)->Fill(var.Pull);
            }
            hChi2->Fill(chi2);
            hProb->Fill(prob);
            hIter->Fill(iterations);

            if (prob < .1) {
                shift_right(comb);
                continue;
            }
            steps->Fill("probability", 1);

            if (chi2 < min_chi2) {
                min_chi2 = chi2;
                best_comb = i;
            }

            // create the fitted final state particles
            //TLorentzVector proton_fit = FitParticle::Make(final_state[3], ParticleTypeDatabase::Proton.Mass());
            TLorentzVector eta_fit(0., 0., 0., 0.);
            eta_fit += FitParticle::Make(final_state[0], ParticleTypeDatabase::Photon.Mass());
            eta_fit += FitParticle::Make(final_state[1], ParticleTypeDatabase::Photon.Mass());
            eta_fit += FitParticle::Make(final_state[2], ParticleTypeDatabase::Photon.Mass());
            etaIM_fit->Fill(eta_fit.M());



            shift_right(comb);
        }

    }

    if (best_comb >= tracks.size() || !isfinite(min_chi2))
        return;
    steps->Fill("best comb", 1);

    // restore combinations with best chi2
    //for (size_t i = 0; i < best_comb; i++)
    while (best_comb-- > 0)
        shift_right(comb);

    proton = Particle(ParticleTypeDatabase::Proton, comb.back());
    eta.SetXYZT(0,0,0,0);
    for (size_t i = 0; i < comb.size()-1; i++)
        eta += Particle(ParticleTypeDatabase::Photon, comb.at(i));
    etaIM_cand->Fill(eta.M());
    protonVeto->Fill(comb.back()->VetoEnergy());
    pTheta->Fill(comb.back()->Theta()*TMath::RadToDeg());
    // at this point a possible eta Dalitz candidate was found, work only with eta final state
    comb.pop_back();

    sort(comb.begin(), comb.end(),
         [] (const TrackPtr& a, const TrackPtr& b) {
            return a->VetoEnergy() > b->VetoEnergy();
         });

    const TrackPtr& l1 = comb.at(0);
    const TrackPtr& l2 = comb.at(1);
    // suppress conversion decays
    if (l1->CentralVeto() == l2->CentralVeto())
        return;
    steps->Fill("distinct PID", 1);
    const double eeIM = (Particle(ParticleTypeDatabase::eMinus, l1) + Particle(ParticleTypeDatabase::eMinus, l2)).M();
    // suppress pi0
    if (eeIM > 130.)
        return;
    steps->Fill("above #pi^{0}", 1);

    etaIM_final->Fill(eta.M());
    hCopl_final->Fill(abs(eta.Phi() - proton.Phi())*TMath::RadToDeg() - 180.);
    for (const TrackPtr& t : comb)
        if (t->VetoEnergy()) {
            eegPID->Fill(t->VetoEnergy(), t->CentralVeto());
            //h_eegPID->Fill(t->VetoEnergy()*sin(t->Theta()), t->CentralVeto());
        }

}

void ant::analysis::etaDalitz_CalibPID::Finish()
{
    double factor = -prompt_window.Length() / (random_window1.Length() + random_window2.Length());
    diff = prompt;
    diff.AddScaled(random, factor);
}

void ant::analysis::etaDalitz_CalibPID::ShowResult()
{
    canvas c("etaDalitz_CalibPID: Overview");
    c << drawoption("colz") << pid
      << padoption::set(padoption_t::Legend)
      << padoption::unset(padoption_t::Legend)
      << steps << eegPID << cb_esum << endc;

    canvas c_pulls("etaDalitz_CalibPID: Pulls");
    c_pulls << padoption::set(padoption_t::LogY);
    for (auto& p : pulls_diff)
        c_pulls << p.second;
    c_pulls << endc;

    canvas c_fitter("etaDalitz_CalibPID: Fitter");
    c_fitter <<  hChi2 <<  hProb << hIter
             <<  etaIM <<  etaIM_fit << etaIM_final << endc;

    // draw all random subtracted histograms
    //diff.Draw();
}

void ant::analysis::etaDalitz_CalibPID::sort_particles(particle_vector& particles)
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

const vector<Particle> ant::analysis::etaDalitz_CalibPID::get_MC_true_particles(const ParticleList& particles)
{
    if (!true_particles.empty())
        return true_particles;

    Particle proton_true(ParticleTypeDatabase::Neutron, 0, 0, 0);
    for (const auto& p : particles) {
        if (p->Type() == ParticleTypeDatabase::Photon)
            true_particles.push_back(Particle(*p));
        else if (p->Type() == ParticleTypeDatabase::eCharged)
            true_particles.push_back(Particle(*p));
        else if (p->Type() == ParticleTypeDatabase::Proton)
            proton_true = Particle(*p);
    }
    if (proton_true.Type() == ParticleTypeDatabase::Neutron || true_particles.size() != nFinalState-1)
        steps->Fill("wrong MC", 1);

    true_particles.push_back(proton_true);
    return true_particles;
}

const Particle ant::analysis::etaDalitz_CalibPID::get_true_particle(const ParticleList& particles, const size_t pos)
{
    if (true_particles.empty())
        return get_MC_true_particles(particles).at(pos);
    else
        return true_particles.at(pos);
}

double ant::analysis::etaDalitz_CalibPID::calculate_energy_weighted_cb_time_average(const TrackList& tracks) const
{
    double time = 0, e_sum = 0;
    for (const auto& track : tracks) {
        if (track->Detector() & detector_t::anyCB) {
            time += track->Time();
            e_sum += track->ClusterEnergy();
        }
    }
    return time/e_sum;
}

void ant::analysis::etaDalitz_CalibPID::GetTracks(const Event& event, TrackList& tracksCB, TrackList& tracksTAPS)
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

void ant::analysis::etaDalitz_CalibPID::GetParticles(const ant::Event& event, particle_vector& particles)
{
    // collect everything as a photon, determine proton via kinematic fit afterwards
    for (const auto& track : event.Reconstructed().Tracks()) {
        // ignore particles below set cluster energy threshold
        if (track->ClusterEnergy() < CLUSTER_TRESH)
            continue;
        if (track->Detector() & ant::detector_t::NaI) {
            particles.push_back(Particle(ParticleTypeDatabase::Photon, track));
            nParticlesCB++;
        } else if (track->Detector() & ant::detector_t::BaF2 || track->Detector() & ant::detector_t::PbWO4) {
            particles.push_back(Particle(ParticleTypeDatabase::Photon, track));
            nParticlesTAPS++;
        }
    }
    nParticles = nParticlesCB + nParticlesTAPS;
}

void ant::analysis::etaDalitz_CalibPID::GetTrueParticles(const ant::Event& event, particle_vector& particles)
{
    for (const auto& p : event.MCTrue().Particles().GetAll())
        if (contains(ParticleTypeDatabase::MCFinalStateTypes(), p->Type()))
            particles.push_back(Particle(*p));
    nParticles = particles.size();
}

TH1* ant::analysis::etaDalitz_CalibPID::read_hist(const char* file, const char* hist_name) const
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

void ant::analysis::etaDalitz_CalibPID::read_file(const char* file, std::vector<double>& values, const int fill_size)
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

void ant::analysis::etaDalitz_CalibPID::set_fit_particle(const Particle& p, FitParticle& part)
{
    vector<double> sigmas;
    get_uncertainties(p, sigmas);
    part.SetFromParticle(p, sigmas);
}

void ant::analysis::etaDalitz_CalibPID::get_uncertainties(const Particle& p, vector<double>& sigmas)
{
    sigmas.clear();
    sigmas.push_back(get_sigma_energy(p));
    sigmas.push_back(get_sigma_theta(p));
    sigmas.push_back(get_sigma_phi(p));
}

double ant::analysis::etaDalitz_CalibPID::get_sigma_energy(const Particle& p) const
{
    if (p.Type() == ParticleTypeDatabase::Proton)
        return 0;//get_sigma(p, proton_energy_uncertainties);
    else
        return get_sigma(p, photon_energy_uncertainties);
}

double ant::analysis::etaDalitz_CalibPID::get_sigma_theta(const Particle& p) const
{
    if (p.Type() == ParticleTypeDatabase::Proton)
        return get_sigma(p, proton_theta_uncertainties);
    else
        return get_sigma(p, photon_theta_uncertainties);
}

double ant::analysis::etaDalitz_CalibPID::get_sigma_phi(const Particle& p) const
{
    if (p.Type() == ParticleTypeDatabase::Proton)
        return get_sigma(p, proton_phi_uncertainties);
    else
        return get_sigma(p, photon_phi_uncertainties);
}

double ant::analysis::etaDalitz_CalibPID::get_sigma(const Particle& p, TH2* const h) const
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
size_t ant::analysis::etaDalitz_CalibPID::get_histogram_neigbours(TH1* const h,
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
double ant::analysis::etaDalitz_CalibPID::local_bin_average(TH1* const h, const int bin) const
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

void ant::analysis::etaDalitz_CalibPID::apply_time_correction(const TrackList& tracksCB, const TrackList& tracksTAPS)
{
    for (const auto track : tracksCB)
        track->SetTime(track->Time() - cb_time_correction.at(track->CentralCrystal()));
}

void ant::analysis::etaDalitz_CalibPID::apply_energy_correction(const TrackList& tracksCB, const TrackList& tracksTAPS)
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

void ant::analysis::etaDalitz_CalibPID::apply_theta_correction(const TrackList& tracksCB, const TrackList& tracksTAPS)
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
