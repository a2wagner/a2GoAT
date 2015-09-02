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

void analysis::SaschaPhysics::FillIM(TH1D *h, const std::vector<analysis::SaschaPhysics::FitParticle>& final_state)
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

ant::analysis::SaschaPhysics::SaschaPhysics(const mev_t energy_scale) :
    Physics("SaschaPhysics"),
    fitter("SaschaPhysics"),
    final_state(nFinalState)
{
    const BinSettings energy_bins(1600,0,energy_scale*1.6);
    const BinSettings tagger_bins(1500,300,1800);
    const BinSettings ntaggerhits_bins(15);
    const BinSettings veto_bins(1000,0,10);
    const BinSettings particle_bins(10,0,10);
    const BinSettings particlecount_bins(16,0,16);
    const BinSettings pull_bins(50,-3,3);
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

    banana = HistFac.makeTH2D("PID Bananas", "Cluster Energy [MeV]", "Veto Energy [MeV]", energy_bins, veto_bins, "pid");
    particles = HistFac.makeTH1D("Identified particles", "Particle Type", "#", particle_bins, "ParticleTypes");
    tagger = HistFac.makeTH1D("Tagger Spectrum", "Photon Beam Energy", "#", tagger_bins, "TaggerSpectrum");
    ntagged = HistFac.makeTH1D("Tagger Hits", "Tagger Hits / event", "#", ntaggerhits_bins, "nTagged");
    cbesum = HistFac.makeTH1D("CB Energy Sum", "E [MeV]", "#", energy_bins, "esum");

    q2_dist_before = HistFac.makeTH1D("q^{2} Distribution before KinFit", "q^{2} [MeV]", "#", q2_bins, "q2_dist_before");
    q2_dist_after = HistFac.makeTH1D("q^{2} Distribution after KinFit", "q^{2} [MeV]", "#", q2_bins, "q2_dist_after");

    // different checks
    lepton_energies = HistFac.makeTH2D("Lepton Energies", "E(lepton1) [MeV]", "E(lepton(2) [MeV]",
                                       energy_bins, energy_bins, "lepton_energies");
    lepton_energies_true = HistFac.makeTH2D("True Lepton Energies", "E(lepton1)_{true} [MeV]", "E(lepton(2)_{true} [MeV]",
                                            energy_bins, energy_bins, "lepton_energies_true");
    photon_energy_vs_opening_angle = HistFac.makeTH2D("Photon energy vs. opening angle", "Opening angle [#circ]",
                                                      "E_{#gamma} [MeV]", theta_bins, energy_bins,
                                                      "photon_energy_vs_opening_angle");
    photon_energy_vs_opening_angle_true = HistFac.makeTH2D("True photon energy vs. opening angle", "Opening angle [#circ]",
                                                           "E_{#gamma} [MeV]", theta_bins, energy_bins,
                                                           "photon_energy_vs_opening_angle_true");
    theta_vs_clusters = HistFac.makeTH2D("Theta vs. #Clusters", "#Clusters", "#vartheta [#circ]",
                                         particlecount_bins, theta_bins, "theta_vs_clusters");
    opening_angle_vs_q2 = HistFac.makeTH2D("Opening Angle Leptons vs. q^{2}", "q^{2} [GeV^{2}]", "Opening angle [#circ]",
                                           q2_bins, theta_bins, "opening_angle_vs_q2");
    opening_angle_vs_E_high = HistFac.makeTH2D("Opening Angle Leptons vs. E_{high}(Lepton)", "E [MeV]", "Opening angle [#circ]",
                                               energy_bins, theta_bins, "opening_angle_vs_E_high");
    opening_angle_vs_E_low = HistFac.makeTH2D("Opening Angle Leptons vs. E_{low}(Lepton)", "E [MeV]", "Opening angle [#circ]",
                                               energy_bins, theta_bins, "opening_angle_vs_E_low");
    dEvE = HistFac.makeTH2D("dE vs. E", "E_{Crystals} [MeV]", "dE_{Veto} [MeV]", energy_bins, veto_bins, "dEvE");
    crystals_vs_ecl_charged = HistFac.makeTH2D("#Crystals vs. Cluster Energy", "Cluster Energy [GeV]", "#Crystals",
                                               q2_bins, count_bins, "crystals_vs_cluster_energy_charged");
    crystals_vs_ecl_uncharged = HistFac.makeTH2D("#Crystals vs. Cluster Energy", "Cluster Energy [GeV]", "#Crystals",
                                                 q2_bins, count_bins, "crystals_vs_cluster_energy_uncharged");
    crystals_vs_ecl_charged_candidates = HistFac.makeTH2D("#Crystals vs. Cluster Energy (Candidates)",
                                                          "Cluster Energy [GeV]", "#Crystals",
                                                          q2_bins, count_bins, "crystals_vs_ecl_charged_candidates");
    energy_vs_momentum_z_balance = HistFac.makeTH2D("Energy vs. p_{z} balance", "p_{z} [GeV]", "Energy [GeV]",
                                                   energy_range_bins, energy_range_bins, "energy_vs_momentum_z_balance");

    opening_angle_leptons = HistFac.makeTH1D("Lepton opening angle", "Opening angle [#circ]", "#",
                                             theta_bins, "opening_angle_leptons");
    opening_angle_leptons_true = HistFac.makeTH1D("Lepton opening angle true", "Opening angle [#circ]", "#",
                                                  theta_bins, "opening_angle_leptons_true");
    energy_lepton1 = HistFac.makeTH1D("Energy 1st lepton", "E [MeV]", "#", energy_bins, "energy_lepton1");
    energy_lepton1_true = HistFac.makeTH1D("True energy 1st lepton", "E_{true} [MeV]", "#", energy_bins, "energy_lepton1_true");
    energy_lepton2 = HistFac.makeTH1D("Energy 2nd lepton", "E [MeV]", "#", energy_bins, "energy_lepton2");
    energy_lepton2_true = HistFac.makeTH1D("True energy 2nd lepton", "E_{true} [MeV]", "#", energy_bins, "energy_lepton2_true");
    energy_photon = HistFac.makeTH1D("Energy photon", "E [MeV]", "#", energy_bins, "energy_photon");
    energy_photon_true = HistFac.makeTH1D("True energy photon", "E_{true} [MeV]", "#", energy_bins, "energy_photon_true");

    proton_energy = HistFac.makeTH1D("Energy proton", "E [MeV]", "#", energy_bins, "proton_energy");
    proton_energy_true = HistFac.makeTH1D("Energy proton true", "E [MeV]", "#", energy_bins, "proton_energy_true");
    proton_energy_fit = HistFac.makeTH1D("Energy proton fitted", "E [MeV]", "#", energy_bins, "proton_energy_fit");
    proton_energy_delta = HistFac.makeTH1D("#DeltaE_{proton} reconstructed - fitted", "E [MeV]", "#", energy_range_bins,
                                           "proton_energy_delta");
    proton_angle_TAPS_expected = HistFac.makeTH1D("Opening Angle reconstr. Cluster_{TAPS} - expected proton",
                                                  "opening angle [#circ]", "#", angle_bins, "proton_angle_TAPS_expected");

    coplanarity = HistFac.makeTH1D("Coplanarity #eta' proton", "Coplanarity [#circ]", "#", phi_bins, "coplanarity");
    missing_mass = HistFac.makeTH1D("Missing Mass Proton", "m_{miss} [MeV]", "#", energy_bins, "missing_mass");

    for( auto& t : ParticleTypeDatabase::DetectableTypes() ) {
        numParticleType[t]= HistFac.makeTH1D("Number of " + t->PrintName(),
                                      "number of " + t->PrintName() + "/ event",
                                      "", particlecount_bins);
    }

    // prepare invM histograms for different q2 ranges
    im_q2.push_back(im_q2_0_50);
    im_q2.push_back(im_q2_50_100);
    im_q2.push_back(im_q2_100_150);
    im_q2.push_back(im_q2_150_200);
    im_q2.push_back(im_q2_200_250);
    im_q2.push_back(im_q2_250_300);
    im_q2.push_back(im_q2_300_350);
    im_q2.push_back(im_q2_350_400);
    im_q2.push_back(im_q2_400_450);
    im_q2.push_back(im_q2_450_500);
    im_q2.push_back(im_q2_500_550);
    im_q2.push_back(im_q2_550_600);
    im_q2.push_back(im_q2_600_650);
    im_q2.push_back(im_q2_650_700);
    im_q2.push_back(im_q2_700_750);
    im_q2.push_back(im_q2_750_800);
    im_q2.push_back(im_q2_800_850);
    im_q2.push_back(im_q2_850_900);
    double start_range = 0.;
    for (auto& h : im_q2) {
        char title[40];
        char name[20];
        sprintf(name, "im_q2_%d_%d", int(start_range), int(start_range + im_q2_mev_steps));
        sprintf(title, "IM %d MeV < q^{2} < %d MeV", int(start_range), int(start_range + im_q2_mev_steps));
        h = HistFac.makeTH1D(title, "IM [MeV]", "#", im_bins, name);
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
    // if the photon originated from (0,0,v_z) instead of origin,
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

    // make fitter histograms
    chisquare   = HistFac.makeTH1D("ChiSquare","#chi^{2}","#",chisquare_bins,"chisquare");
    probability = HistFac.makeTH1D("Probability","Probability","#",probability_bins,"probability");
    iterations = HistFac.makeTH1D("Number of iterations","Iterations","#",iterations_bins,"iterations");

    // create pull histograms
    for (const auto& varname : fitter.VariableNames()) {
        string title(varname);
        size_t pos = title.find("[");
        if (pos != string::npos) {
            const string prop = " " + component.at(atoi(&varname.at(pos+1)));
            title.replace(pos, 3, prop);
        }
        pulls[varname] = HistFac.makeTH1D("Pull " + title, "Pull", "#", pull_bins, "pull_" + varname);
    }

    stringstream fs;
    fs << "e+e-g";
    im_true = HistFac.makeTH1D("IM "+fs.str()+" true","IM","#",im_bins,"im_true");
    im_smeared = HistFac.makeTH1D("IM "+fs.str()+" smeared","IM","#",im_bins,"im_smeared");
    im_fit = HistFac.makeTH1D("IM "+fs.str()+" fit","IM","#",im_bins,"im_fit");

    vertex_z_before =  HistFac.makeTH1D("Vertex Z Before","v_z / cm","#",vertex_bins,"vertex_z_before");
    vertex_z_after =  HistFac.makeTH1D("Vertex Z After","v_z / cm","#",vertex_bins,"vertex_z_after");

    coplanarity_fit = HistFac.makeTH1D("Coplanarity #eta' proton fitted", "Coplanarity [#circ]", "#", phi_bins,
                                        "coplanarity_fit");
    missing_mass_fit = HistFac.makeTH1D("Missing Mass Proton fitted", "m_{miss} [MeV]", "#", energy_bins, "missing_mass_fit");
    energy_vs_momentum_z_balance_fit = HistFac.makeTH2D("Energy vs. p_{z} balance fitted", "p_{z} [GeV]", "Energy [GeV]",
                                                         energy_range_bins, energy_range_bins, "energy_vs_momentum_z_balance_fit");

    // histograms after cuts
    im_cut = HistFac.makeTH1D("IM "+fs.str()+" after cuts","IM","#",im_bins,"im_cut");
    im_fit_cut = HistFac.makeTH1D("IM "+fs.str()+" fit after cuts","IM","#",im_bins,"im_fit_cut");
    q2_dist_cut = HistFac.makeTH1D("q^{2} Distribution after Cuts", "q^{2} [GeV]", "#", q2_bins, "q2_dist_cut");
    q2_dist_fit_cut = HistFac.makeTH1D("q^{2} Distribution fit after Cuts", "q^{2} [GeV]", "#", q2_bins, "q2_dist_fit_cut");
    coplanarity_cut = HistFac.makeTH1D("Coplanarity #eta' proton after Cuts", "Coplanarity [#circ]", "#", phi_bins,
                                       "coplanarity_cut");
    coplanarity_fit_cut = HistFac.makeTH1D("Coplanarity #eta' proton fit after Cuts", "Coplanarity [#circ]", "#", phi_bins,
                                           "coplanarity_fit_cut");
    missing_mass_cut = HistFac.makeTH1D("Missing Mass Proton after Cuts", "m_{miss} [MeV]", "#", energy_bins, "missing_mass_cut");
    missing_mass_fit_cut = HistFac.makeTH1D("Missing Mass Proton fit after Cuts", "m_{miss} [MeV]", "#", energy_bins,
                                    "missing_mass_fit_cut");
    proton_angle_TAPS_expected_cut = HistFac.makeTH1D("Opening Angle reconstr. Cluster_{TAPS} - expected proton after Cuts",
                                                      "opening angle [#circ]", "#", angle_bins, "proton_angle_TAPS_expected_cut");
    energy_vs_momentum_z_balance_cut = HistFac.makeTH2D("Energy vs. p_{z} balance after Ctus", "p_{z} [GeV]", "Energy [GeV]",
                                                        energy_range_bins, energy_range_bins, "energy_vs_momentum_z_balance_cut");
    energy_vs_momentum_z_balance_fit_cut = HistFac.makeTH2D("Energy vs. p_{z} balance fit after Ctus", "p_{z} [GeV]",
                                                            "Energy [GeV]", energy_range_bins, energy_range_bins,
                                                            "energy_vs_momentum_z_balance_fit_cut");

    dEvE_cut = HistFac.makeTH2D("dE vs. E Cut", "E_{Crystals} [MeV]", "dE_{Veto} [MeV]", energy_bins, veto_bins, "dEvE_cut");
    crystals_vs_ecl_cut = HistFac.makeTH2D("#Crystals vs. Cluster Energy Cut", "Cluster Energy [GeV]", "#Crystals",
                                           q2_bins, count_bins, "crystals_vs_cluster_energy_cut");
    crystals_vs_ecl_charged_cut = HistFac.makeTH2D("#Crystals vs. Cluster Energy Cut", "Cluster Energy [GeV]", "#Crystals",
                                                   q2_bins, count_bins, "crystals_vs_cluster_energy_charged_cut");
    crystals_vs_ecl_uncharged_cut = HistFac.makeTH2D("#Crystals vs. Cluster Energy Cut", "Cluster Energy [GeV]", "#Crystals",
                                                     q2_bins, count_bins, "crystals_vs_cluster_energy_uncharged_cut");

    APLCON::Fit_Settings_t settings = fitter.GetSettings();
    settings.MaxIterations = 8;
    fitter.SetSettings(settings);

    cout.precision(3);
    APLCON::PrintFormatting::Width = 11;
}

void ant::analysis::SaschaPhysics::ProcessEvent(const ant::Event &event)
{
    for (auto& track : event.Reconstructed().Tracks()) {
        banana->Fill(track->ClusterEnergy(), track->VetoEnergy());
    }

    for (auto& particle : event.Reconstructed().Particles().GetAll()) {
        particles->Fill(particle->Type().PrintName().c_str(), 1);
    }

    ntagged->Fill(event.Reconstructed().TaggerHits().size());

    cbesum->Fill(event.Reconstructed().TriggerInfos().CBEenergySum());

    for (auto& t : ParticleTypeDatabase::DetectableTypes()) {
        try {
            numParticleType.at(t)->Fill(event.Reconstructed().Particles().Get(*t).size());
        } catch (...) {}
    }

    TaggerHistList tagger_hits;
    //static size_t count = 0;
    const bool MC = false;
    if (MC)
        tagger_hits = event.MCTrue().TaggerHits();
    else
        tagger_hits = event.Reconstructed().TaggerHits();

    for (const auto& taggerhit : tagger_hits) {
        tagger->Fill(taggerhit->PhotonEnergy());

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

        // check theta distribution in dependence of the number of clusters
        TrackPtr tr;
        bool highVeto = false;
        for (const auto& p : particles) {
            tr = p.Tracks().front();
            if (tr->VetoEnergy() > 3.)
                highVeto = true;
            theta_vs_clusters->Fill(nParticles, p.Theta()*TMath::RadToDeg());
            crystals_vs_ecl_charged_candidates->Fill(tr->ClusterEnergy(), tr->ClusterSize());
            dEvE->Fill(tr->ClusterEnergy(), tr->VetoEnergy());
            if (p.Type().Charged())
                crystals_vs_ecl_charged->Fill(tr->ClusterEnergy(), tr->ClusterSize());
            else
                crystals_vs_ecl_uncharged->Fill(tr->ClusterEnergy(), tr->ClusterSize());
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

        const TLorentzVector target(0., 0., 0., ParticleTypeDatabase::Proton.Mass());
        TLorentzVector balanceP4 = taggerhit->PhotonBeam() + target;
        for (const auto& p : particles)
            balanceP4 -= p;
        energy_vs_momentum_z_balance->Fill(balanceP4.Pz(), balanceP4.E());

        TLorentzVector proton = particles.back();
        TLorentzVector etap(0., 0., 0., 0.);
        for (auto it = particles.begin(); it != particles.end()-1; ++it)
            etap += *it;
        const double copl = abs(etap.Phi() - particles.back().Phi())*TMath::RadToDeg();
        coplanarity->Fill(copl);
        TLorentzVector missingProton = taggerhit->PhotonBeam() + target - etap;
        double missM = missingProton.M();
        missing_mass->Fill(missM);

        proton_energy->Fill(proton.T());
        double openAngle_p_TAPS_expected = proton.Angle(missingProton.Vect())*TMath::RadToDeg();
        proton_angle_TAPS_expected->Fill(openAngle_p_TAPS_expected);

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
                opening_angle_leptons_true->Fill(lepton_open_angle_true);
                photon_energy_vs_opening_angle_true->Fill(lepton_open_angle_true, photon_true.Ek());
                lepton_energies_true->Fill(en_lep1_true, en_lep2_true);
                energy_lepton1_true->Fill(en_lep1_true);
                energy_lepton2_true->Fill(en_lep2_true);
            }
            energy_photon_true->Fill(photon_true.Ek());
            proton_energy_true->Fill(proton_true.Ek());
        }

        double q2_before = (particles[0] + particles[1]).M();
        double lepton_open_angle = particles[0].Angle(particles[1].Vect())*TMath::RadToDeg();
        double en_lep1 = particles[0].T();
        double en_lep2 = particles[1].T();
        if (en_lep1 < en_lep2) {
            en_lep1 = en_lep2;
            en_lep2 = particles[0].T();
        }
        opening_angle_leptons->Fill(lepton_open_angle);
        lepton_energies->Fill(en_lep1, en_lep2);
        energy_lepton1->Fill(en_lep1);
        energy_lepton2->Fill(en_lep2);
        energy_photon->Fill(particles[2].E());
        photon_energy_vs_opening_angle->Fill(lepton_open_angle, particles[2].E());
        opening_angle_vs_q2->Fill(q2_before, lepton_open_angle);
        opening_angle_vs_E_high->Fill(en_lep1 > en_lep2 ? en_lep1 : en_lep2, lepton_open_angle);
        opening_angle_vs_E_low->Fill(en_lep1 < en_lep2 ? en_lep1 : en_lep2, lepton_open_angle);
        q2_dist_before->Fill(q2_before);

        // set proton energy sigma to zero to indicate it's unmeasured
        final_state[3].Ek_Sigma = 0;

//        std::transform(myv1.begin(), myv1.end(), myv1.begin(), [](double d) -> double { return d * 3; });
//        for_each(begin(myv1), end(myv1), [](double& a) { a *= 3; });

        beam.SetFromVector(taggerhit->PhotonBeam());
        // set beam energy sigma to 2 MeV
        beam.Ek_Sigma = 2.;

        FillIM(im_true, final_state);

//        // smear the MC true data
//        proton.Smear();
//        for(auto& photon : photons)
//            photon.Smear();
//        beam.Smear();

//        FillIM(im_smeared, photons);

        // Cut on missing mass of the proton
        if (missM < 900. && missM > 990.)
            continue;
        // Cut on Veto energy
        if (highVeto)
            continue;
        // Cut on opening angle between TAPS cluster and expected proton
        // (mainly for suppressing bad events which will be unnecassarily fitted, affects only a few background events)
        if (openAngle_p_TAPS_expected > 5.)
            continue;

        // let APLCON do the work
        const APLCON::Result_t& result = fitter.DoFit();

        //cout << result << endl;

        if (result.Status != APLCON::Result_Status_t::Success) {
            //cout << result << endl;
            continue;
        }

        for (const auto& it_map : result.Variables) {
            const string& varname = it_map.first;
            const APLCON::Result_Variable_t& var = it_map.second;
            pulls.at(varname)->Fill(var.Pull);
        }
        chisquare->Fill(result.ChiSquare);
        probability->Fill(result.Probability);
        iterations->Fill(result.NIterations);

        if (includeVertexFit) {
            vertex_z_after->Fill(result.Variables.at("v_z").Value.After);
            vertex_z_before->Fill(result.Variables.at("v_z").Value.Before);
        }

        FillIM(im_fit, final_state);

        double q2_after = (FitParticle::Make(final_state[0], ParticleTypeDatabase::eMinus.Mass())
                + FitParticle::Make(final_state[1], ParticleTypeDatabase::eMinus.Mass())).M();
        q2_dist_after->Fill(q2_after);

        TLorentzVector balanceP4_fit = target + FitParticle::Make(beam, ParticleTypeDatabase::Photon.Mass());
        // create the fitted final state particles
        TLorentzVector proton_fit = FitParticle::Make(final_state[3], ParticleTypeDatabase::Proton.Mass());
        TLorentzVector etap_fit(0., 0., 0., 0.);
        etap_fit += FitParticle::Make(final_state[0], ParticleTypeDatabase::eMinus.Mass());
        etap_fit += FitParticle::Make(final_state[1], ParticleTypeDatabase::eMinus.Mass());
        etap_fit += FitParticle::Make(final_state[2], ParticleTypeDatabase::Photon.Mass());
        balanceP4_fit -= etap_fit;
        balanceP4_fit -= proton_fit;
        energy_vs_momentum_z_balance_fit->Fill(balanceP4_fit.Pz(), balanceP4_fit.E());

        const double copl_fit = abs(etap_fit.Phi() - particles.back().Phi())*TMath::RadToDeg();
        coplanarity_fit->Fill(copl_fit);
        TLorentzVector missingProtonFit = target + FitParticle::Make(beam, ParticleTypeDatabase::Photon.Mass()) - etap_fit;
        missing_mass_fit->Fill(missingProtonFit.M());

        proton_energy_fit->Fill(proton_fit.T());
        proton_energy_delta->Fill(proton.T() - proton_fit.T());

        // Cut on chi^2
        if (result.ChiSquare > 10.)
            return;

        // fill the invM histograms for different q2 ranges
        if (q2_after < im_q2.size()*im_q2_mev_steps)
            im_q2.at(static_cast<int>(q2_after/im_q2_mev_steps))->Fill(etap_fit.M());

        im_cut->Fill(etap.M());
        FillIM(im_fit_cut, final_state);
        q2_dist_cut->Fill(q2_before);
        q2_dist_fit_cut->Fill(q2_after);
        coplanarity_cut->Fill(copl);
        missing_mass_cut->Fill(missM);
        missing_mass_fit_cut->Fill(missingProtonFit.M());
        coplanarity_fit_cut->Fill(copl_fit);
        proton_angle_TAPS_expected_cut->Fill(openAngle_p_TAPS_expected);
        energy_vs_momentum_z_balance_cut->Fill(balanceP4.Pz(), balanceP4.E());
        energy_vs_momentum_z_balance_fit_cut->Fill(balanceP4_fit.Pz(), balanceP4_fit.E());

        for (const auto& p : particles) {
            tr = p.Tracks().front();
            crystals_vs_ecl_cut->Fill(tr->ClusterEnergy(), tr->ClusterSize());
            dEvE_cut->Fill(tr->ClusterEnergy(), tr->VetoEnergy());
            if (p.Type().Charged())
                crystals_vs_ecl_charged_cut->Fill(tr->ClusterEnergy(), tr->ClusterSize());
            else
                crystals_vs_ecl_uncharged_cut->Fill(tr->ClusterEnergy(), tr->ClusterSize());
        }
    }
}

void ant::analysis::SaschaPhysics::Finish()
{
}

void ant::analysis::SaschaPhysics::ShowResult()
{
    canvas c("SaschaPhysics: Overview");
    c << drawoption("colz") << banana
      << padoption::set(padoption_t::Legend)
      << particles
      << padoption::unset(padoption_t::Legend)
      << tagger << ntagged << cbesum << endc;

    canvas c_pulls("SaschaPhysics: Pulls");
    c_pulls << padoption::set(padoption_t::LogY);
    for (auto& p : pulls)
        c_pulls << p.second;
    c_pulls << endc;

    canvas c_fitter("SaschaPhysics: Fitter");
    c_fitter << chisquare << probability << iterations
             << im_true << im_smeared << im_fit
             << vertex_z_before << vertex_z_after << endc;
}

void ant::analysis::SaschaPhysics::GetParticles(const ant::Event& event, particle_vector& particles)
{
    // example of how to collect particles in a user defined way
    for (const auto& track : event.Reconstructed().Tracks()) {
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

/*    // sort out particles with an energy of less than 10 MeV
    for (particle_it it = particles->begin(); it != particles->end();)
        if (it->E() < 10.)
            it = particles->erase(it);
        else
            ++it;*/
}

void ant::analysis::SaschaPhysics::GetTrueParticles(const ant::Event& event, particle_vector& particles)
{
    for (const auto& p : event.MCTrue().Particles().GetAll())
        if (contains(ParticleTypeDatabase::MCFinalStateTypes(), p->Type()))
            particles.push_back(Particle(*p));
    nParticles = particles.size();
}
