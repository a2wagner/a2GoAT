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

    banana = HistFac.makeTH2D("PID Bananas", "Cluster Energy [MeV]", "Veto Energy [MeV]", energy_bins, veto_bins, "pid");
    particles = HistFac.makeTH1D("Identified particles", "Particle Type", "#", particle_bins, "ParticleTypes");
    tagger = HistFac.makeTH1D("Tagger Spectrum", "Photon Beam Energy", "#", tagger_bins, "TaggerSpectrum");
    ntagged = HistFac.makeTH1D("Tagger Hits", "Tagger Hits / event", "#", ntaggerhits_bins, "nTagged");
    cbesum = HistFac.makeTH1D("CB Energy Sum", "E [MeV]", "#", energy_bins, "esum");

    for( auto& t : ParticleTypeDatabase::DetectableTypes() ) {
        numParticleType[t]= HistFac.makeTH1D("Number of " + t->PrintName(),
                                      "number of " + t->PrintName() + "/ event",
                                      "", particlecount_bins);
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
    chisquare   = HistFac.makeTH1D("ChiSqare","ChiSquare","#",chisquare_bins,"chisquare");
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

    APLCON::Fit_Settings_t settings = fitter.GetSettings();
    settings.MaxIterations = 50;
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

        //std::copy(particles.begin(), particles.end(), final_state);
//        std::transform(particles.begin(), particles.end(), final_state.begin(),
//                       [](Particle& p) -> double { return FitParticle().SetFromVector(p); });
        auto it = final_state.begin();
        for (const auto& p : particles)
            (it++)->SetFromVector(p);

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
                particles.push_back(Particle(ParticleTypeDatabase::eMinus,
                                             track->ClusterEnergy(),
                                             track->Theta(), track->Phi()));
            else
                particles.push_back(Particle(ParticleTypeDatabase::Photon,
                                             track->ClusterEnergy(),
                                             track->Theta(), track->Phi()));
        } else if (track->Detector() & ant::detector_t::BaF2 || track->Detector() & ant::detector_t::PbWO4) {
            nParticlesTAPS++;
            //if (track->VetoEnergy() > 0.)  // Veto entry? --> charged
            if (track->Detector() & detector_t::Veto)
                particles.push_back(Particle(ParticleTypeDatabase::Proton,
                                             track->ClusterEnergy(),
                                             track->Theta(), track->Phi()));
            else
                particles.push_back(Particle(ParticleTypeDatabase::Photon,
                                             track->ClusterEnergy(),
                                             track->Theta(), track->Phi()));
        }
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
