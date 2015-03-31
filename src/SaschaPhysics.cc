#include "SaschaPhysics.h"

SaschaPhysics::SaschaPhysics() : MC(false), pluto(nullptr), geant(nullptr)
{
	n_particles = new TH1I("n_part", "n_particles", 16, 0, 15);
	invM_eeg = new TH1D("invM_eeg", "inv_mass_eeg", 3000, 0, 1500);
	invM_2e = new TH1D("invM_2e", "inv_mass_2e", 2000, 0, 1000);
	missM = new TH1D("missM", "miss_mass", 600, 500, 1100);
	expectProt = new TH1D("expectProt", "expected_proton", 500, 700, 1200);

	balanceE = new TH1D("balanceE", "balanceE", 1200, -200, 200);
	balancePx = new TH1D("balancePx", "balancePx", 1200, -200, 200);
	balancePy = new TH1D("balancePy", "balancePy", 1200, -200, 200);
	balancePz = new TH1D("balancePz", "balancePz", 1200, -200, 200);

	MCbeamDeltaE = new TH1D("MCbeamDeltaE", "MCbeamDeltaE", 2000, -500, 500);
	MCbeamDeltaP = new TH1D("MCbeamDeltaP", "MCbeamDeltaP", 2000, -500, 500);
	MClepton1DeltaE = new TH1D("MClepton1DeltaE", "MClepton1DeltaE", 2000, -500, 500);
	MClepton1DeltaP = new TH1D("MClepton1DeltaP", "MClepton1DeltaP", 2000, -500, 500);
	MClepton2DeltaE = new TH1D("MClepton2DeltaE", "MClepton2DeltaE", 2000, -500, 500);
	MClepton2DeltaP = new TH1D("MClepton2DeltaP", "MClepton2DeltaP", 2000, -500, 500);
	MCphotonDeltaE = new TH1D("MCphotonDeltaE", "MCphotonDeltaE", 2000, -500, 500);
	MCphotonDeltaP = new TH1D("MCphotonDeltaP", "MCphotonDeltaP", 2000, -500, 500);
	MCprotonDeltaE = new TH1D("MCprotonDeltaE", "MCprotonDeltaE", 2000, -500, 500);
	MCprotonDeltaP = new TH1D("MCprotonDeltaP", "MCprotonDeltaP", 2000, -500, 500);

	etap_fit_steps = new TH2D("etap_fit_steps", "etap_fit_steps", 1000, 700, 1200, 50, 0, 49);
	coplanarity = new TH2D("coplanarity", "coplanarity", 720, 0, 360, 50, 0, 49);

	missM_etap_vs_pTheta = new TH2D("missM_etap_vs_pTheta", "missM_etap_vs_pTheta", 1000, 700, 1200, 240, 0, 120);
	etap_invM_gg_vs_beamE = new TH2D("etap_invM_gg_vs_beamE", "etap_invM_gg_vs_beamE", 1000, 700, 1200, 360, 1440, 1620);
	etap_p_kinE_vs_pTheta_true = new TH2D("etap_p_kinE_vs_pTheta_true", "etap_p_kinE_vs_pTheta_true", 1400, 0, 700, 120, 0, 40);
	etap_invM_gg_smeared_vs_invM_gg_fit = new TH2D("etap_invM_gg_smeared_vs_invM_gg_fit", "etap_invM_gg_smeared_vs_invM_gg_fit", 1000, 700, 1200, 1000, 700, 1200);

	lepton1PullE = new TH1D("lepton1PullE", "lepton1_pullE", 5000, -25, 25);
	lepton1PullTheta = new TH1D("lepton1PullTheta", "lepton1_pullTheta", 5000, -25, 25);
	lepton1PullPhi = new TH1D("lepton1PullPhi", "lepton1_pullPhi", 5000, -25, 25);
	lepton2PullE = new TH1D("lepton2PullE", "lepton2_pullE", 5000, -25, 25);
	lepton2PullTheta = new TH1D("lepton2PullTheta", "lepton2_pullTheta", 5000, -25, 25);
	lepton2PullPhi = new TH1D("lepton2PullPhi", "lepton2_pullPhi", 5000, -25, 25);
	photonPullE = new TH1D("photonPullE", "photon_pullE", 5000, -25, 25);
	photonPullTheta = new TH1D("photonPullTheta", "photon_pullTheta", 5000, -25, 25);
	photonPullPhi = new TH1D("photonPullPhi", "photon_pullPhi", 5000, -25, 25);
	protonPullE = new TH1D("protonPullE", "proton_pullE", 5000, -25, 25);
	protonPullTheta = new TH1D("protonPullTheta", "proton_pullTheta", 5000, -25, 25);
	protonPullPhi = new TH1D("protonPullPhi", "proton_pullPhi", 5000, -25, 25);

	chisq_hist = new TH1D("chisq", "chi_square", 1000, 0, 100);
	prob_hist = new TH1D("prob", "probability", 1000, 0, 1);
	etap_smeared = new TH1D("etap_smeared", "etap_mass_smeared", 1000, 700, 1200);
	etap_fitted = new TH1D("etap_fitted", "etap_mass_fitted", 1000, 700, 1200);
	copl_smeared = new TH1D("copl_smeared", "coplanarity_smeared", 720, 0, 360);
	copl_fitted = new TH1D("copl_fitted", "coplanarity_fitted", 720, 0, 360);
	nIter_converged = new TH1I("n_iter_converged", "n_iter_converged", 52, 0, 51);
	nIter_all = new TH1I("n_iter_all", "n_iter_all", 52, 0, 51);
	fitStatus = new TH1I("fit_status", "fit_status", 21, -10, 10);

	q2 = new TH1D("q2", "q2", 2000, 0, 1000);
	q2_cut = new TH1D("q2_cut", "q2_cut", 2000, 0, 1000);
}

SaschaPhysics::~SaschaPhysics()
{
}

Bool_t SaschaPhysics::Init(const char* configfile)
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;
	if(configfile) SetConfigFile(configfile);

	if(!InitBackgroundCuts()) return kFALSE;
	if(!InitTargetMass()) return kFALSE;
	if(!InitTaggerChannelCuts()) return kFALSE;
	if(!InitTaggerScalers()) return kFALSE;
	cout << "--------------------------------------------------" << endl;

	/* Initialise class for kinematic fitting */
	kinFit = KinFit();

	return kTRUE;
}

Bool_t SaschaPhysics::Start()
{
	if(!IsGoATFile())
	{
		cout << "ERROR: Input File is not a GoAT file." << endl;
		return kFALSE;
	}
	SetAsPhysicsFile();


	/* check if the file contains MC data */
	if (!GetScalers()->IsOpenForInput())
		MC = true;
	else
		MC = false;

	if (MC) {
		/* get the MC trees and check if they really exist because you'll still get vectors etc. from non-existent trees! */
		//TODO: FIX THIS?!
		geant = GetGeant();
		pluto = GetPluto();
		if(!(geant_tree = geant->IsOpenForInput()))
			cout << "No Geant Tree in the current file!" << endl;
		if(!(pluto_tree = pluto->IsOpenForInput()))
			cout << "No Pluto Tree in the current file!" << endl;

		if (!pluto_tree && !geant_tree) {
			cout << "No MC trees found in the current file, don't expect any reasonable results if you're proceeding!" << endl;
			exit(EXIT_FAILURE);
		}
	}


	TraverseValidEvents();

	return kTRUE;
}

void SaschaPhysics::ProcessEvent()
{
	/* declare some variables used for kinematic fitting */
	bool dbg = false;
	static int count = 0;//if(count++>5)exit(0);
	double lepton1_pullE, lepton1_pullTheta, lepton1_pullPhi,
			lepton2_pullE, lepton2_pullTheta, lepton2_pullPhi,
			photon_pullE, photon_pullTheta, photon_pullPhi,
			proton_pullE, proton_pullTheta, proton_pullPhi;
	double chisq, prob, etap_mass_smeared, etap_mass_fitted,
			coplanarity_smeared, coplanarity_fitted;
	int n_iter, fit_status;

	/* Gather all particle information, depending on which MC trees are available */
	bool MeV = false;  // use MeV or GeV?
	//TODO fit converges significantly slower (3 times or even more iteration steps) for MeV ?
	nParticles = 0, nParticlesCB = 0, nParticlesTAPS = 0;
	particle_vector particles;
	particle_t beam(photon);
	//TODO add a true flag for MC data to the config file
	bool trueMC = false;
	if (MC && trueMC)
		GetTrueParticles(&particles, &beam);
	else
		GetParticles(&particles, &beam);
	n_particles->Fill(particles.size());
	if (!beam.E())  // does this event contain a beam particle? skip it if this is not the case
		return;
	particle_t target(TLorentzVector(0., 0., 0., MeV ? MASS_PROTON : MASS_PROTON/1000.), proton);
	if (MC && MeV)
		beam.GeV2MeV();
	nParticles = particles.size();
	if (MC && !trueMC && !MeV) {
		for (auto& p : particles)
			p.MeV2GeV();
		//beam.MeV2GeV();
		beam.p4.SetE(beam.E()/1000.);
		beam.p4.SetPz(beam.Pz()/1000.);
	}

	if (dbg) {
		printf("\nEvent %d:\n", ++count);
		for (auto& p : particles) {
			printf("Particle id: %2d (%s), ", p.id, p.is_charged() ? "charged" : "uncharged");
			p.Print();
		}
		printf("Beam: ");
		beam.Print();
	}

	for (particle_t p : particles)
		if (!p.is_final_state())
			throw noFS;


	TLorentzVector tmpState(0., 0., 0., 0.);
	unsigned int nCharged = 0, nNeutral = 0;
	unsigned int charged[N_FINAL_STATE];
	unsigned int neutral[N_FINAL_STATE];
	TLorentzVector missProton;
	// variables to store values only for cuts
	Double_t mMiss, protExpect;

	if (nParticles == N_FINAL_STATE) {  // proceed if there are 3 particles in total (2g and proton)
//		std::cout << "before:" << std::endl;
//		std::for_each(particles.begin(), particles.end(), [](const particle_t& p){ std::cout << p.id << "\n"; });
		// swap the proton and the last vector entry if the proton is not the last particle
		for (auto it = particles.begin(); it != particles.end()-1; ++it)
			if (it->id == proton) {
				std::iter_swap(it, particles.end()-1);
				break;
			}
//		std::cout << "after:" << std::endl;
//		std::for_each(particles.begin(), particles.end(), [](const particle_t& p){ std::cout << p.id << "\n"; });

		int i = 0;
		for (auto& p : particles) {
			if (MeV)
				p.GeV2MeV();
			if (p.is_charged())
				charged[nCharged++] = i;
			else
				neutral[nNeutral++] = i;
			if (p.id != proton)
				tmpState += particles[i];
			i++;
		}

		/* end event selection when we don't have three charged particle (entry in PID or Veto) and a neutral one as well as one proton*/
		std::vector<particle_id> ids;
		std::transform(particles.begin(), particles.end(), back_inserter(ids), [](const particle_t& p){ return p.id; });
		int n_protons = std::count(ids.begin(), ids.end(), proton);
		if (nCharged != 3 || nNeutral != 1 || n_protons != 1)
			return;

//		std::cout << "particle ids:" << std::endl;
//		std::copy(ids.begin(), ids.end(), std::ostream_iterator<particle_id>(std::cout, "\n"));
//std::cout << "Match signature pattern " << ++count << std::endl;
		double protEnergyReconstr = particles[charged[0]].E() - particles[charged[0]].M();
		double invM_2eg = MeV ? tmpState.M() : tmpState.M()*1000.;
		double invM_2charged = (particles[charged[0]] + particles[charged[1]]).M();
		invM_eeg->Fill(invM_2eg);
		invM_2e->Fill(MeV ? invM_2charged : invM_2charged*1000.);

		/* at this point we have one charged and two neutral particles */
		missProton = target + beam - tmpState;
		mMiss = missProton.M();
		protExpect = missProton.E() - missProton.M();
		TLorentzVector balance = target + beam - tmpState - particles[charged[0]];
		missM->Fill(MeV ? mMiss : mMiss*1000.);
		expectProt->Fill(MeV ? protExpect : protExpect*1000.);
		balanceE->Fill(MeV ? balance.E() : balance.E()*1000.);
		balancePx->Fill(MeV ? balance.Px() : balance.Px()*1000.);
		balancePy->Fill(MeV ? balance.Py() : balance.Py()*1000.);
		balancePz->Fill(MeV ? balance.Pz() : balance.Pz()*1000.);

		/* calculate differences between true MC data and detector smeared data */
		double mcBeamDeltaE, mcBeamDeltaP, mcPhotonDeltaE, mcPhotonDeltaP,
				mcLepton1DeltaE, mcLepton1DeltaP, mcLepton2DeltaE, mcLepton2DeltaP,
				mcProtonDeltaE, mcProtonDeltaP;
		particle_t trueBeam;
		particle_vector trueParts;
		TLorentzVector delta;
		match matches;
		GetTrueParticles(&trueParts, &trueBeam);

		MatchMC(&particles, &trueParts, &matches);

		bool mismatch = false;
		std::vector<unsigned int> values;
		std::transform(matches.begin(), matches.end(), back_inserter(values), [](match::value_type& val){ return val.second; });
		for (unsigned int i = 0; i < matches.size(); i++) {
			int sum = std::count(values.begin(), values.end(), i);
			if (sum > 1) {
				mismatch = true;
				if (MC && !trueMC && dbg) {
					printf("WARNING: MC particle mismatch, more than one occurrence of value %d!\n", i);
					std::copy(values.begin(), values.end(), std::ostream_iterator<unsigned int>(std::cout, "\n"));
				}
			}
		}

		if (!mismatch) {
			delta = trueBeam - beam;
			mcBeamDeltaE = delta.E();
			mcBeamDeltaP = delta.P();
			delta = trueParts[matches[charged[0]]] - particles[charged[0]];
			mcLepton1DeltaE = delta.E();
			mcLepton1DeltaP = delta.P();
			delta = trueParts[matches[charged[1]]] - particles[charged[1]];
			mcLepton2DeltaE = delta.E();
			mcLepton2DeltaP = delta.P();
			delta = trueParts[matches[neutral[0]]] - particles[neutral[0]];
			mcPhotonDeltaE = delta.E();
			mcPhotonDeltaP = delta.P();
			delta = trueParts[matches[charged[2]]] - particles[charged[2]];
			mcProtonDeltaE = delta.E();
			mcProtonDeltaP = delta.P();

			MCbeamDeltaE->Fill(MeV ? mcBeamDeltaE : mcBeamDeltaE*1000.);
			MCbeamDeltaP->Fill(MeV ? mcBeamDeltaP : mcBeamDeltaP*1000.);
			MClepton1DeltaE->Fill(MeV ? mcLepton1DeltaE : mcLepton1DeltaE*1000.);
			MClepton1DeltaP->Fill(MeV ? mcLepton1DeltaP : mcLepton1DeltaP*1000.);
			MClepton2DeltaE->Fill(MeV ? mcLepton2DeltaE : mcLepton2DeltaE*1000.);
			MClepton2DeltaP->Fill(MeV ? mcLepton2DeltaP : mcLepton2DeltaP*1000.);
			MCphotonDeltaE->Fill(MeV ? mcPhotonDeltaE : mcPhotonDeltaE*1000.);
			MCphotonDeltaP->Fill(MeV ? mcPhotonDeltaP : mcPhotonDeltaP*1000.);
			MCprotonDeltaE->Fill(MeV ? mcProtonDeltaE : mcProtonDeltaE*1000.);
			MCprotonDeltaP->Fill(MeV ? mcProtonDeltaP : mcProtonDeltaP*1000.);
		}

		/* prepare kinematic fit */
		int err;
		int ndf;

		/* reset the kinematic fit for a new event */
		kinFit.reset();

		TVector3 vecLepton1 = particles[charged[0]].Vect();
		TVector3 vecLepton2 = particles[charged[1]].Vect();
		TVector3 vecPhoton = particles[neutral[0]].Vect();
		TVector3 vecProton = particles[charged[2]].Vect();

		// temp checks
		/*TLorentzVector balance = particles[neutral[0]].p4 + particles[neutral[1]].p4 + particles[charged[0]].p4 - beam.p4 - target.p4;
		printf("Proton: %f, Pi0: %f, BeamE: %f, BeamM: %f\n", particles[charged[0]].M(), (particles[neutral[0]].p4 + particles[neutral[1]].p4).M(), beam.E(), beam.M());
		printf("Energy-Balance: %f, Momentum: (%f, %f, %f)\n", balance.E(), balance.Px(), balance.Py(), balance.Pz());
		printf("MissM: %f, invM_2g: %f\n", mMiss, invM_2neutral);
		printf("Ekin expected Proton: %f, true: %f\n", protExpect, particles[charged[0]].E() - particles[charged[0]].M());
		printf("ProtE = %f, ProtP = %f\n", particles[charged[0]].E(), particles[charged[0]].p4.P());
		printf("Beam: "); beam.Print();
		printf("Photon1: "); particles[neutral[0]].Print();
		printf("Photon2: "); particles[neutral[1]].Print();
		printf("Proton: "); particles[charged[0]].Print();*/

		TLorentzVector fitLepton1;
		TLorentzVector fitLepton2;
		TLorentzVector fitPhoton;
		TLorentzVector fitProton;

		/**
		 * Use absolute or relative smearing?
		 * For absolute smearing use the sigmaX values defined above and comment the ones for each particle
		 * Relative smearing is calculated for each particle using the following data:
		 *   CB: sigma_theta = 2.5°; sigma_phi = sigma_theta/sin(theta); sigma_E = E^.36*.02*E
		 *   TAPS: sigma_theta/phi = 1°; sigma_E = (.018+.008*sqrt(E))*E
		 */
		TMatrixD covLepton1;
		TMatrixD covLepton2;
		TMatrixD covPhoton;
		TMatrixD covProton;

		TRandom3 rand(0);
		unsigned int nPar = 3;
		double sigmaP, sigmaTheta, sigmaPhi;
		int rows = nPar;  // number of rows equal to number of cols
		double errors[nPar];
		double factor = 1.;  // factor to increase/decrease sigma^2 values used for kinematic fitting
		// absolute smearing
		/*sigmaTheta = .05, sigmaPhi = .05;
		if (MeV)
			sigmaP = 20.;
		else
			sigmaP = .02;
		Double_t errors[] = {sigmaP*sigmaP*factor, sigmaTheta*sigmaTheta*factor, sigmaPhi*sigmaPhi*factor};*/
		// relative smearing; assume sigmaP ~ sigmaE
		/* 1st lepton */
		int currPart = charged[0];
		sigmaTheta = sigma_theta(&particles[currPart].p4);
		sigmaPhi = sigma_phi(&particles[currPart].p4);
		sigmaP = sigma_E(&particles[currPart].p4);
		if (trueMC)
			vecLepton1.SetMagThetaPhi(
				rand.Gaus(vecLepton1.Mag(), sigmaP),
				rand.Gaus(vecLepton1.Theta(), sigmaTheta),
				rand.Gaus(vecLepton1.Phi(), sigmaPhi)
			);
		kinFit.fill_array(errors, {sigmaP*sigmaP*factor, sigmaTheta*sigmaTheta*factor, sigmaPhi*sigmaPhi*factor});
		//Double_t errors[] = {1., .01, .01};
		//kinFit.sigmaEThetaPhi(particles[currPart], errors);
		if (kinFit.fillSquareMatrixDiagonal(&covLepton1, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		/* 2nd lepton */
		currPart = charged[1];
		sigmaTheta = sigma_theta(&particles[currPart].p4);
		sigmaPhi = sigma_phi(&particles[currPart].p4);
		sigmaP = sigma_E(&particles[currPart].p4);
		if (trueMC)
			vecLepton2.SetMagThetaPhi(
				rand.Gaus(vecLepton2.Mag(), sigmaP),
				rand.Gaus(vecLepton2.Theta(), sigmaTheta),
				rand.Gaus(vecLepton2.Phi(), sigmaPhi)
			);
		kinFit.fill_array(errors, {sigmaP*sigmaP*factor, sigmaTheta*sigmaTheta*factor, sigmaPhi*sigmaPhi*factor});
		//kinFit.sigmaEThetaPhi(particles[currPart], errors);
		if (kinFit.fillSquareMatrixDiagonal(&covLepton2, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		/* photon */
		currPart = neutral[0];
		sigmaTheta = sigma_theta(&particles[currPart].p4);
		sigmaPhi = sigma_phi(&particles[currPart].p4);
		sigmaP = sigma_E(&particles[currPart].p4);
		if (trueMC)
			vecPhoton.SetMagThetaPhi(
				rand.Gaus(vecPhoton.Mag(), sigmaP),
				rand.Gaus(vecPhoton.Theta(), sigmaTheta),
				rand.Gaus(vecPhoton.Phi(), sigmaPhi)
			);
		kinFit.fill_array(errors, {sigmaP*sigmaP*factor, sigmaTheta*sigmaTheta*factor, sigmaPhi*sigmaPhi*factor});
		//kinFit.sigmaEThetaPhi(particles[currPart], errors);
		if (kinFit.fillSquareMatrixDiagonal(&covPhoton, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		/* proton */
		currPart = charged[2];
		sigmaTheta = sigma_theta(&particles[currPart].p4);
		sigmaPhi = sigma_phi(&particles[currPart].p4);
		sigmaP = sigma_E(&particles[currPart].p4);
		if (trueMC)
			vecProton.SetMagThetaPhi(
				rand.Gaus(vecProton.Mag(), sigmaP),
				rand.Gaus(vecProton.Theta(), sigmaTheta),
				rand.Gaus(vecProton.Phi(), sigmaPhi)
			);
		kinFit.fill_array(errors, {sigmaP*sigmaP*factor, sigmaTheta*sigmaTheta*factor, sigmaPhi*sigmaPhi*factor});
		//kinFit.sigmaEThetaPhi(particles[currPart], errors);
		if (kinFit.fillSquareMatrixDiagonal(&covProton, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");

		TFitParticlePThetaPhi e1("lepton1", "lepton1", &vecLepton1, 0., &covLepton1);
		TFitParticlePThetaPhi e2("lepton2", "lepton2", &vecLepton2, 0., &covLepton2);
		TFitParticlePThetaPhi ph("photon", "photon", &vecPhoton, 0., &covPhoton);
		TFitParticlePThetaPhi pr("proton", "proton", &vecProton, MeV ? MASS_PROTON : MASS_PROTON/1000., &covProton);

		TVector3 vecBeam = beam.Vect();
		TVector3 vecTarget = target.Vect();
		TMatrixD covBeam;
		TMatrixD covTarget;
		//kinFit.sigmaEThetaPhi(Tagged[0], errors);
		errors[0] = .01; errors[1] = .0001; errors[2] = .0001;
		if (MeV)
			errors[0] = 10.;
		if (kinFit.fillSquareMatrixDiagonal(&covBeam, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		covTarget.Zero();
		covTarget.ResizeTo(3, 3);
		TFitParticlePThetaPhi bm("beam", "beam", &vecBeam, 0., &covBeam);
		TFitParticlePThetaPhi trgt("target", "target", &vecTarget, MeV ? MASS_PROTON : MASS_PROTON/1000., &covTarget);

		// energy and momentum constraints have to be defined separately for each component
		// components can be accessed via enum TFitConstraintEp::component
		Double_t constraint = beam.E() + target.E();
		TFitConstraintEp energyConservation("energyConstr", "Energy conservation constraint", 0, TFitConstraintEp::E, constraint);
		energyConservation.addParticles1(&e1, &e2, &ph, &pr);
		constraint = beam.Px();
		TFitConstraintEp pxConservation("pxConstr", "Px conservation constraint", 0, TFitConstraintEp::pX, constraint);
		pxConservation.addParticles1(&e1, &e2, &ph, &pr);
		constraint = beam.Py();
		TFitConstraintEp pyConservation("pyConstr", "Py conservation constraint", 0, TFitConstraintEp::pY, constraint);
		pyConservation.addParticles1(&e1, &e2, &ph, &pr);
		constraint = beam.Pz();
		TFitConstraintEp pzConservation("pzConstr", "Pz conservation constraint", 0, TFitConstraintEp::pZ, constraint);
		pzConservation.addParticles1(&e1, &e2, &ph, &pr);
		TFitConstraintM missMassConstrProton("missMassConstrProton", "missing mass constraint proton", 0, 0, MeV ? MASS_PROTON : MASS_PROTON/1000.);
		missMassConstrProton.addParticles1(&bm, &trgt);
		missMassConstrProton.addParticles2(&e1, &e2, &ph);
		TFitConstraintM massConstrProton("massConstr_proton", "mass constraint proton", 0, 0, MeV ? MASS_PROTON : MASS_PROTON/1000.);
		massConstrProton.addParticle1(&pr);
		// constraint for pi0 mass
		TFitConstraintM massConstrPi0("massConstr_pi0", "mass constraint pi0", 0, 0, MeV ? MASS_PIZERO : MASS_PIZERO/1000.);
		massConstrPi0.addParticles1(&e1, &e2, &ph);

		//TKinFitter fit;
		kinFit.addMeasParticle(&e1);
		kinFit.addMeasParticle(&e2);
		kinFit.addMeasParticle(&ph);
		kinFit.addMeasParticle(&pr);
		if (MC && !trueMC)  // exclude proton energy from fit for detector smeared (and real) data due to punch-throughs
			kinFit.setParamUnmeas(&pr, 0);  // proton energy unmeasured
		kinFit.addMeasParticle(&bm);
		//kinFit.setParamUnmeas(&bm, 1);  // beam theta unmeasured
		//kinFit.setParamUnmeas(&bm, 2);  // beam phi unmeasured
		//kinFit.addMeasParticle(&trgt);

		//kinFit.addConstraint(&massConstrProton);
		kinFit.addConstraint(&energyConservation);
		kinFit.addConstraint(&pxConservation);
		kinFit.addConstraint(&pyConservation);
		kinFit.addConstraint(&pzConservation);
		//kinFit.addConstraint(&missMassConstrProton);  --> leads to a funny beam energy - eta' mass correlation...
		//kinFit.addConstraint(&massConstrPi0);

		/* get the intermediate steps from the fitter */
		//vector<TAbsFitParticle> *temp_results = new vector<TAbsFitParticle>;
		vvP4 temp_results;
		/*vector<TAbsFitParticle*> track_particles;
		track_particles.push_back(&ph1);
		track_particles.push_back(&ph2);
		track_particles.push_back(&pr);*/

		//kinFit.get_intermediate_steps(&temp_results, track_particles);
		kinFit.get_intermediate_steps(&temp_results, {&e1, &e2, &ph, &pr});

		kinFit.setMaxNbIter(50);  // number of maximal iterations
		kinFit.setMaxDeltaS(5e-5);  // max Delta chi2
		kinFit.setMaxF(MeV ? .1 : 1e-4);  // max sum of constraints
		// set verbosity level
		if (dbg)
			kinFit.setVerbosity(3);
		else
			kinFit.setVerbosity(0);
		kinFit.fit();

		// save fit_status and n_iter even if fit diverged
		fit_status = kinFit.getStatus();
		n_iter = kinFit.getNbIter();
		fitStatus->Fill(fit_status);
		nIter_all->Fill(n_iter);
		/* skip if the fit did not converge */
		if (fit_status != 0)  // Status: -1: "NO FIT PERFORMED", 10: "RUNNING", 0: "CONVERGED", 1: "NOT CONVERGED", -10: "ABORTED"; should only return 0 or 1
			return;

		fitLepton1 = (*e1.getCurr4Vec());
		fitLepton2 = (*e2.getCurr4Vec());
		fitPhoton = (*ph.getCurr4Vec());
		fitProton = (*pr.getCurr4Vec());

		/*for (outer_it step = temp_results.begin(); step != temp_results.end(); ++step)
			for (inner_it it = step->begin(); it != step->end(); ++it)
				it->Print();*/

		etap_mass_smeared = MeV ? (*e1.getIni4Vec() + *e2.getIni4Vec() + *ph.getIni4Vec()).M() : (*e1.getIni4Vec() + *e2.getIni4Vec() + *ph.getIni4Vec()).M()*1000.;
		etap_smeared->Fill(etap_mass_smeared);
		coplanarity_smeared = abs((*e1.getIni4Vec() + *e2.getIni4Vec() + *ph.getIni4Vec()).Phi() - (*pr.getIni4Vec()).Phi())*R2D;
		copl_smeared->Fill(coplanarity_smeared);
		int iterStep = 0;
		coplanarity->Fill(coplanarity_smeared, iterStep);
		etap_fit_steps->Fill(etap_mass_smeared, iterStep++);
		for (outer_it step = temp_results.begin(); step != temp_results.end(); ++step) {
			coplanarity_fitted = abs((step->at(0) + step->at(1) + step->at(2)).Phi() - step->at(3).Phi())*R2D;
			coplanarity->Fill(coplanarity_fitted, iterStep);
			etap_mass_fitted = MeV ? (step->at(0) + step->at(1) + step->at(2)).M() : (step->at(0) + step->at(1) + step->at(2)).M()*1000.;
			etap_fit_steps->Fill(etap_mass_fitted, iterStep++);
		}

		// get the pulls from the fit
		TMatrixD pullsLepton1 = *(e1.getPull());
		TMatrixD pullsLepton2 = *(e2.getPull());
		TMatrixD pullsPhoton = *(ph.getPull());
		TMatrixD pullsProton = *(pr.getPull());

		lepton1_pullE = pullsLepton1(0,0);
		lepton1_pullTheta = pullsLepton1(1,0);
		lepton1_pullPhi = pullsLepton1(2,0);
		lepton2_pullE = pullsLepton2(0,0);
		lepton2_pullTheta = pullsLepton2(1,0);
		lepton2_pullPhi = pullsLepton2(2,0);
		photon_pullE = pullsPhoton(0,0);
		photon_pullTheta = pullsPhoton(1,0);
		photon_pullPhi = pullsPhoton(2,0);
		proton_pullE = pullsProton(0,0);
		proton_pullTheta = pullsProton(1,0);
		proton_pullPhi = pullsProton(2,0);

		/*printf("Photon1: pullE: %f, pullTheta: %f, pullPhi: %f\n", photon1_pullE, photon1_pullTheta, photon1_pullPhi);
		printf("Photon2: pullE: %f, pullTheta: %f, pullPhi: %f\n", photon2_pullE, photon2_pullTheta, photon2_pullPhi);
		printf("Proton: pullE: %f, pullTheta: %f, pullPhi: %f\n", proton_pullE, proton_pullTheta, proton_pullPhi);*/

		// skip plotting pulls which are not a number
		if (isfinite(lepton1_pullE))
			lepton1PullE->Fill(lepton1_pullE);
		if (isfinite(lepton1_pullTheta))
			lepton1PullTheta->Fill(lepton1_pullTheta);
		if (isfinite(lepton1_pullPhi))
			lepton1PullPhi->Fill(lepton1_pullPhi);
		if (isfinite(lepton2_pullE))
			lepton2PullE->Fill(lepton2_pullE);
		if (isfinite(lepton2_pullTheta))
			lepton2PullTheta->Fill(lepton2_pullTheta);
		if (isfinite(lepton2_pullPhi))
			lepton2PullPhi->Fill(lepton2_pullPhi);
		if (isfinite(photon_pullE))
			photonPullE->Fill(photon_pullE);
		if (isfinite(photon_pullTheta))
			photonPullTheta->Fill(photon_pullTheta);
		if (isfinite(photon_pullPhi))
			photonPullPhi->Fill(photon_pullPhi);
		if (isfinite(proton_pullE))
			protonPullE->Fill(proton_pullE);
		if (isfinite(proton_pullTheta))
			protonPullTheta->Fill(proton_pullTheta);
		if (isfinite(proton_pullPhi))
			protonPullPhi->Fill(proton_pullPhi);

		//std::cout << "Fit result: " << std::endl;
		//kinFit.print();
		chisq = kinFit.getChi2();
		ndf = kinFit.getNdf();
		prob = kinFit.getProbability();
		chisq_hist->Fill(chisq);
		prob_hist->Fill(prob);
		nIter_converged->Fill(n_iter);
		etap_fitted->Fill(etap_mass_fitted);
		//coplanarity_fitted = abs((fitPhoton1 + fitPhoton2).Phi() - fitProton.Phi())*R2D;
		copl_fitted->Fill(coplanarity_fitted);
		double q2_fit = (fitLepton1 + fitLepton2).M();
		q2->Fill(MeV ? q2_fit : q2_fit*1000.);
		if (prob > .1) {  // perform a cut on the probability
			// eta'
			missM_etap_vs_pTheta->Fill(MeV ? (beam + target - fitProton).M() : (beam + target - fitProton).M()*1000., fitProton.Theta()*R2D);
			etap_invM_gg_vs_beamE->Fill(etap_mass_fitted, MeV ? beam.E() : beam.E()*1000.);
			etap_p_kinE_vs_pTheta_true->Fill((fitProton.E() - fitProton.M())*1000., particles[charged[0]].Theta()*R2D);
			etap_invM_gg_smeared_vs_invM_gg_fit->Fill(etap_mass_smeared, etap_mass_fitted);
			q2_cut->Fill(MeV ? q2_fit : q2_fit*1000.);
		}

		//delete temp_results;

	}

}

void SaschaPhysics::GetParticles(particle_vector* particles, particle_t* beam)
{static int count = 0;
	GetParticles(particles);
	try {
		GetBeam(beam);
	} catch (beam_photon_mismatch_exception& e) {
		//cout << "Caught beam mismatch exception " << ++count << endl;
		//GetTrueBeam(beam);
		beam = nullptr;
	}
}
//        DETECTOR_NONE = 0,
//        DETECTOR_NaI = 1,
//        DETECTOR_PID = 2,
//        DETECTOR_MWPC = 4,
//        DETECTOR_BaF2 = 8,
//        DETECTOR_PbWO4 = 16,
//        DETECTOR_Veto = 32,
void SaschaPhysics::GetParticles(particle_vector* particles)
{  // example of how to collect particles in a user defined way
	GTreeTrack* tracks = GetTracks();
	for (unsigned int i = 0; i < tracks->GetNTracks(); i++) {
		if (tracks->GetDetectors(i) & GTreeTrack::DETECTOR_NaI) {
			nParticlesCB++;
			if (tracks->GetVetoEnergy(i) > 0.)  // PID entry? --> charged
				particles->push_back(particle_t(tracks->GetVector(i, MASS_ELECTRON), electron));
			else
				particles->push_back(particle_t(tracks->GetVector(i), photon));
		} else {  //TAPS: DETECTOR_BaF2, DETECTOR_PbWO4
			nParticlesTAPS++;
			if (tracks->GetVetoEnergy(i) > 0.)  // Veto entry? --> charged
				particles->push_back(particle_t(tracks->GetVector(i, MASS_PROTON), proton));
			else
				particles->push_back(particle_t(tracks->GetVector(i), photon));
		}
	}
	// sort out particles with an energy of less than 10 MeV
	for (particle_it it = particles->begin(); it != particles->end();)
		if (it->E() < 10.)
			it = particles->erase(it);
		else
			++it;
}

void SaschaPhysics::GetBeam(particle_t* beam)
{
	if (GetTagger()->GetNTagged() != 1)
		throw beamMismatch;

	beam->p4 = GetTagger()->GetVector(0);
}

void SaschaPhysics::GetTrueParticles(particle_vector* particles, particle_t* beam)
{
	GetTrueParticles(particles);
	GetTrueBeam(beam);
}

void SaschaPhysics::GetTrueParticles(particle_vector* particles)
{
	if (pluto_tree) {  // in case we have a Pluto tree
		for (const PParticle* p : pluto->GetFinalState())
			particles->push_back(particle_t(p->Vect4(), static_cast<particle_id>(p->ID())));
		//beam->p4 = pluto->GetTrueBeam();
	} else {  // if not, use Geant tree information
		//TODO: maybe do this directly in the GTreeA2Geant class?
		for (unsigned int i = 0; i < geant->GetNTrueParticles(); i++)
			particles->push_back(particle_t(geant->GetTrueVector(i), static_cast<particle_id>(geant->GetTrueID(i))));
		for (particle_it it = particles->begin(); it != particles->end();)
			if (!it->is_final_state())
				it = particles->erase(it);
			else
				++it;
		//beam->p4 = geant->GetBeam();
	}
}

void SaschaPhysics::GetTrueBeam(particle_t* beam)
{
	if (pluto_tree)
		beam->p4 = pluto->GetTrueBeam();
	else
		beam->p4 = geant->GetBeam();
}

void SaschaPhysics::MatchMC(particle_vector* particles, particle_vector* trueParts, match* matches)
{
	double angle, tmp;
	double maxDeltaE = 1.;
	for (unsigned int i = 0; i < trueParts->size(); i++) {
		angle = R2D;
		unsigned int pos;
		// try to find matching particles by their angle
		for (unsigned int j = 0; j < particles->size(); j++)
			if ((tmp = trueParts->at(i).p4.Angle(particles->at(j).Vect())) < angle
					&& abs(trueParts->at(i).E() - particles->at(j).E())/trueParts->at(i).E() < maxDeltaE ) {
				angle = tmp;
				pos = j;
			}
		matches->insert(match_pair(i, pos));
	}
}

void SaschaPhysics::ProcessScalerRead()
{
}

Bool_t SaschaPhysics::Write()
{
	// Write the fit related TH1s
	GTreeManager::Write(n_particles);
	GTreeManager::Write(invM_eeg);
	GTreeManager::Write(invM_2e);
	GTreeManager::Write(missM);
	GTreeManager::Write(expectProt);

	GTreeManager::Write(balanceE);
	GTreeManager::Write(balancePx);
	GTreeManager::Write(balancePy);
	GTreeManager::Write(balancePz);

	GTreeManager::Write(MCbeamDeltaE);
	GTreeManager::Write(MCbeamDeltaP);
	GTreeManager::Write(MClepton1DeltaE);
	GTreeManager::Write(MClepton1DeltaP);
	GTreeManager::Write(MClepton2DeltaE);
	GTreeManager::Write(MClepton2DeltaP);
	GTreeManager::Write(MCphotonDeltaE);
	GTreeManager::Write(MCphotonDeltaP);
	GTreeManager::Write(MCprotonDeltaE);
	GTreeManager::Write(MCprotonDeltaP);

	GTreeManager::Write(etap_fit_steps);
	GTreeManager::Write(coplanarity);

	GTreeManager::Write(missM_etap_vs_pTheta);
	GTreeManager::Write(etap_invM_gg_vs_beamE);
	GTreeManager::Write(etap_p_kinE_vs_pTheta_true);
	GTreeManager::Write(etap_invM_gg_smeared_vs_invM_gg_fit);

	GTreeManager::Write(lepton1PullE);
	GTreeManager::Write(lepton1PullTheta);
	GTreeManager::Write(lepton1PullPhi);
	GTreeManager::Write(lepton2PullE);
	GTreeManager::Write(lepton2PullTheta);
	GTreeManager::Write(lepton2PullPhi);
	GTreeManager::Write(photonPullE);
	GTreeManager::Write(photonPullTheta);
	GTreeManager::Write(photonPullPhi);
	GTreeManager::Write(protonPullE);
	GTreeManager::Write(protonPullTheta);
	GTreeManager::Write(protonPullPhi);

	GTreeManager::Write(chisq_hist);
	GTreeManager::Write(prob_hist);
	GTreeManager::Write(etap_smeared);
	GTreeManager::Write(etap_fitted);
	GTreeManager::Write(copl_smeared);
	GTreeManager::Write(copl_fitted);
	GTreeManager::Write(nIter_converged);
	GTreeManager::Write(nIter_all);
	GTreeManager::Write(fitStatus);

	GTreeManager::Write(q2);
	GTreeManager::Write(q2_cut);

	// Write all GH1's easily
	GTreeManager::Write();
}

double SaschaPhysics::sigma_E(const TLorentzVector* const p)
{
	double E = p->E() - p->M();
	if (p->Theta() < .349)  // TAPS
		return (.018 + .008*sqrt(E))*E;
	else  // CB
		return pow(E, .36)*.02*E;
}

double SaschaPhysics::sigma_theta(const TLorentzVector* const p)
{
	if (p->Theta() < .349)  // TAPS
		return .017;
	else  // CB
		return .044;
}

double SaschaPhysics::sigma_phi(const TLorentzVector* const p)
{
	if (p->Theta() < .349)  // TAPS
		return .017;
	else  // CB
		return .044/sin(p->Theta());
}
