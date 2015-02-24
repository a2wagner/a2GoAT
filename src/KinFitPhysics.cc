#include "KinFitPhysics.h"

KinFitPhysics::KinFitPhysics()
{
	time 	= new GH1("time", 	"time", 	1400, -700, 700);
	time_cut 	= new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

	time_2g 	= new GH1("time_2g",	"time_2g", 	1400, -700, 700);
	time_2g_cut = new GH1("time_2g_cut","time_2g_cut", 	1400, -700, 700);

	IM 		= new GH1("IM", 	"IM", 		400,   0, 400);
	IM_2g 	= new GH1("IM_2g", 	"IM_2g", 	400,   0, 400);

	MM		= new GH1("MM", 	"MM", 	 	400,   800, 1200);
	MM_2g	= new GH1("MM_2g", 	"MM_2g", 	400,   800, 1200);

	TaggerAccScal = new TH1D("TaggerAccScal","TaggerAccScal",352,0,352);

	// own histograms
	n_particles = new TH1I("n_part", "n_particles", 11, 0, 10);
	invM_2g = new TH1D("invM_2g", "inv_mass_2g", 200, 100, 200);
	missM = new TH1D("missM", "miss_mass", 600, 500, 1100);
	expectProt = new TH1D("expectProt", "expected_proton", 500, 700, 1200);

	pi0_fit_steps = new TH2D("pi0_fit_steps", "pi0_fit_steps", 600, 0, 300, 50, 0, 49);
	coplanarity = new TH2D("coplanarity", "coplanarity", 720, 0, 360, 50, 0, 49);

	missM_pi0_vs_pTheta = new TH2D("missM_pi0_vs_pTheta", "missM_pi0_vs_pTheta", 600, 0, 300, 240, 0, 120);
	invM_gg_vs_beamE = new TH2D("invM_gg_vs_beamE", "invM_gg_vs_beamE", 600, 0, 300, 360, 1440, 1620);
	p_kinE_vs_pTheta_true = new TH2D("p_kinE_vs_pTheta_true", "p_kinE_vs_pTheta_true", 24000, 0, 1200, 240, 0, 120);
	invM_gg_smeared_vs_invM_gg_fit = new TH2D("invM_gg_smeared_vs_invM_gg_fit", "invM_gg_smeared_vs_invM_gg_fit", 600, 0, 300, 600, 0, 300);

	photon1PullE = new TH1D("photon1PullE", "photon1_pullE", 1000, -5, 5);
	photon1PullTheta = new TH1D("photon1PullTheta", "photon1_pullTheta", 1000, -5, 5);
	photon1PullPhi = new TH1D("photon1PullPhi", "photon1_pullPhi", 1000, -5, 5);
	photon2PullE = new TH1D("photon2PullE", "photon2_pullE", 1000, -5, 5);
	photon2PullTheta = new TH1D("photon2PullTheta", "photon2_pullTheta", 1000, -5, 5);
	photon2PullPhi = new TH1D("photon2PullPhi", "photon2_pullPhi", 1000, -5, 5);
	protonPullE = new TH1D("protonPullE", "proton_pullE", 1000, -5, 5);
	protonPullTheta = new TH1D("protonPullTheta", "proton_pullTheta", 1000, -5, 5);
	protonPullPhi = new TH1D("protonPullPhi", "proton_pullPhi", 1000, -5, 5);

	chisq_hist = new TH1D("chisq", "chi_square", 1000, 0, 100);
	prob_hist = new TH1D("prob", "probability", 1000, 0, 1);
	pi0_smeared = new TH1D("pi0_smeared", "pi0_mass_smeared", 600, 0, 300);
	pi0_fitted = new TH1D("pi0_fitted", "pi0_mass_fitted", 600, 0, 300);
	copl_smeared = new TH1D("copl_smeared", "coplanarity_smeared", 720, 0, 360);
	copl_fitted = new TH1D("copl_fitted", "coplanarity_fitted", 720, 0, 360);
	nIter_converged = new TH1I("n_iter_converged", "n_iter_converged", 52, 0, 51);
	nIter_all = new TH1I("n_iter_all", "n_iter_all", 52, 0, 51);
	fitStatus = new TH1I("fit_status", "fit_status", 21, -10, 10);

	photon1P_Dsmeared = new TH1D("photon1P_Dsmeared", "photon1P_smeared-true", 2000, -1000, 1000);
	photon2P_Dsmeared = new TH1D("photon2P_Dsmeared", "photon2P_smeared-true", 2000, -1000, 1000);
	protonP_Dsmeared = new TH1D("protonP_Dsmeared", "protonP_smeared-true", 2000, -1000, 1000);
	photon1P_Dfitted = new TH1D("photon1P_Dfitted", "photon1P_fitted-true", 2000, -1000, 1000);
	photon2P_Dfitted = new TH1D("photon2P_Dfitted", "photon2P_fitted-true", 2000, -1000, 1000);
	protonP_Dfitted = new TH1D("protonP_Dfitted", "protonP_fitted-true", 2000, -1000, 1000);
}

KinFitPhysics::~KinFitPhysics()
{
}

Bool_t KinFitPhysics::Init(const char* configfile)
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

Bool_t KinFitPhysics::Start()
{
	if(!IsGoATFile())
	{
		cout << "ERROR: Input File is not a GoAT file." << endl;
		return kFALSE;
	}
	SetAsPhysicsFile();


	/* first get the MC trees and check if they really exist because you'll still get vectors etc. from non-existent trees! */
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


	TraverseValidEvents();

	return kTRUE;
}

void KinFitPhysics::ProcessEvent()
{
	// fill time diff (tagger - pi0), all pi0
	FillTime(*GetNeutralPions(),time);
	FillTimeCut(*GetNeutralPions(),time_cut);

	// fill missing mass, all pi0
	FillMissingMass(*GetNeutralPions(),MM);

	// fill invariant mass, all pi0
	FillMass(*GetNeutralPions(),IM);

	// Some neutral decays
	for (Int_t i = 0; i < GetNeutralPions()->GetNParticles(); i++)
	{
		// Fill MM for 2 photon decay
		if ((GetNeutralPions()->GetNSubParticles(i) == 2) & (GetNeutralPions()->GetNSubPhotons(i) == 2))
		{
		// fill time diff (tagger - pi0), this pi0
		FillTime(*GetNeutralPions(),i,time_2g);
		FillTimeCut(*GetNeutralPions(),i,time_2g_cut);

		// fill missing mass, this pi0
				FillMissingMass(*GetNeutralPions(),i,MM_2g);

		// fill invariant mass, this pi0
			FillMass(*GetNeutralPions(),i,IM_2g);
		}

	}


	/* declare some variables used for kinematic fitting */
	bool dbg = false;
	static int count = 0;
	double photon1_pullE, photon1_pullTheta, photon1_pullPhi,
			photon2_pullE, photon2_pullTheta, photon2_pullPhi,
			proton_pullE, proton_pullTheta, proton_pullPhi;
	double chisq, prob, pi0_mass_smeared, pi0_mass_fitted,
			coplanarity_smeared, coplanarity_fitted;
	int n_iter, fit_status;


	/* Gather all particle information, depending on which MC trees are available */
	bool MeV = false;  // use MeV or GeV?
	//TODO fit converges significantly slower (3 times or even more iteration steps) for MeV ?
	nParticles = 0, nParticlesCB = 0, nParticlesTAPS = 0;
	particle_vector particles;
	particle_t trueBeam(photon);
	GetTrueParticles(&particles, &trueBeam);
	particle_t trueTarget(TLorentzVector(0., 0., 0., MeV ? MASS_PROTON : MASS_PROTON/1000.), proton);
	if (MeV)
		trueBeam.GeV2MeV();
	nParticles = particles.size();

	if (dbg) {
		printf("\nEvent %d:\n", ++count);
		for (auto& p : particles) {
			printf("Particle id: %2d (%s), ", p.id, p.is_charged() ? "charged" : "uncharged");
			p.Print();
		}
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

		int i = 0;
		for (auto& p : particles) {
			if (MeV)
				p.GeV2MeV();
			if (p.is_charged())
				charged[nCharged++] = i;
			else {
				neutral[nNeutral++] = i;
				tmpState += particles[i];
			}
			i++;
		}

		/* end event selection when we don't have one charged particle (entry in PID or Veto) and two neutral ones */
		if (nCharged != 1 && nNeutral != 2)
			return;

		double protEnergyReconstr = particles[charged[0]].E() - particles[charged[0]].M();
		double invM_2neutral = MeV ? tmpState.M() : tmpState.M()*1000.;
		invM_2g->Fill(invM_2neutral);

		/* At this point we have one charged and two neutral particles. Start to examine all tagged photons for prompt and random windows. */
		//for (unsigned int i = 0; i < nBeamPhotons; i++) {
			//missProton = fP4target[0] + Tagged[i].GetP4() - tmpState;
			missProton = trueTarget + trueBeam - tmpState;
			mMiss = missProton.M();
			protExpect = missProton.E() - missProton.M();
			missM->Fill(MeV ? mMiss : mMiss*1000.);
			expectProt->Fill(MeV ? protExpect : protExpect*1000.);
		//}

		/* prepare kinematic fit */
		int err;
		int ndf;

		kinFit.reset();

		TVector3 photon1 = particles[neutral[0]].Vect();
		TVector3 photon2 = particles[neutral[1]].Vect();
		TVector3 proton = particles[charged[0]].Vect();

		// temp checks
		/*TLorentzVector balance = particles[neutral[0]].p4 + particles[neutral[1]].p4 + particles[charged[0]].p4 - trueBeam.p4 - trueTarget.p4;
		printf("Proton: %f, Pi0: %f, BeamE: %f, BeamM: %f\n", particles[charged[0]].M(), (particles[neutral[0]].p4 + particles[neutral[1]].p4).M(), trueBeam.E(), trueBeam.M());
		printf("Energy-Balance: %f, Momentum: (%f, %f, %f)\n", balance.E(), balance.Px(), balance.Py(), balance.Pz());
		printf("MissM: %f, invM_2g: %f\n", mMiss, invM_2neutral);
		printf("Ekin expected Proton: %f, true: %f\n", protExpect, particles[charged[0]].E() - particles[charged[0]].M());
		printf("ProtE = %f, ProtP = %f\n", particles[charged[0]].E(), particles[charged[0]].p4.P());
		printf("Beam: "); trueBeam.Print();
		printf("Photon1: "); particles[neutral[0]].Print();
		printf("Photon2: "); particles[neutral[1]].Print();
		printf("Proton: "); particles[charged[0]].Print();
		printf("MASS_PROTON = %f; MASS_PIZERO = %f\n", MASS_PROTON, MASS_PIZERO);*/

		TLorentzVector fitPhoton1;
		TLorentzVector fitPhoton2;
		TLorentzVector fitProton;

		/**
		 * Use absolute or relative smearing?
		 * For absolute smearing use the sigmaX values defined above and comment the ones for each particle
		 * Relative smearing is calculated for each particle using the following data:
		 *   CB: sigma_theta = 2.5°; sigma_phi = sigma_theta/sin(theta); sigma_E = E^.36*.02*E
		 *   TAPS: sigma_theta/phi = 1°; sigma_E = (.018+.008*sqrt(E))*E
		 */
		TMatrixD covPhoton1;
		TMatrixD covPhoton2;
		TMatrixD covProton;

		TRandom3 rand(0);
		unsigned int nPar = 3;
		double sigmaP, sigmaTheta, sigmaPhi;
		int rows = nPar;  // number of rows equal to number of cols
		double errors[nPar];
		double factor = 1.;//1.2;  // factor to increase/decrease sigma^2 values used for kinematic fitting
		// absolute smearing
		/*sigmaTheta = .05, sigmaPhi = .05;
		if (MeV)
			sigmaP = 20.;
		else
			sigmaP = .02;
		Double_t errors[] = {sigmaP*sigmaP*factor, sigmaTheta*sigmaTheta*factor, sigmaPhi*sigmaPhi*factor};*/
		// relative smearing; assume sigmaP ~ sigmaE
		/* 1st photon */
		int currPart = neutral[0];
		sigmaTheta = sigma_theta(&particles[currPart].p4);
		sigmaPhi = sigma_phi(&particles[currPart].p4);
		sigmaP = sigma_E(&particles[currPart].p4);
		photon1.SetMagThetaPhi(rand.Gaus(photon1.Mag(), sigmaP), rand.Gaus(photon1.Theta(), sigmaTheta), rand.Gaus(photon1.Phi(), sigmaPhi));
		kinFit.fill_array(errors, {sigmaP*sigmaP*factor, sigmaTheta*sigmaTheta*factor, sigmaPhi*sigmaPhi*factor});
		//Double_t errors[] = {1., .01, .01};
		//kinFit.sigmaEThetaPhi(particles[currPart], errors);
		if (kinFit.fillSquareMatrixDiagonal(&covPhoton1, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		/* 2nd photon */
		currPart = neutral[1];
		sigmaTheta = sigma_theta(&particles[currPart].p4);
		sigmaPhi = sigma_phi(&particles[currPart].p4);
		sigmaP = sigma_E(&particles[currPart].p4);
		photon2.SetMagThetaPhi(rand.Gaus(photon2.Mag(), sigmaP), rand.Gaus(photon2.Theta(), sigmaTheta), rand.Gaus(photon2.Phi(), sigmaPhi));
		kinFit.fill_array(errors, {sigmaP*sigmaP*factor, sigmaTheta*sigmaTheta*factor, sigmaPhi*sigmaPhi*factor});
		//kinFit.sigmaEThetaPhi(particles[currPart], errors);
		if (kinFit.fillSquareMatrixDiagonal(&covPhoton2, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		/* proton */
		currPart = charged[0];
		sigmaTheta = sigma_theta(&particles[currPart].p4);
		sigmaPhi = sigma_phi(&particles[currPart].p4);
		sigmaP = sigma_E(&particles[currPart].p4);
		proton.SetMagThetaPhi(rand.Gaus(proton.Mag(), sigmaP), rand.Gaus(proton.Theta(), sigmaTheta), rand.Gaus(proton.Phi(), sigmaPhi));
		kinFit.fill_array(errors, {sigmaP*sigmaP*factor, sigmaTheta*sigmaTheta*factor, sigmaPhi*sigmaPhi*factor});
		//kinFit.sigmaEThetaPhi(particles[currPart], errors);
		if (kinFit.fillSquareMatrixDiagonal(&covProton, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");

//		//printf("Particle 1: sigmaE = %f, sigmaPhi = %f, sigmaTheta = %f\n", particles[0].GetSigmaE(), particles[0].GetSigmaPhi(), particles[0].GetSigmaTheta());
//		particles[0].Print();
//		covPhoton1.Print();
//		covPhoton2.Print();
//		covProton.Print();

		TFitParticlePThetaPhi ph1("neutral1", "neutral1", &photon1, 0., &covPhoton1);
		TFitParticlePThetaPhi ph2("neutral2", "neutral2", &photon2, 0., &covPhoton2);
		TFitParticlePThetaPhi pr("charged1", "charged1", &proton, MeV ? MASS_PROTON : MASS_PROTON/1000., &covProton);

		TVector3 beam = trueBeam.Vect();
		TVector3 target = trueTarget.Vect();
		TMatrixD covBeam;
		TMatrixD covTarget;
		//kinFit.sigmaEThetaPhi(Tagged[0], errors);
		errors[0] = .0001; errors[1] = .0001; errors[2] = .0001;
		if (MeV)
			errors[0] = 10.;
		if (kinFit.fillSquareMatrixDiagonal(&covBeam, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		covTarget.Zero();
		covTarget.ResizeTo(3, 3);
		TFitParticlePThetaPhi bm("beam", "beam", &beam, 0., &covBeam);
		TFitParticlePThetaPhi trgt("target", "target", &target, MeV ? MASS_PROTON : MASS_PROTON/1000., &covTarget);

		// energy and momentum constraints have to be defined separately for each component
		// components can be accessed via enum TFitConstraintEp::component
		Double_t constraint = trueBeam.E() + trueTarget.E();
		TFitConstraintEp energyConservation("energyConstr", "Energy conservation constraint", 0, TFitConstraintEp::E, constraint);
		energyConservation.addParticles1(&ph1, &ph2, &pr);
		constraint = trueBeam.Px();
		TFitConstraintEp pxConservation("pxConstr", "Px conservation constraint", 0, TFitConstraintEp::pX, constraint);
		pxConservation.addParticles1(&ph1, &ph2, &pr);
		constraint = trueBeam.Py();
		TFitConstraintEp pyConservation("pyConstr", "Py conservation constraint", 0, TFitConstraintEp::pY, constraint);
		pyConservation.addParticles1(&ph1, &ph2, &pr);
		constraint = trueBeam.Pz();
		TFitConstraintEp pzConservation("pzConstr", "Pz conservation constraint", 0, TFitConstraintEp::pZ, constraint);
		pzConservation.addParticles1(&ph1, &ph2, &pr);
		TFitConstraintM massConstrProton("massConstr_proton", "mass constraint proton", 0, 0, MeV ? MASS_PROTON : MASS_PROTON/1000.);
		massConstrProton.addParticle1(&pr);
		// constraint for pi0 mass
		TFitConstraintM massConstrPi0("massConstr_pi0", "mass constraint pi0", 0, 0, MeV ? MASS_PIZERO : MASS_PIZERO/1000.);
		massConstrPi0.addParticles1(&ph1, &ph2);

		//TKinFitter fit;
		kinFit.addMeasParticle(&ph1);
		kinFit.addMeasParticle(&ph2);
		kinFit.addMeasParticle(&pr);
		//kinFit.setParamUnmeas(&pr, 0);  // proton energy unmeasured
		//kinFit.addMeasParticle(&bm);
		//kinFit.addMeasParticle(&trgt);

		//kinFit.addConstraint(&massConstrProton);
		kinFit.addConstraint(&energyConservation);
		kinFit.addConstraint(&pxConservation);
		kinFit.addConstraint(&pyConservation);
		kinFit.addConstraint(&pzConservation);
		//kinFit.addConstraint(&massConstrPi0);

		// get the intermediate steps from the fitter
		//vector<TAbsFitParticle> *temp_results = new vector<TAbsFitParticle>;
		vvP4 temp_results;
		/*vector<TAbsFitParticle*> track_particles;
		track_particles.push_back(&ph1);
		track_particles.push_back(&ph2);
		track_particles.push_back(&pr);*/

		//kinFit.get_intermediate_steps(&temp_results, track_particles);
		kinFit.get_intermediate_steps(&temp_results, {&ph1, &ph2, &pr});

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

		fitPhoton1 = (*ph1.getCurr4Vec());
		fitPhoton2 = (*ph2.getCurr4Vec());
		fitProton = (*pr.getCurr4Vec());

		photon1P_Dsmeared->Fill(((*ph1.getIni4Vec()).P() - particles[neutral[0]].P())*1000.);
		photon2P_Dsmeared->Fill(((*ph2.getIni4Vec()).P() - particles[neutral[1]].P())*1000.);
		protonP_Dsmeared->Fill(((*pr.getIni4Vec()).P() - particles[charged[0]].P())*1000.);
		photon1P_Dfitted->Fill((fitPhoton1.P() - particles[neutral[0]].P())*1000.);
		photon2P_Dfitted->Fill((fitPhoton2.P() - particles[neutral[1]].P())*1000.);
		protonP_Dfitted->Fill((fitProton.P() - particles[charged[0]].P())*1000.);

		/*for (outer_it step = temp_results.begin(); step != temp_results.end(); ++step)
			for (inner_it it = step->begin(); it != step->end(); ++it)
				it->Print();*/

		pi0_mass_smeared = MeV ? (*ph1.getIni4Vec() + *ph2.getIni4Vec()).M() : (*ph1.getIni4Vec() + *ph2.getIni4Vec()).M()*1000.;
		pi0_smeared->Fill(pi0_mass_smeared);
		coplanarity_smeared = abs((*ph1.getIni4Vec() + *ph2.getIni4Vec()).Phi() - (*pr.getIni4Vec()).Phi())*R2D;
		copl_smeared->Fill(coplanarity_smeared);
		int iterStep = 0;
		coplanarity->Fill(coplanarity_smeared, iterStep);
		pi0_fit_steps->Fill(pi0_mass_smeared, iterStep++);
		for (outer_it step = temp_results.begin(); step != temp_results.end(); ++step) {
			coplanarity_fitted = abs((step->at(0) + step->at(1)).Phi() - step->at(2).Phi())*R2D;
			coplanarity->Fill(coplanarity_fitted, iterStep);
			pi0_mass_fitted = MeV ? (step->at(0) + step->at(1)).M() : (step->at(0) + step->at(1)).M()*1000.;
			pi0_fit_steps->Fill(pi0_mass_fitted, iterStep++);
		}

		// get the pulls from the fit
		TMatrixD pullsPhoton1 = *(ph1.getPull());
		TMatrixD pullsPhoton2 = *(ph2.getPull());
		TMatrixD pullsProton = *(pr.getPull());

		photon1_pullE = pullsPhoton1(0,0);
		photon1_pullTheta = pullsPhoton1(1,0);
		photon1_pullPhi = pullsPhoton1(2,0);
		photon2_pullE = pullsPhoton2(0,0);
		photon2_pullTheta = pullsPhoton2(1,0);
		photon2_pullPhi = pullsPhoton2(2,0);
		proton_pullE = pullsProton(0,0);
		proton_pullTheta = pullsProton(1,0);
		proton_pullPhi = pullsProton(2,0);

//		printf("Photon1: pullE: %f, pullTheta: %f, pullPhi: %f\n", photon1_pullE, photon1_pullTheta, photon1_pullPhi);
//		printf("Photon2: pullE: %f, pullTheta: %f, pullPhi: %f\n", photon2_pullE, photon2_pullTheta, photon2_pullPhi);
//		printf("Proton: pullE: %f, pullTheta: %f, pullPhi: %f\n", proton_pullE, proton_pullTheta, proton_pullPhi);

		// skip plotting pulls which are not a number
		if (isfinite(photon1_pullE))
			photon1PullE->Fill(photon1_pullE);
		if (isfinite(photon1_pullTheta))
			photon1PullTheta->Fill(photon1_pullTheta);
		if (isfinite(photon1_pullPhi))
			photon1PullPhi->Fill(photon1_pullPhi);
		if (isfinite(photon2_pullE))
			photon2PullE->Fill(photon2_pullE);
		if (isfinite(photon2_pullTheta))
			photon2PullTheta->Fill(photon2_pullTheta);
		if (isfinite(photon2_pullPhi))
			photon2PullPhi->Fill(photon2_pullPhi);
		if (isfinite(proton_pullE))
			protonPullE->Fill(proton_pullE);
		if (isfinite(proton_pullTheta))
			protonPullTheta->Fill(proton_pullTheta);
		if (isfinite(proton_pullPhi))
			protonPullPhi->Fill(proton_pullPhi);

		//std::cout << "Fit result: " << std::endl;
		//kinFit.print();
		ndf = kinFit.getNDF();
		chisq = kinFit.getS();
		//prob = kinFit.getProbability();
		prob = TMath::Prob(chisq, ndf);
		chisq_hist->Fill(chisq);
		prob_hist->Fill(prob);
		nIter_converged->Fill(n_iter);
		//std::cout << "\nProbability: " << prob << "\tchi^2: " << chisq << std::endl;
		//pi0_mass_fitted = (fitPhoton1 + fitPhoton2).M();
		pi0_fitted->Fill(pi0_mass_fitted);
		//coplanarity_fitted = abs((fitPhoton1 + fitPhoton2).Phi() - fitProton.Phi())*R2D;
		copl_fitted->Fill(coplanarity_fitted);
		if (prob > .1) {
			missM_pi0_vs_pTheta->Fill(MeV ? (trueBeam + trueTarget - fitProton).M() : (trueBeam + trueTarget - fitProton).M()*1000., fitProton.Theta()*R2D);
			invM_gg_vs_beamE->Fill(pi0_mass_fitted, MeV ? trueBeam.E() : trueBeam.E()*1000.);
			p_kinE_vs_pTheta_true->Fill((fitProton.E() - fitProton.M())*1000., particles[charged[0]].Theta()*R2D);
			invM_gg_smeared_vs_invM_gg_fit->Fill(pi0_mass_smeared, pi0_mass_fitted);
		}
//if(count++==5)exit(0);
		//delete temp_results;

	}

}

void KinFitPhysics::GetParticles()
{
	//TODO
}

void KinFitPhysics::GetTrueParticles(particle_vector* particles, particle_t* beam)
{
	if (pluto_tree) {  // in case we have a Pluto tree
		for (const PParticle* p : pluto->GetFinalState())
			particles->push_back(particle_t(p->Vect4(), static_cast<particle_id>(p->ID())));
		beam->p4 = pluto->GetTrueBeam();
	} else {  // if not, use Geant tree information
		//TODO: maybe do this directly in the GTreeA2Geant?
		for (unsigned int i = 0; i < geant->GetNTrueParticles(); i++)
			particles->push_back(particle_t(geant->GetTrueVector(i), static_cast<particle_id>(geant->GetTrueID(i))));
		for (particle_it it = particles->begin(); it != particles->end();)
			if (!it->is_final_state())
				it = particles->erase(it);
			else
				++it;
		beam->p4 = geant->GetBeam();
	}
}

void KinFitPhysics::ProcessScalerRead()
{
	// Fill Tagger Scalers
	FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t KinFitPhysics::Write()
{
	// Write some TH1s
	GTreeManager::Write(TaggerAccScal);

	// Write the fit related TH1s
	GTreeManager::Write(n_particles);
	GTreeManager::Write(invM_2g);
	GTreeManager::Write(missM);
	GTreeManager::Write(expectProt);

	GTreeManager::Write(pi0_fit_steps);
	GTreeManager::Write(coplanarity);

	GTreeManager::Write(missM_pi0_vs_pTheta);
	GTreeManager::Write(invM_gg_vs_beamE);
	GTreeManager::Write(p_kinE_vs_pTheta_true);
	GTreeManager::Write(invM_gg_smeared_vs_invM_gg_fit);

	GTreeManager::Write(photon1PullE);
	GTreeManager::Write(photon1PullTheta);
	GTreeManager::Write(photon1PullPhi);
	GTreeManager::Write(photon2PullE);
	GTreeManager::Write(photon2PullTheta);
	GTreeManager::Write(photon2PullPhi);
	GTreeManager::Write(protonPullE);
	GTreeManager::Write(protonPullTheta);
	GTreeManager::Write(protonPullPhi);

	GTreeManager::Write(chisq_hist);
	GTreeManager::Write(prob_hist);
	GTreeManager::Write(pi0_smeared);
	GTreeManager::Write(pi0_fitted);
	GTreeManager::Write(copl_smeared);
	GTreeManager::Write(copl_fitted);
	GTreeManager::Write(nIter_converged);
	GTreeManager::Write(nIter_all);
	GTreeManager::Write(fitStatus);

	GTreeManager::Write(photon1P_Dsmeared);
	GTreeManager::Write(photon2P_Dsmeared);
	GTreeManager::Write(protonP_Dsmeared);
	GTreeManager::Write(photon1P_Dfitted);
	GTreeManager::Write(photon2P_Dfitted);
	GTreeManager::Write(protonP_Dfitted);

	// Write all GH1's easily
	GTreeManager::Write();
}

double KinFitPhysics::sigma_E(const TLorentzVector* const p)
{
	double E = p->E() - p->M();
	if (p->Theta() < .349)  // TAPS
		return (.018 + .008*sqrt(E))*E;
	else  // CB
		return pow(E, .36)*.02*E;
}

double KinFitPhysics::sigma_theta(const TLorentzVector* const p)
{
	if (p->Theta() < .349)  // TAPS
		return .017;
	else  // CB
		return .044;
}

double KinFitPhysics::sigma_phi(const TLorentzVector* const p)
{
	if (p->Theta() < .349)  // TAPS
		return .017;
	else  // CB
		return .044/sin(p->Theta());
}
