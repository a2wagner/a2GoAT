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

	photon1PullE = new TH1D("photon1PullE", "photon1_pullE", 100, -2, 2);
	photon1PullTheta = new TH1D("photon1PullTheta", "photon1_pullTheta", 100, -2, 2);
	photon1PullPhi = new TH1D("photon1PullPhi", "photon1_pullPhi", 100, -2, 2);
	photon2PullE = new TH1D("photon2PullE", "photon2_pullE", 100, -2, 2);
	photon2PullTheta = new TH1D("photon2PullTheta", "photon2_pullTheta", 100, -2, 2);
	photon2PullPhi = new TH1D("photon2PullPhi", "photon2_pullPhi", 100, -2, 2);
	protonPullE = new TH1D("protonPullE", "proton_pullE", 100, -2, 2);
	protonPullTheta = new TH1D("protonPullTheta", "proton_pullTheta", 100, -2, 2);
	protonPullPhi = new TH1D("protonPullPhi", "proton_pullPhi", 100, -2, 2);

	chisq_hist = new TH1D("chisq", "chi_square", 200, 0, 100);
	prob_hist = new TH1D("prob", "probability", 200, 0, 1);
	pi0_fitted = new TH1D("pi0_fitted", "pi0_mass_fitted", 200, 100, 200);
	nIter = new TH1I("n_iter", "n_iter", 51, 0, 50);
	fitStatus = new TH1I("fit_status", "fit_status", 21, -10, 10);
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




	bool dbg = true;
	static int count = 0;
	double photon1_pullE, photon1_pullTheta, photon1_pullPhi,
			photon2_pullE, photon2_pullTheta, photon2_pullPhi,
			proton_pullE, proton_pullTheta, proton_pullPhi;
	double chisq, prob, pi0_mass_fitted;
	int n_iter, fit_status;



	// first get the MC trees and check if they really exist because you'll still get vectors etc. from non-existent trees!
	//TODO: FIX THIS!
	bool pluto_tree, geant_tree;
	GTreeA2Geant* geant = GetGeant();
	GTreePluto* pluto = GetPluto();
	if(!(geant_tree = geant->IsOpenForInput()))
		cout << "No Geant Tree in the current file!" << endl;
	if(!(pluto_tree = pluto->IsOpenForInput()))
		cout << "No Pluto Tree in the current file!" << endl;

	if (!pluto_tree && !geant_tree) {
		cout << "No MC trees found in the current file, don't expect any reasonable results if you're proceeding!" << endl;
		exit(EXIT_FAILURE);
	}

	/* Gather all particle information, depending on which MC trees are available */
	particle_vector particles;
	TLorentzVector trueBeam;
	unsigned int nParticles = 0, nParticlesCB = 0, nParticlesTAPS = 0;
	//if (pluto_tree) {  // in case we have a Pluto tree
	if (!pluto_tree) {  //TODO: for checks, maybe wrong units, GeV instead of MeV Pluto <--> Geant ?
		for (const PParticle* p : pluto->GetFinalState())
			particles.push_back(particle_t(p->Vect4(), static_cast<particle_id>(p->ID())));
		trueBeam = pluto->GetMCTrue(0)->Vect4();
	} else {  // if not, use Geant tree information
		//TODO: maybe do this directly in the GTreeA2Geant?
		for (unsigned int i = 0; i < geant->GetNTrueParticles(); i++)
			particles.push_back(particle_t(geant->GetTrueVector(i), static_cast<particle_id>(geant->GetTrueID(i))));
		for (particle_it it = particles.begin(); it != particles.end();)
			if (!it->is_final_state())
				it = particles.erase(it);
			else
				++it;
		trueBeam = geant->GetBeam();
	}
	nParticles = particles.size();
	TLorentzVector trueTarget = TLorentzVector(0., 0., 0., MASS_PROTON);

	if (dbg) {
		printf("\nEvent %d:\n", ++count);
		for (auto& p : particles) {
			printf("Particle id: %2d (%s), ", p.id, p.is_charged() ? "charged" : "uncharged");
			p.Print();
		}
	}

	for (particle_t p : particles)
		if (!p.is_final_state())
			printf("no final state!\n");  //TODO: an der Stelle vlt. abbrechen?


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
			if (p.is_charged())
				charged[nCharged++] = i;
			else {
				neutral[nNeutral++] = i;
				tmpState += particles[i].p4;
			}
			i++;
		}
		/* end event selection when we don't have one charged particle (entry in PID or Veto) and two neutral ones */
		if (nCharged != 1 && nNeutral != 2)
			return;

		double protEnergyReconstr = particles[charged[0]].E() - particles[charged[0]].M();
		double invM_2neutral = tmpState.M();
		invM_2g->Fill(invM_2neutral);

		/* At this point we have one charged and two neutral particles. Start to examine all tagged photons for prompt and random windows. */
		//for (unsigned int i = 0; i < nBeamPhotons; i++) {
			//missProton = fP4target[0] + Tagged[i].GetP4() - tmpState;
			missProton = trueTarget + trueBeam - tmpState;
			mMiss = missProton.M();
			protExpect = missProton.E() - missProton.M();
			missM->Fill(mMiss);
			expectProt->Fill(protExpect);
		//}

		/* prepare kinematic fit */
		int err;
		int ndf;

		kinFit.reset();

		TVector3 photon1 = particles[neutral[0]].Vect();
		TVector3 photon2 = particles[neutral[1]].Vect();
		TVector3 proton = particles[charged[0]].Vect();

		//TRandom3 rand(0);
		double sigmaPt = .02, sigmaTheta = .05, sigmaPhi = .05;
		//photon1.SetPtThetaPhi(rand.Gaus(photon1.Pt(), sigmaPt), rand.Gaus(photon1.Theta(), sigmaTheta), rand.Gaus(photon1.Phi(), sigmaPhi));
		//photon2.SetPtThetaPhi(rand.Gaus(photon2.Pt(), sigmaPt), rand.Gaus(photon2.Theta(), sigmaTheta), rand.Gaus(photon2.Phi(), sigmaPhi));
		//proton.SetPtThetaPhi(rand.Gaus(proton.Pt(), sigmaPt), rand.Gaus(proton.Theta(), sigmaTheta), rand.Gaus(proton.Phi(), sigmaPhi));
		TMatrixD covPhoton1;
		TMatrixD covPhoton2;
		TMatrixD covProton;

		TLorentzVector fitPhoton1;
		TLorentzVector fitPhoton2;
		TLorentzVector fitProton;

		int rows = 3;  // number of rows equal to number of cols
		//Double_t errors[3];
		//Double_t errors[] = {50, .5, .5};
		Double_t errors[] = {sigmaPt*sigmaPt*1.2, sigmaTheta*sigmaTheta*1.2, sigmaPhi*sigmaPhi*1.2};
		//Double_t errors[] = {1., .01, .01};
		int currPart = neutral[0];
		//kinFit.sigmaEThetaPhi(particles[currPart], errors);
		if (kinFit.fillSquareMatrixDiagonal(&covPhoton1, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		currPart = neutral[1];
		//kinFit.sigmaEThetaPhi(particles[currPart], errors);
		if (kinFit.fillSquareMatrixDiagonal(&covPhoton2, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		currPart = charged[0];
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
		TFitParticlePThetaPhi pr("charged1", "charged1", &proton, MASS_PROTON, &covProton);

		TVector3 beam = trueBeam.Vect();
		TVector3 target = trueTarget.Vect();
		TMatrixD covBeam;
		TMatrixD covTarget;
		//kinFit.sigmaEThetaPhi(Tagged[0], errors);
		errors[0] = .01; errors[1] = .0001; errors[2] = .0001;
		if (kinFit.fillSquareMatrixDiagonal(&covBeam, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		covTarget.Zero();
		covTarget.ResizeTo(3, 3);
		TFitParticlePThetaPhi bm("beam", "beam", &beam, 0., &covBeam);
		TFitParticlePThetaPhi trgt("target", "target", &target, MASS_PROTON, &covTarget);

		// energy and momentum constraints have to be defined separately for each component. components can be accessed via enum TFitConstraintEp::component
		TFitConstraintEp energyConservation("energyConstr", "Energy conservation constraint", 0, TFitConstraintEp::E, 0.);
		energyConservation.addParticles1(&bm, &trgt);
		energyConservation.addParticles2(&ph1, &ph2, &pr);
		TFitConstraintEp pxConservation("pxConstr", "Px conservation constraint", 0, TFitConstraintEp::pX, 0.);
		pxConservation.addParticles1(&bm, &trgt);
		pxConservation.addParticles2(&ph1, &ph2, &pr);
		TFitConstraintEp pyConservation("pyConstr", "Py conservation constraint", 0, TFitConstraintEp::pY, 0.);
		pyConservation.addParticles1(&bm, &trgt);
		pyConservation.addParticles2(&ph1, &ph2, &pr);
		TFitConstraintEp pzConservation("pzConstr", "Pz conservation constraint", 0, TFitConstraintEp::pZ, 0.);
		pzConservation.addParticles1(&bm, &trgt);
		pzConservation.addParticles2(&ph1, &ph2, &pr);
		TFitConstraintM massConstrProton("massConstr_proton", "mass constraint proton", 0, 0, MASS_PROTON);
		massConstrProton.addParticle1(&pr);
		// constraint for pi0 mass
		TFitConstraintM massConstrPi0("massConstr_pi0", "mass constraint pi0", 0, 0, MASS_PIZERO);
		massConstrPi0.addParticles1(&ph1, &ph2);

		//TKinFitter fit;
		kinFit.addMeasParticle(&ph1);
		kinFit.addMeasParticle(&ph2);
		kinFit.addMeasParticle(&pr);
		kinFit.setParamUnmeas(&pr, 0);  // proton energy unmeasured
		kinFit.addMeasParticle(&bm);
		kinFit.addUnmeasParticle(&trgt);

		//kinFit.addConstraint(&massConstrProton);
		kinFit.addConstraint(&energyConservation);
		kinFit.addConstraint(&pxConservation);
		kinFit.addConstraint(&pyConservation);
		kinFit.addConstraint(&pzConservation);
		kinFit.addConstraint(&massConstrPi0);

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
		kinFit.setMaxF(1e-4);  // max sum of constraints
		// set verbosity level
		if (dbg)
			kinFit.setVerbosity(3);
		else
			kinFit.setVerbosity(0);
		kinFit.fit();

		fitPhoton1 = (*ph1.getCurr4Vec());
		fitPhoton2 = (*ph2.getCurr4Vec());
		fitProton = (*pr.getCurr4Vec());

//        for (outer_it step = temp_results.begin(); step != temp_results.end(); ++step)
//            for (inner_it it = step->begin(); it != step->end(); ++it)
//                it->Print();

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

		std::cout << "Fit result: " << std::endl;
		//kinFit.print();
		ndf = kinFit.getNDF();
		chisq = kinFit.getS();
		prob = TMath::Prob(chisq, ndf);
		chisq_hist->Fill(chisq);
		prob_hist->Fill(prob);
		std::cout << "\nProbability: " << prob << "\tchi^2: " << chisq << std::endl;
		//pi0_mass = (TLorentzVector(photon1, photon1.Pt()) + TLorentzVector(photon2, photon2.Pt())).M();
		pi0_mass_fitted = (fitPhoton1 + fitPhoton2).M();
		pi0_fitted->Fill(pi0_mass_fitted);
		n_iter = kinFit.getNbIter();
		nIter->Fill(n_iter);
		fit_status = kinFit.getStatus();  // Status: -1: "NO FIT PERFORMED", 10: "RUNNING", 0: "CONVERGED", 1: "NOT CONVERGED", -10: "ABORTED"; should only return 0 or 1
		fitStatus->Fill(fit_status);

		//delete temp_results;

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
	GTreeManager::Write(pi0_fitted);
	GTreeManager::Write(nIter);
	GTreeManager::Write(fitStatus);

	// Write all GH1's easily
	GTreeManager::Write();
}
