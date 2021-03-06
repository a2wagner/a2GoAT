#ifndef __GParticleReconstruction_h__
#define __GParticleReconstruction_h__

#include <TCutG.h>


#include "GDataChecks.h"

#define	pdg_rootino 0

class	GParticleReconstruction : public GDataChecks
{
public:
    enum ReconstructionType
    {
        ReconstructionType_AllPhotons,
        ReconstructionType_AllProtons,
        ReconstructionType_dEoverE
    };
    enum dEoverE_Type
    {
        dEoverE_Cut_None     = 0,
        dEoverE_Cut_Proton   = 1,
        dEoverE_Cut_PiPlus   = 2,
        dEoverE_Cut_Electron = 4
    };

private:
	std::string config;

    ReconstructionType                 CB_type;
    dEoverE_Type    CB_dEoverE_type;
    ReconstructionType                 TAPS_type;
    dEoverE_Type    TAPS_dEoverE_type;

    char 		cutfilename[256];
    char 		cutname[256];

    TFile* 		CutFile;
    TCutG* 		Cut;
    TCutG* 		Cut_CB_proton;
    TCutG* 		Cut_CB_pion;
    TCutG*		Cut_CB_electron;
    TCutG* 		Cut_TAPS_proton;
    TCutG* 		Cut_TAPS_pion;
    TCutG*		Cut_TAPS_electron;

    Double_t	charged_theta_min;
    Double_t	charged_theta_max;

    Int_t		Cut_CB_proton_active;
    Int_t		Cut_TAPS_proton_active;
    Int_t 		Cut_proton_active;

    Int_t		Cut_CB_pion_active;
    Int_t		Cut_TAPS_pion_active;
    Int_t 		Cut_pion_active;

    Int_t		Cut_CB_electron_active;
    Int_t		Cut_TAPS_electron_active;
    Int_t 		Cut_electron_active;

    Int_t* 		Identified;
    Int_t* 		Charge;
    
    Bool_t 		charge_ignore_PID;
    Bool_t 		charge_ignore_WC0;
    Bool_t 		charge_ignore_WC1;
    Bool_t 		charge_ignore_VETO;    

    Double_t    CBTimeCut[2];
    Double_t    TAPSTimeCut[2];

    Bool_t      DoScalerCorrection;
    Bool_t      DoTrigger;
    Double_t    E_Sum;
    Int_t       multiplicity;

    TCutG*	OpenCutFile(Char_t* filename, Char_t* cutname);
    Bool_t  Trigger();

protected:
            Bool_t  ProcessEventWithoutFilling();
    virtual void    ProcessEvent();
    virtual Bool_t  Start();

public:
    GParticleReconstruction();
    virtual ~GParticleReconstruction();

    Bool_t	Init();
    void    SetCBTimeCut(const Double_t min, const Double_t max)    {CBTimeCut[0]=min; CBTimeCut[1]=max;}
    void    SetTAPSTimeCut(const Double_t min, const Double_t max)  {TAPSTimeCut[0]=min; TAPSTimeCut[1]=max;}
    void    SetScalerCorrection(const Bool_t value)                 {DoScalerCorrection = value;}
    void    SetTrigger(const Double_t esum, const Int_t mult)       {DoTrigger = kTRUE; E_Sum = esum; multiplicity = mult;}
    void    SetCBType(const ReconstructionType type, const dEoverE_Type dEoverE_type = dEoverE_Cut_None)    {CB_type = type; CB_dEoverE_type = dEoverE_type;}
    void    SetTAPSType(const ReconstructionType type, const dEoverE_Type dEoverE_type = dEoverE_Cut_None)  {TAPS_type = type; TAPS_dEoverE_type = dEoverE_type;}
    void    SetThetaRangeChargedParticles(const Double_t min, const Double_t max)  {charged_theta_min = min; charged_theta_max = max;}

    //void	CheckNeutrality();
    //void 	PhotonReconstruction();
    void 	ChargedReconstructionCB(const Int_t index);
    void 	ChargedReconstructionTAPS(const Int_t index);
    void 	MesonReconstruction();
    //void	AddParticle(Int_t pdg_code, Int_t nindex, Int_t index_list[]);
    //void	AddParticle(Int_t pdg_code, Int_t i)                    {Int_t index_list[1]; index_list[0] = i; AddParticle(pdg_code, 1, index_list);}
};

#endif
