#include "GTreeScaler.h"

#include <TLeaf.h>


GTreeScaler::GTreeScaler(GTreeManager *Manager)    :
    GTree(Manager, TString("scaler"), kTRUE),
    EventNumber(0),
    EventID(0),
    NScaler(0)
{
    for(int i=0; i<GTreeScaler_MAX; i++)
        scalers[i] = 0;
}

GTreeScaler::~GTreeScaler()
{

}

void    GTreeScaler::SetBranchAdresses()
{
    tree_in->SetBranchAddress("eventNumber", &EventNumber);
    tree_in->SetBranchAddress("eventID", &EventID);
    NScaler = tree_in->GetLeaf("scalers")->GetLen();
    if(NScaler<=GTreeScaler_MAX)
        tree_in->SetBranchAddress("scalers", scalers);
}

void    GTreeScaler::SetBranches()
{
    tree_out->Branch("eventNumber", &EventNumber, "eventNumber/I");
    tree_out->Branch("eventID", &EventID, "eventID/I");
    Char_t str[256];
    sprintf(str, "scalers[%d]/i", NScaler);
    tree_out->Branch("scalers", scalers, str);
}


UInt_t  GTreeScaler::GetScalerEntry(const Int_t event_number)
{
    if(!IsOpenForInput())
    {
        if(!OpenForInput())
        {
            std::cout << "Can not open Scaler tree in input file." << std::endl;
            return 0;
        }
    }

    for(Int_t i=0; i<GetNEntries(); i++)
    {
        GetEntry(i);
        if(event_number<EventNumber)
            return i;
    }
}


void    GTreeScaler::SetNScaler(const Int_t num)
{
    NScaler = num;
    if(NScaler>GTreeScaler_MAX)
    {
        std::cout << "#ERROR# GTreeScaler::SetNScaler(const Int_t num): Can not handle " << num << " Scalers! Set to max (" << GTreeScaler_MAX << ")." << std::endl;
        NScaler = GTreeScaler_MAX;
    }
}

void    GTreeScaler::Print() const
{
    GTree::Print();

    std::cout << "GTreeScaler: EventNumber->" << EventNumber << " EventID->" << EventID << std::endl;
    std::cout << "             NScaler->" << NScaler << std::endl;
    for(int i=0; i<NScaler; i++)
        std::cout << "Scaler " << i << ": " << scalers[i] << std::endl;
}

void    GTreeScaler::CloneValidEntries()
{
    if(!IsOpenForInput())
    {
        if(!OpenForInput())
        {
            std::cout << "Can not open " << GetName() << " in input file." << std::endl;
            return;
        }
    }
    if(!IsOpenForOutput())
    {
        if(!OpenForOutput())
        {
            std::cout << "Can not create " << GetName() << " in output file." << std::endl;
            return;
        }
    }

    for(int i=1; i<GetNEntries(); i++)
    {
        GetEntry(i);
        Fill();
    }

}
