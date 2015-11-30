#include "Track.h"

using namespace std;
using namespace ant;

ostream& Track::Print(ostream &stream) const
{
    stream << "Track "
           << " ClusterEnergy=" << ClusterEnergy()
           << " Theta=" << Theta()
           << " Phi=" << Phi()
           << " Time=" << Time()
           << " ClusterSize=" << ClusterSize()
           << " Detectors=" << Detector()
           << " VetoEnergy=" << VetoEnergy()
           << " MWPC0Energy=" << MWPC0Energy()
           << " MWPC1Energy=" << MWPC1Energy()
           << " ShortEnergy=" << ShortEnergy();
    return stream;
}

void Track::SetClusterEnergy(const mev_t e)
{
    clusterEnergy = e;
}

void Track::SetTheta(const radian_t t)
{
    theta = t;
}

void Track::SetTime(const ns_t t)
{
    time = t;
}
