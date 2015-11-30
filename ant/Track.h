#ifndef TRACK_H
#define TRACK_H

#include "base/printable.h"
#include "base/types.h"
#include "Detector.h"

#include <ostream>
#include <memory>
#include <vector>

namespace ant {

/**
 * @brief The Track class
 * Representation of GoAT information, with emphasis on
 * physical particle information, not on detector information!
 */
class Track: public ant::printable_traits
{
private:
    mev_t clusterEnergy;
    radian_t theta;
    radian_t phi;
    ns_t time;
    clustersize_t clusterSize;
    element_index_t centralCrystal;
    element_index_t centralVeto;
    detector_t detector;
    mev_t vetoEnergy;
    mev_t _MWPC0Energy;
    mev_t _MWPC1Energy;
    mev_t shortEnergy;
public:
    Track(const mev_t& _clusterEnergy,
          const radian_t& _theta,
          const radian_t& _phi,
          const ns_t& _time,
          const clustersize_t& _clusterSize,
          const element_index_t& _centralCrystal,
          const element_index_t& _centralVeto,
          const detector_t& _detector,
          const mev_t& _vetoEnergy,
          const mev_t& _MWPC0Energy,
          const mev_t& _MWPC1Energy,
          const mev_t& _shortEnergy
          ) :
        clusterEnergy(_clusterEnergy),
        theta(_theta),
        phi(_phi),
        time(_time),
        clusterSize(_clusterSize),
        centralCrystal(_centralCrystal),
        centralVeto(_centralVeto),
        detector(_detector),
        vetoEnergy(_vetoEnergy),
        _MWPC0Energy(_MWPC0Energy),
        _MWPC1Energy(_MWPC1Energy),
        shortEnergy(_shortEnergy)
    {}


    mev_t ClusterEnergy() const { return clusterEnergy; }
    radian_t Theta() const { return theta; }
    radian_t Phi() const { return phi; }
    ns_t Time() const { return time; }
    clustersize_t ClusterSize() const { return clusterSize; }
    element_index_t CentralCrystal() const { return centralCrystal; }
    element_index_t CentralVeto() const { return centralVeto; }
    detector_t Detector() const { return detector; }
    mev_t VetoEnergy() const { return vetoEnergy; }
    mev_t MWPC0Energy() const { return _MWPC0Energy; }
    mev_t MWPC1Energy() const { return _MWPC1Energy; }
    mev_t ShortEnergy() const { return shortEnergy; }

    virtual std::ostream &Print(std::ostream &stream) const;

    void SetClusterEnergy(const mev_t);
    void SetTheta(const radian_t);
    void SetTime(const ns_t);
};

using TrackPtr  = std::shared_ptr<Track>;
using TrackList = std::vector<TrackPtr>;

}

#endif // TRACK_H
