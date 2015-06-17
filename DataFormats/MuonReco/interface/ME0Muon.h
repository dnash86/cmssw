#ifndef MuonReco_ME0Muon_h
#define MuonReco_ME0Muon_h
/** \class reco::ME0Muon ME0Muon.h DataFormats/MuonReco/interface/ME0Muon.h
 *  
 * \author David Nash NEU
 *
 *
 */
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"

#include <DataFormats/GEMRecHit/interface/ME0SegmentCollection.h>

namespace reco {
 
  class ME0Muon {
  public:
    ME0Muon();
    //ME0Muon( const TrackRef & t, const ME0Segment & s) { innerTrack_ = t; me0Segment_ = s;}
    ME0Muon( const TrackRef & t, const ME0Segment & s, const int v, const double c) { innerTrack_ = t; me0Segment_ = s; me0segid_=v; trackCharge_ = c;}
    virtual ~ME0Muon(){}     
    
    /// reference to Track reconstructed in the tracker only
    TrackRef innerTrack() const { return innerTrack_; }
    TrackRef track() const { return innerTrack(); }
    /// set reference to Track
    void setInnerTrack( const TrackRef & t ) { innerTrack_ = t; }
    void setTrack( const TrackRef & t ) { setInnerTrack(t); }
    /// set reference to our new ME0Segment type
    void setME0Segment( const ME0Segment & s ) { me0Segment_ = s; }

    const ME0Segment& me0segment() const { return me0Segment_; }
    
    //Added for testing
    void setme0segid( const int v){me0segid_=v;}
    int me0segid() const {return me0segid_;}

    /// a bunch of useful accessors
    const int& charge() const { return innerTrack_.get()->charge(); }
    /// polar angle  
    const double& theta() const { return innerTrack_.get()->theta(); }
    /// momentum vector magnitude
    const double& p() const { return innerTrack_.get()->p(); }
    /// track transverse momentum
    const double& pt() const { return innerTrack_.get()->pt(); }
    /// x coordinate of momentum vector
    const double& px() const { return innerTrack_.get()->px(); }
    /// y coordinate of momentum vector
    const double& py() const { return innerTrack_.get()->py(); }
    /// z coordinate of momentum vector
    const double& pz() const { return innerTrack_.get()->pz(); }
    /// azimuthal angle of momentum vector
    const double& phi() const { return innerTrack_.get()->phi(); }
    /// pseudorapidity of momentum vector
    const double& eta() const { return innerTrack_.get()->eta(); }

    const GlobalPoint& globalTrackPosAtSurface() const { return globalTrackPosAtSurface_; }
    const GlobalVector& globalTrackMomAtSurface() const { return globalTrackMomAtSurface_; }
    int trackCharge() const { return trackCharge_; }
    const AlgebraicSymMatrix66& trackCov() const { return trackCov_; }

    void setGlobalTrackPosAtSurface(const GlobalPoint globalTrackPosAtSurface) { globalTrackPosAtSurface_ = globalTrackPosAtSurface; }
    void setGlobalTrackMomAtSurface(const GlobalVector globalTrackMomAtSurface) { globalTrackMomAtSurface_ = globalTrackMomAtSurface; }
    void setTrackCharge(const int trackCharge) { trackCharge_ = trackCharge; }
    void setTrackCov(const AlgebraicSymMatrix66 trackCov) { trackCov_ = trackCov; }
     
  private:
    /// reference to Track reconstructed in the tracker only
    TrackRef innerTrack_;
    ME0Segment me0Segment_;
    int me0segid_;

    GlobalPoint globalTrackPosAtSurface_;
    GlobalVector globalTrackMomAtSurface_;
    int trackCharge_;
    AlgebraicSymMatrix66 trackCov_;

    //double xpull_,ypull_,xdiff_,ydiff_,phidirdiff_;
  };

}


#endif


