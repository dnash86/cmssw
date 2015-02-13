/** \file ME0SegmentProducer.cc
 *
 */

#include <RecoLocalMuon/GEMRecHit/src/ME0SegmentProducer.h>
#include <RecoLocalMuon/GEMRecHit/src/ME0SegmentBuilder.h>

#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 

#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include <DataFormats/GEMRecHit/interface/ME0RecHitCollection.h>
#include <DataFormats/GEMRecHit/interface/ME0SegmentCollection.h>
#include <DataFormats/GEMRecHit/interface/ME0Segment.h>

ME0SegmentProducer::ME0SegmentProducer(const edm::ParameterSet& pas) : iev(0) {
	
    inputObjectsTag = pas.getParameter<edm::InputTag>("me0RecHitLabel");
    segmentBuilder_ = new ME0SegmentBuilder(pas); // pass on the PS
    LogDebug("ME0Segment|ME0") << "Constructor of ME0SegmentProducer"<<std::endl;
  	// register what this produces
    produces<ME0SegmentCollection>();
}

ME0SegmentProducer::~ME0SegmentProducer() {

    LogDebug("ME0Segment|ME0") << "deleting ME0SegmentBuilder after " << iev << " events w/csc data.";
    delete segmentBuilder_;
}

void ME0SegmentProducer::produce(edm::Event& ev, const edm::EventSetup& setup) {

    LogDebug("ME0Segment|ME0") << "start producing segments for " << ++iev << "th event with csc data";
	
    // find the geometry (& conditions?) for this event & cache it in the builder
    edm::ESHandle<ME0Geometry> me0g;
    setup.get<MuonGeometryRecord>().get(me0g);
    const ME0Geometry* mgeom = &*me0g;
    segmentBuilder_->setGeometry(mgeom);
  

    // get the collection of ME0RecHit
    LogDebug("ME0Segment|ME0") << "Getting me0rechits by label: "<<inputObjectsTag<<std::endl;
    edm::Handle<ME0RecHitCollection> me0RecHits;
    ev.getByLabel(inputObjectsTag, me0RecHits);  

    // create empty collection of Segments
    LogDebug("ME0Segment|ME0") << "Creating empty me0segmentcollection"<<std::endl;
    std::auto_ptr<ME0SegmentCollection> oc( new ME0SegmentCollection );

  	// fill the collection
    LogDebug("ME0Segment|ME0") << "Calling the builder"<<std::endl;
    segmentBuilder_->build(me0RecHits.product(), *oc); //@@ FILL oc

    // put collection in event
    ev.put(oc);
}
