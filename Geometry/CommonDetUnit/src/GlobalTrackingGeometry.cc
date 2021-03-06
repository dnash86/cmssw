/** \file GlobalTrackingGeometry.cc
 *
 *  \author M. Sani
 */

#include <Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h>
#include <FWCore/Utilities/interface/Exception.h>
#include <memory>

GlobalTrackingGeometry::GlobalTrackingGeometry(std::vector<const TrackingGeometry*>& geos)
    : theGeometries(geos),
    theDetTypes(nullptr), theDetUnits(nullptr), theDets(nullptr), theDetUnitIds(nullptr), theDetIds(nullptr)
{}

GlobalTrackingGeometry::~GlobalTrackingGeometry()
{
    delete theDetTypes;
    theDetTypes = nullptr;
    delete theDetUnits;
    theDetUnits = nullptr;
    delete theDets;
    theDets = nullptr;
    delete theDetUnitIds;
    theDetUnitIds = nullptr;
    delete theDetIds;
    theDetIds = nullptr;
}

const GeomDetUnit* GlobalTrackingGeometry::idToDetUnit(DetId id) const {
    
    const TrackingGeometry* tg = slaveGeometry(id);
    
    if (tg != 0) {
      return tg->idToDetUnit(id);
    } else {
      return 0;
    }
}


const GeomDet* GlobalTrackingGeometry::idToDet(DetId id) const{
  
    const TrackingGeometry* tg = slaveGeometry(id);
    
    if (tg != 0) {
        return tg->idToDet(id);
    } else {
      return 0;
    }
}

const TrackingGeometry* GlobalTrackingGeometry::slaveGeometry(DetId id) const {  
  
    int idx = id.det()-1;
    if (id.det() == DetId::Muon) {
        
        idx+=id.subdetId()-1;
    }

    if (theGeometries[idx]==0) throw cms::Exception("NoGeometry") << "No Tracking Geometry is available for DetId " << id.rawId() << std::endl;

    return theGeometries[idx];
}

const TrackingGeometry::DetTypeContainer&
GlobalTrackingGeometry::detTypes( void ) const
{    
   if (!theDetTypes.load(std::memory_order_acquire)) {
       std::unique_ptr<DetTypeContainer> ptr{new DetTypeContainer()};
       for( auto geom = theGeometries.cbegin(), geomEnd = theGeometries.cend(); geom != geomEnd; ++geom )
       {
        if( *geom == 0 ) continue;
        DetTypeContainer detTypes(( *geom )->detTypes());
        if( detTypes.size() + ptr->size() < ptr->capacity()) ptr->resize( detTypes.size() + ptr->size());
        for( auto detType = detTypes.cbegin(), detTypeEnd = detTypes.cend(); detType != detTypeEnd; ++detType )
          ptr->push_back( *detType );
       }
       DetTypeContainer* expect = nullptr;
       if(theDetTypes.compare_exchange_strong(expect, ptr.get(), std::memory_order_acq_rel)) {
           ptr.release();
       }
   }
   return *theDetTypes.load(std::memory_order_acquire);
}

const TrackingGeometry::DetUnitContainer&
GlobalTrackingGeometry::detUnits( void ) const
{
   if (!theDetUnits.load(std::memory_order_acquire)) {
       std::unique_ptr<DetUnitContainer> ptr{new DetUnitContainer()};
       for( auto geom = theGeometries.cbegin(), geomEnd = theGeometries.cend(); geom != geomEnd; ++geom )
       {
        if( *geom == 0 ) continue;
        DetUnitContainer detUnits(( *geom )->detUnits());
        if( detUnits.size() + ptr->size() < ptr->capacity()) ptr->resize( detUnits.size() + ptr->size());
        for( auto detUnit = detUnits.cbegin(), detUnitEnd = detUnits.cend(); detUnit != detUnitEnd; ++detUnit )
          ptr->push_back( *detUnit );
       }
       DetUnitContainer* expect = nullptr;
       if(theDetUnits.compare_exchange_strong(expect, ptr.get(), std::memory_order_acq_rel)) {
           ptr.release();
       }
   }
   return *theDetUnits.load(std::memory_order_acquire);
}

const TrackingGeometry::DetContainer&
GlobalTrackingGeometry::dets( void ) const
{
   if (!theDets.load(std::memory_order_acquire)) {
       std::unique_ptr<DetContainer> ptr{new DetContainer()};
       for( auto geom = theGeometries.cbegin(), geomEnd = theGeometries.cend(); geom != geomEnd; ++geom )
       {
        if( *geom == 0 ) continue;
        DetContainer dets(( *geom )->dets());
        if( dets.size() + ptr->size() < ptr->capacity()) ptr->resize( dets.size() + ptr->size());
        for( auto det = dets.cbegin(), detEnd = dets.cend(); det != detEnd; ++det )
          ptr->push_back( *det );
       }
       DetContainer* expect = nullptr;
       if(theDets.compare_exchange_strong(expect, ptr.get(), std::memory_order_acq_rel)) {
           ptr.release();
       }
   }
   return *theDets.load(std::memory_order_acquire);
}

const TrackingGeometry::DetIdContainer&
GlobalTrackingGeometry::detUnitIds( void ) const
{
   if (!theDetUnitIds.load(std::memory_order_acquire)) {
       std::unique_ptr<DetIdContainer> ptr{new DetIdContainer()};
       for( auto geom = theGeometries.cbegin(), geomEnd = theGeometries.cend(); geom != geomEnd; ++geom )
       {
        if( *geom == 0 ) continue;
        DetIdContainer detUnitIds(( *geom )->detUnitIds());
        if( detUnitIds.size() + ptr->size() < ptr->capacity()) ptr->resize( detUnitIds.size() + ptr->size());
        for( auto detUnitId = detUnitIds.cbegin(), detUnitIdEnd = detUnitIds.cend(); detUnitId != detUnitIdEnd; ++detUnitId )
          ptr->push_back( *detUnitId );
       }
       DetIdContainer* expect = nullptr;
       if(theDetUnitIds.compare_exchange_strong(expect, ptr.get(), std::memory_order_acq_rel)) {
           ptr.release();
       }
   }
   return *theDetUnitIds.load(std::memory_order_acquire);
}

const TrackingGeometry::DetIdContainer&
GlobalTrackingGeometry::detIds( void ) const
{
   if (!theDetIds.load(std::memory_order_acquire)) {
       std::unique_ptr<DetIdContainer> ptr{new DetIdContainer()};
       for( auto geom = theGeometries.cbegin(), geomEnd = theGeometries.cend(); geom != geomEnd; ++geom )
       {
        if( *geom == 0 ) continue;
        DetIdContainer detIds(( *geom )->detIds());
        if( detIds.size() + ptr->size() < ptr->capacity()) ptr->resize( detIds.size() + ptr->size());
        for( auto detId = detIds.cbegin(), detIdEnd = detIds.cend(); detId != detIdEnd; ++detId )
          ptr->push_back( *detId );
       }
       DetIdContainer* expect = nullptr;
       if(theDetIds.compare_exchange_strong(expect, ptr.get(), std::memory_order_acq_rel)) {
           ptr.release();
       }
   }
   return *theDetIds.load(std::memory_order_acquire);
}
