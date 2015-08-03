#include "TrackerParametersESModule.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/PTrackerParametersRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerParametersFromDD.h"
#include "CondFormats/GeometryObjects/interface/PTrackerParameters.h"

TrackerParametersESModule::TrackerParametersESModule( const edm::ParameterSet& )
{
  edm::LogInfo("TRACKER") << "TrackerParametersESModule::TrackerParametersESModule";

  setWhatProduced(this);
}

TrackerParametersESModule::~TrackerParametersESModule()
{ 
}

void
TrackerParametersESModule::fillDescriptions( edm::ConfigurationDescriptions & descriptions ) 
{
  edm::ParameterSetDescription desc;
  descriptions.add( "trackerParameters", desc );
}

TrackerParametersESModule::ReturnType
TrackerParametersESModule::produce( const PTrackerParametersRcd& iRecord )
{
  //edm::LogInfo("TrackerParametersESModule")
  std::cout <<  "TrackerParametersESModule::produce(const PTrackerParametersRcd& iRecord)" << std::endl;
  edm::ESTransientHandle<DDCompactView> cpv;
  std::cout <<  "TrackerParametersESModule::produce(const PTrackerParametersRcd& iRecord)" << std::endl;
  iRecord.getRecord<IdealGeometryRecord>().get( cpv );
  std::cout <<  "TrackerParametersESModule::produce(const PTrackerParametersRcd& iRecord)" << std::endl;
  PTrackerParameters* ptp = new PTrackerParameters();
  std::cout <<  "TrackerParametersESModule::produce(const PTrackerParametersRcd& iRecord)" << std::endl;
  TrackerParametersFromDD builder;
  std::cout <<  "TrackerParametersESModule::produce(const PTrackerParametersRcd& iRecord)" << std::endl;
  builder.build( &(*cpv), *ptp );
  std::cout <<  "TrackerParametersESModule::produce(const PTrackerParametersRcd& iRecord)" << std::endl;
  return ReturnType( ptp ) ;
}

DEFINE_FWK_EVENTSETUP_MODULE( TrackerParametersESModule);
