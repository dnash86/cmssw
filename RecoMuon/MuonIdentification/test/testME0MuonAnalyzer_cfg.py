import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")


#process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')

process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')


process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")

process.load( "DQMServices/Core/DQMStore_cfg" )

process.load('Validation.RecoMuon.associators_cff')
process.load('Validation.RecoMuon.selectors_cff')
process.load('Validation.RecoMuon.MuonTrackValidator_cfi')
process.load('Validation.RecoMuon.RecoMuonValidator_cfi')



process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)



process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:///somewhere/simevent.root') ##/somewhere/simevent.root" }


)

from Validation.RecoMuon.selectors_cff import *
from Validation.RecoMuon.associators_cff import *
# Configurations for MuonTrackValidators
import Validation.RecoMuon.MuonTrackValidator_cfi


# Configurations for RecoMuonValidators
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from Validation.RecoMuon.RecoMuonValidator_cfi import *

#import SimGeneral.MixingModule.mixNoPU_cfi
from SimMuon.MCTruth.MuonAssociatorByHitsESProducer_NoSimHits_cfi import *
from SimMuon.MCTruth.MuonAssociatorByHits_cfi import muonAssociatorByHitsCommonParameters

#process.TrackAssociatorByChi2ESProducer = Validation.RecoMuon.associators_cff.TrackAssociatorByChi2ESProducer.clone(chi2cut = 100.0,ComponentName = 'TrackAssociatorByChi2')     ####Uncomment this for association by chi2

###  Thsi is for association by hits
import SimMuon.MCTruth.MuonAssociatorByHitsESProducer_cfi

process.muonAssociatorByHits = SimMuon.MCTruth.MuonAssociatorByHitsESProducer_cfi.muonAssociatorByHitsESProducer.clone(ComponentName = 'muonAssociatorByHits',
 #tpTag = 'mix:MergedTrackTruth',
 UseTracker = True, 
 UseMuon = False,
 EfficiencyCut_track = cms.double(0.0),   
 PurityCut_track = cms.double(0.0)       
 )
###Note: not yet configurable for ME0 ID, set ID here: Validation/RecoMuon/plugins/ME0MuonTrackCollProducer.cc
#----------ME0Muon Collection Production for association by chi2
process.me0muon = cms.EDProducer("ME0MuonTrackCollProducer",     
                         me0MuonTag = cms.InputTag("me0SegmentMatching"),
                         selectionTags = cms.vstring('All'),
                         doTimeWindow=cms.untracked.bool(False),
                         TimeWindowCut=cms.untracked.double(12.5)  ## 18 +/- this parameter
                         )
#--------------------
process.me0muonColl_seq = cms.Sequence(
    process.me0muon
    )

process.recoMuonValidation = cms.Sequence( me0muonColl_seq   )

process.Test = cms.EDAnalyzer("ME0MuonAnalyzer",
                              

                              HistoFolder = cms.string('OUTPUTTEMPLATE'),  #This will make a folder with plots
                              HistoFile = cms.string('OUTPUTTEMPLATE.root'),  #The root file where plots also go

                              ME0MuonSelectionType = cms.string('Loose'),   ##Not actually implemented yet in analyzer, must choose ME0 type inside analyzer
                              FakeRatePtCut = cms.double(5.0),  ##pT cut
                              MatchingWindowDelR = cms.double (.15),  ##matching window for DelR matching
                              RejectEndcapMuons = cms.bool(False),  
                              UseAssociators = cms.bool(True),   ##Must be true for association by hits/chi2 to work
                              associators = cms.vstring('muonAssociatorByHits'),  ##Parameters for assocByHits can be modified above (Efficiency & Purity)

                              label = cms.VInputTag('me0muon'),
                              
                              
)

process.p = cms.Path(process.recoMuonValidation*process.Test)



process.PoolSource.fileNames = [

FILETEMPLATE

]
