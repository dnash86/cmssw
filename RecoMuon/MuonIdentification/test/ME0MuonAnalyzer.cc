#include "FWCore/Framework/interface/Event.h"

#include <FWCore/PluginManager/interface/ModuleDef.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 
#include "FWCore/Utilities/interface/InputTag.h"

#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include <DataFormats/MuonReco/interface/EmulatedME0Segment.h>
#include <DataFormats/MuonReco/interface/EmulatedME0SegmentCollection.h>

#include <DataFormats/MuonReco/interface/ME0Muon.h>
#include <DataFormats/MuonReco/interface/ME0MuonCollection.h>

// #include "CLHEP/Matrix/SymMatrix.h"
// #include "CLHEP/Matrix/Matrix.h"
// #include "CLHEP/Vector/ThreeVector.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TLorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

//#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
//#include "TRandom3.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/Math/interface/deltaPhi.h"
//#include <deltaR.h>
//#include <DataFormats/GEMRecHit/interface/ME0SegmentCollection.h>

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"


//Associator for chi2: Including header files
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

#include "DataFormats/MuonReco/interface/Muon.h"

#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimTracker/TrackAssociation/plugins/CosmicParametersDefinerForTPESProducer.h"

#include "CommonTools/CandAlgos/interface/GenParticleCustomSelector.h"
//#include "CommonTools/CandAlgos/interface/TrackingParticleCustomSelector.h"

//#include "RecoMuon/MuonIdentification/plugins/ME0MuonSelector.cc"

#include "Fit/FitResult.h"
#include "TF1.h" 


#include "TMath.h"
#include "TLorentzVector.h"

#include "TH1.h" 
#include <TH2.h>
#include "TFile.h"
#include <TProfile.h>
#include "TStyle.h"
#include <TCanvas.h>
#include <TLatex.h>
//#include "CMSStyle.C"
//#include "tdrstyle.C"
//#include "lumi.C"
//#include <DataFormats/GEMRecHit/interface/ME0SegmentCollection.h>

//#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
//#include <Geometry/GEMGeometry/interface/ME0EtaPartition.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
//#include <DataFormats/MuonDetId/interface/ME0DetId.h>


#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToLocal.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCartesian.h"
#include "TGraph.h"

#include <sstream>    

#include <iostream>
#include <fstream>
#include <sys/stat.h>

class ME0MuonAnalyzer : public edm::EDAnalyzer {
public:
  explicit ME0MuonAnalyzer(const edm::ParameterSet&);
  ~ME0MuonAnalyzer();
  FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& , 
			     int , const AlgebraicSymMatrix66& ,
			     const MagneticField* );

  FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& , 
			     int , const AlgebraicSymMatrix55& ,
			     const MagneticField* );
    void getFromFTS(const FreeTrajectoryState& ,
		  GlobalVector& , GlobalVector& , 
		  int& , AlgebraicSymMatrix66& );


  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void endJob();
  //virtual void beginJob(const edm::EventSetup&);
  //void beginJob();
  void beginRun(edm::Run const&, edm::EventSetup const&);
  void endRun(edm::Run const&, edm::EventSetup const&);


  //For Track Association



  //protected:
  
  private:
  //Associator for chi2: objects
  //edm::InputTag associatormap;
  bool UseAssociators;
  bool RejectEndcapMuons;
  //const TrackAssociatorByChi2* associatorByChi2;
  //const TrackAssociatorByHits* associatorByHits;
  std::vector<std::string> associators;
  std::vector<const TrackAssociatorBase*> associator;
  std::vector<edm::InputTag> label;
  //GenParticleCustomSelector gpSelector;	
  //TrackingParticleCustomSelector gpSelector;	
  //std::string parametersDefiner;


  TString histoFolder;
  TString me0MuonSelector;
  TFile* histoFile; 
  TH1F *Candidate_Eta;  TH1F *Mass_h; 
  TH1F *Segment_Eta;    TH1F *Segment_Phi;    TH1F *Segment_R;  TH2F *Segment_Pos;  
  TH1F *Rechit_Eta;    TH1F *Rechit_Phi;    TH1F *Rechit_R;  TH2F *Rechit_Pos;  
  TH1F *GenMuon_Phi;    TH1F *GenMuon_R;  TH2F *GenMuon_Pos;  
  TH1F *Track_Eta; TH1F *Track_Pt;  TH1F *ME0Muon_Eta; TH1F *ME0Muon_Pt;  TH1F *CheckME0Muon_Eta; 
  TH1F *ME0Muon_Cuts_Eta_5_10;   TH1F *ME0Muon_Cuts_Eta_9_11; TH1F *ME0Muon_Cuts_Eta_10_50; TH1F *ME0Muon_Cuts_Eta_50_100; TH1F *ME0Muon_Cuts_Eta_100; 
  TH1F *UnmatchedME0Muon_Eta; TH1F *UnmatchedME0Muon_Pt;    TH1F *UnmatchedME0Muon_Window_Pt;    TH1F *Chi2UnmatchedME0Muon_Eta; 
  TH1F *UnmatchedME0Muon_Cuts_Eta_5_10;    TH1F *UnmatchedME0Muon_Cuts_Eta_9_11;  TH1F *UnmatchedME0Muon_Cuts_Eta_10_50;  TH1F *UnmatchedME0Muon_Cuts_Eta_50_100;  TH1F *UnmatchedME0Muon_Cuts_Eta_100;
  TH1F *TracksPerSegment_h;  TH2F *TracksPerSegment_s;  TProfile *TracksPerSegment_p;
  TH2F *ClosestDelR_s; TProfile *ClosestDelR_p;
  TH2F *PtDiff_s; TProfile *PtDiff_p; TH1F *PtDiff_h; TH1F *QOverPtDiff_h; TH1F *PtDiff_rms; TH1F *PtDiff_gaus_narrow; TH1F *PtDiff_gaus_wide;
  TH2F *StandalonePtDiff_s; TProfile *StandalonePtDiff_p; TH1F *StandalonePtDiff_h; TH1F *StandaloneQOverPtDiff_h; TH1F *StandalonePtDiff_rms; TH1F *StandalonePtDiff_gaus_narrow; TH1F *StandalonePtDiff_gaus_wide;
  TH1F *PtDiff_gaus_5_10;  TH1F *PtDiff_gaus_10_50;  TH1F *PtDiff_gaus_50_100; TH1F *PtDiff_gaus_100;
  TH1F *StandalonePtDiff_gaus;
  TH1F *VertexDiff_h;
  TH2F *PDiff_s; TProfile *PDiff_p; TH1F *PDiff_h;
  TH2F *PtDiff_s_5_10;    TH2F *PtDiff_s_10_50;    TH2F *PtDiff_s_50_100;    TH2F *PtDiff_s_100;
  TH1F *FakeTracksPerSegment_h;  TH2F *FakeTracksPerSegment_s;  TProfile *FakeTracksPerSegment_p;
  TH1F *FakeTracksPerAssociatedSegment_h;  TH2F *FakeTracksPerAssociatedSegment_s;  TProfile *FakeTracksPerAssociatedSegment_p;
  TH1F *GenMuon_Eta; TH1F *GenMuon_Pt;   TH1F *MatchedME0Muon_Eta; TH1F *MatchedME0Muon_Pt; TH1F *Chi2MatchedME0Muon_Eta; TH1F *Chi2MatchedME0Muon_Pt; 
  TH1F *TPMuon_Eta;
  TH1F *MatchedME0Muon_Eta_5_10;    TH1F *MatchedME0Muon_Eta_9_11;  TH1F *MatchedME0Muon_Eta_10_50;  TH1F *MatchedME0Muon_Eta_50_100;  TH1F *MatchedME0Muon_Eta_100;
  TH1F *Chi2MatchedME0Muon_Eta_5_10;   TH1F *Chi2MatchedME0Muon_Eta_9_11; TH1F *Chi2MatchedME0Muon_Eta_10_50;  TH1F *Chi2MatchedME0Muon_Eta_50_100;  TH1F *Chi2MatchedME0Muon_Eta_100;
  TH1F *GenMuon_Eta_5_10;   TH1F *GenMuon_Eta_9_11;  TH1F *GenMuon_Eta_10_50;  TH1F *GenMuon_Eta_50_100;  TH1F *GenMuon_Eta_100;
  TH1F *MuonRecoEff_Eta;  TH1F *MuonRecoEff_Pt;   TH1F *Chi2MuonRecoEff_Eta;  
  TH1F *MuonRecoEff_Eta_5_10;   TH1F *MuonRecoEff_Eta_9_11;  TH1F *MuonRecoEff_Eta_10_50;  TH1F *MuonRecoEff_Eta_50_100;  TH1F *MuonRecoEff_Eta_100;
  TH1F *Chi2MuonRecoEff_Eta_5_10;    TH1F *Chi2MuonRecoEff_Eta_9_11;  TH1F *Chi2MuonRecoEff_Eta_10_50;  TH1F *Chi2MuonRecoEff_Eta_50_100;  TH1F *Chi2MuonRecoEff_Eta_100;
  TH1F *FakeRate_Eta;  TH1F *FakeRate_Pt;  TH1F *FakeRate_Eta_PerEvent;    TH1F *Chi2FakeRate_Eta;  

  TH1F *Chi2FakeRate_WideBinning_Eta;  
  TH1F *Chi2FakeRate_WidestBinning_Eta;  
  TH1F *FakeRate_WideBinning_Eta;
  TH1F *FakeRate_WidestBinning_Eta;
  TH1F *UnmatchedME0Muon_Cuts_WideBinning_Eta;
  TH1F *UnmatchedME0Muon_Cuts_WidestBinning_Eta;
  TH1F *ME0Muon_Cuts_WideBinning_Eta; 
  TH1F *ME0Muon_Cuts_WidestBinning_Eta;
  TH1F *Chi2UnmatchedME0Muon_WideBinning_Eta; 
  TH1F *Chi2UnmatchedME0Muon_WidestBinning_Eta; 
  TH1F *TPMuon_WideBinning_Eta;
  TH1F *TPMuon_WidestBinning_Eta;
  TH1F *GenMuon_WideBinning_Eta;
  TH1F *GenMuon_WidestBinning_Eta;
  TH1F *MatchedME0Muon_WideBinning_Eta;
  TH1F *MatchedME0Muon_WidestBinning_Eta;
  TH1F *Chi2MatchedME0Muon_WideBinning_Eta;
  TH1F *Chi2MatchedME0Muon_WidestBinning_Eta;
  TH1F *MuonRecoEff_WideBinning_Eta;
  TH1F *MuonRecoEff_WidestBinning_Eta;
  TH1F *Chi2MuonRecoEff_WideBinning_Eta;  
  TH1F *Chi2MuonRecoEff_WidestBinning_Eta;  


  TH1F *FakeRate_Eta_5_10;    TH1F *FakeRate_Eta_9_11;  TH1F *FakeRate_Eta_10_50;  TH1F *FakeRate_Eta_50_100;  TH1F *FakeRate_Eta_100;
  TH1F *MuonAllTracksEff_Eta;  TH1F *MuonAllTracksEff_Pt;
  TH1F *MuonUnmatchedTracksEff_Eta;  TH1F *MuonUnmatchedTracksEff_Pt; TH1F *FractionMatched_Eta;

  TH1F *StandaloneMuonRecoEff_Eta;   TH1F *StandaloneMuonRecoEff_WideBinning_Eta;   TH1F *StandaloneMuonRecoEff_WidestBinning_Eta;
  TH1F *UnmatchedME0Muon_Cuts_Eta;TH1F *ME0Muon_Cuts_Eta;
  TH1F *StandaloneMatchedME0Muon_Eta;    TH1F *StandaloneMatchedME0Muon_WideBinning_Eta;    TH1F *StandaloneMatchedME0Muon_WidestBinning_Eta;
  TH1F *DelR_Segment_GenMuon;

  TH1F *SegPosDirPhiDiff_True_h;    TH1F *SegPosDirEtaDiff_True_h;     TH1F *SegPosDirPhiDiff_All_h;    TH1F *SegPosDirEtaDiff_All_h;   
  TH1F *SegTrackDirPhiDiff_True_h;    TH1F *SegTrackDirEtaDiff_True_h;     TH1F *SegTrackDirPhiDiff_All_h;    TH1F *SegTrackDirEtaDiff_All_h;   TH1F *SegTrackDirPhiPull_True_h;   TH1F *SegTrackDirPhiPull_All_h;   

  TH1F *SegGenDirPhiDiff_True_h;    TH1F *SegGenDirEtaDiff_True_h;     TH1F *SegGenDirPhiDiff_All_h;    TH1F *SegGenDirEtaDiff_All_h;   TH1F *SegGenDirPhiPull_True_h;   TH1F *SegGenDirPhiPull_All_h;   

  TH1F *XDiff_h;   TH1F *YDiff_h;   TH1F *XPull_h;   TH1F *YPull_h;


  TH1F *DelR_Window_Under5; TH1F  *Pt_Window_Under5;
  TH1F *DelR_Track_Window_Under5; TH1F  *Pt_Track_Window_Under5;  TH1F  *Pt_Track_Window;
  TH1F *DelR_Track_Window_Failed_Under5; TH1F  *Pt_Track_Window_Failed_Under5;  TH1F  *Pt_Track_Window_Failed;

  TH1F *FailedTrack_Window_XPull;    TH1F *FailedTrack_Window_YPull;    TH1F *FailedTrack_Window_PhiDiff;
  TH1F *FailedTrack_Window_XDiff;    TH1F *FailedTrack_Window_YDiff;    

  TH1F *NormChi2_h;    TH1F *NormChi2Prob_h; TH2F *NormChi2VsHits_h;	TH2F *chi2_vs_eta_h;  TH1F *AssociatedChi2_h;  TH1F *AssociatedChi2_Prob_h;

  TH1F *PreMatch_TP_R;   TH1F *PostMatch_TP_R;  TH1F *PostMatch_BX0_TP_R;

  double  FakeRatePtCut, MatchingWindowDelR;

  double Nevents;

  
//Removing this
};

ME0MuonAnalyzer::ME0MuonAnalyzer(const edm::ParameterSet& iConfig) 
{
  std::cout<<"Contructor"<<std::endl;
  histoFile = new TFile(iConfig.getParameter<std::string>("HistoFile").c_str(), "recreate");
  histoFolder = iConfig.getParameter<std::string>("HistoFolder").c_str();
  me0MuonSelector = iConfig.getParameter<std::string>("ME0MuonSelectionType").c_str();
  RejectEndcapMuons = iConfig.getParameter< bool >("RejectEndcapMuons");
  UseAssociators = iConfig.getParameter< bool >("UseAssociators");

  FakeRatePtCut   = iConfig.getParameter<double>("FakeRatePtCut");
  MatchingWindowDelR   = iConfig.getParameter<double>("MatchingWindowDelR");

  //Associator for chi2: getting parametters
  //associatormap = iConfig.getParameter< edm::InputTag >("associatormap");
  //UseAssociators = iConfig.getParameter< bool >("UseAssociators");
  associators = iConfig.getParameter< std::vector<std::string> >("associators");

  label = iConfig.getParameter< std::vector<edm::InputTag> >("label");

  // gpSelector = GenParticleCustomSelector(iConfig.getParameter<double>("ptMinGP"),
  // 					 iConfig.getParameter<double>("minRapidityGP"),
  // 					 iConfig.getParameter<double>("maxRapidityGP"),
  // 					 iConfig.getParameter<double>("tipGP"),
  // 					 iConfig.getParameter<double>("lipGP"),
  // 					 iConfig.getParameter<bool>("chargedOnlyGP"),
  // 					 iConfig.getParameter<int>("statusGP"),
  // 					 iConfig.getParameter<std::vector<int> >("pdgIdGP"));

    // gpSelector = TrackingParticleCustomSelector(iConfig.getParameter<double>("ptMinGP"),
  // 					 iConfig.getParameter<double>("minRapidityGP"),
  // 					 iConfig.getParameter<double>("maxRapidityGP"),
  // 					 iConfig.getParameter<double>("tipGP"),
  // 					 iConfig.getParameter<double>("lipGP"),
  // 					 iConfig.getParameter<bool>("chargedOnlyGP"),
  // 					 iConfig.getParameter<int>("statusGP"),
  // 					 iConfig.getParameter<std::vector<int> >("pdgIdGP"));
  // //
    //parametersDefiner =iConfig.getParameter<std::string>("parametersDefiner");

  std::cout<<"Contructor end"<<std::endl;
}



//void ME0MuonAnalyzer::beginJob(const edm::EventSetup& iSetup)
void ME0MuonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {

//void ME0MuonAnalyzer::beginJob()
//{

  std::cout<<"At start of begin run"<<std::endl;  
  mkdir(histoFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  Candidate_Eta = new TH1F("Candidate_Eta"      , "Candidate #eta"   , 4, 2.0, 2.8 );

  Track_Eta = new TH1F("Track_Eta"      , "Track #eta"   , 4, 2.0, 2.8 );
  Track_Pt = new TH1F("Track_Pt"      , "Muon p_{T}"   , 120,0 , 120. );

  Segment_Eta = new TH1F("Segment_Eta"      , "Segment #eta"   , 4, 2.0, 2.8 );
  Segment_Phi = new TH1F("Segment_Phi"      , "Segment #phi"   , 60, -3, 3. );
  Segment_R = new TH1F("Segment_R"      , "Segment r"   , 30, 0, 150 );
  Segment_Pos = new TH2F("Segment_Pos"      , "Segment x,y"   ,100,-100.,100., 100,-100.,100. );

  Rechit_Eta = new TH1F("Rechit_Eta"      , "Rechit #eta"   , 4, 2.0, 2.8 );
  Rechit_Phi = new TH1F("Rechit_Phi"      , "Rechit #phi"   , 60, -3, 3. );
  Rechit_R = new TH1F("Rechit_R"      , "Rechit r"   , 30, 0, 150 );
  Rechit_Pos = new TH2F("Rechit_Pos"      , "Rechit x,y"   ,100,-100.,100., 100,-100.,100. );

  //  GenMuon_Eta = new TH1F("GenMuon_Eta"      , "GenMuon #eta"   , 4, 2.0, 2.8 );
  GenMuon_Phi = new TH1F("GenMuon_Phi"      , "GenMuon #phi"   , 60, -3, 3. );
  GenMuon_R = new TH1F("GenMuon_R"      , "GenMuon r"   , 30, 0, 150 );
  GenMuon_Pos = new TH2F("GenMuon_Pos"      , "GenMuon x,y"   ,100,-100.,100., 100,-100.,100. );

  ME0Muon_Eta = new TH1F("ME0Muon_Eta"      , "Muon #eta"   , 4, 2.0, 2.8 );
  ME0Muon_Cuts_Eta_5_10 = new TH1F("ME0Muon_Cuts_Eta_5_10"      , "Muon #eta"   , 4, 2.0, 2.8 );
  ME0Muon_Cuts_Eta_9_11 = new TH1F("ME0Muon_Cuts_Eta_9_11"      , "Muon #eta"   , 4, 2.0, 2.8 );
  ME0Muon_Cuts_Eta_10_50 = new TH1F("ME0Muon_Cuts_Eta_10_50"      , "Muon #eta"   , 4, 2.0, 2.8 );
  ME0Muon_Cuts_Eta_50_100 = new TH1F("ME0Muon_Cuts_Eta_50_100"      , "Muon #eta"   , 4, 2.0, 2.8 );
  ME0Muon_Cuts_Eta_100 = new TH1F("ME0Muon_Cuts_Eta_100"      , "Muon #eta"   , 4, 2.0, 2.8 );

  CheckME0Muon_Eta = new TH1F("CheckME0Muon_Eta"      , "Muon #eta"   , 4, 2.0, 2.8 );
  ME0Muon_Pt = new TH1F("ME0Muon_Pt"      , "Muon p_{T}"   , 120,0 , 120. );

  GenMuon_Eta = new TH1F("GenMuon_Eta"      , "Muon #eta"   , 4, 2.0, 2.8 );
  GenMuon_Eta_5_10 = new TH1F("GenMuon_Eta_5_10"      , "Muon #eta"   , 4, 2.0, 2.8 );
  GenMuon_Eta_9_11 = new TH1F("GenMuon_Eta_9_11"      , "Muon #eta"   , 4, 2.0, 2.8 );
  GenMuon_Eta_10_50 = new TH1F("GenMuon_Eta_10_50"      , "Muon #eta"   , 4, 2.0, 2.8 );
  GenMuon_Eta_50_100 = new TH1F("GenMuon_Eta_50_100"      , "Muon #eta"   , 4, 2.0, 2.8 );
  GenMuon_Eta_100 = new TH1F("GenMuon_Eta_100"      , "Muon #eta"   , 4, 2.0, 2.8 );

  TPMuon_Eta = new TH1F("TPMuon_Eta"      , "Muon #eta"   , 4, 2.0, 2.8 );

  GenMuon_Pt = new TH1F("GenMuon_Pt"      , "Muon p_{T}"   , 120,0 , 120. );

  MatchedME0Muon_Eta = new TH1F("MatchedME0Muon_Eta"      , "Muon #eta"   , 4, 2.0, 2.8 );
  StandaloneMatchedME0Muon_Eta = new TH1F("StandaloneMatchedME0Muon_Eta"      , "Muon #eta"   , 4, 2.0, 2.8 );
  StandaloneMatchedME0Muon_WideBinning_Eta = new TH1F("StandaloneMatchedME0Muon_WideBinning_Eta"      , "Muon #eta"   , 8, 2.0, 2.8 );
  StandaloneMatchedME0Muon_WidestBinning_Eta = new TH1F("StandaloneMatchedME0Muon_WidestBinning_Eta"      , "Muon #eta"   , 20, 2.0, 3.0 );
  MatchedME0Muon_Eta_5_10 = new TH1F("MatchedME0Muon_Eta_5_10"      , "Muon #eta"   , 4, 2.0, 2.8 );
  MatchedME0Muon_Eta_9_11 = new TH1F("MatchedME0Muon_Eta_9_11"      , "Muon #eta"   , 4, 2.0, 2.8 );
  MatchedME0Muon_Eta_10_50 = new TH1F("MatchedME0Muon_Eta_10_50"      , "Muon #eta"   , 4, 2.0, 2.8 );
  MatchedME0Muon_Eta_50_100 = new TH1F("MatchedME0Muon_Eta_50_100"      , "Muon #eta"   , 4, 2.0, 2.8 );
  MatchedME0Muon_Eta_100 = new TH1F("MatchedME0Muon_Eta_100"      , "Muon #eta"   , 4, 2.0, 2.8 );


  Chi2MatchedME0Muon_Eta_5_10 = new TH1F("Chi2MatchedME0Muon_Eta_5_10"      , "Muon #eta"   , 4, 2.0, 2.8 );
  Chi2MatchedME0Muon_Eta_9_11 = new TH1F("Chi2MatchedME0Muon_Eta_9_11"      , "Muon #eta"   , 4, 2.0, 2.8 );
  Chi2MatchedME0Muon_Eta_10_50 = new TH1F("Chi2MatchedME0Muon_Eta_10_50"      , "Muon #eta"   , 4, 2.0, 2.8 );
  Chi2MatchedME0Muon_Eta_50_100 = new TH1F("Chi2MatchedME0Muon_Eta_50_100"      , "Muon #eta"   , 4, 2.0, 2.8 );
  Chi2MatchedME0Muon_Eta_100 = new TH1F("Chi2MatchedME0Muon_Eta_100"      , "Muon #eta"   , 4, 2.0, 2.8 );

  MatchedME0Muon_Pt = new TH1F("MatchedME0Muon_Pt"      , "Muon p_{T}"   , 40,0 , 20 );

  Chi2MatchedME0Muon_Eta = new TH1F("Chi2MatchedME0Muon_Eta"      , "Muon #eta"   , 4, 2.0, 2.8 );
  Chi2MatchedME0Muon_Pt = new TH1F("Chi2MatchedME0Muon_Pt"      , "Muon p_{T}"   , 40,0 , 20 );

  Chi2UnmatchedME0Muon_Eta = new TH1F("Chi2UnmatchedME0Muon_Eta"      , "Muon #eta"   , 4, 2.0, 2.8 );

  UnmatchedME0Muon_Eta = new TH1F("UnmatchedME0Muon_Eta"      , "Muon #eta"   , 4, 2.0, 2.8 );
  UnmatchedME0Muon_Cuts_Eta_5_10 = new TH1F("UnmatchedME0Muon_Cuts_Eta_5_10"      , "Muon #eta"   , 4, 2.0, 2.8 );
  UnmatchedME0Muon_Cuts_Eta_9_11 = new TH1F("UnmatchedME0Muon_Cuts_Eta_9_11"      , "Muon #eta"   , 4, 2.0, 2.8 );
  UnmatchedME0Muon_Cuts_Eta_10_50 = new TH1F("UnmatchedME0Muon_Cuts_Eta_10_50"      , "Muon #eta"   , 4, 2.0, 2.8 );
  UnmatchedME0Muon_Cuts_Eta_50_100 = new TH1F("UnmatchedME0Muon_Cuts_Eta_50_100"      , "Muon #eta"   , 4, 2.0, 2.8 );
  UnmatchedME0Muon_Cuts_Eta_100 = new TH1F("UnmatchedME0Muon_Cuts_Eta_100"      , "Muon #eta"   , 4, 2.0, 2.8 );

  UnmatchedME0Muon_Pt = new TH1F("UnmatchedME0Muon_Pt"      , "Muon p_{T}"   , 500,0 , 50 );
  UnmatchedME0Muon_Window_Pt = new TH1F("UnmatchedME0Muon_Window_Pt"      , "Muon p_{T}"   , 500,0 , 50 );

  UnmatchedME0Muon_Cuts_Eta = new TH1F("UnmatchedME0Muon_Cuts_Eta"      , "Muon #eta"   , 4, 2.0, 2.8 );
  ME0Muon_Cuts_Eta = new TH1F("ME0Muon_Cuts_Eta"      , "Muon #eta"   , 4, 2.0, 2.8 );

  Mass_h = new TH1F("Mass_h"      , "Mass"   , 100, 0., 200 );

  MuonRecoEff_Eta = new TH1F("MuonRecoEff_Eta"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );

  MuonRecoEff_Eta_5_10 = new TH1F("MuonRecoEff_Eta_5_10"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );
  MuonRecoEff_Eta_9_11 = new TH1F("MuonRecoEff_Eta_9_11"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );
  MuonRecoEff_Eta_10_50 = new TH1F("MuonRecoEff_Eta_10_50"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );
  MuonRecoEff_Eta_50_100 = new TH1F("MuonRecoEff_Eta_50_100"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );
  MuonRecoEff_Eta_100 = new TH1F("MuonRecoEff_Eta_100"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );
  Chi2MuonRecoEff_Eta = new TH1F("Chi2MuonRecoEff_Eta"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );
  Chi2MuonRecoEff_Eta_5_10 = new TH1F("Chi2MuonRecoEff_Eta_5_10"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );
  Chi2MuonRecoEff_Eta_9_11 = new TH1F("Chi2MuonRecoEff_Eta_9_11"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );
  Chi2MuonRecoEff_Eta_10_50 = new TH1F("Chi2MuonRecoEff_Eta_10_50"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );
  Chi2MuonRecoEff_Eta_50_100 = new TH1F("Chi2MuonRecoEff_Eta_50_100"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );
  Chi2MuonRecoEff_Eta_100 = new TH1F("Chi2MuonRecoEff_Eta_100"      , "Fraction of ME0Muons matched to gen muons"   ,4, 2.0, 2.8  );

  MuonRecoEff_Pt = new TH1F("MuonRecoEff_Pt"      , "Fraction of ME0Muons matched to gen muons"   ,8, 0,40  );

  StandaloneMuonRecoEff_Eta = new TH1F("StandaloneMuonRecoEff_Eta"      , "Fraction of Standalone Muons matched to gen muons"   ,4, 2.0, 2.8  );
  StandaloneMuonRecoEff_WideBinning_Eta = new TH1F("StandaloneMuonRecoEff_WideBinning_Eta"      , "Fraction of Standalone Muons matched to gen muons"   ,8, 2.0, 2.8  );
  StandaloneMuonRecoEff_WidestBinning_Eta = new TH1F("StandaloneMuonRecoEff_WidestBinning_Eta"      , "Fraction of Standalone Muons matched to gen muons"   ,20, 2.0, 3.0  );

  FakeRate_Eta = new TH1F("FakeRate_Eta"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,4, 2.0, 2.8  );
  FakeRate_Eta_5_10 = new TH1F("FakeRate_Eta_5_10"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,4, 2.0, 2.8  );
  FakeRate_Eta_9_11 = new TH1F("FakeRate_Eta_9_11"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,4, 2.0, 2.8  );
  FakeRate_Eta_10_50 = new TH1F("FakeRate_Eta_10_50"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,4, 2.0, 2.8  );
  FakeRate_Eta_50_100 = new TH1F("FakeRate_Eta_50_100"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,4, 2.0, 2.8  );
  FakeRate_Eta_100 = new TH1F("FakeRate_Eta_100"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,4, 2.0, 2.8  );

  Chi2FakeRate_Eta = new TH1F("Chi2FakeRate_Eta"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,4, 2.0, 2.8  );

  FakeRate_Eta_PerEvent = new TH1F("FakeRate_Eta_PerEvent"      , "PU140, unmatched ME0Muons/all ME0Muons normalized by N_{events}"   ,4, 2.0, 2.8  );
  FakeRate_Pt = new TH1F("FakeRate_Pt"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,8, 0,40  );

  MuonAllTracksEff_Eta = new TH1F("MuonAllTracksEff_Eta"      , "All ME0Muons over all tracks"   ,4, 2.0, 2.8  );
  MuonAllTracksEff_Pt = new TH1F("MuonAllTracksEff_Pt"      , "All ME0Muons over all tracks"   ,8, 0,40  );

  MuonUnmatchedTracksEff_Eta = new TH1F("MuonUnmatchedTracksEff_Eta"      , "Unmatched ME0Muons over all ME0Muons"   ,4, 2.0, 2.8  );
  MuonUnmatchedTracksEff_Pt = new TH1F("MuonUnmatchedTracksEff_Pt"      , "Unmatched ME0Muons over all ME0Muons"   ,8, 0,40  );

  TracksPerSegment_h = new TH1F("TracksPerSegment_h", "Number of tracks", 60,0.,60.);
  TracksPerSegment_s = new TH2F("TracksPerSegment_s" , "Tracks per segment vs |#eta|", 4, 2.0, 2.8, 60,0.,60.);
  TracksPerSegment_p = new TProfile("TracksPerSegment_p" , "Tracks per segment vs |#eta|", 4, 2.0, 2.8, 0.,60.);

  FakeTracksPerSegment_h = new TH1F("FakeTracksPerSegment_h", "Number of fake tracks", 60,0.,60.);
  FakeTracksPerSegment_s = new TH2F("FakeTracksPerSegment_s" , "Fake tracks per segment", 10, 2.4, 4.0, 100,0.,60.);
  FakeTracksPerSegment_p = new TProfile("FakeTracksPerSegment_p" , "Average N_{tracks}/segment not matched to genmuons", 10, 2.4, 4.0, 0.,60.);

  FakeTracksPerAssociatedSegment_h = new TH1F("FakeTracksPerAssociatedSegment_h", "Number of fake tracks", 60,0.,60.);
  FakeTracksPerAssociatedSegment_s = new TH2F("FakeTracksPerAssociatedSegment_s" , "Fake tracks per segment", 10, 2.4, 4.0, 100,0.,60.);
  FakeTracksPerAssociatedSegment_p = new TProfile("FakeTracksPerAssociatedSegment_p" , "Average N_{tracks}/segment not matched to genmuons", 10, 2.4, 4.0, 0.,60.);

  ClosestDelR_s = new TH2F("ClosestDelR_s" , "#Delta R", 4, 2.0, 2.8, 15,0.,0.15);
  ClosestDelR_p = new TProfile("ClosestDelR_p" , "#Delta R", 4, 2.0, 2.8, 0.,0.15);

  DelR_Window_Under5 = new TH1F("DelR_Window_Under5","#Delta R", 15, 0,0.15  );
  Pt_Window_Under5 = new TH1F("Pt_Window_Under5","pt",500, 0,50  );

  DelR_Track_Window_Under5 = new TH1F("DelR_Track_Window_Under5","#Delta R", 15, 0,0.15  );
  Pt_Track_Window_Under5 = new TH1F("Pt_Track_Window_Under5","pt",20, 0,5  );
  Pt_Track_Window = new TH1F("Pt_Track_Window","pt",500, 0,  50);

  DelR_Track_Window_Failed_Under5 = new TH1F("DelR_Track_Window_Failed_Under5","#Delta R", 15, 0,0.15  );
  Pt_Track_Window_Failed_Under5 = new TH1F("Pt_Track_Window_Failed_Under5","pt",20, 0,5  );
  Pt_Track_Window_Failed = new TH1F("Pt_Track_Window_Failed","pt",500, 0,  50);

  FailedTrack_Window_XPull = new TH1F("FailedTrack_Window_XPull", "X Pull failed tracks", 100, 0,20);
  FailedTrack_Window_YPull = new TH1F("FailedTrack_Window_YPull", "Y  Pull failed tracks", 100, 0,20);
  FailedTrack_Window_XDiff = new TH1F("FailedTrack_Window_XDiff", "X Diff failed tracks", 100, 0,20);
  FailedTrack_Window_YDiff = new TH1F("FailedTrack_Window_YDiff", "Y  Diff failed tracks", 100, 0,20);

  FailedTrack_Window_PhiDiff = new TH1F("FailedTrack_Window_PhiDiff", "Phi Dir Diff failed tracks", 100,0 ,2.0);

  DelR_Segment_GenMuon = new TH1F("DelR_Segment_GenMuon", "#Delta R between me0segment and gen muon",200,0,2);
  FractionMatched_Eta = new TH1F("FractionMatched_Eta"      , "Fraction of ME0Muons that end up successfully matched (matched/all)"   ,4, 2.0, 2.8  );

  PtDiff_s = new TH2F("PtDiff_s" , "Relative pt difference", 4, 2.0, 2.8, 200,-1,1.0);

  PtDiff_s_5_10 = new TH2F("PtDiff_s_5_10" , "Relative pt difference", 4, 2.0, 2.8, 200,-1,1.0);
  PtDiff_s_10_50 = new TH2F("PtDiff_s_10_50" , "Relative pt difference", 4, 2.0, 2.8, 200,-1,1.0);
  PtDiff_s_50_100 = new TH2F("PtDiff_s_50_100" , "Relative pt difference", 4, 2.0, 2.8, 200,-1,1.0);
  PtDiff_s_100 = new TH2F("PtDiff_s_100" , "Relative pt difference", 4, 2.0, 2.8, 200,-1,1.0);

  PtDiff_h = new TH1F("PtDiff_h" , "pt resolution", 100,-0.5,0.5);
  QOverPtDiff_h = new TH1F("QOverPtDiff_h" , "q/pt resolution", 100,-0.5,0.5);
  PtDiff_p = new TProfile("PtDiff_p" , "pt resolution vs. #eta", 4, 2.0, 2.8, -1.0,1.0,"s");

  StandalonePtDiff_s = new TH2F("StandalonePtDiff_s" , "Relative pt difference", 4, 2.0, 2.8, 200,-1,1.0);
  StandalonePtDiff_h = new TH1F("StandalonePtDiff_h" , "pt resolution", 100,-0.5,0.5);
  StandaloneQOverPtDiff_h = new TH1F("StandaloneQOverPtDiff_h" , "q/pt resolution", 100,-0.5,0.5);
  StandalonePtDiff_p = new TProfile("StandalonePtDiff_p" , "pt resolution vs. #eta", 4, 2.0, 2.8, -1.0,1.0,"s");

  PtDiff_rms    = new TH1F( "PtDiff_rms",    "RMS", 4, 2.0, 2.8 ); 
  PtDiff_gaus_wide    = new TH1F( "PtDiff_gaus_wide",    "GAUS_WIDE", 4, 2.0, 2.8 ); 
  PtDiff_gaus_narrow    = new TH1F( "PtDiff_gaus_narrow",    "GAUS_NARROW", 4, 2.0, 2.8 ); 

  PtDiff_gaus_5_10    = new TH1F( "PtDiff_gaus_5_10",    "GAUS_WIDE", 4, 2.0, 2.8 ); 
  PtDiff_gaus_10_50    = new TH1F( "PtDiff_gaus_10_50",    "GAUS_WIDE", 4, 2.0, 2.8 ); 
  PtDiff_gaus_50_100    = new TH1F( "PtDiff_gaus_50_100",    "GAUS_WIDE", 4, 2.0, 2.8 ); 
  PtDiff_gaus_100    = new TH1F( "PtDiff_gaus_100",    "GAUS_WIDE", 4, 2.0, 2.8 ); 

  StandalonePtDiff_gaus    = new TH1F( "StandalonePtDiff_gaus",    "GAUS_WIDE", 4, 2.0, 2.8 ); 

  PDiff_s = new TH2F("PDiff_s" , "Relative p difference", 4, 2.0, 2.8, 50,0.,0.5);
  PDiff_h = new TH1F("PDiff_s" , "Relative p difference", 50,0.,0.5);
  PDiff_p = new TProfile("PDiff_p" , "Relative p difference", 4, 2.0, 2.8, 0.,1.0,"s");

  VertexDiff_h = new TH1F("VertexDiff_h", "Difference in vertex Z", 50, 0, 0.2);

  SegPosDirPhiDiff_True_h = new TH1F("SegPosDirPhiDiff_True_h", "#phi Dir. Diff. Real Muons", 50, -2,2);
  SegPosDirEtaDiff_True_h = new TH1F("SegPosDirEtaDiff_True_h", "#eta Dir. Diff. Real Muons", 50, -2,2);

  SegPosDirPhiDiff_All_h = new TH1F("SegPosDirPhiDiff_All_h", "#phi Dir. Diff. All Muons", 50, -3,3);
  SegPosDirEtaDiff_All_h = new TH1F("SegPosDirEtaDiff_All_h", "#eta Dir. Diff. All Muons", 50, -3,3);

  SegTrackDirPhiDiff_True_h = new TH1F("SegTrackDirPhiDiff_True_h", "#phi Dir. Diff. Real Muons", 50, -2,2);
  SegTrackDirEtaDiff_True_h = new TH1F("SegTrackDirEtaDiff_True_h", "#eta Dir. Diff. Real Muons", 50, -2,2);

  SegTrackDirPhiPull_True_h = new TH1F("SegTrackDirPhiPull_True_h", "#phi Dir. Pull. Real Muons", 50, -3,3);
  SegTrackDirPhiPull_All_h = new TH1F("SegTrackDirPhiPull_True_h", "#phi Dir. Pull. All Muons", 50, -3,3);

  SegTrackDirPhiDiff_All_h = new TH1F("SegTrackDirPhiDiff_All_h", "#phi Dir. Diff. All Muons", 50, -3,3);
  SegTrackDirEtaDiff_All_h = new TH1F("SegTrackDirEtaDiff_All_h", "#eta Dir. Diff. All Muons", 50, -3,3);

  SegGenDirPhiDiff_True_h = new TH1F("SegGenDirPhiDiff_True_h", "#phi Dir. Diff. Real Muons", 50, -2,2);
  SegGenDirEtaDiff_True_h = new TH1F("SegGenDirEtaDiff_True_h", "#eta Dir. Diff. Real Muons", 50, -2,2);

  SegGenDirPhiPull_True_h = new TH1F("SegGenDirPhiPull_True_h", "#phi Dir. Pull. Real Muons", 50, -3,3);
  SegGenDirPhiPull_All_h = new TH1F("SegGenDirPhiPull_True_h", "#phi Dir. Pull. All Muons", 50, -3,3);

  SegGenDirPhiDiff_All_h = new TH1F("SegGenDirPhiDiff_All_h", "#phi Dir. Diff. All Muons", 50, -3,3);
  SegGenDirEtaDiff_All_h = new TH1F("SegGenDirEtaDiff_All_h", "#eta Dir. Diff. All Muons", 50, -3,3);


  PreMatch_TP_R = new TH1F("PreMatch_TP_R", "r distance from TP pre match to beamline", 100, 0, 10);
  PostMatch_TP_R = new TH1F("PostMatch_TP_R", "r distance from TP post match to beamline", 200, 0, 20);
  PostMatch_BX0_TP_R = new TH1F("PostMatch_BX0_TP_R", "r distance from TP post match to beamline", 200, 0, 20);


  XDiff_h = new TH1F("XDiff_h", "X Diff", 100, -10.0, 10.0 );
  YDiff_h = new TH1F("YDiff_h", "Y Diff", 100, -50.0, 50.0 ); 
  XPull_h = new TH1F("XPull_h", "X Pull", 100, -5.0, 5.0 );
  YPull_h = new TH1F("YPull_h", "Y Pull", 40, -50.0, 50.0 );

  MuonRecoEff_WideBinning_Eta = new TH1F("MuonRecoEff_WideBinning_Eta"      , "Fraction of ME0Muons matched to gen muons"   ,8, 2.0, 2.8  );
  MuonRecoEff_WidestBinning_Eta = new TH1F("MuonRecoEff_WidestBinning_Eta"      , "Fraction of ME0Muons matched to gen muons"   ,20, 2.0, 3.0  );
  Chi2MuonRecoEff_WideBinning_Eta = new TH1F("Chi2MuonRecoEff_WideBinning_Eta"      , "Fraction of ME0Muons matched to gen muons"   ,8, 2.0, 2.8  );
  Chi2MuonRecoEff_WidestBinning_Eta = new TH1F("Chi2MuonRecoEff_WidestBinning_Eta"      , "Fraction of ME0Muons matched to gen muons"   ,20, 2.0, 3.0  );
  Chi2FakeRate_WideBinning_Eta = new TH1F("Chi2FakeRate_WideBinning_Eta"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,8, 2.0, 2.8  );
  Chi2FakeRate_WidestBinning_Eta = new TH1F("Chi2FakeRate_WidestBinning_Eta"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,20, 2.0, 3.0  );
  FakeRate_WideBinning_Eta = new TH1F("FakeRate_WideBinning_Eta"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,8, 2.0, 2.8  );
  FakeRate_WidestBinning_Eta = new TH1F("FakeRate_WidestBinning_Eta"      , "PU140, unmatched ME0Muons/all ME0Muons"   ,20, 2.0, 3.0  );

  UnmatchedME0Muon_Cuts_WideBinning_Eta = new TH1F("UnmatchedME0Muon_Cuts_WideBinning_Eta"      , "Muon #eta"   , 8, 2.0, 2.8 );
  UnmatchedME0Muon_Cuts_WidestBinning_Eta = new TH1F("UnmatchedME0Muon_Cuts_WidestBinning_Eta"      , "Muon #eta"   , 20, 2.0, 3.0 );
  ME0Muon_Cuts_WideBinning_Eta = new TH1F("ME0Muon_Cuts_WideBinning_Eta"      , "Muon #eta"   , 8, 2.0, 2.8 );
  ME0Muon_Cuts_WidestBinning_Eta = new TH1F("ME0Muon_Cuts_WidestBinning_Eta"      , "Muon #eta"   , 20, 2.0, 3.0 );
  Chi2UnmatchedME0Muon_WideBinning_Eta = new TH1F("Chi2UnmatchedME0Muon_WideBinning_Eta"      , "Muon #eta"   , 8, 2.0, 2.8 );
  Chi2UnmatchedME0Muon_WidestBinning_Eta = new TH1F("Chi2UnmatchedME0Muon_WidestBinning_Eta"      , "Muon #eta"   , 20, 2.0, 3.0 );
  TPMuon_WideBinning_Eta = new TH1F("TPMuon_WideBinning_Eta"      , "Muon #eta"   , 8, 2.0, 2.8 );
  TPMuon_WidestBinning_Eta = new TH1F("TPMuon_WidestBinning_Eta"      , "Muon #eta"   , 20, 2.0, 3.0 );
  GenMuon_WideBinning_Eta = new TH1F("GenMuon_WideBinning_Eta"      , "Muon #eta"   , 8, 2.0, 2.8 );
  GenMuon_WidestBinning_Eta = new TH1F("GenMuon_WidestBinning_Eta"      , "Muon #eta"   , 20, 2.0, 3.0 );
  MatchedME0Muon_WideBinning_Eta = new TH1F("MatchedME0Muon_WideBinning_Eta"      , "Muon #eta"   , 8, 2.0, 2.8 );
  MatchedME0Muon_WidestBinning_Eta = new TH1F("MatchedME0Muon_WidestBinning_Eta"      , "Muon #eta"   , 20, 2.0, 3.0 );
  Chi2MatchedME0Muon_WideBinning_Eta = new TH1F("Chi2MatchedME0Muon_WideBinning_Eta"      , "Muon #eta"   , 8, 2.0, 2.8 );
  Chi2MatchedME0Muon_WidestBinning_Eta = new TH1F("Chi2MatchedME0Muon_WidestBinning_Eta"      , "Muon #eta"   , 20, 2.0, 3.0 );

 
  AssociatedChi2_h = new TH1F("AssociatedChi2_h","Associated #chi^{2}",50,0,50);
  AssociatedChi2_Prob_h = new TH1F("AssociatedChi2_h","Associated #chi^{2}",50,0,1);
  NormChi2_h = new TH1F("NormChi2_h","normalized #chi^{2}", 200, 0, 20);
  NormChi2Prob_h = new TH1F("NormChi2Prob_h","normalized #chi^{2} probability", 100, 0, 1);
  NormChi2VsHits_h = new TH2F("NormChi2VsHits_h","#chi^{2} vs nhits",25,0,25,100,0,10);
  chi2_vs_eta_h = new TH2F("chi2_vs_eta_h","#chi^{2} vs #eta",4, 2.0, 2.8 , 200, 0, 20);

  Nevents=0;
  std::cout<<"HERE NOW, about to check if get associator"<<std::endl;
  if (UseAssociators) {
    std::cout<<"Getting the associator"<<std::endl;
    edm::ESHandle<TrackAssociatorBase> theAssociator;
    for (unsigned int w=0;w<associators.size();w++) {
      iSetup.get<TrackAssociatorRecord>().get(associators[w],theAssociator);
      associator.push_back( theAssociator.product() );
    }
  }
  std::cout<<"HERE NOW"<<std::endl;

}


ME0MuonAnalyzer::~ME0MuonAnalyzer(){}

void
ME0MuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{

  std::cout<<"ANALYZER"<<std::endl;
  
  using namespace edm;


  using namespace reco;

  Handle<GenParticleCollection> genParticles;

  iEvent.getByLabel<GenParticleCollection>("genParticles", genParticles);

  unsigned int gensize=genParticles->size();




  Nevents++;

  std::cout<<"About to get muons:"<<std::endl;
  edm::Handle<std::vector<Muon> > muons;
  //iEvent.getByLabel("Muons", muons);
  iEvent.getByLabel("muons", muons);
  
  std::cout<<"Have muons, about to start"<<std::endl;
  for(unsigned int i=0; i<gensize; ++i) {
    const reco::GenParticle& CurrentParticle=(*genParticles)[i];
    if ( (CurrentParticle.status()==1) && ( (CurrentParticle.pdgId()==13)  || (CurrentParticle.pdgId()==-13) ) ){  

      double LowestDelR = 9999;
      double thisDelR = 9999;
      
      std::vector<double> ReferenceTrackPt;

      if (CurrentParticle.pt() > FakeRatePtCut ) {
	GenMuon_Eta->Fill(fabs(CurrentParticle.eta()));
	GenMuon_WideBinning_Eta->Fill(fabs(CurrentParticle.eta()));
	GenMuon_WidestBinning_Eta->Fill(fabs(CurrentParticle.eta()));
      }
      //double VertexDiff=-1,PtDiff=-1,QOverPtDiff=-1,PDiff=-1,MatchedEta=-1;
      double PtDiff=-1,QOverPtDiff=-1,MatchedEta=-1;

      std::cout<<"Size = "<<muons->size()<<std::endl;
      for (std::vector<Muon>::const_iterator thisMuon = muons->begin();
  	   thisMuon != muons->end(); ++thisMuon){

  	if (thisMuon->isME0Muon()){
	  //if (thisMuon->isTrackerMuon()){
	  
  	  TrackRef tkRef = thisMuon->innerTrack();
	  if (tkRef.isNull() ) continue;
  	  thisDelR = reco::deltaR(CurrentParticle,*tkRef);
  	  if (CurrentParticle.pt() > FakeRatePtCut ) {
  	    if (thisDelR < MatchingWindowDelR ){
  	      if (thisDelR < LowestDelR){
  		LowestDelR = thisDelR;

  		MatchedEta=fabs(CurrentParticle.eta());
  		//VertexDiff = fabs(tkRef->vz()-CurrentParticle.vz());
  		QOverPtDiff = ( (tkRef->charge() /tkRef->pt()) - (CurrentParticle.charge()/CurrentParticle.pt() ) )/  (CurrentParticle.charge()/CurrentParticle.pt() );
  		PtDiff = (tkRef->pt() - CurrentParticle.pt())/CurrentParticle.pt();
  		//PDiff = (tkRef->p() - CurrentParticle.p())/CurrentParticle.p();
  	      }
  	    }
  	  }
  	}
      }
      StandaloneMatchedME0Muon_Eta->Fill(MatchedEta);
      StandaloneMatchedME0Muon_WideBinning_Eta->Fill(MatchedEta);
      StandaloneMatchedME0Muon_WidestBinning_Eta->Fill(MatchedEta);

      //StandaloneVertexDiff_h->Fill(VertexDiff);
      StandalonePtDiff_h->Fill(PtDiff);	
      StandaloneQOverPtDiff_h->Fill(QOverPtDiff);
      StandalonePtDiff_s->Fill(CurrentParticle.eta(),PtDiff);

    }
  }

  // for (std::vector<ME0Muon>::const_iterator thisMuon = OurMuons->begin();
  //      thisMuon != OurMuons->end(); ++thisMuon){
  //   TrackRef tkRef = thisMuon->innerTrack();
  //   ME0Muon_Eta->Fill(tkRef->eta());
  //   if ( (TMath::Abs(tkRef->eta()) > 2.4) && (TMath::Abs(tkRef->eta()) < 4.0) ) ME0Muon_Pt->Fill(tkRef->pt());
  // }
  
  // std::vector<double> SegmentEta, SegmentPhi, SegmentR, SegmentX, SegmentY;
  // // std::vector<const ME0Segment*> Ids;
  // // std::vector<const ME0Segment*> Ids_NonGenMuons;
  // // std::vector<const ME0Segment*> UniqueIdList;
  //  std::vector<int> Ids;
  //  std::vector<int> Ids_NonGenMuons;
  //  std::vector<int> UniqueIdList;
  //  int TrackID=0;

  //  //std::cout<<"Doing some propagation"<<std::endl;
  //  int MuID = 0;
  //  for (std::vector<ME0Muon>::const_iterator thisMuon = OurMuons->begin();
  // 	thisMuon != OurMuons->end(); ++thisMuon){
  //   if (!muon::isGoodMuon(me0Geom, *thisMuon, muon::Tight)) continue;
  //   TrackRef tkRef = thisMuon->innerTrack();
  //   //ME0Segment segRef = thisMuon->me0segment();
  //   //const ME0Segment* SegId = segRef->get();

  //   ME0Segment Seg = thisMuon->me0segment();
  //   ME0DetId id =Seg.me0DetId();
  //   auto roll = me0Geom->etaPartition(id); 
  //   auto GlobVect(roll->toGlobal(Seg.localPosition()));

  //   int SegId=thisMuon->me0segid();

  //   //std::cout<<SegId<<std::endl;


  //   //For a direction study...
  //   // if ( (tkRef->pt() > 5.0) && (IsMatched[MuID]) ){
  //   //   SegPosDirPhiDiff_True_h->Fill(roll->toGlobal(Seg.localPosition()).phi()-roll->toGlobal(Seg.localDirection()).phi() );
  //   //   SegPosDirEtaDiff_True_h->Fill(roll->toGlobal(Seg.localPosition()).eta()-roll->toGlobal(Seg.localDirection()).eta() );
  //   // }

  //   // SegPosDirPhiDiff_All_h->Fill(roll->toGlobal(Seg.localPosition()).phi()-roll->toGlobal(Seg.localDirection()).phi() );
  //   // SegPosDirEtaDiff_All_h->Fill(roll->toGlobal(Seg.localPosition()).eta()-roll->toGlobal(Seg.localDirection()).eta() );

  //   // // For another direction study...
  //   // float zSign  = tkRef->pz()/fabs(tkRef->pz());

  //   // //float zValue = 560. * zSign;
  //   // float zValue = 526.75 * zSign;
  //   // Plane *plane = new Plane(Surface::PositionType(0,0,zValue),Surface::RotationType());
  //   // //Getting the initial variables for propagation
  //   // int chargeReco = tkRef->charge(); 
  //   // GlobalVector p3reco, r3reco;

  //   // p3reco = GlobalVector(tkRef->outerPx(), tkRef->outerPy(), tkRef->outerPz());
  //   // r3reco = GlobalVector(tkRef->outerX(), tkRef->outerY(), tkRef->outerZ());

  //   // AlgebraicSymMatrix66 covReco;
  //   // //This is to fill the cov matrix correctly
  //   // AlgebraicSymMatrix55 covReco_curv;
  //   // covReco_curv = tkRef->outerStateCovariance();
  //   // FreeTrajectoryState initrecostate = getFTS(p3reco, r3reco, chargeReco, covReco_curv, &*bField);
  //   // getFromFTS(initrecostate, p3reco, r3reco, chargeReco, covReco);

  //   // //Now we propagate and get the propagated variables from the propagated state
  //   // SteppingHelixStateInfo startrecostate(initrecostate);
  //   // SteppingHelixStateInfo lastrecostate;

  //   // const SteppingHelixPropagator* ThisshProp = 
  //   //   dynamic_cast<const SteppingHelixPropagator*>(&*shProp);
	
  //   // lastrecostate = ThisshProp->propagate(startrecostate, *plane);
	
  //   // FreeTrajectoryState finalrecostate;
  //   // lastrecostate.getFreeState(finalrecostate);
      
  //   // AlgebraicSymMatrix66 covFinalReco;
  //   // GlobalVector p3FinalReco_glob, r3FinalReco_globv;
  //   // getFromFTS(finalrecostate, p3FinalReco_glob, r3FinalReco_globv, chargeReco, covFinalReco);
  //   // GlobalPoint r3FinalReco_glob(r3FinalReco_globv.x(),r3FinalReco_globv.y(),r3FinalReco_globv.z());

  //   // double DirectionPull, DirectionPullNum, DirectionPullDenom;

  //   // //Computing the sigma for the track direction
  //   // Double_t mag_track = p3FinalReco_glob.perp();
  //   // //Double_t phi_track = p3FinalReco_glob.phi();

  //   // //Double_t dmagdx_track = p3FinalReco_glob.x()/mag_track;
  //   // //Double_t dmagdy_track = p3FinalReco_glob.y()/mag_track;
  //   // Double_t dphidx_track = -p3FinalReco_glob.y()/(mag_track*mag_track);
  //   // Double_t dphidy_track = p3FinalReco_glob.x()/(mag_track*mag_track);
  //   // Double_t sigmaphi_track = sqrt( dphidx_track*dphidx_track*covFinalReco(3,3)+
  //   // 				dphidy_track*dphidy_track*covFinalReco(4,4)+
  //   // 				dphidx_track*dphidy_track*2*covFinalReco(3,4) );

  //   // DirectionPullNum = p3FinalReco_glob.phi()-roll->toGlobal(Seg.localDirection()).phi();
  //   // DirectionPullDenom = sqrt( pow(roll->toGlobal(Seg.localPosition()).phi(),2) + pow(sigmaphi_track,2) );
  //   // DirectionPull = DirectionPullNum / DirectionPullDenom;

  //   // if ( (tkRef->pt() > 5.0)&& (IsMatched[MuID]) ){
  //   //   SegTrackDirPhiDiff_True_h->Fill(p3FinalReco_glob.phi()-roll->toGlobal(Seg.localDirection()).phi() );
  //   //   SegTrackDirEtaDiff_True_h->Fill(p3FinalReco_glob.eta()-roll->toGlobal(Seg.localDirection()).eta() );
  //   //   SegTrackDirPhiPull_True_h->Fill(DirectionPull);
  //   // }
  //   // SegTrackDirPhiDiff_All_h->Fill(p3FinalReco_glob.phi()-roll->toGlobal(Seg.localDirection()).phi() );
  //   // SegTrackDirPhiPull_All_h->Fill(DirectionPull);
    
  //   // SegTrackDirEtaDiff_All_h->Fill(p3FinalReco_glob.eta()-roll->toGlobal(Seg.localDirection()).eta() );


  //   // LocalPoint r3FinalReco = roll->toLocal(r3FinalReco_glob);
  //   // LocalVector p3FinalReco=roll->toLocal(p3FinalReco_glob);
  //   // LocalTrajectoryParameters ltp(r3FinalReco,p3FinalReco,chargeReco);
  //   // JacobianCartesianToLocal jctl(roll->surface(),ltp);
  //   // AlgebraicMatrix56 jacobGlbToLoc = jctl.jacobian(); 

  //   // AlgebraicMatrix55 Ctmp =  (jacobGlbToLoc * covFinalReco) * ROOT::Math::Transpose(jacobGlbToLoc); 
  //   // AlgebraicSymMatrix55 C;  // I couldn't find any other way, so I resort to the brute force
  //   // for(int i=0; i<5; ++i) {
  //   //   for(int j=0; j<5; ++j) {
  //   // 	C[i][j] = Ctmp[i][j]; 

  //   //   }
  //   // }  

  //   // LocalPoint thisPosition(Seg.localPosition());

  //   // Double_t sigmax = sqrt(C[3][3]+Seg.localPositionError().xx() );      
  //   // Double_t sigmay = sqrt(C[4][4]+Seg.localPositionError().yy() );

  //   // XPull_h->Fill((thisPosition.x()-r3FinalReco.x())/sigmax);
  //   // YPull_h->Fill((thisPosition.y()-r3FinalReco.y())/sigmay);
    
  //   // XDiff_h->Fill((thisPosition.x()-r3FinalReco.x()));
  //   // YDiff_h->Fill((thisPosition.y()-r3FinalReco.y()));

  //   // //std::cout<<"AM HERE"<<std::endl;
  //   // if ( (tkRef->pt() > FakeRatePtCut)&& (IsMatched[MuID]) ){
      


  //   //   //std::cout<<"thisPosition = "<<thisPosition<<std::endl;
  //   //   //std::cout<<"r3FinalReco = "<<r3FinalReco<<std::endl;
  //   // }

  //   //End Direction studies


  //   bool IsNew = true;
  //   for (unsigned int i =0; i < Ids.size(); i++){
  //     if (SegId == Ids[i]) IsNew=false;
  //   }

  //   if (IsNew) {
  //     UniqueIdList.push_back(SegId);
  //     //std::cout<<"New SegId = "<<SegId<<std::endl;
  //     //std::cout<<GlobVect<<std::endl;
  //     SegmentEta.push_back(GlobVect.eta());
  //     SegmentPhi.push_back(GlobVect.phi());
  //     SegmentR.push_back(GlobVect.perp());
  //     SegmentX.push_back(GlobVect.x());
  //     SegmentY.push_back(GlobVect.y());
  //   }
  //   Ids.push_back(SegId);
  //   if (!IsMatched[TrackID]) Ids_NonGenMuons.push_back(SegId);

  //   ME0Muon_Eta->Fill(fabs(tkRef->eta()));


  //   if ((tkRef->pt() > FakeRatePtCut) && (TMath::Abs(tkRef->eta()) < 2.8)){
  //     ME0Muon_Cuts_Eta->Fill(fabs(tkRef->eta()));
  //     ME0Muon_Cuts_WideBinning_Eta->Fill(fabs(tkRef->eta()));
  //     ME0Muon_Cuts_WidestBinning_Eta->Fill(fabs(tkRef->eta()));
  //     if ( (tkRef->pt() > 5.0) && (tkRef->pt() <= 10.0) )  	ME0Muon_Cuts_Eta_5_10->Fill(fabs(tkRef->eta()));
  //     if ( (tkRef->pt() > 9.0) && (tkRef->pt() <= 11.0) )  	ME0Muon_Cuts_Eta_9_11->Fill(fabs(tkRef->eta()));
  //     if ( (tkRef->pt() > 10.0) && (tkRef->pt() <= 20.0) )	ME0Muon_Cuts_Eta_10_50->Fill(fabs(tkRef->eta()));
  //     if ( (tkRef->pt() > 20.0) && (tkRef->pt() <= 40.0) )	ME0Muon_Cuts_Eta_50_100->Fill(fabs(tkRef->eta()));
  //     if ( tkRef->pt() > 40.0) 		ME0Muon_Cuts_Eta_100->Fill(fabs(tkRef->eta()));
  //   }
    


  //   if ( (TMath::Abs(tkRef->eta()) > 2.0) && (TMath::Abs(tkRef->eta()) < 2.8) ) ME0Muon_Pt->Fill(tkRef->pt());

  //   TrackID++;
  //   MuID++;
  // }
  
  //  //std::cout<<UniqueIdList.size()<<" unique segments per event"<<std::endl;
  // for (unsigned int i = 0; i < UniqueIdList.size(); i++){
  //   int Num_Total=0, Num_Fake = 0, Num_Fake_Associated = 0;
  //   for (unsigned int j = 0; j < Ids.size(); j++){
  //     if (Ids[j] == UniqueIdList[i]) Num_Total++;
  //   }

  //   for (unsigned int j = 0; j < Ids_NonGenMuons.size(); j++){
  //     if (Ids_NonGenMuons[j] == UniqueIdList[i]) Num_Fake++;
  //     bool AssociatedWithMatchedSegment = false;
  //     for (unsigned int isegid=0;isegid < MatchedSegIds.size();isegid++){
  // 	if (MatchedSegIds[isegid]==Ids_NonGenMuons[j]) AssociatedWithMatchedSegment=true;
  //     }
  //     if (AssociatedWithMatchedSegment) Num_Fake_Associated++;
  //   }

  //   TracksPerSegment_h->Fill((double)Num_Total);
  //   TracksPerSegment_s->Fill(SegmentEta[i], (double)Num_Total);
  //   TracksPerSegment_p->Fill(SegmentEta[i], (double)Num_Total);

  //   FakeTracksPerSegment_h->Fill((double)Num_Fake);
  //   FakeTracksPerSegment_s->Fill(SegmentEta[i], (double)Num_Fake);
  //   FakeTracksPerSegment_p->Fill(SegmentEta[i], (double)Num_Fake);

  //   FakeTracksPerAssociatedSegment_h->Fill((double)Num_Fake_Associated);
  //   FakeTracksPerAssociatedSegment_s->Fill(SegmentEta[i], (double)Num_Fake_Associated);
  //   FakeTracksPerAssociatedSegment_p->Fill(SegmentEta[i], (double)Num_Fake_Associated);

  //   // if (SegmentEta[i] > 2.4){
  //   //   Segment_Eta->Fill(SegmentEta[i]);
  //   //   Segment_Phi->Fill(SegmentPhi[i]);
  //   //   Segment_R->Fill(SegmentR[i]);
  //   //   Segment_Pos->Fill(SegmentX[i],SegmentY[i]);
  //   // }
  // }

  // //================  For Segment Plotting
  // for (auto thisSegment = OurSegments->begin(); thisSegment != OurSegments->end(); 
  //      ++thisSegment){
  //   ME0DetId id = thisSegment->me0DetId();
  //   //std::cout<<"ME0DetId =  "<<id<<std::endl;
  //   auto roll = me0Geom->etaPartition(id); 
  //   auto GlobVect(roll->toGlobal(thisSegment->localPosition()));
  //   Segment_Eta->Fill(fabs(GlobVect.eta()));
  //   Segment_Phi->Fill(GlobVect.phi());
  //   Segment_R->Fill(GlobVect.perp());
  //   Segment_Pos->Fill(GlobVect.x(),GlobVect.y());



  //   auto theseRecHits = thisSegment->specificRecHits();
  //   //std::cout <<"ME0 Ensemble Det Id "<<id<<"  Number of RecHits "<<theseRecHits.size()<<std::endl;
    
  //   for (auto thisRecHit = theseRecHits.begin(); thisRecHit!= theseRecHits.end(); thisRecHit++){
  //     auto me0id = thisRecHit->me0Id();
  //     auto rollForRechit = me0Geom->etaPartition(me0id);
      
  //     auto thisRecHitGlobalPoint = rollForRechit->toGlobal(thisRecHit->localPosition()); 
      
  //     Rechit_Eta->Fill(fabs(thisRecHitGlobalPoint.eta()));
  //     Rechit_Phi->Fill(thisRecHitGlobalPoint.phi());
  //     Rechit_R->Fill(thisRecHitGlobalPoint.perp());
  //     Rechit_Pos->Fill(thisRecHitGlobalPoint.x(),thisRecHitGlobalPoint.y());
      
  //   }
  // }
  // //==================

  //   for(unsigned int i=0; i<gensize; ++i) {
  //     const reco::GenParticle& CurrentParticle=(*genParticles)[i];
  //     if ( (CurrentParticle.status()==1) && ( (CurrentParticle.pdgId()==13)  || (CurrentParticle.pdgId()==-13) ) ){  
  // 	double SmallestDelR = 999.;
  // 	for (auto thisSegment = OurSegments->begin(); thisSegment != OurSegments->end(); 
  // 	     ++thisSegment){
  // 	  ME0DetId id = thisSegment->me0DetId();
  // 	  //std::cout<<"ME0DetId =  "<<id<<std::endl;
  // 	  auto roll = me0Geom->etaPartition(id); 
  // 	  auto GlobVect(roll->toGlobal(thisSegment->localPosition()));
  // 	  Segment_Eta->Fill(fabs(GlobVect.eta()));
  // 	  Segment_Phi->Fill(GlobVect.phi());
  // 	  Segment_R->Fill(GlobVect.perp());
  // 	  Segment_Pos->Fill(GlobVect.x(),GlobVect.y());
	  
  // 	  if (reco::deltaR(CurrentParticle,GlobVect) < SmallestDelR) SmallestDelR = reco::deltaR(CurrentParticle,GlobVect);

  // 	}
  // 	if ((fabs(CurrentParticle.eta()) < 2.0 ) ||(fabs(CurrentParticle.eta()) > 2.8 )) continue;
  // 	DelR_Segment_GenMuon->Fill(SmallestDelR);
  //     }
  //   }
  

  // //std::cout<<recosize<<std::endl;
  // for (std::vector<RecoChargedCandidate>::const_iterator thisCandidate = OurCandidates->begin();
  //      thisCandidate != OurCandidates->end(); ++thisCandidate){
  //   TLorentzVector CandidateVector;
  //   CandidateVector.SetPtEtaPhiM(thisCandidate->pt(),thisCandidate->eta(),thisCandidate->phi(),0);
  //   //std::cout<<"On a muon"<<std::endl;
  //   //std::cout<<thisCandidate->eta()<<std::endl;
  //   Candidate_Eta->Fill(fabs(thisCandidate->eta()));
  // }

  // if (OurCandidates->size() == 2){
  //   TLorentzVector CandidateVector1,CandidateVector2;
  //   CandidateVector1.SetPtEtaPhiM((*OurCandidates)[0].pt(),(*OurCandidates)[0].eta(),(*OurCandidates)[0].phi(),0);
  //   CandidateVector2.SetPtEtaPhiM((*OurCandidates)[1].pt(),(*OurCandidates)[1].eta(),(*OurCandidates)[1].phi(),0);
  //   Double_t Mass = (CandidateVector1+CandidateVector2).M();
  //   Mass_h->Fill(Mass);
  // }
  
  
}

//void ME0MuonAnalyzer::endJob() 
void ME0MuonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) 

{

  //std::cout<<"Nevents = "<<Nevents<<std::endl;
  //TString cmsText     = "CMS Prelim.";
  //TString cmsText     = "#splitline{CMS PhaseII Simulation}{Prelim}";
  TString cmsText     = "CMS PhaseII Simulation Prelim.";

  TString lumiText = "PU 140, 14 TeV";
  float cmsTextFont   = 61;  // default is helvetic-bold

 
  //float extraTextFont = 52;  // default is helvetica-italics
  float lumiTextSize     = 0.05;

  float lumiTextOffset   = 0.2;
  float cmsTextSize      = 0.05;
  //float cmsTextOffset    = 0.1;  // only used in outOfFrame version

  // float relPosX    = 0.045;
  // float relPosY    = 0.035;
  // float relExtraDY = 1.2;

  // //ratio of "CMS" and extra text size
  // float extraOverCmsTextSize  = 0.76;



  histoFile->cd();


  TCanvas *c1 = new TCanvas("c1", "canvas" );


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //setTDRStyle();

  gStyle->SetOptStat(1);     
  //XPull_h->Fit("gaus","","",-1.,1.);
  XPull_h->Draw(); 
  XPull_h->GetXaxis()->SetTitle("Local pulls: X");
  XPull_h->GetXaxis()->SetTitleSize(0.05);
  c1->Print("PullX.png");

  //YPull_h->Fit("gaus");
  YPull_h->Draw(); 
  YPull_h->GetXaxis()->SetTitle("Local pulls: Y");
  YPull_h->GetXaxis()->SetTitleSize(0.05);
  c1->Print("PullY.png");

  gStyle->SetOptStat(1);
  //  XDiff_h->Fit("gaus","","",-1.,1.);
  XDiff_h->Draw(); 
  XDiff_h->GetXaxis()->SetTitle("Local residuals : X");
  XDiff_h->GetXaxis()->SetTitleSize(0.05);
  c1->Print("DiffX.png");

  //  YDiff_h->Fit("gaus");
  YDiff_h->Draw(); 
  YDiff_h->GetXaxis()->SetTitle("Local residuals : Y");
  YDiff_h->GetXaxis()->SetTitleSize(0.05);
  c1->Print("DiffY.png");

  gStyle->SetOptStat(0);

  SegPosDirPhiDiff_True_h->Write();   SegPosDirPhiDiff_True_h->Draw();  c1->Print(histoFolder+"/SegPosDirPhiDiff_True_h.png");
  SegPosDirEtaDiff_True_h->Write();   SegPosDirEtaDiff_True_h->Draw();  c1->Print(histoFolder+"/SegPosDirEtaDiff_True_h.png");

  c1->SetLogy();
  SegPosDirPhiDiff_All_h->Write();   SegPosDirPhiDiff_All_h->Draw();  c1->Print(histoFolder+"/SegPosDirPhiDiff_All_h.png");
  c1->SetLogy();
  SegPosDirEtaDiff_All_h->Write();   SegPosDirEtaDiff_All_h->Draw();  c1->Print(histoFolder+"/SegPosDirEtaDiff_All_h.png");

  SegTrackDirPhiDiff_True_h->Write();   SegTrackDirPhiDiff_True_h->Draw();  c1->Print(histoFolder+"/SegTrackDirPhiDiff_True_h.png");
  SegTrackDirEtaDiff_True_h->Write();   SegTrackDirEtaDiff_True_h->Draw();  c1->Print(histoFolder+"/SegTrackDirEtaDiff_True_h.png");

  SegTrackDirPhiPull_True_h->Write();   SegTrackDirPhiPull_True_h->Draw();  c1->Print(histoFolder+"/SegTrackDirPhiPull_True_h.png");
  SegTrackDirPhiPull_All_h->Write();   SegTrackDirPhiPull_All_h->Draw();  c1->Print(histoFolder+"/SegTrackDirPhiPull_All_h.png");

  c1->SetLogy();
  SegTrackDirPhiDiff_All_h->Write();   SegTrackDirPhiDiff_All_h->Draw();  c1->Print(histoFolder+"/SegTrackDirPhiDiff_All_h.png");
  c1->SetLogy();
  SegTrackDirEtaDiff_All_h->Write();   SegTrackDirEtaDiff_All_h->Draw();  c1->Print(histoFolder+"/SegTrackDirEtaDiff_All_h.png");


  SegGenDirPhiDiff_True_h->Write();   SegGenDirPhiDiff_True_h->Draw();  c1->Print(histoFolder+"/SegGenDirPhiDiff_True_h.png");
  SegGenDirEtaDiff_True_h->Write();   SegGenDirEtaDiff_True_h->Draw();  c1->Print(histoFolder+"/SegGenDirEtaDiff_True_h.png");

  SegGenDirPhiPull_True_h->Write();   SegGenDirPhiPull_True_h->Draw();  c1->Print(histoFolder+"/SegGenDirPhiPull_True_h.png");
  SegGenDirPhiPull_All_h->Write();   SegGenDirPhiPull_All_h->Draw();  c1->Print(histoFolder+"/SegGenDirPhiPull_All_h.png");

  c1->SetLogy();
  SegGenDirPhiDiff_All_h->Write();   SegGenDirPhiDiff_All_h->Draw();  c1->Print(histoFolder+"/SegGenDirPhiDiff_All_h.png");
  c1->SetLogy();
  SegGenDirEtaDiff_All_h->Write();   SegGenDirEtaDiff_All_h->Draw();  c1->Print(histoFolder+"/SegGenDirEtaDiff_All_h.png");

  Candidate_Eta->Write();   Candidate_Eta->Draw();  c1->Print(histoFolder+"/Candidate_Eta.png");
  Track_Eta->Write();   Track_Eta->Draw();  c1->Print(histoFolder+"/Track_Eta.png");
  Track_Pt->Write();   Track_Pt->Draw();  c1->Print(histoFolder+"/Track_Pt.png");

  Segment_Eta->GetXaxis()->SetTitle("me0segment |#eta|");
  Segment_Eta->GetYaxis()->SetTitle(" \# of Segments");
  Segment_Eta->Write();   Segment_Eta->Draw();  
  //GenMuon_Eta->SetLineColor(2);GenMuon_Eta->Draw("SAME"); 
  c1->Print(histoFolder+"/Segment_Eta.png");

  Segment_Phi->Write();   Segment_Phi->Draw();  c1->Print(histoFolder+"/Segment_Phi.png");
  Segment_R->Write();   Segment_R->Draw();  c1->Print(histoFolder+"/Segment_R.png");
  Segment_Pos->Write();   Segment_Pos->Draw();  c1->Print(histoFolder+"/Segment_Pos.png");

  Rechit_Eta->Write();   Rechit_Eta->Draw();   c1->Print(histoFolder+"/Rechit_Eta.png");
  Rechit_Phi->Write();   Rechit_Phi->Draw();  c1->Print(histoFolder+"/Rechit_Phi.png");
  Rechit_R->Write();   Rechit_R->Draw();  c1->Print(histoFolder+"/Rechit_R.png");
  Rechit_Pos->Write();   Rechit_Pos->Draw();  c1->Print(histoFolder+"/Rechit_Pos.png");

  ME0Muon_Eta->Write();   ME0Muon_Eta->Draw();  
  ME0Muon_Eta->GetXaxis()->SetTitle("ME0Muon |#eta|");
  ME0Muon_Eta->GetXaxis()->SetTitleSize(0.05);
  c1->Print(histoFolder+"/ME0Muon_Eta.png");

  CheckME0Muon_Eta->Write();   CheckME0Muon_Eta->Draw();  
  CheckME0Muon_Eta->GetXaxis()->SetTitle("CheckME0Muon |#eta|");
  CheckME0Muon_Eta->GetXaxis()->SetTitleSize(0.05);
  c1->Print(histoFolder+"/CheckME0Muon_Eta.png");

  ME0Muon_Cuts_Eta->Write();   ME0Muon_Cuts_Eta->Draw();  c1->Print(histoFolder+"/ME0Muon_Cuts_Eta.png");
  ME0Muon_Cuts_WidestBinning_Eta->Write();   ME0Muon_Cuts_WidestBinning_Eta->Draw();  c1->Print(histoFolder+"/ME0Muon_Cuts_WidestBinning_Eta.png");
  ME0Muon_Cuts_WideBinning_Eta->Write();   ME0Muon_Cuts_WideBinning_Eta->Draw();  c1->Print(histoFolder+"/ME0Muon_Cuts_WideBinning_Eta.png");
  //c1->SetLogy();
  ME0Muon_Pt->Write();   ME0Muon_Pt->Draw();  
  ME0Muon_Pt->GetXaxis()->SetTitle("ME0Muon p_{T}");
  ME0Muon_Pt->GetXaxis()->SetTitleSize(0.05);
  c1->Print(histoFolder+"/ME0Muon_Pt.png");

  GenMuon_Eta->Write();   GenMuon_Eta->Draw();  c1->Print(histoFolder+"/GenMuon_Eta.png");
  GenMuon_WideBinning_Eta->Write();   GenMuon_WideBinning_Eta->Draw();  c1->Print(histoFolder+"/GenMuon_WideBinning_Eta.png");
  GenMuon_WidestBinning_Eta->Write();   GenMuon_WidestBinning_Eta->Draw();  c1->Print(histoFolder+"/GenMuon_WidestBinning_Eta.png");

  TPMuon_Eta->Write();   TPMuon_Eta->Draw();  c1->Print(histoFolder+"/TPMuon_Eta.png");
  TPMuon_WideBinning_Eta->Write();   TPMuon_WideBinning_Eta->Draw();  c1->Print(histoFolder+"/TPMuon_WideBinning_Eta.png");
  TPMuon_WidestBinning_Eta->Write();   TPMuon_WidestBinning_Eta->Draw();  c1->Print(histoFolder+"/TPMuon_WidestBinning_Eta.png");

  GenMuon_Pt->Write();   GenMuon_Pt->Draw();  c1->Print(histoFolder+"/GenMuon_Pt.png");

  MatchedME0Muon_Eta->Write();   MatchedME0Muon_Eta->Draw();  c1->Print(histoFolder+"/MatchedME0Muon_Eta.png");
  MatchedME0Muon_WideBinning_Eta->Write();   MatchedME0Muon_WideBinning_Eta->Draw();  c1->Print(histoFolder+"/MatchedME0Muon_WideBinning_Eta.png");
  MatchedME0Muon_WidestBinning_Eta->Write();   MatchedME0Muon_WidestBinning_Eta->Draw();  c1->Print(histoFolder+"/MatchedME0Muon_WidestBinning_Eta.png");

  StandaloneMatchedME0Muon_Eta->Write();   StandaloneMatchedME0Muon_Eta->Draw();  c1->Print(histoFolder+"/StandaloneMatchedME0Muon_Eta.png");
  StandaloneMatchedME0Muon_WideBinning_Eta->Write();   StandaloneMatchedME0Muon_WideBinning_Eta->Draw();  c1->Print(histoFolder+"/StandaloneMatchedME0Muon_WideBinning_Eta.png");
  StandaloneMatchedME0Muon_WidestBinning_Eta->Write();   StandaloneMatchedME0Muon_WidestBinning_Eta->Draw();  c1->Print(histoFolder+"/StandaloneMatchedME0Muon_WidestBinning_Eta.png");

  Chi2MatchedME0Muon_Eta->Write();   Chi2MatchedME0Muon_Eta->Draw();  c1->Print(histoFolder+"/Chi2MatchedME0Muon_Eta.png");
  Chi2MatchedME0Muon_WideBinning_Eta->Write();   Chi2MatchedME0Muon_WideBinning_Eta->Draw();  c1->Print(histoFolder+"/Chi2MatchedME0Muon_WideBinning_Eta.png");
  Chi2MatchedME0Muon_WidestBinning_Eta->Write();   Chi2MatchedME0Muon_WidestBinning_Eta->Draw();  c1->Print(histoFolder+"/Chi2MatchedME0Muon_WidestBinning_Eta.png");
  Chi2UnmatchedME0Muon_Eta->Write();   Chi2UnmatchedME0Muon_Eta->Draw();  c1->Print(histoFolder+"/Chi2UnmatchedME0Muon_Eta.png");
  Chi2UnmatchedME0Muon_WideBinning_Eta->Write();   Chi2UnmatchedME0Muon_WideBinning_Eta->Draw();  c1->Print(histoFolder+"/Chi2UnmatchedME0Muon_WideBinning_Eta.png");
  Chi2UnmatchedME0Muon_WidestBinning_Eta->Write();   Chi2UnmatchedME0Muon_WidestBinning_Eta->Draw();  c1->Print(histoFolder+"/Chi2UnmatchedME0Muon_WidestBinning_Eta.png");

  gStyle->SetOptStat(1);
  MatchedME0Muon_Pt->GetXaxis()->SetTitle("ME0Muon p_{T}");
  //MatchedME0Muon_Pt->GetYaxis()->SetTitle(" \# of Se");

  MatchedME0Muon_Pt->Write();   MatchedME0Muon_Pt->Draw();  c1->Print(histoFolder+"/MatchedME0Muon_Pt.png");
  gStyle->SetOptStat(0);

  UnmatchedME0Muon_Eta->Write();   UnmatchedME0Muon_Eta->Draw();  c1->Print(histoFolder+"/UnmatchedME0Muon_Eta.png");
  UnmatchedME0Muon_Cuts_Eta->Write();   UnmatchedME0Muon_Cuts_Eta->Draw();  c1->Print(histoFolder+"/UnmatchedME0Muon_Cuts_Eta.png");
  UnmatchedME0Muon_Cuts_WideBinning_Eta->Write();   UnmatchedME0Muon_Cuts_WideBinning_Eta->Draw();  c1->Print(histoFolder+"/UnmatchedME0Muon_Cuts_WideBinning_Eta.png");
  UnmatchedME0Muon_Cuts_WidestBinning_Eta->Write();   UnmatchedME0Muon_Cuts_WidestBinning_Eta->Draw();  c1->Print(histoFolder+"/UnmatchedME0Muon_Cuts_WidestBinning_Eta.png");


  //gStyle->SetOptStat('oue');
  c1->SetLogy();
  UnmatchedME0Muon_Pt->Write();   UnmatchedME0Muon_Pt->Draw();  c1->Print(histoFolder+"/UnmatchedME0Muon_Pt.png");
  UnmatchedME0Muon_Window_Pt->Write();   UnmatchedME0Muon_Window_Pt->Draw();  c1->Print(histoFolder+"/UnmatchedME0Muon_Window_Pt.png");
  gStyle->SetOptStat(0);

  FailedTrack_Window_XPull->Write();   FailedTrack_Window_XPull->Draw();  c1->Print(histoFolder+"/FailedTrack_Window_XPull.png");
  FailedTrack_Window_YPull->Write();   FailedTrack_Window_YPull->Draw();  c1->Print(histoFolder+"/FailedTrack_Window_YPull.png");
  FailedTrack_Window_XDiff->Write();   FailedTrack_Window_XDiff->Draw();  c1->Print(histoFolder+"/FailedTrack_Window_XDiff.png");
  FailedTrack_Window_YDiff->Write();   FailedTrack_Window_YDiff->Draw();  c1->Print(histoFolder+"/FailedTrack_Window_YDiff.png");
  FailedTrack_Window_PhiDiff->Write();   FailedTrack_Window_PhiDiff->Draw();  c1->Print(histoFolder+"/FailedTrack_Window_PhiDiff.png");

  c1->SetLogy(0);
  TH1F *UnmatchedME0Muon_Cuts_Eta_PerEvent;
  UnmatchedME0Muon_Cuts_Eta_PerEvent = new TH1F("UnmatchedME0Muon_Cuts_Eta_PerEvent"      , "Muon |#eta|"   , 4, 2.0, 2.8 );
  //UnmatchedME0Muon_Cuts_Eta_PerEvent->Sumw2();
  for (int i=1; i<=UnmatchedME0Muon_Cuts_Eta_PerEvent->GetNbinsX(); ++i){
    UnmatchedME0Muon_Cuts_Eta_PerEvent->SetBinContent(i,(UnmatchedME0Muon_Cuts_Eta->GetBinContent(i)));
  }
  UnmatchedME0Muon_Cuts_Eta_PerEvent->Scale(1/Nevents);

  UnmatchedME0Muon_Cuts_Eta_PerEvent->GetXaxis()->SetTitle("ME0Muon |#eta|");
  UnmatchedME0Muon_Cuts_Eta_PerEvent->GetXaxis()->SetTitleSize(0.05);
  
  UnmatchedME0Muon_Cuts_Eta_PerEvent->GetYaxis()->SetTitle("Average \# ME0Muons per event");
  UnmatchedME0Muon_Cuts_Eta_PerEvent->GetYaxis()->SetTitleSize(0.05);

  UnmatchedME0Muon_Cuts_Eta_PerEvent->Write();   UnmatchedME0Muon_Cuts_Eta_PerEvent->Draw();  c1->Print(histoFolder+"/UnmatchedME0Muon_Cuts_Eta_PerEvent.png");

  TH1F *ME0Muon_Cuts_Eta_PerEvent;
  ME0Muon_Cuts_Eta_PerEvent = new TH1F("ME0Muon_Cuts_Eta_PerEvent"      , "Muon |#eta|"   , 4, 2.0, 2.8 );

  for (int i=1; i<=ME0Muon_Cuts_Eta_PerEvent->GetNbinsX(); ++i){
    ME0Muon_Cuts_Eta_PerEvent->SetBinContent(i,(ME0Muon_Cuts_Eta->GetBinContent(i)));
  }
  ME0Muon_Cuts_Eta_PerEvent->Scale(1/Nevents);

  ME0Muon_Cuts_Eta_PerEvent->Write();   ME0Muon_Cuts_Eta_PerEvent->Draw();  c1->Print(histoFolder+"/ME0Muon_Cuts_Eta_PerEvent.png");

  

  Mass_h->Write();   Mass_h->Draw();  c1->Print(histoFolder+"/Mass_h.png");
  TracksPerSegment_s->SetMarkerStyle(1);
  TracksPerSegment_s->SetMarkerSize(3.0);
  TracksPerSegment_s->Write();     TracksPerSegment_s->Draw();  c1->Print(histoFolder+"/TracksPerSegment_s.png");

  TracksPerSegment_h->Write();     TracksPerSegment_h->Draw();  c1->Print(histoFolder+"/TracksPerSegment_h.png");

  TracksPerSegment_p->GetXaxis()->SetTitle("Gen Muon #eta");
  TracksPerSegment_p->GetYaxis()->SetTitle("Average N_{Tracks} per segment");
  TracksPerSegment_p->Write();     TracksPerSegment_p->Draw();  c1->Print(histoFolder+"/TracksPerSegment_p.png");

  ClosestDelR_s->SetMarkerStyle(1);
  ClosestDelR_s->SetMarkerSize(3.0);
  ClosestDelR_s->Write();     ClosestDelR_s->Draw();  c1->Print(histoFolder+"/ClosestDelR_s.png");

  DelR_Window_Under5->Write();     DelR_Window_Under5->Draw();     c1->Print(histoFolder+"/DelR_Window_Under5.png");
  Pt_Window_Under5->Write();    Pt_Window_Under5->Draw();    c1->Print(histoFolder+"/Pt_Window_Under5.png");

  DelR_Track_Window_Under5->Write();     DelR_Track_Window_Under5->Draw();     c1->Print(histoFolder+"/DelR_Track_Window_Under5.png");
  Pt_Track_Window_Under5->Write();    Pt_Track_Window_Under5->Draw();    c1->Print(histoFolder+"/Pt_Track_Window_Under5.png");
  c1->SetLogy(1);
  Pt_Track_Window->Write();    Pt_Track_Window->Draw();    c1->Print(histoFolder+"/Pt_Track_Window.png");
  c1->SetLogy(0);

  DelR_Track_Window_Failed_Under5->Write();     DelR_Track_Window_Failed_Under5->Draw();     c1->Print(histoFolder+"/DelR_Track_Window_Failed_Under5.png");
  Pt_Track_Window_Failed_Under5->Write();    Pt_Track_Window_Failed_Under5->Draw();    c1->Print(histoFolder+"/Pt_Track_Window_Failed_Under5.png");
  c1->SetLogy(1);
  Pt_Track_Window_Failed->Write();    Pt_Track_Window_Failed->Draw();    c1->Print(histoFolder+"/Pt_Track_Window_Failed.png");
  c1->SetLogy(0);

  DelR_Segment_GenMuon->Write();   DelR_Segment_GenMuon->Draw();  c1->Print(histoFolder+"/DelR_Segment_GenMuon.png");

  ClosestDelR_p->GetXaxis()->SetTitle("Gen Muon #eta");
  ClosestDelR_p->GetYaxis()->SetTitle("Average closest #Delta R track");
  std::cout<<"  ClosestDelR_p values:"<<std::endl;
  for (int i=1; i<=ClosestDelR_p->GetNbinsX(); ++i){
    std::cout<<2.4+(double)i*((4.0-2.4)/40.)<<","<<ClosestDelR_p->GetBinContent(i)<<std::endl;
  }
  ClosestDelR_p->Write();     ClosestDelR_p->Draw();  c1->Print(histoFolder+"/ClosestDelR_p.png");

  FakeTracksPerSegment_s->SetMarkerStyle(1);
  FakeTracksPerSegment_s->SetMarkerSize(3.0);
  FakeTracksPerSegment_s->Write();     FakeTracksPerSegment_s->Draw();  c1->Print(histoFolder+"/FakeTracksPerSegment_s.png");

  FakeTracksPerSegment_h->Write();     FakeTracksPerSegment_h->Draw();  c1->Print(histoFolder+"/FakeTracksPerSegment_h.png");

  FakeTracksPerSegment_p->GetXaxis()->SetTitle("Gen Muon #eta");
  FakeTracksPerSegment_p->GetYaxis()->SetTitle("Average N_{Tracks} per segment");
  FakeTracksPerSegment_p->Write();     FakeTracksPerSegment_p->Draw();  c1->Print(histoFolder+"/FakeTracksPerSegment_p.png");

  FakeTracksPerAssociatedSegment_s->SetMarkerStyle(1);
  FakeTracksPerAssociatedSegment_s->SetMarkerSize(3.0);
  FakeTracksPerAssociatedSegment_s->Write();     FakeTracksPerAssociatedSegment_s->Draw();  c1->Print(histoFolder+"/FakeTracksPerAssociatedSegment_s.png");

  FakeTracksPerAssociatedSegment_h->Write();     FakeTracksPerAssociatedSegment_h->Draw();  c1->Print(histoFolder+"/FakeTracksPerAssociatedSegment_h.png");

  FakeTracksPerAssociatedSegment_p->GetXaxis()->SetTitle("Gen Muon #eta");
  FakeTracksPerAssociatedSegment_p->GetYaxis()->SetTitle("Average N_{Tracks} per segment");
  FakeTracksPerAssociatedSegment_p->Write();     FakeTracksPerAssociatedSegment_p->Draw();  c1->Print(histoFolder+"/FakeTracksPerAssociatedSegment_p.png");

  PreMatch_TP_R->Write(); PreMatch_TP_R->Draw();  c1->Print(histoFolder+"/PreMatch_TP_R.png");
  PostMatch_TP_R->Write(); PostMatch_TP_R->Draw();  c1->Print(histoFolder+"/PostMatch_TP_R.png");
  PostMatch_BX0_TP_R->Write(); PostMatch_BX0_TP_R->Draw();  c1->Print(histoFolder+"/PostMatch_BX0_TP_R.png");

  GenMuon_Eta->Sumw2();  MatchedME0Muon_Eta->Sumw2();  Chi2MatchedME0Muon_Eta->Sumw2();   Chi2UnmatchedME0Muon_Eta->Sumw2();TPMuon_Eta->Sumw2();
  GenMuon_Pt->Sumw2();  MatchedME0Muon_Pt->Sumw2();
  StandaloneMatchedME0Muon_Eta->Sumw2();
  StandaloneMatchedME0Muon_WideBinning_Eta->Sumw2();
  StandaloneMatchedME0Muon_WidestBinning_Eta->Sumw2();
  

  Track_Eta->Sumw2();  ME0Muon_Eta->Sumw2();
  Track_Pt->Sumw2();  ME0Muon_Pt->Sumw2();

  UnmatchedME0Muon_Eta->Sumw2();
  UnmatchedME0Muon_Pt->Sumw2();
  
  UnmatchedME0Muon_Cuts_Eta->Sumw2();    ME0Muon_Cuts_Eta->Sumw2();

  ME0Muon_Cuts_Eta_5_10->Sumw2();  ME0Muon_Cuts_Eta_9_11->Sumw2();  ME0Muon_Cuts_Eta_10_50->Sumw2();  ME0Muon_Cuts_Eta_50_100->Sumw2();  ME0Muon_Cuts_Eta_100->Sumw2();
  UnmatchedME0Muon_Cuts_Eta_5_10->Sumw2();    UnmatchedME0Muon_Cuts_Eta_9_11->Sumw2();  UnmatchedME0Muon_Cuts_Eta_10_50->Sumw2();  UnmatchedME0Muon_Cuts_Eta_50_100->Sumw2();  UnmatchedME0Muon_Cuts_Eta_100->Sumw2();
  GenMuon_Eta_5_10->Sumw2();    GenMuon_Eta_9_11->Sumw2();  GenMuon_Eta_10_50->Sumw2();  GenMuon_Eta_50_100->Sumw2();  GenMuon_Eta_100->Sumw2();
  MatchedME0Muon_Eta_5_10->Sumw2();   MatchedME0Muon_Eta_9_11->Sumw2();  MatchedME0Muon_Eta_10_50->Sumw2();  MatchedME0Muon_Eta_50_100->Sumw2();  MatchedME0Muon_Eta_100->Sumw2();

  Chi2MatchedME0Muon_Eta_5_10->Sumw2();   Chi2MatchedME0Muon_Eta_9_11->Sumw2();  Chi2MatchedME0Muon_Eta_10_50->Sumw2();  Chi2MatchedME0Muon_Eta_50_100->Sumw2();  Chi2MatchedME0Muon_Eta_100->Sumw2();

  UnmatchedME0Muon_Cuts_WideBinning_Eta->Sumw2();
  UnmatchedME0Muon_Cuts_WidestBinning_Eta->Sumw2();
  GenMuon_WideBinning_Eta->Sumw2();
  GenMuon_WidestBinning_Eta->Sumw2();
  TPMuon_WideBinning_Eta->Sumw2();
  TPMuon_WidestBinning_Eta->Sumw2();
  MatchedME0Muon_WideBinning_Eta->Sumw2();
  MatchedME0Muon_WidestBinning_Eta->Sumw2();
  Chi2MatchedME0Muon_WideBinning_Eta->Sumw2();
  Chi2MatchedME0Muon_WidestBinning_Eta->Sumw2();
  ME0Muon_Cuts_WideBinning_Eta->Sumw2();
  ME0Muon_Cuts_WidestBinning_Eta->Sumw2();
  Chi2UnmatchedME0Muon_WideBinning_Eta->Sumw2();
  Chi2UnmatchedME0Muon_WidestBinning_Eta->Sumw2();


  //Captions/labels
  std::stringstream PtCutString;

  PtCutString<<"#splitline{DY }{Reco Track p_{T} > "<<FakeRatePtCut<<" GeV}";
  const std::string& ptmp = PtCutString.str();
  const char* pcstr = ptmp.c_str();


  TLatex* txt =new TLatex;
  //txt->SetTextAlign(12);
  //txt->SetTextFont(42);
  txt->SetNDC();
  //txt->SetTextSize(0.05);
  txt->SetTextFont(132);
  txt->SetTextSize(0.05);


  float t = c1->GetTopMargin();
  float r = c1->GetRightMargin();

  TLatex* latex1 = new TLatex;
  latex1->SetNDC();
  latex1->SetTextAngle(0);
  latex1->SetTextColor(kBlack);    


  latex1->SetTextFont(42);
  latex1->SetTextAlign(31); 
  //latex1->SetTextSize(lumiTextSize*t);    
  latex1->SetTextSize(lumiTextSize);    
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  TLatex* latex = new TLatex;
  latex->SetTextFont(cmsTextFont);
  latex->SetNDC();
  latex->SetTextSize(cmsTextSize);
  //latex->SetTextAlign(align_);

  //End captions/labels
  std::cout<<"GenMuon_Eta =  "<<GenMuon_Eta->Integral()<<std::endl;
  std::cout<<"MatchedME0Muon_Eta =  "<<MatchedME0Muon_Eta->Integral()<<std::endl;

  MuonRecoEff_Eta->Divide(MatchedME0Muon_Eta, GenMuon_Eta, 1, 1, "B");
  MuonRecoEff_Eta->GetXaxis()->SetTitle("Gen Muon |#eta|");
  MuonRecoEff_Eta->GetXaxis()->SetTitleSize(0.05);
  MuonRecoEff_Eta->GetYaxis()->SetTitle("ME0Muon Efficiency");
  MuonRecoEff_Eta->GetYaxis()->SetTitleSize(0.05);
  //MuonRecoEff_Eta->SetMinimum(MuonRecoEff_Eta->GetMinimum()-0.1);
  MuonRecoEff_Eta->SetMinimum(0);
  //MuonRecoEff_Eta->SetMaximum(MuonRecoEff_Eta->GetMaximum()+0.1);
  MuonRecoEff_Eta->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  MuonRecoEff_Eta->Write();   MuonRecoEff_Eta->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestMuonRecoEff_Eta.png");
  c1->Print(histoFolder+"/MuonRecoEff_Eta.png");


  MuonRecoEff_WideBinning_Eta->Divide(MatchedME0Muon_WideBinning_Eta, GenMuon_WideBinning_Eta, 1, 1, "B");
  MuonRecoEff_WideBinning_Eta->GetXaxis()->SetTitle("Gen Muon |#eta|");
  MuonRecoEff_WideBinning_Eta->GetXaxis()->SetTitleSize(0.05);
  MuonRecoEff_WideBinning_Eta->GetYaxis()->SetTitle("ME0Muon Efficiency");
  MuonRecoEff_WideBinning_Eta->GetYaxis()->SetTitleSize(0.05);
  //MuonRecoEff_WideBinning_Eta->SetMinimum(MuonRecoEff_WideBinning_Eta->GetMinimum()-0.1);
  MuonRecoEff_WideBinning_Eta->SetMinimum(0);
  //MuonRecoEff_WideBinning_Eta->SetMaximum(MuonRecoEff_WideBinning_Eta->GetMaximum()+0.1);
  MuonRecoEff_WideBinning_Eta->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  MuonRecoEff_WideBinning_Eta->Write();   MuonRecoEff_WideBinning_Eta->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestMuonRecoEff_WideBinning_Eta.png");
  c1->Print(histoFolder+"/MuonRecoEff_WideBinning_Eta.png");


  MuonRecoEff_WidestBinning_Eta->Divide(MatchedME0Muon_WidestBinning_Eta, GenMuon_WidestBinning_Eta, 1, 1, "B");
  MuonRecoEff_WidestBinning_Eta->GetXaxis()->SetTitle("Gen Muon |#eta|");
  MuonRecoEff_WidestBinning_Eta->GetXaxis()->SetTitleSize(0.05);
  MuonRecoEff_WidestBinning_Eta->GetYaxis()->SetTitle("ME0Muon Efficiency");
  MuonRecoEff_WidestBinning_Eta->GetYaxis()->SetTitleSize(0.05);
  //MuonRecoEff_WidestBinning_Eta->SetMinimum(MuonRecoEff_WidestBinning_Eta->GetMinimum()-0.1);
  MuonRecoEff_WidestBinning_Eta->SetMinimum(0);
  //MuonRecoEff_WidestBinning_Eta->SetMaximum(MuonRecoEff_WidestBinning_Eta->GetMaximum()+0.1);
  MuonRecoEff_WidestBinning_Eta->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  MuonRecoEff_WidestBinning_Eta->Write();   MuonRecoEff_WidestBinning_Eta->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestMuonRecoEff_WidestBinning_Eta.png");
  c1->Print(histoFolder+"/MuonRecoEff_WidestBinning_Eta.png");


  StandaloneMuonRecoEff_Eta->Divide(StandaloneMatchedME0Muon_Eta, GenMuon_Eta, 1, 1, "B");
  StandaloneMuonRecoEff_Eta->GetXaxis()->SetTitle("Gen Muon |#eta|");
  StandaloneMuonRecoEff_Eta->GetXaxis()->SetTitleSize(0.05);
  StandaloneMuonRecoEff_Eta->GetYaxis()->SetTitle("Standalone Muon Efficiency");
  StandaloneMuonRecoEff_Eta->GetYaxis()->SetTitleSize(0.05);
  //StandaloneMuonRecoEff_Eta->SetMinimum(StandaloneMuonRecoEff_Eta->GetMinimum()-0.1);
  StandaloneMuonRecoEff_Eta->SetMinimum(0);
  //StandaloneMuonRecoEff_Eta->SetMaximum(StandaloneMuonRecoEff_Eta->GetMaximum()+0.1);
  StandaloneMuonRecoEff_Eta->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  StandaloneMuonRecoEff_Eta->Write();   StandaloneMuonRecoEff_Eta->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestStandaloneMuonRecoEff_Eta.png");
  c1->Print(histoFolder+"/StandaloneMuonRecoEff_Eta.png");


  StandaloneMuonRecoEff_WideBinning_Eta->Divide(StandaloneMatchedME0Muon_WideBinning_Eta, GenMuon_WideBinning_Eta, 1, 1, "B");
  StandaloneMuonRecoEff_WideBinning_Eta->GetXaxis()->SetTitle("Gen Muon |#eta|");
  StandaloneMuonRecoEff_WideBinning_Eta->GetXaxis()->SetTitleSize(0.05);
  StandaloneMuonRecoEff_WideBinning_Eta->GetYaxis()->SetTitle("Standalone Muon Efficiency");
  StandaloneMuonRecoEff_WideBinning_Eta->GetYaxis()->SetTitleSize(0.05);
  //StandaloneMuonRecoEff_WideBinning_Eta->SetMinimum(StandaloneMuonRecoEff_WideBinning_Eta->GetMinimum()-0.1);
  StandaloneMuonRecoEff_WideBinning_Eta->SetMinimum(0);
  //StandaloneMuonRecoEff_WideBinning_Eta->SetMaximum(StandaloneMuonRecoEff_WideBinning_Eta->GetMaximum()+0.1);
  StandaloneMuonRecoEff_WideBinning_Eta->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  StandaloneMuonRecoEff_WideBinning_Eta->Write();   StandaloneMuonRecoEff_WideBinning_Eta->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestStandaloneMuonRecoEff_WideBinning_Eta.png");
  c1->Print(histoFolder+"/StandaloneMuonRecoEff_WideBinning_Eta.png");


  StandaloneMuonRecoEff_WidestBinning_Eta->Divide(StandaloneMatchedME0Muon_WidestBinning_Eta, GenMuon_WidestBinning_Eta, 1, 1, "B");
  StandaloneMuonRecoEff_WidestBinning_Eta->GetXaxis()->SetTitle("Gen Muon |#eta|");
  StandaloneMuonRecoEff_WidestBinning_Eta->GetXaxis()->SetTitleSize(0.05);
  StandaloneMuonRecoEff_WidestBinning_Eta->GetYaxis()->SetTitle("Standalone Muon Efficiency");
  StandaloneMuonRecoEff_WidestBinning_Eta->GetYaxis()->SetTitleSize(0.05);
  //StandaloneMuonRecoEff_WidestBinning_Eta->SetMinimum(StandaloneMuonRecoEff_WidestBinning_Eta->GetMinimum()-0.1);
  StandaloneMuonRecoEff_WidestBinning_Eta->SetMinimum(0);
  //StandaloneMuonRecoEff_WidestBinning_Eta->SetMaximum(StandaloneMuonRecoEff_WidestBinning_Eta->GetMaximum()+0.1);
  StandaloneMuonRecoEff_WidestBinning_Eta->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  StandaloneMuonRecoEff_WidestBinning_Eta->Write();   StandaloneMuonRecoEff_WidestBinning_Eta->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestStandaloneMuonRecoEff_WidestBinning_Eta.png");
  c1->Print(histoFolder+"/StandaloneMuonRecoEff_WidestBinning_Eta.png");


  MuonRecoEff_Eta_5_10->Divide(MatchedME0Muon_Eta_5_10, GenMuon_Eta_5_10, 1, 1, "B");
  MuonRecoEff_Eta_5_10->GetXaxis()->SetTitle("Gen Muon |#eta|");
  MuonRecoEff_Eta_5_10->GetXaxis()->SetTitleSize(0.05);
  MuonRecoEff_Eta_5_10->GetYaxis()->SetTitle("ME0Muon Efficiency");
  MuonRecoEff_Eta_5_10->GetYaxis()->SetTitleSize(0.05);
  //MuonRecoEff_Eta_5_10->SetMinimum(MuonRecoEff_Eta_5_10->GetMinimum()-0.1);
  MuonRecoEff_Eta_5_10->SetMinimum(0);
  //MuonRecoEff_Eta_5_10->SetMaximum(MuonRecoEff_Eta_5_10->GetMaximum()+0.1);
  MuonRecoEff_Eta_5_10->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  MuonRecoEff_Eta_5_10->Write();   MuonRecoEff_Eta_5_10->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestMuonRecoEff_Eta_5_10.png");
  c1->Print(histoFolder+"/MuonRecoEff_Eta_5_10.png");

  MuonRecoEff_Eta_9_11->Divide(MatchedME0Muon_Eta_9_11, GenMuon_Eta_9_11, 1, 1, "B");
  MuonRecoEff_Eta_9_11->GetXaxis()->SetTitle("Gen Muon |#eta|");
  MuonRecoEff_Eta_9_11->GetXaxis()->SetTitleSize(0.05);
  MuonRecoEff_Eta_9_11->GetYaxis()->SetTitle("ME0Muon Efficiency");
  MuonRecoEff_Eta_9_11->GetYaxis()->SetTitleSize(0.05);
  //MuonRecoEff_Eta_9_11->SetMinimum(MuonRecoEff_Eta_9_11->GetMinimum()-0.1);
  MuonRecoEff_Eta_9_11->SetMinimum(0);
  //MuonRecoEff_Eta_9_11->SetMaximum(MuonRecoEff_Eta_9_11->GetMaximum()+0.1);
  MuonRecoEff_Eta_9_11->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  MuonRecoEff_Eta_9_11->Write();   MuonRecoEff_Eta_9_11->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestMuonRecoEff_Eta_9_11.png");
  c1->Print(histoFolder+"/MuonRecoEff_Eta_9_11.png");

  MuonRecoEff_Eta_10_50->Divide(MatchedME0Muon_Eta_10_50, GenMuon_Eta_10_50, 1, 1, "B");
  MuonRecoEff_Eta_10_50->GetXaxis()->SetTitle("Gen Muon |#eta|");
  MuonRecoEff_Eta_10_50->GetXaxis()->SetTitleSize(0.05);
  MuonRecoEff_Eta_10_50->GetYaxis()->SetTitle("ME0Muon Efficiency");
  MuonRecoEff_Eta_10_50->GetYaxis()->SetTitleSize(0.05);
  //MuonRecoEff_Eta_10_50->SetMinimum(MuonRecoEff_Eta_10_50->GetMinimum()-0.1);
  MuonRecoEff_Eta_10_50->SetMinimum(0);
  //MuonRecoEff_Eta_10_50->SetMaximum(MuonRecoEff_Eta_10_50->GetMaximum()+0.1);
  MuonRecoEff_Eta_10_50->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  MuonRecoEff_Eta_10_50->Write();   MuonRecoEff_Eta_10_50->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestMuonRecoEff_Eta_10_50.png");
  c1->Print(histoFolder+"/MuonRecoEff_Eta_10_50.png");


  MuonRecoEff_Eta_50_100->Divide(MatchedME0Muon_Eta_50_100, GenMuon_Eta_50_100, 1, 1, "B");
  MuonRecoEff_Eta_50_100->GetXaxis()->SetTitle("Gen Muon |#eta|");
  MuonRecoEff_Eta_50_100->GetXaxis()->SetTitleSize(0.05);
  MuonRecoEff_Eta_50_100->GetYaxis()->SetTitle("ME0Muon Efficiency");
  MuonRecoEff_Eta_50_100->GetYaxis()->SetTitleSize(0.05);
  //MuonRecoEff_Eta_50_100->SetMinimum(MuonRecoEff_Eta_50_100->GetMinimum()-0.1);
  MuonRecoEff_Eta_50_100->SetMinimum(0);
  //MuonRecoEff_Eta_50_100->SetMaximum(MuonRecoEff_Eta_50_100->GetMaximum()+0.1);
  MuonRecoEff_Eta_50_100->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  MuonRecoEff_Eta_50_100->Write();   MuonRecoEff_Eta_50_100->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestMuonRecoEff_Eta_50_100.png");
  c1->Print(histoFolder+"/MuonRecoEff_Eta_50_100.png");


  MuonRecoEff_Eta_100->Divide(MatchedME0Muon_Eta_100, GenMuon_Eta_100, 1, 1, "B");
  MuonRecoEff_Eta_100->GetXaxis()->SetTitle("Gen Muon |#eta|");
  MuonRecoEff_Eta_100->GetXaxis()->SetTitleSize(0.05);
  MuonRecoEff_Eta_100->GetYaxis()->SetTitle("ME0Muon Efficiency");
  MuonRecoEff_Eta_100->GetYaxis()->SetTitleSize(0.05);
  //MuonRecoEff_Eta_100->SetMinimum(MuonRecoEff_Eta_100->GetMinimum()-0.1);
  MuonRecoEff_Eta_100->SetMinimum(0);
  //MuonRecoEff_Eta_100->SetMaximum(MuonRecoEff_Eta_100->GetMaximum()+0.1);
  MuonRecoEff_Eta_100->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  MuonRecoEff_Eta_100->Write();   MuonRecoEff_Eta_100->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestMuonRecoEff_Eta_100.png");
  c1->Print(histoFolder+"/MuonRecoEff_Eta_100.png");





  Chi2MuonRecoEff_Eta->Divide(Chi2MatchedME0Muon_Eta, TPMuon_Eta, 1, 1, "B");
  std::cout<<"TPMuon_Eta =  "<<TPMuon_Eta->Integral()<<std::endl;
  std::cout<<"Chi2MatchedME0Muon_Eta =  "<<Chi2MatchedME0Muon_Eta->Integral()<<std::endl;
  Chi2MuonRecoEff_Eta->GetXaxis()->SetTitle("Gen Muon |#eta|");
  Chi2MuonRecoEff_Eta->GetXaxis()->SetTitleSize(0.05);
  
  Chi2MuonRecoEff_Eta->GetYaxis()->SetTitle("ME0Muon Efficiency");
  Chi2MuonRecoEff_Eta->GetYaxis()->SetTitleSize(0.05);
  Chi2MuonRecoEff_Eta->SetMinimum(0);
  Chi2MuonRecoEff_Eta->SetMaximum(1.2);

  Chi2MuonRecoEff_Eta->Write();   Chi2MuonRecoEff_Eta->Draw();  

  txt->DrawLatex(0.15,0.2,pcstr);

   
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //latex->SetTextAlign(align_);
  latex->DrawLatex(0.4, 0.85, cmsText);

  c1->Print(histoFolder+"/Chi2MuonRecoEff_Eta.png");

  Chi2MuonRecoEff_WideBinning_Eta->Divide(Chi2MatchedME0Muon_WideBinning_Eta, TPMuon_WideBinning_Eta, 1, 1, "B");
  std::cout<<"TPMuon_WideBinning_Eta =  "<<TPMuon_WideBinning_Eta->Integral()<<std::endl;
  std::cout<<"Chi2MatchedME0Muon_WideBinning_Eta =  "<<Chi2MatchedME0Muon_WideBinning_Eta->Integral()<<std::endl;
  Chi2MuonRecoEff_WideBinning_Eta->GetXaxis()->SetTitle("Gen Muon |#eta|");
  Chi2MuonRecoEff_WideBinning_Eta->GetXaxis()->SetTitleSize(0.05);
  
  Chi2MuonRecoEff_WideBinning_Eta->GetYaxis()->SetTitle("ME0Muon Efficiency");
  Chi2MuonRecoEff_WideBinning_Eta->GetYaxis()->SetTitleSize(0.05);
  Chi2MuonRecoEff_WideBinning_Eta->SetMinimum(0);
  Chi2MuonRecoEff_WideBinning_Eta->SetMaximum(1.2);

  Chi2MuonRecoEff_WideBinning_Eta->Write();   Chi2MuonRecoEff_WideBinning_Eta->Draw();  

  txt->DrawLatex(0.15,0.2,pcstr);

   
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //latex->SetTextAlign(align_);
  latex->DrawLatex(0.4, 0.85, cmsText);

  c1->Print(histoFolder+"/Chi2MuonRecoEff_WideBinning_Eta.png");

  Chi2MuonRecoEff_WidestBinning_Eta->Divide(Chi2MatchedME0Muon_WidestBinning_Eta, TPMuon_WidestBinning_Eta, 1, 1, "B");
  std::cout<<"TPMuon_WidestBinning_Eta =  "<<TPMuon_WidestBinning_Eta->Integral()<<std::endl;
  std::cout<<"Chi2MatchedME0Muon_WidestBinning_Eta =  "<<Chi2MatchedME0Muon_WidestBinning_Eta->Integral()<<std::endl;
  Chi2MuonRecoEff_WidestBinning_Eta->GetXaxis()->SetTitle("Gen Muon |#eta|");
  Chi2MuonRecoEff_WidestBinning_Eta->GetXaxis()->SetTitleSize(0.05);
  
  Chi2MuonRecoEff_WidestBinning_Eta->GetYaxis()->SetTitle("ME0Muon Efficiency");
  Chi2MuonRecoEff_WidestBinning_Eta->GetYaxis()->SetTitleSize(0.05);
  Chi2MuonRecoEff_WidestBinning_Eta->SetMinimum(0);
  Chi2MuonRecoEff_WidestBinning_Eta->SetMaximum(1.2);

  Chi2MuonRecoEff_WidestBinning_Eta->Write();   Chi2MuonRecoEff_WidestBinning_Eta->Draw();  

  txt->DrawLatex(0.15,0.2,pcstr);

   
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //latex->SetTextAlign(align_);
  latex->DrawLatex(0.4, 0.85, cmsText);

  c1->Print(histoFolder+"/Chi2MuonRecoEff_WidestBinning_Eta.png");

  std::cout<<"Here0"<<std::endl;

  Chi2MuonRecoEff_Eta_5_10->Divide(Chi2MatchedME0Muon_Eta_5_10, GenMuon_Eta_5_10, 1, 1, "B");
  std::cout<<"Here0"<<std::endl;
  Chi2MuonRecoEff_Eta_5_10->GetXaxis()->SetTitle("Gen Muon |#eta|");
  Chi2MuonRecoEff_Eta_5_10->GetXaxis()->SetTitleSize(0.05);
  Chi2MuonRecoEff_Eta_5_10->GetYaxis()->SetTitle("ME0Muon Efficiency");
  Chi2MuonRecoEff_Eta_5_10->GetYaxis()->SetTitleSize(0.05);
  //Chi2MuonRecoEff_Eta_5_10->SetMinimum(Chi2MuonRecoEff_Eta_5_10->GetMinimum()-0.1);
  Chi2MuonRecoEff_Eta_5_10->SetMinimum(0);
  //Chi2MuonRecoEff_Eta_5_10->SetMaximum(Chi2MuonRecoEff_Eta_5_10->GetMaximum()+0.1);
  Chi2MuonRecoEff_Eta_5_10->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  Chi2MuonRecoEff_Eta_5_10->Write();   Chi2MuonRecoEff_Eta_5_10->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestChi2MuonRecoEff_Eta_5_10.png");
  c1->Print(histoFolder+"/Chi2MuonRecoEff_Eta_5_10.png");


  Chi2MuonRecoEff_Eta_9_11->Divide(Chi2MatchedME0Muon_Eta_9_11, GenMuon_Eta_9_11, 1, 1, "B");
  std::cout<<"Here0"<<std::endl;
  Chi2MuonRecoEff_Eta_9_11->GetXaxis()->SetTitle("Gen Muon |#eta|");
  Chi2MuonRecoEff_Eta_9_11->GetXaxis()->SetTitleSize(0.05);
  Chi2MuonRecoEff_Eta_9_11->GetYaxis()->SetTitle("ME0Muon Efficiency");
  Chi2MuonRecoEff_Eta_9_11->GetYaxis()->SetTitleSize(0.05);
  //Chi2MuonRecoEff_Eta_9_11->SetMinimum(Chi2MuonRecoEff_Eta_9_11->GetMinimum()-0.1);
  Chi2MuonRecoEff_Eta_9_11->SetMinimum(0);
  //Chi2MuonRecoEff_Eta_9_11->SetMaximum(Chi2MuonRecoEff_Eta_9_11->GetMaximum()+0.1);
  Chi2MuonRecoEff_Eta_9_11->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  Chi2MuonRecoEff_Eta_9_11->Write();   Chi2MuonRecoEff_Eta_9_11->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestChi2MuonRecoEff_Eta_9_11.png");
  c1->Print(histoFolder+"/Chi2MuonRecoEff_Eta_9_11.png");

  std::cout<<"Here"<<std::endl;

  Chi2MuonRecoEff_Eta_10_50->Divide(Chi2MatchedME0Muon_Eta_10_50, GenMuon_Eta_10_50, 1, 1, "B");
  Chi2MuonRecoEff_Eta_10_50->GetXaxis()->SetTitle("Gen Muon |#eta|");
  Chi2MuonRecoEff_Eta_10_50->GetXaxis()->SetTitleSize(0.05);
  Chi2MuonRecoEff_Eta_10_50->GetYaxis()->SetTitle("ME0Muon Efficiency");
  Chi2MuonRecoEff_Eta_10_50->GetYaxis()->SetTitleSize(0.05);
  //Chi2MuonRecoEff_Eta_10_50->SetMinimum(Chi2MuonRecoEff_Eta_10_50->GetMinimum()-0.1);
  Chi2MuonRecoEff_Eta_10_50->SetMinimum(0);
  //Chi2MuonRecoEff_Eta_10_50->SetMaximum(Chi2MuonRecoEff_Eta_10_50->GetMaximum()+0.1);
  Chi2MuonRecoEff_Eta_10_50->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  Chi2MuonRecoEff_Eta_10_50->Write();   Chi2MuonRecoEff_Eta_10_50->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestChi2MuonRecoEff_Eta_10_50.png");
  c1->Print(histoFolder+"/Chi2MuonRecoEff_Eta_10_50.png");


  Chi2MuonRecoEff_Eta_50_100->Divide(Chi2MatchedME0Muon_Eta_50_100, GenMuon_Eta_50_100, 1, 1, "B");
  Chi2MuonRecoEff_Eta_50_100->GetXaxis()->SetTitle("Gen Muon |#eta|");
  Chi2MuonRecoEff_Eta_50_100->GetXaxis()->SetTitleSize(0.05);
  Chi2MuonRecoEff_Eta_50_100->GetYaxis()->SetTitle("ME0Muon Efficiency");
  Chi2MuonRecoEff_Eta_50_100->GetYaxis()->SetTitleSize(0.05);
  //Chi2MuonRecoEff_Eta_50_100->SetMinimum(Chi2MuonRecoEff_Eta_50_100->GetMinimum()-0.1);
  Chi2MuonRecoEff_Eta_50_100->SetMinimum(0);
  //Chi2MuonRecoEff_Eta_50_100->SetMaximum(Chi2MuonRecoEff_Eta_50_100->GetMaximum()+0.1);
  Chi2MuonRecoEff_Eta_50_100->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  Chi2MuonRecoEff_Eta_50_100->Write();   Chi2MuonRecoEff_Eta_50_100->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestChi2MuonRecoEff_Eta_50_100.png");
  c1->Print(histoFolder+"/Chi2MuonRecoEff_Eta_50_100.png");


  Chi2MuonRecoEff_Eta_100->Divide(Chi2MatchedME0Muon_Eta_100, GenMuon_Eta_100, 1, 1, "B");
  Chi2MuonRecoEff_Eta_100->GetXaxis()->SetTitle("Gen Muon |#eta|");
  Chi2MuonRecoEff_Eta_100->GetXaxis()->SetTitleSize(0.05);
  Chi2MuonRecoEff_Eta_100->GetYaxis()->SetTitle("ME0Muon Efficiency");
  Chi2MuonRecoEff_Eta_100->GetYaxis()->SetTitleSize(0.05);
  //Chi2MuonRecoEff_Eta_100->SetMinimum(Chi2MuonRecoEff_Eta_100->GetMinimum()-0.1);
  Chi2MuonRecoEff_Eta_100->SetMinimum(0);
  //Chi2MuonRecoEff_Eta_100->SetMaximum(Chi2MuonRecoEff_Eta_100->GetMaximum()+0.1);
  Chi2MuonRecoEff_Eta_100->SetMaximum(1.2);
  //CMS_lumi( c1, 7, 11 );
  Chi2MuonRecoEff_Eta_100->Write();   Chi2MuonRecoEff_Eta_100->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);

  //c1->SaveAs("TestChi2MuonRecoEff_Eta_100.png");
  c1->Print(histoFolder+"/Chi2MuonRecoEff_Eta_100.png");

  std::cout<<"  MuonRecoEff_Eta values:"<<std::endl;
  //MuonRecoEff_Eta->Sumw2();
  for (int i=1; i<=MuonRecoEff_Eta->GetNbinsX(); ++i){
    std::cout<<2.4+(double)i*((4.0-2.4)/40.)<<","<<MuonRecoEff_Eta->GetBinContent(i)<<","<<MuonRecoEff_Eta->GetBinError(i)<<std::endl;
    
  }
  

  // MuonRecoEff_Pt->Divide(MatchedME0Muon_Pt, GenMuon_Pt, 1, 1, "B");
  // MuonRecoEff_Pt->GetXaxis()->SetTitle("Gen Muon p_{T}");
  // MuonRecoEff_Pt->GetYaxis()->SetTitle("Matching Efficiency");
  // MuonRecoEff_Pt->SetMinimum(.85);
  // MuonRecoEff_Pt->Write();   MuonRecoEff_Pt->Draw();  c1->Print(histoFolder+"/MuonRecoEff_Pt.png");

  std::cout<<"UnmatchedME0Muon_Eta =  "<<UnmatchedME0Muon_Cuts_Eta->Integral()<<std::endl;
  std::cout<<"ME0Muon_Eta =  "<<ME0Muon_Cuts_Eta->Integral()<<std::endl;


  FakeRate_Eta->Divide(UnmatchedME0Muon_Cuts_Eta, ME0Muon_Cuts_Eta, 1, 1, "B");
  FakeRate_Eta->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  FakeRate_Eta->GetXaxis()->SetTitleSize(0.05);
  FakeRate_Eta->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  FakeRate_Eta->GetYaxis()->SetTitleSize(0.05);
  //FakeRate_Eta->SetMinimum(FakeRate_Eta->GetMinimum()-0.1);
  FakeRate_Eta->SetMinimum(0);
  //FakeRate_Eta->SetMaximum(FakeRate_Eta->GetMaximum()+0.1);
  FakeRate_Eta->SetMaximum(1.2);
  FakeRate_Eta->Write();   FakeRate_Eta->Draw();  

  txt->DrawLatex(0.15,0.4,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestFakeRate_Eta.png');
  c1->Print(histoFolder+"/FakeRate_Eta.png");

  FakeRate_WideBinning_Eta->Divide(UnmatchedME0Muon_Cuts_WideBinning_Eta, ME0Muon_Cuts_WideBinning_Eta, 1, 1, "B");
  FakeRate_WideBinning_Eta->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  FakeRate_WideBinning_Eta->GetXaxis()->SetTitleSize(0.05);
  FakeRate_WideBinning_Eta->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  FakeRate_WideBinning_Eta->GetYaxis()->SetTitleSize(0.05);
  //FakeRate_WideBinning_Eta->SetMinimum(FakeRate_WideBinning_Eta->GetMinimum()-0.1);
  FakeRate_WideBinning_Eta->SetMinimum(0);
  //FakeRate_WideBinning_Eta->SetMaximum(FakeRate_WideBinning_Eta->GetMaximum()+0.1);
  FakeRate_WideBinning_Eta->SetMaximum(1.2);
  FakeRate_WideBinning_Eta->Write();   FakeRate_WideBinning_Eta->Draw();  

  txt->DrawLatex(0.15,0.4,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestFakeRate_WideBinning_Eta.png');
  c1->Print(histoFolder+"/FakeRate_WideBinning_Eta.png");


  FakeRate_WidestBinning_Eta->Divide(UnmatchedME0Muon_Cuts_WidestBinning_Eta, ME0Muon_Cuts_WidestBinning_Eta, 1, 1, "B");
  FakeRate_WidestBinning_Eta->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  FakeRate_WidestBinning_Eta->GetXaxis()->SetTitleSize(0.05);
  FakeRate_WidestBinning_Eta->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  FakeRate_WidestBinning_Eta->GetYaxis()->SetTitleSize(0.05);
  //FakeRate_WidestBinning_Eta->SetMinimum(FakeRate_WidestBinning_Eta->GetMinimum()-0.1);
  FakeRate_WidestBinning_Eta->SetMinimum(0);
  //FakeRate_WidestBinning_Eta->SetMaximum(FakeRate_WidestBinning_Eta->GetMaximum()+0.1);
  FakeRate_WidestBinning_Eta->SetMaximum(1.2);
  FakeRate_WidestBinning_Eta->Write();   FakeRate_WidestBinning_Eta->Draw();  

  txt->DrawLatex(0.15,0.4,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestFakeRate_WidestBinning_Eta.png');
  c1->Print(histoFolder+"/FakeRate_WidestBinning_Eta.png");


  FakeRate_Eta_5_10->Divide(UnmatchedME0Muon_Cuts_Eta_5_10, ME0Muon_Cuts_Eta_5_10, 1, 1, "B");
  FakeRate_Eta_5_10->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  FakeRate_Eta_5_10->GetXaxis()->SetTitleSize(0.05);
  FakeRate_Eta_5_10->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  FakeRate_Eta_5_10->GetYaxis()->SetTitleSize(0.05);
  //FakeRate_Eta_5_10->SetMinimum(FakeRate_Eta_5_10->GetMinimum()-0.1);
  FakeRate_Eta_5_10->SetMinimum(0);
  //FakeRate_Eta_5_10->SetMaximum(FakeRate_Eta_5_10->GetMaximum()+0.1);
  FakeRate_Eta_5_10->SetMaximum(1.2);
  FakeRate_Eta_5_10->Write();   FakeRate_Eta_5_10->Draw();  

  txt->DrawLatex(0.15,0.4,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestFakeRate_Eta_5_10.png');
  c1->Print(histoFolder+"/FakeRate_Eta_5_10.png");


  FakeRate_Eta_9_11->Divide(UnmatchedME0Muon_Cuts_Eta_9_11, ME0Muon_Cuts_Eta_9_11, 1, 1, "B");
  FakeRate_Eta_9_11->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  FakeRate_Eta_9_11->GetXaxis()->SetTitleSize(0.05);
  FakeRate_Eta_9_11->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  FakeRate_Eta_9_11->GetYaxis()->SetTitleSize(0.05);
  //FakeRate_Eta_9_11->SetMinimum(FakeRate_Eta_9_11->GetMinimum()-0.1);
  FakeRate_Eta_9_11->SetMinimum(0);
  //FakeRate_Eta_9_11->SetMaximum(FakeRate_Eta_9_11->GetMaximum()+0.1);
  FakeRate_Eta_9_11->SetMaximum(1.2);
  FakeRate_Eta_9_11->Write();   FakeRate_Eta_9_11->Draw();  

  txt->DrawLatex(0.15,0.4,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestFakeRate_Eta_9_11.png');
  c1->Print(histoFolder+"/FakeRate_Eta_9_11.png");



  FakeRate_Eta_10_50->Divide(UnmatchedME0Muon_Cuts_Eta_10_50, ME0Muon_Cuts_Eta_10_50, 1, 1, "B");
  FakeRate_Eta_10_50->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  FakeRate_Eta_10_50->GetXaxis()->SetTitleSize(0.05);
  FakeRate_Eta_10_50->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  FakeRate_Eta_10_50->GetYaxis()->SetTitleSize(0.05);
  //FakeRate_Eta_10_50->SetMinimum(FakeRate_Eta_10_50->GetMinimum()-0.1);
  FakeRate_Eta_10_50->SetMinimum(0);
  //FakeRate_Eta_10_50->SetMaximum(FakeRate_Eta_10_50->GetMaximum()+0.1);
  FakeRate_Eta_10_50->SetMaximum(1.2);
  FakeRate_Eta_10_50->Write();   FakeRate_Eta_10_50->Draw();  

  txt->DrawLatex(0.15,0.4,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestFakeRate_Eta_10_50.png');
  c1->Print(histoFolder+"/FakeRate_Eta_10_50.png");



  FakeRate_Eta_50_100->Divide(UnmatchedME0Muon_Cuts_Eta_50_100, ME0Muon_Cuts_Eta_50_100, 1, 1, "B");
  FakeRate_Eta_50_100->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  FakeRate_Eta_50_100->GetXaxis()->SetTitleSize(0.05);
  FakeRate_Eta_50_100->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  FakeRate_Eta_50_100->GetYaxis()->SetTitleSize(0.05);
  //FakeRate_Eta_50_100->SetMinimum(FakeRate_Eta_50_100->GetMinimum()-0.1);
  FakeRate_Eta_50_100->SetMinimum(0);
  //FakeRate_Eta_50_100->SetMaximum(FakeRate_Eta_50_100->GetMaximum()+0.1);
  FakeRate_Eta_50_100->SetMaximum(1.2);
  FakeRate_Eta_50_100->Write();   FakeRate_Eta_50_100->Draw();  

  txt->DrawLatex(0.15,0.4,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestFakeRate_Eta_50_100.png');
  c1->Print(histoFolder+"/FakeRate_Eta_50_100.png");



  FakeRate_Eta_100->Divide(UnmatchedME0Muon_Cuts_Eta_100, ME0Muon_Cuts_Eta_100, 1, 1, "B");
  FakeRate_Eta_100->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  FakeRate_Eta_100->GetXaxis()->SetTitleSize(0.05);
  FakeRate_Eta_100->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  FakeRate_Eta_100->GetYaxis()->SetTitleSize(0.05);
  //FakeRate_Eta_100->SetMinimum(FakeRate_Eta_100->GetMinimum()-0.1);
  FakeRate_Eta_100->SetMinimum(0);
  //FakeRate_Eta_100->SetMaximum(FakeRate_Eta_100->GetMaximum()+0.1);
  FakeRate_Eta_100->SetMaximum(1.2);
  FakeRate_Eta_100->Write();   FakeRate_Eta_100->Draw();  

  txt->DrawLatex(0.15,0.4,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestFakeRate_Eta_100.png');
  c1->Print(histoFolder+"/FakeRate_Eta_100.png");



  Chi2FakeRate_Eta->Divide(Chi2UnmatchedME0Muon_Eta, ME0Muon_Cuts_Eta, 1, 1, "B");
  std::cout<<"Chi2UnmatchedME0Muon_Eta =  "<<Chi2UnmatchedME0Muon_Eta->Integral()<<std::endl;
  std::cout<<"UnmatchedME0Muon_Eta =  "<<UnmatchedME0Muon_Eta->Integral()<<std::endl;
  std::cout<<"ME0Muon_Eta =  "<<ME0Muon_Cuts_Eta->Integral()<<std::endl;
  std::cout<<"ME0Muon_Eta without cuts =  "<<ME0Muon_Eta->Integral()<<std::endl;

  Chi2FakeRate_Eta->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  Chi2FakeRate_Eta->GetXaxis()->SetTitleSize(0.05);
  Chi2FakeRate_Eta->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  Chi2FakeRate_Eta->GetYaxis()->SetTitleSize(0.05);
  Chi2FakeRate_Eta->SetMinimum(0);
  Chi2FakeRate_Eta->SetMaximum(1.2);
  Chi2FakeRate_Eta->Write();   Chi2FakeRate_Eta->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestChi2FakeRate_Eta.png');
  c1->Print(histoFolder+"/Chi2FakeRate_Eta.png");



  Chi2FakeRate_WideBinning_Eta->Divide(Chi2UnmatchedME0Muon_WideBinning_Eta, ME0Muon_Cuts_WideBinning_Eta, 1, 1, "B");
  //std::cout<<"Chi2UnmatchedME0Muon_WideBinning_Eta =  "<<Chi2UnmatchedME0Muon_WideBinning_Eta->Integral()<<std::endl;
  //std::cout<<"UnmatchedME0Muon_WideBinning_Eta =  "<<UnmatchedME0Muon_WideBinning_Eta->Integral()<<std::endl;
  //  std::cout<<"ME0Muon_WideBinning_Eta =  "<<ME0Muon_Cuts_WideBinning_Eta->Integral()<<std::endl;
  //std::cout<<"ME0Muon_WideBinning_Eta without cuts =  "<<ME0Muon_WideBinning_Eta->Integral()<<std::endl;

  Chi2FakeRate_WideBinning_Eta->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  Chi2FakeRate_WideBinning_Eta->GetXaxis()->SetTitleSize(0.05);
  Chi2FakeRate_WideBinning_Eta->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  Chi2FakeRate_WideBinning_Eta->GetYaxis()->SetTitleSize(0.05);
  Chi2FakeRate_WideBinning_Eta->SetMinimum(0);
  Chi2FakeRate_WideBinning_Eta->SetMaximum(1.2);
  Chi2FakeRate_WideBinning_Eta->Write();   Chi2FakeRate_WideBinning_Eta->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestChi2FakeRate_WideBinning_Eta.png');
  c1->Print(histoFolder+"/Chi2FakeRate_WideBinning_Eta.png");


  Chi2FakeRate_WidestBinning_Eta->Divide(Chi2UnmatchedME0Muon_WidestBinning_Eta, ME0Muon_Cuts_WidestBinning_Eta, 1, 1, "B");
  //std::cout<<"Chi2UnmatchedME0Muon_WidestBinning_Eta =  "<<Chi2UnmatchedME0Muon_WidestBinning_Eta->Integral()<<std::endl;
  //std::cout<<"UnmatchedME0Muon_WidestBinning_Eta =  "<<UnmatchedME0Muon_WidestBinning_Eta->Integral()<<std::endl;
  //std::cout<<"ME0Muon_WidestBinning_Eta =  "<<ME0Muon_Cuts_WidestBinning_Eta->Integral()<<std::endl;
  //std::cout<<"ME0Muon_WidestBinning_Eta without cuts =  "<<ME0Muon_WidestBinning_Eta->Integral()<<std::endl;

  Chi2FakeRate_WidestBinning_Eta->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  Chi2FakeRate_WidestBinning_Eta->GetXaxis()->SetTitleSize(0.05);
  Chi2FakeRate_WidestBinning_Eta->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  Chi2FakeRate_WidestBinning_Eta->GetYaxis()->SetTitleSize(0.05);
  Chi2FakeRate_WidestBinning_Eta->SetMinimum(0);
  Chi2FakeRate_WidestBinning_Eta->SetMaximum(1.2);
  Chi2FakeRate_WidestBinning_Eta->Write();   Chi2FakeRate_WidestBinning_Eta->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestChi2FakeRate_WidestBinning_Eta.png');
  c1->Print(histoFolder+"/Chi2FakeRate_WidestBinning_Eta.png");


  //Fake Rate per Event:

  FakeRate_Eta_PerEvent->Divide(UnmatchedME0Muon_Cuts_Eta_PerEvent, ME0Muon_Cuts_Eta_PerEvent, 1, 1, "B");
  std::cout<<"UnmatchedME0Muon_Eta_PerEvent =  "<<UnmatchedME0Muon_Cuts_Eta_PerEvent->Integral()<<std::endl;
  std::cout<<"ME0Muon_Eta_PerEvent =  "<<ME0Muon_Cuts_Eta_PerEvent->Integral()<<std::endl;

  FakeRate_Eta_PerEvent->GetXaxis()->SetTitle("Reconstructed track |#eta|");
  FakeRate_Eta_PerEvent->GetXaxis()->SetTitleSize(0.05);
  FakeRate_Eta_PerEvent->GetYaxis()->SetTitle("ME0 Muon Fake Rate");
  FakeRate_Eta_PerEvent->GetYaxis()->SetTitleSize(0.05);
  FakeRate_Eta_PerEvent->SetMinimum(0);
  FakeRate_Eta_PerEvent->SetMaximum(1.2);
  FakeRate_Eta_PerEvent->Write();   FakeRate_Eta_PerEvent->Draw();  
  txt->DrawLatex(0.15,0.2,pcstr);
  latex->DrawLatex(0.4, 0.85, cmsText);
  latex1->DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  //c1->SaveAs('TestFakeRate_Eta_PerEvent.png');
  c1->Print(histoFolder+"/FakeRate_Eta_PerEvent.png");

  std::cout<<"  FakeRate_Eta values:"<<std::endl;
  for (int i=1; i<=FakeRate_Eta->GetNbinsX(); ++i){
    std::cout<<2.4+(double)i*((4.0-2.4)/40.)<<","<<FakeRate_Eta->GetBinContent(i)<<std::endl;
  }
  

  // FakeRate_Pt->Divide(MatchedME0Muon_Pt, GenMuon_Pt, 1, 1, "B");
  // FakeRate_Pt->GetXaxis()->SetTitle("Gen Muon p_{T}");
  // FakeRate_Pt->GetYaxis()->SetTitle("Matching Efficiency");
  // FakeRate_Pt->SetMinimum(.85);
  // FakeRate_Pt->Write();   FakeRate_Pt->Draw();  c1->Print(histoFolder+"/FakeRate_Pt.png");

  MuonAllTracksEff_Eta->Divide(ME0Muon_Eta, Track_Eta, 1, 1, "B");
  MuonAllTracksEff_Eta->Write();   MuonAllTracksEff_Eta->Draw();  c1->Print(histoFolder+"/MuonAllTracksEff_Eta.png");

  // MuonAllTracksEff_Pt->Divide(ME0Muon_Pt, Track_Pt, 1, 1, "B");
  // MuonAllTracksEff_Pt->Write();   MuonAllTracksEff_Pt->Draw();  c1->Print(histoFolder+"/MuonAllTracksEff_Pt.png");

  // MuonUnmatchedTracksEff_Eta->Divide(UnmatchedME0Muon_Eta, ME0Muon_Eta, 1, 1, "B");
  // MuonUnmatchedTracksEff_Eta->Write();   Candidate_Eta->Draw();  c1->Print(histoFolder+"/Candidate_Eta.png");

  // MuonUnmatchedTracksEff_Pt->Divide(UnmatchedME0Muon_Pt, ME0Muon_Pt, 1, 1, "B");
  // MuonUnmatchedTracksEff_Pt->Write();   Candidate_Eta->Draw();  c1->Print(histoFolder+"/Candidate_Eta.png");
  FractionMatched_Eta->Divide(MatchedME0Muon_Eta, ME0Muon_Eta, 1, 1, "B");
  FractionMatched_Eta->GetXaxis()->SetTitle("Gen Muon |#eta|");
  FractionMatched_Eta->GetYaxis()->SetTitle("Matched/All ME0Muons");
  FractionMatched_Eta->Write();   FractionMatched_Eta->Draw();  c1->Print(histoFolder+"/FractionMatched_Eta.png");

  gStyle->SetOptStat(1);
  PtDiff_h->GetXaxis()->SetTitle("(pt track-ptgen)/ptgen");
  PtDiff_h->Write();     PtDiff_h->Draw();  c1->Print(histoFolder+"/PtDiff_h.png");

  QOverPtDiff_h->GetXaxis()->SetTitle("(q/pt track-q/ptgen)/(q/ptgen)");
  QOverPtDiff_h->Write();     QOverPtDiff_h->Draw();  c1->Print(histoFolder+"/QOverPtDiff_h.png");

  gStyle->SetOptStat(0);
  PtDiff_s->SetMarkerStyle(1);
  PtDiff_s->SetMarkerSize(3.0);
  PtDiff_s->GetXaxis()->SetTitle("Gen Muon #eta");
  PtDiff_s->GetYaxis()->SetTitle("|(pt track-ptgen)/ptgen|");
  PtDiff_s->Write();     PtDiff_s->Draw();  c1->Print(histoFolder+"/PtDiff_s.png");
  
  for(int i=1; i<=PtDiff_p->GetNbinsX(); ++i) {
    PtDiff_rms->SetBinContent(i, PtDiff_p->GetBinError(i)); 
    std::cout<<"Pt_rms = "<<PtDiff_p->GetBinError(i)<<std::endl;
  }


  TH1D *test;
  std::cout<<"Total integral is "<<PtDiff_s->Integral()<<std::endl;
    test= new TH1D("test"   , "pt resolution"   , 200, -1.0, 1.0 );      
  for(Int_t i=1; i<=PtDiff_s->GetNbinsX(); ++i) {


    std::stringstream tempstore;
    tempstore<<i;
    const std::string& thistemp = tempstore.str();
    test->Draw();


    PtDiff_s->ProjectionY("test",i,i,"");
    std::cout<<"Bin = "<<PtDiff_s->GetBinContent(i)<<std::endl;
    std::cout<<"Integral = "<<test->Integral()<<std::endl;
    if (test->Integral() < 1.0) continue;
    std::cout<<"Running some gaussian fits"<<std::endl;

    // TF1 *gaus_narrow = new TF1("gaus_narrow","gaus",-.1,.1);
    // test->Fit(gaus_narrow,"R");

    // Double_t n0  = gaus_narrow->GetParameter(0);
    // Double_t n1  = gaus_narrow->GetParameter(1);
    // Double_t n2  = gaus_narrow->GetParameter(2);

    // //Double_t e_n0  = gaus_narrow->GetParameterError(0);
    // //Double_t e_n1  = gaus_narrow->GetParameterError(1);
    // Double_t e_n2  = gaus_narrow->GetParError(2);

    // std::cout<<n0<<", "<<n1<<", "<<n2<<std::endl;

    //TF1 *gaus_wide = new TF1("gaus_wide","gaus",-.2,.2);
    TF1 *gaus_wide = new TF1("gaus_wide","gaus",-1.,1.);
    std::cout<<"About to fit"<<std::endl;
    test->Fit(gaus_wide,"R");

    std::cout<<"Getting values"<<std::endl;
    Double_t w2  = gaus_wide->GetParameter(2);

    Double_t e_w2  = gaus_wide->GetParError(2);

    std::cout<<"Got values"<<std::endl;
    // PtDiff_gaus_narrow->SetBinContent(i, n2); 
    // PtDiff_gaus_narrow->SetBinError(i, e_n2); 
    PtDiff_gaus_wide->SetBinContent(i, w2); 
    PtDiff_gaus_wide->SetBinError(i, e_w2); 
    
    test->Write();
    TString FileName = "Bin"+thistemp+"Fit.png";
    c1->Print(histoFolder+"/"+FileName);

    //test->Draw();
    //delete test;
    
    //continue;


    delete test;
    test= new TH1D("test"   , "pt resolution"   , 200, -1.0, 1.0 );  
    test->Draw();
    // Redoing for pt 5 to 10
    std::cout<<"About to project"<<std::endl;
    PtDiff_s_5_10->ProjectionY("test",i,i,"");
    std::cout<<"About to check, "<<std::endl;
    std::cout<<test->Integral()<<std::endl;
    if (test->Integral() < 1.0) continue;

    std::cout<<"Running the 5-10 fit"<<std::endl;
    TF1 *gaus_5_10 = new TF1("gaus_5_10","gaus",-.2,.2);
    test->Fit(gaus_5_10,"R");

     w2  = gaus_5_10->GetParameter(2);
     e_w2  = gaus_5_10->GetParError(2);

    PtDiff_gaus_5_10->SetBinContent(i, w2); 
    PtDiff_gaus_5_10->SetBinError(i, e_w2); 

    test->Draw();
    FileName = "Bin"+thistemp+"Fit_5_10.png";
    c1->Print(histoFolder+"/"+FileName);

    delete test;
    test= new TH1D("test"   , "pt resolution"   , 200, -1.0, 1.0 );  
    test->Draw();
    // Redoing for pt 10 to 20
    PtDiff_s_10_50->ProjectionY("test",i,i,"");
    if (test->Integral() < 1.0) continue;

    TF1 *gaus_10_50 = new TF1("gaus_10_50","gaus",-.2,.2);
    test->Fit(gaus_10_50,"R");

     w2  = gaus_10_50->GetParameter(2);
     e_w2  = gaus_10_50->GetParError(2);

    PtDiff_gaus_10_50->SetBinContent(i, w2); 
    PtDiff_gaus_10_50->SetBinError(i, e_w2); 

    test->Draw();
    FileName = "Bin"+thistemp+"Fit_10_50.png";
    c1->Print(histoFolder+"/"+FileName);

    delete test;

    test= new TH1D("test"   , "pt resolution"   , 200, -1.0, 1.0 );  
    test->Draw();
    // Redoing for pt 20 to 40
    PtDiff_s_50_100->ProjectionY("test",i,i,"");
    if (test->Integral() < 1.0) continue;

    TF1 *gaus_50_100 = new TF1("gaus_50_100","gaus",-.2,.2);
    test->Fit(gaus_50_100,"R");

     w2  = gaus_50_100->GetParameter(2);
     e_w2  = gaus_50_100->GetParError(2);

    PtDiff_gaus_50_100->SetBinContent(i, w2); 
    PtDiff_gaus_50_100->SetBinError(i, e_w2); 

    test->Draw();
    FileName = "Bin"+thistemp+"Fit_50_100.png";
    c1->Print(histoFolder+"/"+FileName);

    delete test;

    test= new TH1D("test"   , "pt resolution"   , 200, -1.0, 1.0 );  
    test->Draw();
    // Redoing for pt 40+
    PtDiff_s_100->ProjectionY("test",i,i,"");
    if (test->Integral() < 1.0) continue;

    TF1 *gaus_100 = new TF1("gaus_100","gaus",-.2,.2);
    test->Fit(gaus_100,"R");

     w2  = gaus_100->GetParameter(2);
     e_w2  = gaus_100->GetParError(2);

    PtDiff_gaus_100->SetBinContent(i, w2); 
    PtDiff_gaus_100->SetBinError(i, e_w2); 

    test->Draw();
    FileName = "Bin"+thistemp+"Fit_100.png";
    c1->Print(histoFolder+"/"+FileName);

    delete test;

    test= new TH1D("test"   , "pt resolution"   , 200, -1.0, 1.0 );  
    test->Draw();
    // Redoing for pt 40+
    StandalonePtDiff_s->ProjectionY("test",i,i,"");
    if (test->Integral() < 1.0) continue;

    TF1 *Standalonegaus = new TF1("Standalonegaus","gaus",-.2,.2);
    test->Fit(Standalonegaus,"R");

     w2  = gaus_100->GetParameter(2);
     e_w2  = gaus_100->GetParError(2);

     StandalonePtDiff_gaus->SetBinContent(i, w2); 
     StandalonePtDiff_gaus->SetBinError(i, e_w2); 

     test->Draw();
     FileName = "Bin"+thistemp+"StandaloneFit.png";
     c1->Print(histoFolder+"/"+FileName);
     
     delete test;
     //test->Clear();
  }
  // PtDiff_gaus_narrow->SetMarkerStyle(22); 
  // PtDiff_gaus_narrow->SetMarkerSize(1.2); 
  // PtDiff_gaus_narrow->SetMarkerColor(kRed); 
  // //PtDiff_gaus_narrow->SetLineColor(kRed); 
  
  // //PtDiff_gaus_narrow->Draw("PL"); 

  // PtDiff_gaus_narrow->GetXaxis()->SetTitle("Gen Muon #eta");
  // PtDiff_gaus_narrow->GetYaxis()->SetTitle("Gaussian width of (pt track-ptgen)/ptgen");
  // PtDiff_gaus_narrow->Write();     PtDiff_gaus_narrow->Draw("PE");

  PtDiff_gaus_wide->SetMarkerStyle(22); 
  PtDiff_gaus_wide->SetMarkerSize(1.2); 
  PtDiff_gaus_wide->SetMarkerColor(kBlue); 
  //PtDiff_gaus_wide->SetLineColor(kRed); 
  
  //PtDiff_gaus_wide->Draw("PL"); 

  PtDiff_gaus_wide->GetXaxis()->SetTitle("Gen Muon |#eta|");
  PtDiff_gaus_wide->GetYaxis()->SetTitle("Gaussian width of (pt track-ptgen)/ptgen");
  PtDiff_gaus_wide->GetYaxis()->SetTitleSize(.04);
  PtDiff_gaus_wide->Write();     PtDiff_gaus_wide->Draw("PE");  c1->Print(histoFolder+"/PtDiff_gaus.png");

  PtDiff_gaus_5_10->SetMarkerStyle(22); 
  PtDiff_gaus_5_10->SetMarkerSize(1.2); 
  PtDiff_gaus_5_10->SetMarkerColor(kBlue); 
  //PtDiff_gaus_5_10->SetLineColor(kRed); 
  
  //PtDiff_gaus_5_10->Draw("PL"); 

  PtDiff_gaus_5_10->GetXaxis()->SetTitle("Gen Muon |#eta|");
  PtDiff_gaus_5_10->GetYaxis()->SetTitle("Gaussian width of (pt track-ptgen)/ptgen");
  PtDiff_gaus_5_10->Write();     PtDiff_gaus_5_10->Draw("PE");  c1->Print(histoFolder+"/PtDiff_gaus_5_10.png");

  PtDiff_gaus_10_50->SetMarkerStyle(22); 
  PtDiff_gaus_10_50->SetMarkerSize(1.2); 
  PtDiff_gaus_10_50->SetMarkerColor(kBlue); 
  //PtDiff_gaus_10_50->SetLineColor(kRed); 
  
  //PtDiff_gaus_10_50->Draw("PL"); 

  PtDiff_gaus_10_50->GetXaxis()->SetTitle("Gen Muon |#eta|");
  PtDiff_gaus_10_50->GetYaxis()->SetTitle("Gaussian width of (pt track-ptgen)/ptgen");
  PtDiff_gaus_10_50->Write();     PtDiff_gaus_10_50->Draw("PE");  c1->Print(histoFolder+"/PtDiff_gaus_10_50.png");

  PtDiff_gaus_50_100->SetMarkerStyle(22); 
  PtDiff_gaus_50_100->SetMarkerSize(1.2); 
  PtDiff_gaus_50_100->SetMarkerColor(kBlue); 
  //PtDiff_gaus_50_100->SetLineColor(kRed); 
  
  //PtDiff_gaus_50_100->Draw("PL"); 

  PtDiff_gaus_50_100->GetXaxis()->SetTitle("Gen Muon |#eta|");
  PtDiff_gaus_50_100->GetYaxis()->SetTitle("Gaussian width of (pt track-ptgen)/ptgen");
  PtDiff_gaus_50_100->Write();     PtDiff_gaus_50_100->Draw("PE");  c1->Print(histoFolder+"/PtDiff_gaus_50_100.png");

  PtDiff_gaus_100->SetMarkerStyle(22); 
  PtDiff_gaus_100->SetMarkerSize(1.2); 
  PtDiff_gaus_100->SetMarkerColor(kBlue); 
  //PtDiff_gaus_100->SetLineColor(kRed); 
  
  //PtDiff_gaus_100->Draw("PL"); 

  PtDiff_gaus_100->GetXaxis()->SetTitle("Gen Muon |#eta|");
  PtDiff_gaus_100->GetYaxis()->SetTitle("Gaussian width of (pt track-ptgen)/ptgen");
  PtDiff_gaus_100->Write();     PtDiff_gaus_100->Draw("PE");  c1->Print(histoFolder+"/PtDiff_gaus_100.png");

  StandalonePtDiff_gaus->SetMarkerStyle(22); 
  StandalonePtDiff_gaus->SetMarkerSize(1.2); 
  StandalonePtDiff_gaus->SetMarkerColor(kBlue); 
  //StandalonePtDiff_gaus->SetLineColor(kRed); 
  
  //StandalonePtDiff_gaus->Draw("PL"); 

  StandalonePtDiff_gaus->GetXaxis()->SetTitle("Gen Muon |#eta|");
  StandalonePtDiff_gaus->GetYaxis()->SetTitle("Gaussian width of (pt track-ptgen)/ptgen");
  StandalonePtDiff_gaus->Write();     StandalonePtDiff_gaus->Draw("PE");  c1->Print(histoFolder+"/StandalonePtDiff_gaus.png");


  PtDiff_p->SetMarkerStyle(1);
  PtDiff_p->SetMarkerSize(3.0);
  PtDiff_p->GetXaxis()->SetTitle("Gen Muon |#eta|");
  PtDiff_p->GetYaxis()->SetTitle("Average (pt track-ptgen)/ptgen");
  PtDiff_p->Write();     PtDiff_p->Draw();  c1->Print(histoFolder+"/PtDiff_p.png");

  //PtDiff_rms->SetMarkerStyle(1);
  //PtDiff_rms->SetMarkerSize(3.0);

  PtDiff_rms->SetMarkerStyle(22); 
  PtDiff_rms->SetMarkerSize(1.2); 
  PtDiff_rms->SetMarkerColor(kBlue); 
  //PtDiff_rms->SetLineColor(kRed); 
  
  //PtDiff_rms->Draw("PL"); 

  PtDiff_rms->GetXaxis()->SetTitle("Gen Muon |#eta|");
  PtDiff_rms->GetYaxis()->SetTitle("RMS of (pt track-ptgen)/ptgen");
  PtDiff_rms->Write();     PtDiff_rms->Draw("P");  c1->Print(histoFolder+"/PtDiff_rms.png");

  PDiff_h->GetXaxis()->SetTitle("|(p track-pgen)/pgen|");
  PDiff_h->Write();     PDiff_h->Draw();  c1->Print(histoFolder+"/PDiff_h.png");

  PDiff_s->SetMarkerStyle(1);
  PDiff_s->SetMarkerSize(3.0);
  PDiff_s->GetXaxis()->SetTitle("Gen Muon |#eta|");
  PDiff_s->GetYaxis()->SetTitle("|(p track-pgen)/pgen|");
  PDiff_s->Write();     PDiff_s->Draw();  c1->Print(histoFolder+"/PDiff_s.png");
  
  PDiff_p->SetMarkerStyle(1);
  PDiff_p->SetMarkerSize(3.0);
  PDiff_p->GetXaxis()->SetTitle("Gen Muon #eta");
  PDiff_p->GetYaxis()->SetTitle("Average |(p track-pgen)/pgen|");
  PDiff_p->Write();     PDiff_p->Draw();  c1->Print(histoFolder+"/PDiff_p.png");

  VertexDiff_h->Write();     VertexDiff_h->Draw();  c1->Print(histoFolder+"/VertexDiff_h.png");



  NormChi2_h->Write();       NormChi2_h->Draw();          c1->Print(histoFolder+"/NormChi2_h.png");
  NormChi2Prob_h->Write();    NormChi2Prob_h->Draw();      c1->Print(histoFolder+"/NormChi2Prob_h.png");
  NormChi2VsHits_h->Write();    NormChi2VsHits_h->Draw(); c1->Print(histoFolder+"/NormChi2VsHits_h.png");
  chi2_vs_eta_h->Write();      chi2_vs_eta_h->Draw();     c1->Print(histoFolder+"/chi2_vs_eta_h.png");

  c1->SetLogy(0);
  AssociatedChi2_h->Write();    AssociatedChi2_h->Draw();   c1->Print(histoFolder+"/AssociatedChi2_h.png");
  AssociatedChi2_Prob_h->Write();    AssociatedChi2_Prob_h->Draw();   c1->Print(histoFolder+"/AssociatedChi2_Prob_h.png");

  //Printing needing values to an output log file
  using namespace std;
  ofstream logout;
  logout.open (histoFolder+"/Log.txt");

  logout<<"Chi 2 Efficiencies and errors:\n";
  for (int i=1; i<=Chi2MuonRecoEff_Eta->GetNbinsX(); ++i){
    logout<<Chi2MuonRecoEff_Eta->GetBinContent(i)<<","<<Chi2MuonRecoEff_Eta->GetBinError(i)<<"\n";
  }    

  logout<<"Efficiencies and errors:\n";
  for (int i=1; i<=MuonRecoEff_Eta->GetNbinsX(); ++i){
    logout<<MuonRecoEff_Eta->GetBinContent(i)<<","<<MuonRecoEff_Eta->GetBinError(i)<<"\n";
  }    

  logout<<"Fake Rate:\n";
  for (int i=1; i<=FakeRate_Eta->GetNbinsX(); ++i){
    logout<<FakeRate_Eta->GetBinContent(i)<<","<<FakeRate_Eta->GetBinError(i)<<"\n";
  }    

  logout<<"Resolution vs eta:\n";
  for (int i=1; i<=PtDiff_gaus_wide->GetNbinsX(); ++i){
    logout<<PtDiff_gaus_wide->GetBinContent(i)<<","<<PtDiff_gaus_wide->GetBinError(i)<<"\n";
  }    


  logout<<"Efficiencies and errors 5_10:\n";
  for (int i=1; i<=MuonRecoEff_Eta_5_10->GetNbinsX(); ++i){
    logout<<MuonRecoEff_Eta_5_10->GetBinContent(i)<<","<<MuonRecoEff_Eta_5_10->GetBinError(i)<<"\n";
  }    


  logout<<"Efficiencies and errors 9_11:\n";
  for (int i=1; i<=MuonRecoEff_Eta_9_11->GetNbinsX(); ++i){
    logout<<MuonRecoEff_Eta_9_11->GetBinContent(i)<<","<<MuonRecoEff_Eta_9_11->GetBinError(i)<<"\n";
  }    


  logout<<"Chi 2 Efficiencies and errors 5_10:\n";
  for (int i=1; i<=Chi2MuonRecoEff_Eta_5_10->GetNbinsX(); ++i){
    logout<<Chi2MuonRecoEff_Eta_5_10->GetBinContent(i)<<","<<Chi2MuonRecoEff_Eta_5_10->GetBinError(i)<<"\n";
  }    

  logout<<"Fake Rate 5_10:\n";
  for (int i=1; i<=FakeRate_Eta_5_10->GetNbinsX(); ++i){
    logout<<FakeRate_Eta_5_10->GetBinContent(i)<<","<<FakeRate_Eta_5_10->GetBinError(i)<<"\n";
  }    

  logout<<"Resolution vs eta 5_10:\n";
  for (int i=1; i<=PtDiff_gaus_5_10->GetNbinsX(); ++i){
    logout<<PtDiff_gaus_5_10->GetBinContent(i)<<","<<PtDiff_gaus_5_10->GetBinError(i)<<"\n";
  }    


  logout<<"Efficiencies and errors 10_50:\n";
  for (int i=1; i<=MuonRecoEff_Eta_10_50->GetNbinsX(); ++i){
    logout<<MuonRecoEff_Eta_10_50->GetBinContent(i)<<","<<MuonRecoEff_Eta_10_50->GetBinError(i)<<"\n";
  }    

  logout<<"Chi 2 Efficiencies and errors 10_50:\n";
  for (int i=1; i<=Chi2MuonRecoEff_Eta_10_50->GetNbinsX(); ++i){
    logout<<Chi2MuonRecoEff_Eta_10_50->GetBinContent(i)<<","<<Chi2MuonRecoEff_Eta_10_50->GetBinError(i)<<"\n";
  }    


  logout<<"Fake Rate 10_50:\n";
  for (int i=1; i<=FakeRate_Eta_10_50->GetNbinsX(); ++i){
    logout<<FakeRate_Eta_10_50->GetBinContent(i)<<","<<FakeRate_Eta_10_50->GetBinError(i)<<"\n";
  }    

  logout<<"Resolution vs eta 10_50:\n";
  for (int i=1; i<=PtDiff_gaus_10_50->GetNbinsX(); ++i){
    logout<<PtDiff_gaus_10_50->GetBinContent(i)<<","<<PtDiff_gaus_10_50->GetBinError(i)<<"\n";
  }    


  logout<<"Efficiencies and errors 50_100:\n";
  for (int i=1; i<=MuonRecoEff_Eta_50_100->GetNbinsX(); ++i){
    logout<<MuonRecoEff_Eta_50_100->GetBinContent(i)<<","<<MuonRecoEff_Eta_50_100->GetBinError(i)<<"\n";
  }    


  logout<<"Chi 2 Efficiencies and errors 50_100:\n";
  for (int i=1; i<=Chi2MuonRecoEff_Eta_50_100->GetNbinsX(); ++i){
    logout<<Chi2MuonRecoEff_Eta_50_100->GetBinContent(i)<<","<<Chi2MuonRecoEff_Eta_50_100->GetBinError(i)<<"\n";
  }    

  logout<<"Fake Rate 50_100:\n";
  for (int i=1; i<=FakeRate_Eta_50_100->GetNbinsX(); ++i){
    logout<<FakeRate_Eta_50_100->GetBinContent(i)<<","<<FakeRate_Eta_50_100->GetBinError(i)<<"\n";
  }    

  logout<<"Resolution vs eta 50_100:\n";
  for (int i=1; i<=PtDiff_gaus_50_100->GetNbinsX(); ++i){
    logout<<PtDiff_gaus_50_100->GetBinContent(i)<<","<<PtDiff_gaus_50_100->GetBinError(i)<<"\n";
  }    


  logout<<"Efficiencies and errors 40:\n";
  for (int i=1; i<=MuonRecoEff_Eta_100->GetNbinsX(); ++i){
    logout<<MuonRecoEff_Eta_100->GetBinContent(i)<<","<<MuonRecoEff_Eta_100->GetBinError(i)<<"\n";
  }    


  logout<<"Chi 2 Efficiencies and errors 40:\n";
  for (int i=1; i<=Chi2MuonRecoEff_Eta_100->GetNbinsX(); ++i){
    logout<<Chi2MuonRecoEff_Eta_100->GetBinContent(i)<<","<<Chi2MuonRecoEff_Eta_100->GetBinError(i)<<"\n";
  }    

  logout<<"Fake Rate 40:\n";
  for (int i=1; i<=FakeRate_Eta_100->GetNbinsX(); ++i){
    logout<<FakeRate_Eta_100->GetBinContent(i)<<","<<FakeRate_Eta_100->GetBinError(i)<<"\n";
  }    

  logout<<"Resolution vs eta 40:\n";
  for (int i=1; i<=PtDiff_gaus_100->GetNbinsX(); ++i){
    logout<<PtDiff_gaus_100->GetBinContent(i)<<","<<PtDiff_gaus_100->GetBinError(i)<<"\n";
  }    


  logout<<"Background yield:\n";
  for (int i=1; i<=UnmatchedME0Muon_Cuts_Eta_PerEvent->GetNbinsX(); ++i){
    logout<<UnmatchedME0Muon_Cuts_Eta_PerEvent->GetBinContent(i)<<","<<UnmatchedME0Muon_Cuts_Eta_PerEvent->GetBinError(i)<<"\n";
  }    

  
  //logout << "Writing this to a file.\n";
  logout.close();

  delete histoFile; histoFile = 0;
}



FreeTrajectoryState
ME0MuonAnalyzer::getFTS(const GlobalVector& p3, const GlobalVector& r3, 
			   int charge, const AlgebraicSymMatrix55& cov,
			   const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CurvilinearTrajectoryError tCov(cov);
  
  return cov.kRows == 5 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState
ME0MuonAnalyzer::getFTS(const GlobalVector& p3, const GlobalVector& r3, 
			   int charge, const AlgebraicSymMatrix66& cov,
			   const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CartesianTrajectoryError tCov(cov);
  
  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

void ME0MuonAnalyzer::getFromFTS(const FreeTrajectoryState& fts,
				    GlobalVector& p3, GlobalVector& r3, 
				    int& charge, AlgebraicSymMatrix66& cov){
  GlobalVector p3GV = fts.momentum();
  GlobalPoint r3GP = fts.position();

  GlobalVector p3T(p3GV.x(), p3GV.y(), p3GV.z());
  GlobalVector r3T(r3GP.x(), r3GP.y(), r3GP.z());
  p3 = p3T;
  r3 = r3T;  //Yikes, was setting this to p3T instead of r3T!?!
  // p3.set(p3GV.x(), p3GV.y(), p3GV.z());
  // r3.set(r3GP.x(), r3GP.y(), r3GP.z());
  
  charge = fts.charge();
  cov = fts.hasError() ? fts.cartesianError().matrix() : AlgebraicSymMatrix66();

}

DEFINE_FWK_MODULE(ME0MuonAnalyzer);
