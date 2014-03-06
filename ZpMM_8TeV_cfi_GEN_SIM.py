# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: ZpMM_8TeV_cfi --conditions auto:upgrade2019 -n 10 --eventcontent FEVTDEBUG --relval 9000,200 -s GEN,SIM --datatier GEN-SIM --beamspot Gauss --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_phase2_BE5DPixel10D --geometry ExtendedPhase2TkBE5DPixel10D --magField 38T_PostLS1 --fileout file:step1.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('ZpMM_8TeV_cfi nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:step1.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(8000.0),
    crossSection = cms.untracked.double(2.497e-05),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'PARP(82)=1.921 ! pt cutoff for multiparton interactions', 
            'PARP(89)=1800. ! sqrts for which PARP82 is set', 
            'PARP(90)=0.227 ! Multiple interactions: rescaling power', 
            'MSTP(95)=6     ! CR (color reconnection parameters)', 
            'PARP(77)=1.016 ! CR', 
            'PARP(78)=0.538 ! CR', 
            'PARP(80)=0.1   ! Prob. colored parton from BBR', 
            'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
            'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 
            'PARP(62)=1.025 ! ISR cutoff', 
            'MSTP(91)=1     ! Gaussian primordial kT', 
            'PARP(93)=10.0  ! primordial kT-max', 
            'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model'),
        processParameters = cms.vstring('MSEL=0             ! User defined processes', 
            "MSUB(141)   = 1    ! ff -> gamma/Z0/Z\'", 
            'MSTP(44) = 3       ! only select the Z process', 
            "PMAS(32,1)  = 3000 ! Z\' mass (GeV)", 
            'CKIN(1)     = -1  ! lower invariant mass cutoff (GeV)', 
            'CKIN(2)     = -1   ! no upper invariant mass cutoff', 
            'PARU(121)=  0.        ! vd', 
            'PARU(122)=  0.506809  ! ad', 
            'PARU(123)=  0.        ! vu', 
            'PARU(124)=  0.506809  ! au', 
            'PARU(125)=  0.        ! ve', 
            'PARU(126)=  0.506809  ! ae', 
            'PARU(127)= -0.253405  ! vnu', 
            'PARU(128)=  0.253405  ! anu', 
            'PARJ(180)=  0.        ! vd', 
            'PARJ(181)=  0.506809  ! ad', 
            'PARJ(182)=  0.        ! vu', 
            'PARJ(183)=  0.506809  ! au', 
            'PARJ(184)=  0.        ! ve', 
            'PARJ(185)=  0.506809  ! ae', 
            'PARJ(186)= -0.253405  ! vnu', 
            'PARJ(187)=  0.253405  ! anu', 
            'PARJ(188)=  0.        ! vd', 
            'PARJ(189)=  0.506809  ! ad', 
            'PARJ(190)=  0.        ! vu', 
            'PARJ(191)=  0.506809  ! au', 
            'PARJ(192)=  0.        ! ve', 
            'PARJ(193)=  0.506809  ! ae', 
            'PARJ(194)= -0.253405  ! vnu', 
            'PARJ(195)=  0.253405  ! anu', 
            'MDME(289,1) = 0    ! d dbar', 
            'MDME(290,1) = 0    ! u ubar', 
            'MDME(291,1) = 0    ! s sbar', 
            'MDME(292,1) = 0    ! c cbar', 
            'MDME(293,1) = 0    ! b bar', 
            'MDME(294,1) = 0    ! t tbar', 
            'MDME(295,1) = -1   ! 4th gen q qbar', 
            'MDME(296,1) = -1   ! 4th gen q qbar', 
            'MDME(297,1) = 0    ! e-     e+', 
            'MDME(298,1) = 0    ! nu_e   nu_ebar', 
            'MDME(299,1) = 1    ! mu-    mu+', 
            'MDME(300,1) = 0    ! nu_mu  nu_mubar', 
            'MDME(301,1) = 0    ! tau    tau', 
            'MDME(302,1) = 0    ! nu_tau nu_taubar', 
            'MDME(303,1) = -1   ! 4th gen l- l+', 
            'MDME(304,1) = -1   ! 4th gen nu nubar', 
            'MDME(305,1) = -1   ! W+ W-', 
            'MDME(306,1) = -1   ! H+ H-', 
            'MDME(307,1) = -1   ! Z0 gamma', 
            'MDME(308,1) = -1   ! Z0 h0', 
            'MDME(309,1) = -1   ! h0 A0', 
            'MDME(310,1) = -1   ! H0 A0'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_phase2_BE5DPixel10D 

#call to customisation function cust_phase2_BE5DPixel10D imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_phase2_BE5DPixel10D(process)

# End of customisation functions
