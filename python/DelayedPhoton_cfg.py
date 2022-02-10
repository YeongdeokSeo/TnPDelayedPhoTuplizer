import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #"file:MyOutputFile.root"
        #'/store/data/Run2018B/EGamma/MINIAOD/PromptReco-v1/000/317/696/00000/FC9C403C-9970-E811-90FE-FA163E1FAB71.root'
        '/store/data/Run2018B/EGamma/MINIAOD/17Sep2018-v1/60000/FC233A84-73D8-054E-B774-CE80543407D8.root'                
    )
                                                                                            )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#TFileService for output 
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("razorNtuple.root"),
    closeFileFast = cms.untracked.bool(True)
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

#------ Declare the correct global tag ------#

#Global Tag for Run2015D
#process.GlobalTag.globaltag = '101X_dataRun2_Prompt_v11'
process.GlobalTag.globaltag = '102X_dataRun2_v12'

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.taggingMode = cms.bool(True)

#------ Electron MVA Setup ------#
# define which IDs we want to produce
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff','RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


#------ Analyzer ------#

process.ntuples = cms.EDAnalyzer('DelayedAnalyzer',
    photons = cms.VInputTag(cms.InputTag("slimmedPhotons")),
    photonHLTFilterNamesFile = cms.string("RazorNtuple/DelayedAnalyzer/data/RazorPhotonHLTFilterNames2018.dat"),

    mets = cms.InputTag("slimmedMETs"),

    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    triggerObjects = cms.InputTag("slimmedPatTrigger"),

    ebRecHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "RECO"),
    eeRecHits = cms.InputTag("reducedEgamma", "reducedEERecHits", "RECO"),
)

#------ run ------#
process.p = cms.Path(process.ntuples)
