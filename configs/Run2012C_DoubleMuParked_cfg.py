import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
process = cms.Process("AOD2NanoAOD")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "INFO"
process.MessageLogger.categories.append("AOD2NanoAOD")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

# Set the maximum number of events to be processed (-1 processes all events)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

# Define files of dataset
files = FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10000_file_index.txt")
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10001_file_index.txt"))
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10002_file_index.txt"))
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10003_file_index.txt"))
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10010_file_index.txt"))
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10011_file_index.txt"))
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10013_file_index.txt"))
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10016_file_index.txt"))
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10018_file_index.txt"))
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10021_file_index.txt"))
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10022_file_index.txt"))
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_10024_file_index.txt"))
files.extend(FileUtils.loadListFromFile("data/CMS_Run2012C_DoubleMuParked_AOD_22Jan2013-v1_20000_file_index.txt"))

process.source = cms.Source(
    "PoolSource", fileNames=cms.untracked.vstring(*files))

# Apply JSON file with lumi mask (needs to be done after the process.source definition)
goodJSON = "data/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt"
myLumis = LumiList.LumiList(filename=goodJSON).getCMSSWString().split(",")
process.source.lumisToProcess = CfgTypes.untracked(
    CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)

# Number of events to be skipped (0 by default)
process.source.skipEvents = cms.untracked.uint32(0)

# Register fileservice for output file
process.aod2nanoaod = cms.EDAnalyzer("AOD2NanoAOD")
process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("Run2012C_DoubleMuParked.root"))

process.p = cms.Path(process.aod2nanoaod)