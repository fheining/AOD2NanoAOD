// -*- C++ -*-
//
// Package:    AOD2NanoAOD
// Class:      AOD2NanoAOD

#include <memory>
#include <vector>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "math.h"

#include "TFile.h"
#include "TTree.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

class AOD2NanoAOD : public edm::EDAnalyzer {
public:
  explicit AOD2NanoAOD(const edm::ParameterSet &);
  ~AOD2NanoAOD();

private:
	virtual void beginJob();
	virtual void analyze(const edm::Event &, const edm::EventSetup &);
	virtual void endJob();
	bool providesGoodLumisection(const edm::Event &iEvent);
	bool isDataFlag_;
	
	TTree *tree;

	// Event information
	Int_t value_run;
	UInt_t value_lumi_block;
	ULong64_t value_event;

	// Vertices
	int value_ve_n;
	float value_ve_x;
	float value_ve_y;
	float value_ve_z;

	// Muons
	const static int max_mu = 100;
	UInt_t value_mu_n;
	float value_mu_pt[max_mu];
	float value_mu_eta[max_mu];
	float value_mu_phi[max_mu];
	float value_mu_mass[max_mu];
	int value_mu_charge[max_mu];
	float value_mu_pfreliso03all[max_mu];
	float value_mu_pfreliso04all[max_mu];
	bool value_mu_tightid[max_mu];
	bool value_mu_softid[max_mu];
	float value_mu_dxy[max_mu];
	float value_mu_dxyErr[max_mu];
	float value_mu_dz[max_mu];
	float value_mu_dzErr[max_mu];

	// Electrons
	const static int max_el = 100;
	UInt_t value_el_n;
	float value_el_pt[max_el];
	float value_el_eta[max_el];
	float value_el_phi[max_el];
	float value_el_mass[max_el];
	int value_el_charge[max_el];
	float value_el_pfreliso03all[max_el];
	float value_el_dxy[max_el];
	float value_el_dxyErr[max_el];
	float value_el_dz[max_el];
	float value_el_dzErr[max_el];

	// Taus
	const static int max_tau = 100;
	UInt_t value_tau_n;
	float value_tau_pt[max_tau];
	float value_tau_eta[max_tau];
	float value_tau_phi[max_tau];
	float value_tau_mass[max_tau];
	int value_tau_charge[max_tau];
	int value_tau_decaymode[max_tau];
	float value_tau_chargediso[max_tau];
	float value_tau_neutraliso[max_tau];
	float value_tau_relisochg[max_tau]; //build additional branches  
	float value_tau_relisoneut[max_tau];
	
	bool value_tau_matched_tau[max_tau];
	bool value_tau_matched_mu[max_tau];
	bool value_tau_matched_el[max_tau];
	
		//Electron rejection
		float value_tau_emFraction[max_tau];	// HCAL/ECAL
		float value_tau_hcalTotOverPLead[max_tau];
		float value_tau_hcalMaxOverPLead[max_tau];
		float value_tau_hcal3x3OverPLead[max_tau];
		float value_tau_ecalStripSumEOverPLead[max_tau];
		float value_tau_bremsRecoveryEOverPLead[max_tau];
		//reco::TrackRef value_electronPreIDTrack[max_tau];
		float value_tau_electronPreIDOutput[max_tau];
		bool value_tau_electronPreIDDecision[max_tau];
		

		//Muon rejection
		bool value_tau_hasMuonReference[max_tau];
		float value_tau_caloComp[max_tau];
		float value_tau_segComp[max_tau];
		bool value_tau_muonDecision[max_tau];



	// Photons
	const static int max_ph = 100;
	UInt_t value_ph_n;
	float value_ph_pt[max_ph];
	float value_ph_eta[max_ph];
	float value_ph_phi[max_ph];
	float value_ph_mass[max_ph];
	int value_ph_charge[max_ph];
	float value_ph_pfreliso03all[max_ph];

	// MET
	float value_met_pt;
	float value_met_phi;
	float value_met_sumet;
	float value_met_significance;
	float value_met_covxx;
	float value_met_covxy;
	float value_met_covyy;

	// Jets
	const static int max_jet = 300;
	UInt_t value_jet_n;
	float value_jet_pt[max_jet];
	float value_jet_eta[max_jet];
	float value_jet_phi[max_jet];
	float value_jet_mass[max_jet];

	// Generator particles

		//Electrons
	const static int max_gen = 100;
		


		UInt_t value_gen_el_n;
		
		std::vector<float> value_gen_el_pt;
                std::vector<float>* pvalue_gen_el_pt = &value_gen_el_pt;

                std::vector<float> value_gen_el_eta;
                std::vector<float>* pvalue_gen_el_eta = &value_gen_el_eta;

                std::vector<float> value_gen_el_phi;
                std::vector<float>* pvalue_gen_el_phi = &value_gen_el_phi;

                std::vector<float> value_gen_el_mass;
                std::vector<float>* pvalue_gen_el_mass = &value_gen_el_mass;
	
		
		//Muons
		UInt_t value_gen_mu_n;
		
		std::vector<float> value_gen_mu_pt;
	      	std::vector<float>* pvalue_gen_mu_pt = &value_gen_mu_pt;		

		std::vector<float> value_gen_mu_eta;
                std::vector<float>* pvalue_gen_mu_eta = &value_gen_mu_eta;

		std::vector<float> value_gen_mu_phi;
                std::vector<float>* pvalue_gen_mu_phi = &value_gen_mu_phi;
	
		std::vector<float> value_gen_mu_mass;
                std::vector<float>* pvalue_gen_mu_mass = &value_gen_mu_mass;

		
		//Taus
	
		UInt_t value_gen_tau_n;	

		std::vector<float> value_gen_tau_pt;
                std::vector<float>* pvalue_gen_tau_pt = &value_gen_tau_pt;

                std::vector<float> value_gen_tau_eta;
                std::vector<float>* pvalue_gen_tau_eta = &value_gen_tau_eta;

                std::vector<float> value_gen_tau_phi;
                std::vector<float>* pvalue_gen_tau_phi = &value_gen_tau_phi;

                std::vector<float> value_gen_tau_mass;
                std::vector<float>* pvalue_gen_tau_mass = &value_gen_tau_mass;
 
};

AOD2NanoAOD::AOD2NanoAOD(const edm::ParameterSet &iConfig):
	isDataFlag_(iConfig.getParameter<bool>("isData"))
{
	edm::Service<TFileService> fs;

	tree = fs->make<TTree>("Events", "Events");

	// Event information
	tree->Branch("run", &value_run);
	tree->Branch("luminosityBlock", &value_lumi_block);
	tree->Branch("event", &value_event);

	// Vertices
	tree->Branch("PV_npvs", &value_ve_n, "PV_npvs/I");
	tree->Branch("PV_x", &value_ve_x, "PV_x/F");
	tree->Branch("PV_y", &value_ve_y, "PV_y/F");
	tree->Branch("PV_z", &value_ve_z, "PV_z/F");

	// Muons
	tree->Branch("nMuon", &value_mu_n, "nMuon/i");
	tree->Branch("Muon_pt", value_mu_pt, "Muon_pt[nMuon]/F");
	tree->Branch("Muon_eta", value_mu_eta, "Muon_eta[nMuon]/F");
	tree->Branch("Muon_phi", value_mu_phi, "Muon_phi[nMuon]/F");
	tree->Branch("Muon_mass", value_mu_mass, "Muon_mass[nMuon]/F");
	tree->Branch("Muon_charge", value_mu_charge, "Muon_charge[nMuon]/I");
	tree->Branch("Muon_pfRelIso03_all", value_mu_pfreliso03all, "Muon_pfRelIso03_all[nMuon]/F");
	tree->Branch("Muon_pfRelIso04_all", value_mu_pfreliso04all, "Muon_pfRelIso04_all[nMuon]/F");
	tree->Branch("Muon_tightId", value_mu_tightid, "Muon_tightId[nMuon]/O");
	tree->Branch("Muon_softId", value_mu_softid, "Muon_softId[nMuon]/O");
	tree->Branch("Muon_dxy", value_mu_dxy, "Muon_dxy[nMuon]/F");
	tree->Branch("Muon_dxyErr", value_mu_dxyErr, "Muon_dxyErr[nMuon]/F");
	tree->Branch("Muon_dz", value_mu_dz, "Muon_dz[nMuon]/F");
	tree->Branch("Muon_dzErr", value_mu_dzErr, "Muon_dzErr[nMuon]/F");

	// Electrons
	tree->Branch("nElectron", &value_el_n, "nElectron/i");
	tree->Branch("Electron_pt", value_el_pt, "Electron_pt[nElectron]/F");
	tree->Branch("Electron_eta", value_el_eta, "Electron_eta[nElectron]/F");
	tree->Branch("Electron_phi", value_el_phi, "Electron_phi[nElectron]/F");
	tree->Branch("Electron_mass", value_el_mass, "Electron_mass[nElectron]/F");
	tree->Branch("Electron_charge", value_el_charge, "Electron_charge[nElectron]/I");
	tree->Branch("Electron_pfRelIso03_all", value_el_pfreliso03all, "Electron_pfRelIso03_all[nElectron]/F");
	tree->Branch("Electron_dxy", value_el_dxy, "Electron_dxy[nElectron]/F");
	tree->Branch("Electron_dxyErr", value_el_dxyErr, "Electron_dxyErr[nElectron]/F");
	tree->Branch("Electron_dz", value_el_dz, "Electron_dz[nElectron]/F");
	tree->Branch("Electron_dzErr", value_el_dzErr, "Electron_dzErr[nElectron]/F");

	// Taus
	tree->Branch("nTau", &value_tau_n, "nTau/i");
	tree->Branch("Tau_pt", value_tau_pt, "Tau_pt[nTau]/F");
	tree->Branch("Tau_eta", value_tau_eta, "Tau_eta[nTau]/F");
	tree->Branch("Tau_phi", value_tau_phi, "Tau_phi[nTau]/F");
	tree->Branch("Tau_mass", value_tau_mass, "Tau_mass[nTau]/F");
	tree->Branch("Tau_charge", value_tau_charge, "Tau_charge[nTau]/I");
	tree->Branch("Tau_decayMode", value_tau_decaymode, "Tau_decayMode[nTau]/I");
	tree->Branch("Tau_chargedIso", value_tau_chargediso, "Tau_chargedIso[nTau]/F");
	tree->Branch("Tau_neutralIso", value_tau_neutraliso, "Tau_neutralIso[nTau]/F");
	tree->Branch("Tau_reliso_charged", value_tau_relisochg, "Tau_reliso_charged[nTau]/F"); 
	tree->Branch("Tau_reliso_gamma", value_tau_relisoneut, "Tau_reliso_gamma[nTau]/F");
		
	tree->Branch("Tau_emFraction", value_tau_emFraction, "Tau_emFraction[nTau]/F");
	tree->Branch("Tau_hcalTotOverPLead", value_tau_hcalTotOverPLead, "Tau_hcalTotOverPLead[nTau]/F");
	tree->Branch("Tau_hcalMaxOverPLead", value_tau_hcalMaxOverPLead, "Tau_hcalMaxOverPLead[nTau]/F");
	tree->Branch("Tau_hcal3x3OverPLead", value_tau_hcal3x3OverPLead, "Tau_hcal3x3OverPLead[nTau]/F");
	tree->Branch("Tau_ecalStripSumEOverPLead", value_tau_ecalStripSumEOverPLead, "Tau_ecalStripSumEOverPLead[nTau]/F");
	tree->Branch("Tau_bremsRecoveryEOverPLead", value_tau_bremsRecoveryEOverPLead, "Tau_bremsRecoveryEOverPLead[nTau]/F");
	tree->Branch("Tau_electronPreIDOutput", value_tau_electronPreIDOutput, "Tau_electronPreIDOutput[nTau]/F");
	tree->Branch("Tau_electronPreIDDecision", value_tau_electronPreIDDecision, "Tau_electronPreIDDecision[nTau]/O");
		
	tree->Branch("Tau_hasMuonReference", value_tau_hasMuonReference, "Tau_hasMuonReference[nTau]/O");
	tree->Branch("Tau_caloComp", value_tau_caloComp, "Tau_caloComp[nTau]/F");
	tree->Branch("Tau_segComp", value_tau_segComp, "Tau_segComp[nTau]/F");
	tree->Branch("Tau_muonDecision", value_tau_muonDecision, "Tau_muonDecision[nTau]/O");

if(isDataFlag_ == false){	
	tree->Branch("Tau_matched_tau", value_tau_matched_tau, "Tau_matched_tau[nTau]/O");
	tree->Branch("Tau_matched_mu", value_tau_matched_mu, "Tau_matched_mu[nTau]/O");
	tree->Branch("Tau_matched_el", value_tau_matched_el, "Tau_matched_el[nTau]/O");	
}





	// Photons
	tree->Branch("nPhoton", &value_ph_n, "nPhoton/i");
	tree->Branch("Photon_pt", value_ph_pt, "Photon_pt[nPhoton]/F");
	tree->Branch("Photon_eta", value_ph_eta, "Photon_eta[nPhoton]/F");
	tree->Branch("Photon_phi", value_ph_phi, "Photon_phi[nPhoton]/F");
	tree->Branch("Photon_mass", value_ph_mass, "Photon_mass[nPhoton]/F");
	tree->Branch("Photon_charge", value_ph_charge, "Photon_charge[nPhoton]/I");
	tree->Branch("Photon_pfRelIso03_all", value_ph_pfreliso03all, "Photon_pfRelIso03_all[nPhoton]/F");

	// MET
	tree->Branch("MET_pt", &value_met_pt, "MET_pt/F");
	tree->Branch("MET_phi", &value_met_phi, "MET_phi/F");
	tree->Branch("MET_sumet", &value_met_sumet, "MET_sumet/F");
	tree->Branch("MET_significance", &value_met_significance, "MET_significance/F");
	tree->Branch("MET_CovXX", &value_met_covxx, "MET_CovXX/F");
	tree->Branch("MET_CovXY", &value_met_covxy, "MET_CovXY/F");
	tree->Branch("MET_CovYY", &value_met_covyy, "MET_CovYY/F");

	// Jets
	tree->Branch("nJet", &value_jet_n, "nJet/i");
	tree->Branch("Jet_pt", value_jet_pt, "Jet_pt[nJet]/F");
	tree->Branch("Jet_eta", value_jet_eta, "Jet_eta[nJet]/F");
	tree->Branch("Jet_phi", value_jet_phi, "Jet_phi[nJet]/F");
	tree->Branch("Jet_mass", value_jet_mass, "Jet_mass[nJet]/F");


	if(isDataFlag_ == false){	
  // Generator particles

		//Electrons  
		  tree->Branch("nGenPart_el", &value_gen_el_n, "nGenPart_el/i");
		  tree->Branch("GenPart_el_pt", "std::vector <float>", &pvalue_gen_el_pt);
		  tree->Branch("GenPart_el_eta", "std::vector <float>", &pvalue_gen_el_eta);
		  tree->Branch("GenPart_el_phi", "std::vector <float>", &pvalue_gen_el_phi);
		  tree->Branch("GenPart_el_mass", "std::vector <float>", &pvalue_gen_el_mass);
		//Muons
	 
		  tree->Branch("nGenPart_mu", &value_gen_mu_n, "nGenPart_mu/i");
		  tree->Branch("GenPart_mu_pt", "std::vector <float>", &pvalue_gen_mu_pt);
		  tree->Branch("GenPart_mu_eta", "std::vector <float>", &pvalue_gen_mu_eta );
		  tree->Branch("GenPart_mu_phi", "std::vector <float>", &pvalue_gen_mu_phi);
		  tree->Branch("GenPart_mu_mass", "std::vector <float>", &pvalue_gen_mu_mass);
		
		//Taus
		  tree->Branch("nGenPart_tau", &value_gen_tau_n, "nGenPart_tau/i");
		  tree->Branch("GenPart_tau_pt", "std::vector <float>", &value_gen_tau_pt);
		  tree->Branch("GenPart_tau_eta", "std::vector <float>", &value_gen_tau_eta);
		  tree->Branch("GenPart_tau_phi", "std::vector <float>", &value_gen_tau_phi);
		  tree->Branch("GenPart_tau_mass", "std::vector <float>", &value_gen_tau_mass);

	}

 
}

AOD2NanoAOD::~AOD2NanoAOD() {}

void AOD2NanoAOD::analyze(const edm::Event &iEvent,
                          const edm::EventSetup &iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;

	pvalue_gen_el_pt->clear();
	pvalue_gen_el_eta->clear();
	pvalue_gen_el_phi->clear();
	pvalue_gen_el_mass->clear();

	pvalue_gen_mu_pt->clear();
	pvalue_gen_mu_eta->clear();
	pvalue_gen_mu_phi->clear();
	pvalue_gen_mu_mass->clear();

	pvalue_gen_tau_pt->clear();
	pvalue_gen_tau_eta->clear();
	pvalue_gen_tau_phi->clear();
	pvalue_gen_tau_mass->clear();

	
  // Event information
  value_run = iEvent.run();
  value_lumi_block = iEvent.luminosityBlock();
  value_event = iEvent.id().event();

  // Vertex
  Handle<VertexCollection> vertices;
  iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
  value_ve_n = vertices->size();
  value_ve_x = vertices->begin()->x();
  value_ve_y = vertices->begin()->y();
  value_ve_z = vertices->begin()->z();
  math::XYZPoint pv(vertices->begin()->position());

  // Muons
  Handle<MuonCollection> muons;
  iEvent.getByLabel(InputTag("muons"), muons);

  value_mu_n = 0;
  const float mu_min_pt = 10;
  for (auto it = muons->begin(); it != muons->end(); it++) {
    if (it->pt() > mu_min_pt) {
      value_mu_pt[value_mu_n] = it->pt();
      value_mu_eta[value_mu_n] = it->eta();
      value_mu_phi[value_mu_n] = it->phi();
      value_mu_charge[value_mu_n] = it->charge();
      value_mu_mass[value_mu_n] = it->mass();
      if (it->isPFMuon() && it->isPFIsolationValid()) {
        auto iso03 = it->pfIsolationR03();
        value_mu_pfreliso03all[value_mu_n] =
            (iso03.sumChargedHadronPt + iso03.sumNeutralHadronEt + iso03.sumPhotonEt)/it->pt(); 
        auto iso04 = it->pfIsolationR04();
        value_mu_pfreliso04all[value_mu_n] =
            (iso04.sumChargedHadronPt + iso04.sumNeutralHadronEt + iso04.sumPhotonEt)/it->pt();
      } else {
        value_mu_pfreliso03all[value_mu_n] = -999;
        value_mu_pfreliso04all[value_mu_n] = -999;
      }
      value_mu_tightid[value_mu_n] = muon::isTightMuon(*it, *vertices->begin());
      value_mu_softid[value_mu_n] = muon::isSoftMuon(*it, *vertices->begin());
      auto trk = it->globalTrack();
      if (trk.isNonnull()) {
        value_mu_dxy[value_mu_n] = trk->dxy(pv);
        value_mu_dz[value_mu_n] = trk->dz(pv);
        value_mu_dxyErr[value_mu_n] = trk->d0Error();
        value_mu_dzErr[value_mu_n] = trk->dzError();
      } else {
        value_mu_dxy[value_mu_n] = -999;
        value_mu_dxyErr[value_mu_n] = -999;
        value_mu_dz[value_mu_n] = -999;
        value_mu_dzErr[value_mu_n] = -999;
      }
      value_mu_n++;
    }
  }

  // Electrons
  Handle<GsfElectronCollection> electrons;
  iEvent.getByLabel(InputTag("gsfElectrons"), electrons);

  value_el_n = 0;
  const float el_min_pt = 10;
  for (auto it = electrons->begin(); it != electrons->end(); it++) {
    if (it->pt() > el_min_pt) {
      value_el_pt[value_el_n] = it->pt();
      value_el_eta[value_el_n] = it->eta();
      value_el_phi[value_el_n] = it->phi();
      value_el_charge[value_el_n] = it->charge();
      value_el_mass[value_el_n] = it->mass();
      if (it->passingPflowPreselection()) {
        auto iso03 = it->pfIsolationVariables();
        value_el_pfreliso03all[value_el_n] =
            (iso03.chargedHadronIso + iso03.neutralHadronIso + iso03.photonIso)/it->pt();
      } else {
        value_el_pfreliso03all[value_el_n] = -999;
      }
      auto trk = it->gsfTrack();
      value_el_dxy[value_el_n] = trk->dxy(pv);
      value_el_dz[value_el_n] = trk->dz(pv);
      value_el_dxyErr[value_el_n] = trk->d0Error();
      value_el_dzErr[value_el_n] = trk->dzError();
      value_el_n++;
    }
  }

  // Taus
  Handle<PFTauCollection> taus;
  iEvent.getByLabel(InputTag("hpsPFTauProducer"), taus);
//passt das so?
  Handle<GenParticleCollection> gens;
  iEvent.getByLabel(InputTag("genParticles"), gens);

  const float tau_min_pt = 15;
  value_tau_n = 0;
for (auto it = taus->begin(); it != taus->end(); it++) {
    if (it->pt() > tau_min_pt) {
	value_tau_pt[value_tau_n] = it->pt();
	value_tau_eta[value_tau_n] = it->eta();
	value_tau_phi[value_tau_n] = it->phi();
	value_tau_charge[value_tau_n] = it->charge();
	value_tau_mass[value_tau_n] = it->mass();
	value_tau_decaymode[value_tau_n] = it->decayMode();
	value_tau_chargediso[value_tau_n] = it->isolationPFChargedHadrCandsPtSum(); //not relative Isolation(?)
	value_tau_neutraliso[value_tau_n] = it->isolationPFGammaCandsEtSum();	//" " only gammas!
	
	//electron rejection
	value_tau_emFraction[value_tau_n] = it->emFraction(); // ECAL/HCAL
	value_tau_hcalTotOverPLead[value_tau_n] = it->hcalTotOverPLead();
	value_tau_hcalMaxOverPLead[value_tau_n] = it->hcalMaxOverPLead();
	value_tau_hcal3x3OverPLead[value_tau_n] = it->hcal3x3OverPLead();
	value_tau_ecalStripSumEOverPLead[value_tau_n] = it->ecalStripSumEOverPLead();
	value_tau_bremsRecoveryEOverPLead[value_tau_n] = it->bremsRecoveryEOverPLead();
	value_tau_electronPreIDOutput[value_tau_n] = it->electronPreIDOutput();
	value_tau_electronPreIDDecision[value_tau_n] = it->electronPreIDDecision();
		
	//muon rejection
	value_tau_hasMuonReference[value_tau_n] = it->hasMuonReference();
	value_tau_caloComp[value_tau_n] = it->caloComp();
	value_tau_segComp[value_tau_n] = it ->segComp();
	value_tau_muonDecision[value_tau_n] = it->muonDecision();

	value_tau_relisochg[value_tau_n] = it->isolationPFChargedHadrCandsPtSum()/it->pt();
	value_tau_relisoneut[value_tau_n] =it->isolationPFGammaCandsEtSum()/it->pt();
	
	//Delta R matching
		if(isDataFlag_ == false){

		Handle<GenParticleCollection> gens;
		iEvent.getByLabel(InputTag("genParticles"), gens);
			float pos_eta = it->eta();
			float pos_phi = it->phi();
			bool matched_tau = false;
			bool matched_el = false;
			bool matched_mu = false;
			float DeltaRsq = 0.09;
			
		

			for(auto genp = gens->begin(); genp != gens->end(); genp++){	

				float DeltaRsq_genpart = (pos_eta - (genp->eta() ))*(pos_eta - (genp->eta() )) +(pos_phi - (genp->phi() ))*(pos_phi - (genp->phi() ));

				if( DeltaRsq_genpart < DeltaRsq ){   //in DeltaR cone of 0.3
						
					

					
					if(std::abs(genp->pdgId()) == 15){
					matched_tau = true;
					matched_el = false;
					matched_mu = false;
					DeltaRsq = DeltaRsq_genpart;
					}

					if(std::abs(genp->pdgId()) == 11){
					matched_el = true;
					matched_tau = false;
					matched_mu = false;
					DeltaRsq = DeltaRsq_genpart;
					}

					if(std::abs(genp->pdgId()) == 13){				
					matched_mu = true;
					matched_tau = false;
					matched_el = false;
					DeltaRsq = DeltaRsq_genpart;			
					} 
				}		
  				
  			}
			value_tau_matched_tau[value_tau_n] = matched_tau;
			value_tau_matched_mu[value_tau_n] = matched_mu;
			value_tau_matched_el[value_tau_n] = matched_el;
		}
	value_tau_n++;
	}	
}
  // Photons
  Handle<PhotonCollection> photons;
  iEvent.getByLabel(InputTag("photons"), photons);

  value_ph_n = 0;
  const float ph_min_pt = 5;
  for (auto it = photons->begin(); it != photons->end(); it++) {
    if (it->pt() > ph_min_pt) {
      value_ph_pt[value_ph_n] = it->pt();
      value_ph_eta[value_ph_n] = it->eta();
      value_ph_phi[value_ph_n] = it->phi();
      value_ph_charge[value_ph_n] = it->charge();
      value_ph_mass[value_ph_n] = it->mass();
      value_ph_pfreliso03all[value_ph_n] = it->ecalRecHitSumEtConeDR03();
      value_ph_n++;
    }
  }

  // MET
  Handle<PFMETCollection> met;
  iEvent.getByLabel(InputTag("pfMet"), met);
  value_met_pt = met->begin()->pt();
  value_met_phi = met->begin()->phi();
  value_met_sumet = met->begin()->sumEt();
  value_met_significance = met->begin()->significance();
  auto cov = met->begin()->getSignificanceMatrix();
  value_met_covxx = cov[0][0];
  value_met_covxy = cov[0][1];
  value_met_covyy = cov[1][1];

  // Jets
  Handle<PFJetCollection> jets;
  iEvent.getByLabel(InputTag("ak5PFJets"), jets);

  const float jet_min_pt = 15;
  value_jet_n = 0;
  for (auto it = jets->begin(); it != jets->end(); it++) {
    if (it->pt() > jet_min_pt) {
      value_jet_pt[value_jet_n] = it->pt();
      value_jet_eta[value_jet_n] = it->eta();
      value_jet_phi[value_jet_n] = it->phi();
      value_jet_mass[value_jet_n] = it->mass();
      value_jet_n++;
    }
  }

  // Generator particles
if (isDataFlag_ == false){	 
const float gen_max_pt = 15000;


  	//Electrons
	value_gen_mu_n = 0;
	value_gen_el_n = 0;
	value_gen_tau_n = 0; 
 for (auto it = gens->begin(); it != gens->end(); it++) {

   	 if (it->status() == 1 && std::abs(it->pdgId())==11 ) {	//find generatorinfo for PDG Electrons
   		if(it->pt() < gen_max_pt){

			pvalue_gen_el_pt->push_back(it->pt());
			std:: cout << "gen_el_pt is: " << it->pt() << std::endl;
			pvalue_gen_el_eta->push_back(it->eta());
			pvalue_gen_el_phi->push_back(it->phi());
			pvalue_gen_el_mass->push_back(it->mass());
 		  	value_gen_el_n++;
		}
	}
	//Muons
	 if (it->status() == 1 && std::abs(it->pdgId()) == 13 ) {   //find generatorinfo for PDG Muons
                if(it->pt() < gen_max_pt ){

			pvalue_gen_mu_pt->push_back(it->pt());
			std:: cout << "gen_mu_pt is: " << it->pt() << std::endl;
			pvalue_gen_mu_eta->push_back(it->eta());
                        pvalue_gen_mu_phi->push_back(it->phi());
                        pvalue_gen_mu_mass->push_back(it->mass());              
                	value_gen_mu_n++;
		}


	}

	
	//Taus
	if (it->status() == 2 && std::abs(it->pdgId()) == 15){
		if(it->pt() < gen_max_pt ){
			
			pvalue_gen_tau_pt->push_back(it->pt());
			std:: cout << "gen_mu_pt is: " << it->pt() << std::endl;
			pvalue_gen_tau_eta->push_back(it->eta());
                        pvalue_gen_tau_phi->push_back(it->phi());
                        pvalue_gen_tau_mass->push_back(it->mass());
                        value_gen_tau_n++;
            	            

		}


   	 }
 }
}


	
  // Fill Tree
tree->Fill();
	 	
}

void AOD2NanoAOD::beginJob() {}

void AOD2NanoAOD::endJob() {}

DEFINE_FWK_MODULE(AOD2NanoAOD);
