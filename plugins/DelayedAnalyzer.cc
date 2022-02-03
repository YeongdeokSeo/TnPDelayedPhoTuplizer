// -*- C++ -*-
//
// Package:    RazorNtuple/DelayedAnalyzer
// Class:      DelayedAnalyzer
//
/**\class DelayedAnalyzer DelayedAnalyzer.cc RazorNtuple/DelayedAnalyzer/plugins/DelayedAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  YeongDeok Seo
//         Created:  Fri, 28 Jan 2022 06:32:18 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <tuple>
#include <string>

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// CMSSW package includes
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RazorNtuple/DelayedAnalyzer/interface/EGammaMvaEleEstimatorCSA14.h"
#include "RazorNtuple/DelayedAnalyzer/interface/ElectronMVAEstimatorRun2NonTrig.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RazorNtuple/DelayedAnalyzer/interface/EGammaMvaPhotonEstimator.h"
#include "RazorNtuple/DelayedAnalyzer/interface/RazorPDFWeightsHelper.h"

//ECAL tools
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

// Geometry
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"


//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

#define OBJECTARRAYSIZE 1000
#define MAX_NPV 200
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


//class DelayedAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
class DelayedAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DelayedAnalyzer(const edm::ParameterSet&);
      ~DelayedAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      virtual void enablePhotonBranches();
      virtual void resetBranches();


   //private:
   protected:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      static const int MAX_PhotonHLTFilters = 100;
      string photonHLTFilterNames[MAX_PhotonHLTFilters];

      std::vector<edm::EDGetTokenT<pat::PhotonCollection>> v_photonsToken_;
      std::vector<edm::InputTag> v_photonsInputTag;
      std::vector<edm::Handle<pat::PhotonCollection>> photons;
      string photonHLTFilterNamesFile_;

      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebRecHitsToken_;
      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeRecHitsToken_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;

      edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebRecHits;
      edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeRecHits;
      edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
      edm::Handle<edm::TriggerResults> triggerBits;
      edm::Handle<pat::METCollection> mets;

      TTree *RazorEvents;

      //Photons
      int nPhotons;
      int nPhotons_overlap;
      float phoE[OBJECTARRAYSIZE];
      float phoPt[OBJECTARRAYSIZE];
      float phoEta[OBJECTARRAYSIZE];
      float phoPhi[OBJECTARRAYSIZE];
      float phoSigmaIetaIeta[OBJECTARRAYSIZE];
      float phoFull5x5SigmaIetaIeta[OBJECTARRAYSIZE];
      float phoR9[OBJECTARRAYSIZE];
      float pho_sminor[OBJECTARRAYSIZE];
      float pho_smajor[OBJECTARRAYSIZE];
      float pho_HoverE[OBJECTARRAYSIZE];
      float pho_pfIsoChargedHadronIso[OBJECTARRAYSIZE];
      float pho_pfIsoNeutralHadronIso[OBJECTARRAYSIZE];
      float pho_pfIsoPhotonIso[OBJECTARRAYSIZE];
      bool  pho_isConversion[OBJECTARRAYSIZE];
      bool  pho_passEleVeto[OBJECTARRAYSIZE];
      float pho_superClusterEnergy[OBJECTARRAYSIZE];
      float pho_superClusterRawEnergy[OBJECTARRAYSIZE];
      float pho_superClusterEta[OBJECTARRAYSIZE];
      float pho_superClusterPhi[OBJECTARRAYSIZE];
      float pho_superClusterX[OBJECTARRAYSIZE];
      float pho_superClusterY[OBJECTARRAYSIZE];
      float pho_superClusterZ[OBJECTARRAYSIZE];
      bool pho_hasPixelSeed[OBJECTARRAYSIZE];
      bool pho_passHLTFilter[OBJECTARRAYSIZE][MAX_PhotonHLTFilters];

      float metType1Px;
      float metType1Py;
      float metType1Pt;
      float metType1Phi;

};


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DelayedAnalyzer::DelayedAnalyzer(const edm::ParameterSet& iConfig)
 :
  v_photonsInputTag(iConfig.getParameter<std::vector<edm::InputTag>>("photons")),
  photonHLTFilterNamesFile_(iConfig.getParameter<string> ("photonHLTFilterNamesFile")),
  ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
  eeRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("eeRecHits"))),
  triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
  

{
  for(unsigned int i=0;i<v_photonsInputTag.size();i++)
  {
    v_photonsToken_.push_back(consumes<pat::PhotonCollection>(v_photonsInputTag[i]));
  }

  //declare the TFileService for output
  edm::Service<TFileService> fs;

  //set up output tree
  RazorEvents = fs->make<TTree>("RazorEvents", "selected miniAOC information");



  //*****************************************************************************************
  //Read in Photon HLT Filters List from config file
  //*****************************************************************************************
  for (int i = 0; i<MAX_PhotonHLTFilters; ++i) photonHLTFilterNames[i] = "";
  ifstream myPhotonHLTFilterFile (edm::FileInPath(photonHLTFilterNamesFile_.c_str()).fullPath().c_str()) ;
  if (myPhotonHLTFilterFile.is_open()) {
    char tmp[1024];
    string line;
    int index;
    string hltfiltername;

    while(myPhotonHLTFilterFile>>line) {
      
      if ( line.empty() || line.substr(0,1) == "#") {
	myPhotonHLTFilterFile.getline(tmp,1024);
	continue;
      }

      index = atoi(line.c_str());
      myPhotonHLTFilterFile >> hltfiltername;
      
      if (index < MAX_PhotonHLTFilters) {
	photonHLTFilterNames[index] = hltfiltername;
      }    
    }    
    myPhotonHLTFilterFile.close();
  } else {
    cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(photonHLTFilterNamesFile_.c_str()).fullPath().c_str() << "\n";
  }
}


DelayedAnalyzer::~DelayedAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void DelayedAnalyzer::enablePhotonBranches(){
  
  cout << "line #255" << endl;
  RazorEvents->Branch("nPhotons", &nPhotons,"nPhotons/I");
  RazorEvents->Branch("nPhotons_overlap", &nPhotons_overlap,"nPhotons_overlap/I");
  RazorEvents->Branch("phoE", phoE,"phoE[nPhotons]/F");
  RazorEvents->Branch("phoPt", phoPt,"phoPt[nPhotons]/F");
  RazorEvents->Branch("phoEta", phoEta,"phoEta[nPhotons]/F");
  RazorEvents->Branch("phoPhi", phoPhi,"phoPhi[nPhotons]/F");
  RazorEvents->Branch("phoSigmaIetaIeta", phoSigmaIetaIeta, "phoSigmaIetaIeta[nPhotons]/F");
  RazorEvents->Branch("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta, "phoFull5x5SigmaIetaIeta[nPhotons]/F");
  RazorEvents->Branch("phoR9", phoR9, "phoR9[nPhotons]/F");
  RazorEvents->Branch("pho_sminor", pho_sminor, "pho_sminor[nPhotons]/F");
  RazorEvents->Branch("pho_smajor", pho_smajor, "pho_smajor[nPhotons]/F");
  RazorEvents->Branch("pho_HoverE", pho_HoverE, "pho_HoverE[nPhotons]/F");
  RazorEvents->Branch("pho_pfIsoChargedHadronIso", pho_pfIsoChargedHadronIso, "pho_pfIsoChargedHadronIso[nPhotons]/F");
  RazorEvents->Branch("pho_pfIsoNeutralHadronIso", pho_pfIsoNeutralHadronIso, "pho_pfIsoNeutralHadronIso[nPhotons]/F");
  RazorEvents->Branch("pho_pfIsoPhotonIso", pho_pfIsoPhotonIso, "pho_pfIsoPhotonIso[nPhotons]/F");
  RazorEvents->Branch("pho_isConversion", pho_isConversion, "pho_isConversion[nPhotons]/O");
  RazorEvents->Branch("pho_passEleVeto", pho_passEleVeto, "pho_passEleVeto[nPhotons]/O");
  RazorEvents->Branch("pho_superClusterEnergy", pho_superClusterEnergy, "pho_superClusterEnergy[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterRawEnergy", pho_superClusterRawEnergy, "pho_superClusterRawEnergy[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterEta", pho_superClusterEta, "pho_superClusterEta[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterPhi", pho_superClusterPhi, "pho_superClusterPhi[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterX", pho_superClusterX, "pho_superClusterX[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterY", pho_superClusterY, "pho_superClusterY[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterZ", pho_superClusterZ, "pho_superClusterZ[nPhotons]/F");
  RazorEvents->Branch("pho_hasPixelSeed", pho_hasPixelSeed, "pho_hasPixelSeed[nPhotons]/O");
  RazorEvents->Branch("pho_passHLTFilter", &pho_passHLTFilter, Form("pho_passHLTFilter[nPhotons][%d]/O",MAX_PhotonHLTFilters));

  RazorEvents->Branch("metType1Px", &metType1Px, "metType1Px/F");
  RazorEvents->Branch("metType1Py", &metType1Py, "metType1Py/F");
  RazorEvents->Branch("metType1Pt", &metType1Pt, "metType1Pt/F");
  RazorEvents->Branch("metType1Phi", &metType1Phi, "metType1Phi/F");
}

void DelayedAnalyzer::resetBranches(){

    nPhotons = 0;
    nPhotons_overlap = 0;

    for(int i = 0; i < OBJECTARRAYSIZE; i++)
    {      
        phoE[i] = 0.0;
        phoPt[i] = 0.0;
        phoEta[i] = 0.0;
        phoPhi[i] = 0.0;
        phoSigmaIetaIeta[i] = -99.0;
        phoFull5x5SigmaIetaIeta[i] = -99.0;
        phoR9[i] = -99.0;
        pho_sminor[i] = -99.0;
        pho_smajor[i] = -99.0;
        pho_HoverE[i] = -99.0;
        pho_pfIsoChargedHadronIso[i] = -99.0;
        pho_pfIsoNeutralHadronIso[i] = -99.0;
        pho_pfIsoPhotonIso[i] = -99.0;
        pho_isConversion[i] = false;
        pho_passEleVeto[i] = false;    
        pho_superClusterEnergy[i] = -99.0;
        pho_superClusterRawEnergy[i] = -99.0;
        pho_superClusterEta[i]    = -99.0;
        pho_superClusterPhi[i]    = -99.0;
        pho_superClusterX[i]      = -99.0;
        pho_superClusterY[i]      = -99.0;
        pho_superClusterZ[i]      = -99.0;
        pho_hasPixelSeed[i] = false;
        for (int q=0;q<MAX_PhotonHLTFilters;q++) pho_passHLTFilter[i][q] = false;
      
    }

    metType1Pt = -99.0;
    metType1Px = -99.0;
    metType1Py = -99.0;
    metType1Phi = -99.0;

}

// ------------ method called for each event  ------------
void
DelayedAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   resetBranches();
   using namespace edm;

   for(unsigned int i=0;i<v_photonsToken_.size();i++)
   {
	 edm::Handle<pat::PhotonCollection> temp_photons;
    	iEvent.getByToken(v_photonsToken_[i],temp_photons); 
	 if(temp_photons->size() < OBJECTARRAYSIZE) photons.push_back(temp_photons);
   //}
   iEvent.getByToken(ebRecHitsToken_,ebRecHits);
   iEvent.getByToken(eeRecHitsToken_,eeRecHits);
   iEvent.getByToken(triggerObjectsToken_, triggerObjects);
   iEvent.getByToken(triggerBitsToken_, triggerBits);
   iEvent.getByToken(metToken_, mets);

   const pat::MET &Met = mets->front();

   metType1Pt = Met.pt();
   metType1Px = Met.px();
   metType1Py = Met.py();
    
   noZS::EcalClusterLazyTools *lazyToolnoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, ebRecHitsToken_, eeRecHitsToken_);

   std::vector<unsigned int> idx_OOTphotonsToSkip;

   cout << "line #352" << endl;

   for(unsigned int ind_photons = 0; ind_photons<photons.size();ind_photons++)
   {

       cout << photons.size() << " : photons.size" << endl;
       unsigned int idx_OOTphotons = 0;

       for (const pat::Photon &pho : *photons[ind_photons]) 
       {
           cout << "B" << endl;
           if(ind_photons > 0) idx_OOTphotons ++;
           bool toBeSkipped = false;
           cout << "C" << endl;
           //for OOT photons, check if this photon is marked as to be skipped
           if (ind_photons > 0)
           {
               cout << "C2" << endl;
               for(unsigned int i=0; i<idx_OOTphotonsToSkip.size(); i++)
               {
                   cout << "D" << endl;
                   cout << idx_OOTphotonsToSkip.size() << " : idx_OOT size" << endl;
                   cout << i << " : i in for loop" << endl;
                   if(idx_OOTphotons == idx_OOTphotonsToSkip[i])
                   {
                       toBeSkipped = true;
                       continue;
                   }
               }	
           } 
           cout << "D2" << endl;
           //for in-time photons, check if it is overlapped with OOT photons
           //if overlap, and pT_inTime > pT_OOT, then remove OOT photon
           //if overlap, and pT_OOT > pT_inTime, then remove in-time photon
           if (ind_photons == 0 && photons.size()>0){
               cout << "D3" << endl;
               float pt_inTime = pho.pt();
               float eta_inTime = pho.eta();
               float phi_inTime = pho.phi();
               unsigned int idx_OOTphotons_beingchecked = 0;
               for(unsigned int ind_photons_OOT = 1; ind_photons_OOT<photons.size();ind_photons_OOT++)
               {
                   cout << "E" << endl;
                   for (const pat::Photon &pho_OOT : *photons[ind_photons_OOT]) 
                   {
                       cout << "F" << endl;
                       idx_OOTphotons_beingchecked ++;
                       float deltaR_inTime_OOT = deltaR(eta_inTime, phi_inTime, pho_OOT.eta(), pho_OOT.phi());
                       if(deltaR_inTime_OOT > 0.3) continue;
                       else
                       {
                           if(pt_inTime < pho_OOT.pt()) toBeSkipped = true; // remove the in time photon
                           else idx_OOTphotonsToSkip.push_back(idx_OOTphotons_beingchecked);
                       }
                   }
               }
           }
           cout << "line #390" << endl;

            toBeSkipped = false;
           if (toBeSkipped) 
           {
               //MET correction: if remove inTime photon, add its pt to met pt, because it should not be considered in the met calculation but it was
               if (ind_photons == 0) 
               {
                   metType1Px = metType1Px + pho.px();
                   metType1Py = metType1Py + pho.py();
                   metType1Pt = sqrt(metType1Px*metType1Px+metType1Py*metType1Py); 
                   TVector3 vec_met_temp(metType1Px, metType1Py, 0);
                   metType1Phi = vec_met_temp.Phi();	
               }

               if (pho.pt() > 15) nPhotons_overlap ++;
               continue;
           }
           //MET correction: if keep OOT photon, subtract its pt to met pt, because it was originally not considered in the met calculation

           if(ind_photons > 0)
           {
               metType1Px = metType1Px - pho.px();
               cout << metType1Px << " : metType1Px" << endl;
               metType1Py = metType1Py - pho.py();
               metType1Pt = sqrt(metType1Px*metType1Px+metType1Py*metType1Py);
               TVector3 vec_met_temp(metType1Px, metType1Py, 0);
               metType1Phi = vec_met_temp.Phi();	
           }

           cout << pho.pt() << " : pho.pt()" << endl;
           if (pho.pt() < 15) continue;
           cout << "line #419" << endl;

           std::vector<float> vCov = lazyToolnoZS->localCovariances( *(pho.superCluster()->seed()) );

           cout << "line #423" << endl;
           //get photon smajor and sminor
           const auto recHits = (pho.isEB() ? ebRecHits.product() : eeRecHits.product());
           cout << "line #426" << endl;
           if(recHits->size() > 0) 
           {
               cout << "line #429" << endl;
               Cluster2ndMoments ph2ndMoments = noZS::EcalClusterTools::cluster2ndMoments( *(pho.superCluster()), *recHits);
               cout << "line #431" << endl;
               cout << ph2ndMoments.sMaj << " : ph2ndMoments.sMaj" << endl;
               cout << nPhotons << " : nPhotons" << endl;
               cout << pho_smajor[nPhotons] << " : pho_smajor[nPhotons]" << endl;
               pho_smajor[nPhotons] = ph2ndMoments.sMaj;
               cout << "line #433" << endl;
               pho_sminor[nPhotons] = ph2ndMoments.sMin;
               cout << "line #435" << endl;
           }
           cout << "line #437" << endl;
           //-------------------------------------------------
           //default photon 4-mometum already vertex corrected
           //-------------------------------------------------
           //phoE[nPhotons] = pho.getCorrectedEnergy(reco::Photon::P4type::ecal_standard);
           //phoE[nPhotons]   = pho.energy();
           //phoPt[nPhotons]  = pho.pt();
           //phoEta[nPhotons] = pho.eta(); //correct this for the vertex
           //phoPhi[nPhotons] = pho.phi(); //correct this for the vertex

           /*std::cout << "phoE: " << pho.energy() << " phoCorr En:" << pho.getCorrectedEnergy(reco::Photon::P4type::regression2) << " un: " 
             << pho.getCorrectedEnergyError(reco::Photon::P4type::regression2) << " " 
             << pho.getCorrectedEnergyError( pho.getCandidateP4type() ) << std::endl;
             */

           //phoSigmaIetaIeta[nPhotons] = pho.see();
           //phoFull5x5SigmaIetaIeta[nPhotons] = pho.full5x5_sigmaIetaIeta();    

           //phoR9[nPhotons] = pho.r9();
           //Use the noZS version of this according to Emanuele
           //phoR9[nPhotons] = pho.full5x5_r9();

           //pho_HoverE[nPhotons] = pho.hadTowOverEm();
           //pho_isConversion[nPhotons] = pho.hasConversionTracks();

           //pho_passEleVeto[nPhotons] = !hasMatchedPromptElectron(pho.superCluster(),electrons, 
           //conversions, beamSpot->position());
           //use this for 2017 dataset and later - originally for synchronization with Myriam (ETH)
           pho_passEleVeto[nPhotons] = pho.passElectronVeto();

           //**********************************************************
           // Fill default miniAOD isolation quantities
           //**********************************************************
           pho_pfIsoChargedHadronIso[nPhotons] = pho.chargedHadronIso();
           pho_pfIsoNeutralHadronIso[nPhotons] = pho.neutralHadronIso();
           pho_pfIsoPhotonIso[nPhotons] = pho.photonIso();



           //-----------------------
           // super cluster position
           //-----------------------  
           pho_superClusterEnergy[nPhotons] = pho.superCluster()->energy();
           pho_superClusterRawEnergy[nPhotons] = pho.superCluster()->rawEnergy();
           pho_superClusterEta[nPhotons]    = pho.superCluster()->eta();
           pho_superClusterPhi[nPhotons]    = pho.superCluster()->phi();
           pho_superClusterX[nPhotons]      = pho.superCluster()->x();
           pho_superClusterY[nPhotons]      = pho.superCluster()->y();
           pho_superClusterZ[nPhotons]      = pho.superCluster()->z();
           pho_hasPixelSeed[nPhotons]       = pho.hasPixelSeed();

           //*************************************************
           //Trigger Object Matching
           //*************************************************
           for (pat::TriggerObjectStandAlone trigObject : *triggerObjects) {
               if (deltaR(trigObject.eta(), trigObject.phi(),pho.eta(),pho.phi()) > 0.3) continue;
               trigObject.unpackFilterLabels(iEvent, *triggerBits); 

               //check all filters
               for ( int q=0; q<MAX_PhotonHLTFilters;q++) {
                   if (trigObject.hasFilterLabel(photonHLTFilterNames[q].c_str())) pho_passHLTFilter[nPhotons][q] = true;
               }
           }
           nPhotons++;
       }
       cout << nPhotons << " : nPhotons" << endl;
   }

}
}


// ------------ method called once each job just before starting event loop  ------------
void
DelayedAnalyzer::beginJob()
{
  enablePhotonBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void
DelayedAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DelayedAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DelayedAnalyzer);
