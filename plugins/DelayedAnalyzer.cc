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
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// CMSSW package includes
#include "DataFormats/PatCandidates/interface/Photon.h"


//ECAL tools
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"


//ROOT includes
#include "TTree.h"
#include "TFile.h"

#define OBJECTARRAYSIZE 100
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

      std::vector<edm::EDGetTokenT<pat::PhotonCollection>> v_photonsToken_;
      std::vector<edm::InputTag> v_photonsInputTag;
      std::vector<edm::Handle<pat::PhotonCollection>> photons;

      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebRecHitsToken_;
      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeRecHitsToken_;

      edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebRecHits;
      edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeRecHits;

      TTree *RazorEvents;

      //Photons
      int nPhotons;
      float pho_sminor[OBJECTARRAYSIZE];
      float pho_smajor[OBJECTARRAYSIZE];

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
  ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
  eeRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("eeRecHits")))
{
  for(unsigned int i=0;i<v_photonsInputTag.size();i++)
  {
    v_photonsToken_.push_back(consumes<pat::PhotonCollection>(v_photonsInputTag[i]));
  }

  //declare the TFileService for output
  edm::Service<TFileService> fs;

  //set up output tree
  RazorEvents = fs->make<TTree>("RazorEvents", "selected miniAOC information");
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
  
  RazorEvents->Branch("nPhotons", &nPhotons,"nPhotons/I");
  RazorEvents->Branch("pho_sminor", pho_sminor, "pho_sminor[nPhotons]/F");
  RazorEvents->Branch("pho_smajor", pho_smajor, "pho_smajor[nPhotons]/F");
}

void DelayedAnalyzer::resetBranches(){

    nPhotons = 0;

    for(int i = 0; i < OBJECTARRAYSIZE; i++)
    {      
        pho_sminor[i] = -99.0;
        pho_smajor[i] = -99.0;
      
    }

}

// ------------ method called for each event  ------------
void
DelayedAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   resetBranches();
   using namespace edm;

   photons.clear();
   //edm::Handle<edm::View<pat::Photon> > photons;
   //iEvent.getByToken(v_photonsToken_, photons);

   for(unsigned int i=0;i<v_photonsToken_.size();i++)
   {
	 edm::Handle<pat::PhotonCollection> temp_photons;
   iEvent.getByToken(v_photonsToken_[i],temp_photons); 
	 if(temp_photons->size() < OBJECTARRAYSIZE) photons.push_back(temp_photons);
   }
   	
   iEvent.getByToken(ebRecHitsToken_,ebRecHits);
   iEvent.getByToken(eeRecHitsToken_,eeRecHits);

    
   noZS::EcalClusterLazyTools *lazyToolnoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, ebRecHitsToken_, eeRecHitsToken_);


   for(unsigned int ind_photons = 0; ind_photons<photons.size();ind_photons++)
   {
       for (const pat::Photon &pho : *photons[ind_photons]) 
       {
         std::vector<float> vCov = lazyToolnoZS->localCovariances( *(pho.superCluster()->seed()) );
         const auto recHits = (pho.isEB() ? ebRecHits.product() : eeRecHits.product());
         if(recHits->size() > 0) 
         {
             Cluster2ndMoments ph2ndMoments = noZS::EcalClusterTools::cluster2ndMoments( *(pho.superCluster()), *recHits);
             pho_smajor[nPhotons] = ph2ndMoments.sMaj;
             pho_sminor[nPhotons] = ph2ndMoments.sMin;
         }
         nPhotons++;
       }

   }


   RazorEvents->Fill();
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
