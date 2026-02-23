#include <memory>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

// for the fit
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"


struct TauAux {
  std::vector<reco::TransientTrack> ttks;
  RefCountedKinematicParticle fitPart;
  RefCountedKinematicVertex fitVtx;
  reco::Vertex refPV;
};

struct TauSortWrapper {
    pat::CompositeCandidate tau;
    TauAux aux;
};

struct CleanedTauWithAux {
    pat::CompositeCandidate tau;
    TauAux aux;
};

struct OneProng {
  const reco::Candidate* cand;
  int type;      // 1: Hadron, 2: Muon, 3: Electron
  float mass;
  int pdgId;
  int isFromBsTau;
  int genMatchId;
  // --- New variables for Lepton details ---
  float dxy = -999.f;
  float dz = -999.f;
  float iso = -999.f;
  float id = -1.f; // id value
};


class TauProducer : public edm::stream::EDProducer<> {
public:
  explicit TauProducer(const edm::ParameterSet&);
private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  void fillGenMatchInfo(const reco::Candidate& recoCand, int pdgSearch, const reco::GenParticleCollection& genParticles, int& sig, int& gid);
  
  // Helper for Collinear Reconstruction (3+3 exact solution)
  struct KinRes { float x; float alpha; reco::Candidate::LorentzVector p_full; bool valid; };
  KinRes calculate_reconstruction(reco::Candidate::LorentzVector p_vis, TVector3 sv, reco::Vertex pv);

  const bool isMC_;
  const edm::EDGetTokenT<pat::JetCollection> jetToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pcToken_;
  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  const edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;
};

TauProducer::TauProducer(const edm::ParameterSet& iConfig)
  : isMC_(iConfig.getParameter<bool>("isMC")),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
    pcToken_(consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    genToken_(isMC_ ? consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))
                    : edm::EDGetTokenT<reco::GenParticleCollection>()),
    ttBuilderToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder")))
{

  produces<pat::CompositeCandidateCollection>("DiTausHad3P");
  produces<pat::CompositeCandidateCollection>("DiTausMu3P");
  produces<pat::CompositeCandidateCollection>("DiTausEle3P");
  produces<pat::CompositeCandidateCollection>("DiTaus3P3P");
}

TauProducer::KinRes TauProducer::calculate_reconstruction(reco::Candidate::LorentzVector p_vis, TVector3 sv, reco::Vertex pv) {
  TVector3 flight(sv.X() - pv.x(), sv.Y() - pv.y(), sv.Z() - pv.z());
  double dist = flight.Mag();
  if (dist < 1e-5) return {0, 0, {}, false};
  
  TVector3 n = flight.Unit();
  TVector3 p_vis_vec(p_vis.px(), p_vis.py(), p_vis.pz());
  double p_vis_mag = p_vis_vec.Mag();
  double p_vis_n = p_vis_vec.Dot(n);
  
  double m_tau = 1.77686;
  double m_vis2 = std::max(0.0, p_vis.energy()*p_vis.energy() - p_vis_mag*p_vis_mag);
  
  double denom = 2.0 * (p_vis.energy() - p_vis_n);
  if (std::abs(denom) < 1e-8) return {0, 0, {}, false};
  
  double p_nu = (m_tau*m_tau - m_vis2) / denom;
  if (p_nu < 0) return {0, 0, {}, false};
  
  double p_full_mag = p_vis_n + p_nu;
  float x = (p_full_mag > 0) ? (float)(p_vis_n / p_full_mag) : -1.0f;
  float alpha = (p_vis_mag > 0) ? (float)std::acos(std::clamp(p_vis_n / p_vis_mag, -1.0, 1.0)) : 0.0f;
  
  TVector3 p_full_vec = p_vis_vec + (n * p_nu);
  reco::Candidate::LorentzVector p_full(p_full_vec.X(), p_full_vec.Y(), p_full_vec.Z(), p_vis.energy() + p_nu);
  
  return {x, alpha, p_full, true};
}


void TauProducer::fillGenMatchInfo(const reco::Candidate& recoCand, int pdgSearch, const reco::GenParticleCollection& genParticles, int& sig, int& gid) {
  float minDR = 0.02;
  const reco::GenParticle* bestGP = nullptr;
  for (const auto& gp : genParticles) {
    if (gp.status() == 1 && std::abs(gp.pdgId()) == pdgSearch && recoCand.charge() == gp.charge()) {
      double dr = reco::deltaR(recoCand.p4(), gp.p4());
      if (dr < minDR) { minDR = dr; bestGP = &gp; }
    }
  }
  if (bestGP) {
    const reco::Candidate* ancestor = bestGP->mother();
    while (ancestor) {
      if (std::abs(ancestor->pdgId()) == 15) {
        bool hasTauDaughter = false;
        for (unsigned int d = 0; d < ancestor->numberOfDaughters(); ++d) 
          if (std::abs(ancestor->daughter(d)->pdgId()) == 15) hasTauDaughter = true;
        if (!hasTauDaughter) gid = ancestor->pdgId();
      }
      if (std::abs(ancestor->pdgId()) == 531) sig = 1;
      ancestor = ancestor->mother();
    }
  }
}


void TauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //  bool debug = true;
  
  
  using namespace edm;
  using namespace reco;

  auto const& Jets = iEvent.get(jetToken_);
  auto const& vertices = iEvent.get(vertexToken_);
  Handle<reco::BeamSpot> beamspotHandle; iEvent.getByToken(beamspotToken_, beamspotHandle);
  Handle<pat::PackedCandidateCollection> allPackedCandidates; iEvent.getByToken(pcToken_, allPackedCandidates);
  Handle<pat::MuonCollection> muons; iEvent.getByToken(muonToken_, muons);
  Handle<pat::ElectronCollection> electrons; iEvent.getByToken(electronToken_, electrons);
  Handle<reco::GenParticleCollection> genParticles; if (isMC_) iEvent.getByToken(genToken_, genParticles);
  auto const& ttBuilder = iSetup.getData(ttBuilderToken_);

  if (vertices.empty()) return;
  const auto& pv = vertices.at(0); 
  
  KalmanVertexFitter kvf(true);
  AdaptiveVertexFitter avf;
  VertexDistance3D vdist;
  KinematicParticleFactoryFromTransientTrack pFactory;
  
  auto outHad3P = std::make_unique<pat::CompositeCandidateCollection>();
  auto outMu3P  = std::make_unique<pat::CompositeCandidateCollection>();
  auto outEle3P = std::make_unique<pat::CompositeCandidateCollection>();
  auto out3P3P  = std::make_unique<pat::CompositeCandidateCollection>();

  //  if(debug) std::cout << "!!!!!!!!!!!!!!!!! NEW EVENT !!!!!!!!!!!!!!!!!!!!!" << std::endl;
  
  int counter_jet = -1;
  for (const auto& jet : Jets) {
    counter_jet += 1;
    float btagScore = jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb") + jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
    if (btagScore < 0.0494) continue; 

    //    if(debug) std::cout << "------------------------ [INFO] jet_counter = " << counter_jet << std::endl;

    std::vector<OneProng> JetLevelOneProngCands;
    
    std::vector<const pat::PackedCandidate*> tracks;
    for (unsigned int d = 0; d < jet.numberOfDaughters(); d++) {

      const auto* part = static_cast<const pat::PackedCandidate*>(jet.daughter(d));
      if (part && part->hasTrackDetails() && part->charge() != 0) {
        if (part->pt() > 0.5 && abs(part->dz()) < 1.0){

	  int absPdg = std::abs(part->pdgId()); // 11:Ele, 13:Mu, 211:Hadron...
	  //	  if(debug) std::cout << "track found: " << d << " " << part->pt() << " " << part->dz() << " " << part->pdgId() << std::endl;

	  //	  const reco::Track* partTrack = part->bestTrack();

	  if(absPdg==13){
	    bool isMuon = false;
	    for (const auto& mu : *muons) {
	      if (mu.pt() > 2.0 && mu.isLooseMuon()) {
		//		const reco::Track* muTrack = mu.bestTrack(); 

		double dR = reco::deltaR(part->eta(), part->phi(), mu.eta(), mu.phi());
                double relativePtDiff = std::abs(part->pt() - mu.pt()) / mu.pt();

                if (dR < 0.002 && relativePtDiff < 0.02) {

		  int sig = 0, gid = 0;
		  if (isMC_ && genParticles.isValid()) fillGenMatchInfo(mu, 13, *genParticles, sig, gid);

		  float muIso = (mu.pfIsolationR04().sumChargedHadronPt + std::max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5 * mu.pfIsolationR04().sumPUPt)) / mu.pt();
		  float mudxy = (mu.bestTrack() != nullptr) ? (float)mu.bestTrack()->dxy(pv.position()) : -999.f;
		  float mudz  = (mu.bestTrack() != nullptr) ? (float)mu.bestTrack()->dz(pv.position()) : -999.f;
		  //		  float muID = mu.isMediumMuon() ? 1.0f : 0.0f;
		  float muScore = mu.softMvaValue();
                  // Update: Added isolation, ID, and IP variables to the struct
                  JetLevelOneProngCands.push_back({ &mu, 2, 0.10566f, 13, sig, gid, mudxy, mudz, muIso, muScore });

		  isMuon = true;
		  break;
		}
	      }
	    }
	    if (isMuon){
	      //	      if(debug) std::cout << "muon found!" << std::endl;
	      continue;
	    }
	  }

	  if(absPdg==11){
	    bool isEle = false;
	    for (const auto& ele : *electrons) {
	      if (ele.pt() > 2.0 && ele.electronID("mvaEleID-Fall17-noIso-V2-wpLoose")) {
		//		const reco::Track* eleTrack = ele.bestTrack();

		double dR = reco::deltaR(part->eta(), part->phi(), ele.eta(), ele.phi());
		double relativePtDiff = std::abs(part->pt() - ele.pt()) / ele.pt();

		if (dR < 0.002 && relativePtDiff < 0.02) {

		  int sig = 0, gid = 0;
		  if (isMC_ && genParticles.isValid()) fillGenMatchInfo(ele, 11, *genParticles, sig, gid);

		  float eleIso = (ele.pfIsolationVariables().sumChargedHadronPt + std::max(0.0, ele.pfIsolationVariables().sumNeutralHadronEt + ele.pfIsolationVariables().sumPhotonEt - 0.5 * ele.pfIsolationVariables().sumPUPt)) / ele.pt();

		  float eledxy = ele.gsfTrack().isNonnull() ? (float)ele.gsfTrack()->dxy(pv.position()) : -999.f;
		  float eledz  = ele.gsfTrack().isNonnull() ? (float)ele.gsfTrack()->dz(pv.position()) : -999.f;
		  //float eleID = ele.electronID("mvaEleID-Fall17-noIso-V2-wp90") ? 1.0f : 0.0f; // Example: wp90 as "tight"
		  //float eleScore = ele.electronID("mvaEleID-Fall17-noIso-V2-value");
		  float eleScore = ele.userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values");
		  
                  // Update: Added isolation, ID, and IP variables to the struct
                  JetLevelOneProngCands.push_back({ &ele, 3, 0.000511f, 11, sig, gid, eledxy, eledz, eleIso, eleScore });

		  isEle = true;
		  break;
		}
	      }
	    }
	    if (isEle){
	      //	      if(debug) std::cout << "electron found!" << std::endl;
	      continue;
	    }
	  }


	  if (absPdg== 211) tracks.push_back(part);
	  
	}
      }
    }
    if (tracks.size() < 3) continue;

    //    if(debug) std::cout << "Number of tracks = " << tracks.size() << std::endl;
    
    std::vector<pat::CompositeCandidate> tempTausInJet;
    std::vector<TauAux> auxInJet;

    for (size_t i = 0; i < tracks.size(); ++i) {
      for (size_t k = i + 1; k < tracks.size(); ++k) {
        for (size_t l = k + 1; l < tracks.size(); ++l) {
          const auto *p1 = tracks[i], *p2 = tracks[k], *p3 = tracks[l];

	  if (abs(p1->charge() + p2->charge() + p3->charge()) != 1) continue;

	  if (std::abs(p1->dz() - p2->dz()) > 0.2 || std::abs(p1->dz() - p3->dz()) > 0.2) continue;

          auto p4_3prong = p1->p4() + p2->p4() + p3->p4();
          if (p4_3prong.M() > 1.7) continue;

          float m_rho_best = -1.0; float min_diff = 999.0;
          const std::vector<const pat::PackedCandidate*> ps = {p1, p2, p3};
          for(int a=0; a<3; ++a){
            for(int b=a+1; b<3; ++b){
              if(ps[a]->charge() + ps[b]->charge() == 0){
                float m_cand = (float)(ps[a]->p4() + ps[b]->p4()).M();
                if(std::abs(m_cand - 0.770) < min_diff){ min_diff = std::abs(m_cand - 0.770); m_rho_best = m_cand; }
              }
            }
          }

          std::vector<reco::TransientTrack> ttks;
          ttks.push_back(ttBuilder.build(p1->bestTrack()));
          ttks.push_back(ttBuilder.build(p2->bestTrack()));
          ttks.push_back(ttBuilder.build(p3->bestTrack()));

          TransientVertex tv = kvf.vertex(ttks);
          if (!tv.isValid()) continue;
          float vtxProb = TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
          if (vtxProb < 0.1) continue;

          reco::Vertex refittedPV = pv;
          std::vector<reco::TransientTrack> pvTracks;

          if (allPackedCandidates.isValid()) {
            for (const auto& cand : *allPackedCandidates) {
              if (cand.charge() == 0 || !cand.hasTrackDetails() || cand.vertexRef().isNull()) continue;
              if (std::abs(cand.vertexRef()->z() - pv.position().z()) > 0.01) continue;
              bool isSig = false;
              for (const auto* p_sig : {p1, p2, p3}) {
                if (reco::deltaR2(cand.eta(), cand.phi(), p_sig->eta(), p_sig->phi()) < 0.0001) { isSig = true; break; }
              }
              if (!isSig && cand.bestTrack()) pvTracks.push_back(ttBuilder.build(cand.bestTrack()));
            }
          }
          float deltaChi2 = -1; float vweight = -1;
          if (pvTracks.size() >= 2) {
            TransientVertex tvPV = avf.vertex(pvTracks, *beamspotHandle);
            if (tvPV.isValid()){ refittedPV = reco::Vertex(tvPV); deltaChi2 = pv.chi2() - refittedPV.chi2(); vweight = (refittedPV.ndof() + 2.) / (2. * pvTracks.size()); }
          }

          float piMass = 0.13957f; float piSigma = 0.0000001f; float piEmpty = 0.0f;
          std::vector<RefCountedKinematicParticle> kinParticles;
          for(auto const& t : ttks) kinParticles.push_back(pFactory.particle(t, piMass, piEmpty, piEmpty, piSigma));

          KinematicParticleVertexFitter kpvFitter;

	  RefCountedKinematicTree tauTree = kpvFitter.fit(kinParticles);
	  
          if (!tauTree->isValid()) continue;
          tauTree->movePointerToTheTop();
          RefCountedKinematicParticle tauFitPart = tauTree->currentParticle();
          RefCountedKinematicVertex tauFitVtx = tauTree->currentDecayVertex();

          reco::Vertex svFit = reco::Vertex(reco::Vertex::Point(tauFitVtx->position()), tauFitVtx->error().matrix(), tauFitVtx->chiSquared(), tauFitVtx->degreesOfFreedom(), 3);
          Measurement1D dist3D = vdist.distance(refittedPV, svFit);
          if (dist3D.significance() < 3.0) continue;


	  //	  if(debug) std::cout << "(i,k,l, vtxProb, p4_3prong.M(), dist3D) = " << i << " " << k << " " << l << " " << vtxProb << " " << p4_3prong.M() << " " << dist3D.significance() << std::endl;
	  
          TwoTrackMinimumDistance md; double maxDoca = -1.0; double minDoca = 99999.9;

	  //	  if(debug) std::cout << "kinParticles.size: " << kinParticles.size() << std::endl;
	  
	  for (size_t m = 0; m < kinParticles.size(); ++m) {
            for (size_t n = m + 1; n < kinParticles.size(); ++n) {

	      if (kinParticles[m] == kinParticles[n]) continue;

	      //	      std::cout << "check:" << kinParticles[m]->currentState().globalMomentum().perp() << " " << kinParticles[n]->currentState().globalMomentum().perp() << std::endl;
	      
              md.calculate(kinParticles[m]->currentState().freeTrajectoryState(), kinParticles[n]->currentState().freeTrajectoryState());
              if (md.distance() > maxDoca) maxDoca = md.distance();
              if (md.distance() < minDoca) minDoca = md.distance();
            }
          }

	  //	  if(debug) std::cout << "maxDoca, minDoca=" << maxDoca << " " << minDoca << std::endl;
	  
	  reco::TransientTrack tauTT = ttBuilder.build(tauFitPart->currentState().freeTrajectoryState());
          std::pair<bool, Measurement1D> cur3DIP = IPTools::absoluteImpactParameter3D(tauTT, refittedPV);
          float pvips = (cur3DIP.first && cur3DIP.second.error() != 0) ? (float)(cur3DIP.second.value() / cur3DIP.second.error()) : -1.0f;
          TVector3 plab(tauFitPart->currentState().globalMomentum().x(), tauFitPart->currentState().globalMomentum().y(), tauFitPart->currentState().globalMomentum().z());
          TVector3 tv3diff(tauFitVtx->position().x() - refittedPV.x(), tauFitVtx->position().y() - refittedPV.y(), tauFitVtx->position().z() - refittedPV.z());
          float alpha_val = (plab.Mag() != 0. && tv3diff.Mag() != 0.) ? (float)(plab.Dot(tv3diff) / (plab.Mag() * tv3diff.Mag())) : -1.0f;

          int nExtra = 0;
          for (const auto& other : tracks) {
            if (other == p1 || other == p2 || other == p3) continue;
            double dx = svFit.x() - other->vertex().x(); double dy = svFit.y() - other->vertex().y(); double dz = svFit.z() - other->vertex().z();
            if (std::sqrt(dx*dx + dy*dy + dz*dz) < 0.1) nExtra++;
          }

          int genMatchId = 0; bool isFromBsTau = false;
          if(isMC_ && genParticles.isValid()){
            std::vector<const reco::Candidate*> currentTriplet = {p1, p2, p3};
            std::vector<const reco::GenParticle*> matchedGenPions;
            bool possibleMatch = true;
            for (const auto* recoTrk : currentTriplet) {
              float minDR = 0.02; const reco::GenParticle* bestMatch = nullptr;
              for (const auto& gp : *genParticles) {
                if (gp.status() != 1 || std::abs(gp.pdgId()) != 211) continue;
                if (recoTrk->charge() != gp.charge()) continue;
                double dr = reco::deltaR(recoTrk->p4(), gp.p4());
                if (dr < minDR) {
                  if (std::abs(recoTrk->pt() - gp.pt()) / gp.pt() > 0.10) continue;
                  bool alreadyUsed = false;
                  for (const auto* usedGP : matchedGenPions) { if (usedGP == &gp) { alreadyUsed = true; break; } }
                  if (alreadyUsed) continue;
                  minDR = dr; bestMatch = &gp;
                }
              }
              if (bestMatch) matchedGenPions.push_back(bestMatch); else { possibleMatch = false; break; }
            }
            if (possibleMatch && matchedGenPions.size() == 3) {
              const reco::Candidate* mom0 = matchedGenPions[0]->mother();
              if (mom0 && mom0 == matchedGenPions[1]->mother() && mom0 == matchedGenPions[2]->mother()) {
                const reco::Candidate* ancestor = mom0;
                while (ancestor) {
                  int pdg = std::abs(ancestor->pdgId());
                  if (pdg == 15) {
                    bool hasTauDaughter = false;
                    for (unsigned int d = 0; d < ancestor->numberOfDaughters(); ++d) { if (std::abs(ancestor->daughter(d)->pdgId()) == 15) { hasTauDaughter = true; break; } }
                    if (!hasTauDaughter) {
                      int nGenChargedHad = 0;
                      for (unsigned int gd = 0; gd < ancestor->numberOfDaughters(); ++gd) { if (ancestor->daughter(gd)->charge() != 0) nGenChargedHad++; }
                      if (nGenChargedHad == 3) genMatchId = ancestor->pdgId();
                    }
                  }
                  if (pdg == 531) isFromBsTau = true;
                  ancestor = ancestor->mother();
                }
              }
            }
          }

          pat::CompositeCandidate myTau;
          myTau.setP4(p4_3prong); 
          myTau.setCharge(p1->charge() + p2->charge() + p3->charge());
          myTau.addDaughter(*p1, "pion1"); myTau.addDaughter(*p2, "pion2"); myTau.addDaughter(*p3, "pion3");
	  
	  for (int idx = 0; idx < 3; ++idx) {
            std::string pName = "pion" + std::to_string(idx + 1);
            const auto* p_cand = static_cast<const pat::PackedCandidate*>(ps[idx]);
            myTau.addUserFloat(pName + "_pt",     (float)p_cand->pt());
            myTau.addUserFloat(pName + "_eta",    (float)p_cand->eta());
            myTau.addUserFloat(pName + "_phi",    (float)p_cand->phi());
            myTau.addUserInt(pName + "_charge", (int)p_cand->charge());
	  }
	  
	  myTau.addUserFloat("vtxProb", vtxProb);
          myTau.addUserFloat("l3dSig", (float)dist3D.significance());
          myTau.addUserFloat("flightLen", (float)dist3D.value());
          myTau.addUserFloat("pvips", pvips);
          myTau.addUserFloat("alpha", alpha_val);
	  myTau.addUserFloat("maxDoca", (float)maxDoca);
	  myTau.addUserFloat("minDoca", (float)minDoca);
          myTau.addUserFloat("fitMass", (float)tauFitPart->currentState().mass());
          myTau.addUserFloat("fitPt", (float)tauFitPart->currentState().globalMomentum().perp());
          myTau.addUserInt("genMatchId", genMatchId);
          myTau.addUserInt("isFromBsTau", isFromBsTau ? 1 : 0);
          myTau.addUserFloat("svX", (float)svFit.x()); myTau.addUserFloat("svY", (float)svFit.y()); myTau.addUserFloat("svZ", (float)svFit.z());
          myTau.addUserFloat("pvX", (float)refittedPV.x()); myTau.addUserFloat("pvY", (float)refittedPV.y()); myTau.addUserFloat("pvZ", (float)refittedPV.z());
          myTau.addUserFloat("origPvX", (float)pv.x()); myTau.addUserFloat("origPvY", (float)pv.y()); myTau.addUserFloat("origPvZ", (float)pv.z());
          myTau.addUserFloat("nExtra", (float)nExtra);
          myTau.addUserFloat("energyFraction", (float)(p4_3prong.energy() / jet.energy()));
          myTau.addUserInt("counter_jet", counter_jet);
          myTau.addUserFloat("jetPt", (float)jet.pt()); myTau.addUserFloat("jetEta", (float)jet.eta());
          myTau.addUserFloat("jetPhi", (float)jet.phi()); myTau.addUserFloat("jetMass", (float)jet.mass());
          myTau.addUserInt("jetNCharged", jet.chargedHadronMultiplicity());
          myTau.addUserInt("jetNNeutral", jet.neutralHadronMultiplicity());
          myTau.addUserFloat("deltaChi2", (float)deltaChi2);
          myTau.addUserFloat("vweight", (float)vweight);
          myTau.addUserFloat("m_rho", m_rho_best);
          
          myTau.addUserInt("trkIdx1", (int)i); myTau.addUserInt("trkIdx2", (int)k); myTau.addUserInt("trkIdx3", (int)l);
          bool isTight = (jet.neutralHadronEnergyFraction() < 0.90 && jet.neutralEmEnergyFraction() < 0.90 && (jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1 && jet.chargedHadronEnergyFraction() > 0 && jet.chargedMultiplicity() > 0);
          myTau.addUserInt("jetId", isTight ? 1 : 0);
          myTau.addUserInt("puId", jet.hasUserInt("pileupJetId:fullId") ? jet.userInt("pileupJetId:fullId") : -999);
          myTau.addUserFloat("puScore", jet.hasUserFloat("pileupJetId:fullDiscriminant") ? (float)jet.userFloat("pileupJetId:fullDiscriminant") : -999.0f);
          myTau.addUserFloat("deepFlavB", (float)jet.bDiscriminator("pfDeepFlavourJetTags:probb"));
          myTau.addUserFloat("deepFlavBB", (float)jet.bDiscriminator("pfDeepFlavourJetTags:probbb"));
          myTau.addUserFloat("deepFlavLep", (float)jet.bDiscriminator("pfDeepFlavourJetTags:problepb"));
          myTau.addUserFloat("jetBTag", (float)(jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")));
          myTau.addUserInt("hadronFlavour", jet.hadronFlavour());

          tempTausInJet.push_back(myTau);
          auxInJet.push_back({ttks, tauFitPart, tauFitVtx, refittedPV});
        }
      }
    }

    std::vector<TauSortWrapper> wrappedTaus;
    for (size_t i = 0; i < tempTausInJet.size(); ++i) {
      wrappedTaus.push_back({tempTausInJet[i], auxInJet[i]});
    }
    
    // Sort the wrapper vector (Both Tau and Aux move together)
    std::sort(wrappedTaus.begin(), wrappedTaus.end(), [](const TauSortWrapper& a, const TauSortWrapper& b) {
							return a.tau.userFloat("vtxProb") > b.tau.userFloat("vtxProb");
						      });
    
    std::vector<CleanedTauWithAux> finalTausWithAux;
    std::set<std::pair<float, float>> usedTrackIDs;

    
    for (size_t iT = 0; iT < wrappedTaus.size(); ++iT) {

      const auto& tau = wrappedTaus[iT].tau;
      bool hasOverlap = false;
      std::vector<std::pair<float, float>> currentTauTrackIDs;

      //      if(debug) std::cout << "Checking Candidate #" << iT << ": vtxProb=" << tau.userFloat("vtxProb") 
      //			  << " | TrkIdx=(" << tau.userInt("trkIdx1") << "," << tau.userInt("trkIdx2") << "," << tau.userInt("trkIdx3") << ")" << std::endl;
      
      for (unsigned int d = 0; d < tau.numberOfDaughters(); ++d) {
        const auto* daug = tau.daughter(d);
        if (!daug) continue;
        std::pair<float, float> trackID = { (float)daug->pt(), (float)daug->eta() };
        if (usedTrackIDs.count(trackID)) {
	  hasOverlap = true;
	  //	  if(debug) std::cout << "  [OVERLAP FOUND] Track Pt=" << trackID.first << " Eta=" << trackID.second << " already in usedTrackIDs" << std::endl;
	  break;
	}
        currentTauTrackIDs.push_back(trackID);
      }
      if (!hasOverlap) {
	//	if(debug) std::cout << "  [CLEANED] No overlap. Adding to cleanedTaus." << std::endl;
	//        cleanedTaus.push_back(tau);
        //tauAuxMap[&cleanedTaus.back()] = auxInJet[iT]; 
	finalTausWithAux.push_back({wrappedTaus[iT].tau, wrappedTaus[iT].aux});
	
        for (const auto& id : currentTauTrackIDs){
	  //	  if(debug) std::cout << "  [INSERT] Registered Track: Pt=" << id.first << " Eta=" << id.second << std::endl;
	  usedTrackIDs.insert(id);
	}
      }
      //else{
	//	if(debug) std::cout << "  [SKIP] Candidate rejected due to overlap." << std::endl;
      //      }
    }

    //    if(debug) std::cout << "cleandTau.size = " << finalTausWithAux.size() << std::endl;

    // tau->3p loop


    // only allow one electron or one muon within jet ...  take highest in pT ... 


    
    
    std::vector<OneProng> filteredCands;
    const OneProng* bestMu = nullptr;
    const OneProng* bestEle = nullptr;
    
    for (const auto& cand : JetLevelOneProngCands) {
      if (std::abs(cand.pdgId) == 13) { // Muon
	if (!bestMu || cand.cand->pt() > bestMu->cand->pt()) bestMu = &cand;
      } 
      else if (std::abs(cand.pdgId) == 11) { // Electron
        if (!bestEle || cand.cand->pt() > bestEle->cand->pt()) bestEle = &cand;
      } 
      else {
        // Keep all other candidates (pions/hadrons)
        filteredCands.push_back(cand);
      }
    }
    
    // Add the leading muon and electron back to the list
    if (bestMu) filteredCands.push_back(*bestMu);
    if (bestEle) filteredCands.push_back(*bestEle);
    
    // Replace the original vector with the filtered one
    JetLevelOneProngCands = filteredCands;

    // cleaning done 
    
    for (size_t i = 0; i < finalTausWithAux.size(); ++i) {
      const auto* tau3p = &finalTausWithAux[i].tau;
      const auto& aux1 = finalTausWithAux[i].aux;
      
      //    for (size_t iT = 0; iT < cleanedTaus.size(); ++iT) {
      //      const auto* tau3p = &cleanedTaus[iT];
      //      const auto& aux1 = tauAuxMap[tau3p];

      int idx1 = tau3p->userInt("trkIdx1");
      int idx2 = tau3p->userInt("trkIdx2");
      int idx3 = tau3p->userInt("trkIdx3");
      
      //      if(debug) std::cout << "THIS IS TAU for building di-Tau!" <<  idx1 << " " << idx2 << " " << idx3 << std::endl;

      std::vector<OneProng> oneProngCands = JetLevelOneProngCands;

      // for tracks ... (need to be done here to exclude tracks that belong to tau->3p candidate 
      for (size_t idx = 0; idx < tracks.size(); ++idx) {
	//	for (const auto* trk : tracks) {
	const auto* trk = tracks[idx];
	
	//        bool isDaughter = false;
	//	std::cout << "considering: " << idx << std::endl;
	
	if ((int)idx == idx1 || (int)idx == idx2 || (int)idx == idx3) continue;
	//	  for (unsigned int d = 0; d < 3; ++d) if (trk == tau3p->daughter(d)) isDaughter = true;
	//	std::cout << "  THIS TRACK is USED: Pt=" << trk->pt() << " Eta=" << trk->phi() << std::endl;
	
	//        if (!isDaughter) {
	int sig = 0, gid = 0;
	if (isMC_ && genParticles.isValid()) fillGenMatchInfo(*trk, 211, *genParticles, sig, gid);
	oneProngCands.push_back({ trk, 1, 0.13957f, 211, sig, gid, -999.f, -999.f, -1.f, -1.f});

      }
      
      for (const auto& op : oneProngCands) {
	
	//	if(debug) std::cout << "***********************"<< op.type << "1p+3p started" << std::endl;

	if ((tau3p->charge() + op.cand->charge()) != 0) continue;

	reco::TransientTrack tt_1p;
        if (op.type == 1) tt_1p = ttBuilder.build(static_cast<const pat::PackedCandidate*>(op.cand)->bestTrack());
        else if (op.type == 2) tt_1p = ttBuilder.build(static_cast<const pat::Muon*>(op.cand)->bestTrack());
        else if (op.type == 3) tt_1p = ttBuilder.build(static_cast<const pat::Electron*>(op.cand)->bestTrack());
        if (!tt_1p.isValid()) continue;


	
        std::vector<RefCountedKinematicParticle> kin3; float piS = 0.0000001f;
        for(auto const& t : aux1.ttks) kin3.push_back(pFactory.particle(t, 0.13957f, 0.0f, 0.0f, piS));
        kin3.push_back(pFactory.particle(tt_1p, op.mass, 0.0f, 0.0f, piS));
        if (kin3.size() < 2){
	  //	  if(debug) std::cout << op.type << " [WARNING] kin3.size()" << kin3.size() << "is less than 2 !!!" << std::endl;
	  continue;
	}

	//	std::cout << op.type << " kin3.size()" << kin3.size() <<std::endl;

	//// CHECK


	// --- DEBUG: Physical and Index Overlap Check (1p vs 3p daughters) ---
	//	std::cout << "[DEBUG_MATCH] Checking One-Prong candidate against Three-Prong daughters" << std::endl;
	
	// 1. Get One-Prong kinematics from the TransientTrack
	double opPt = tt_1p.track().pt();
	double opEta = tt_1p.track().eta();
	double opPhi = tt_1p.track().phi();
	
	//	std::cout << "  [One-Prong] pt: " << opPt << " | eta: " << opEta << " | phi: " << opPhi << std::endl;
	
	// 2. Get Three-Prong indices saved in UserInts of the tau3p object
	//	int tpIdx1 = tau3p->userInt("trkIdx1");
	//	int tpIdx2 = tau3p->userInt("trkIdx2");
	//	int tpIdx3 = tau3p->userInt("trkIdx3");
	
	//	std::cout << "  [Three-Prong Indices] " << tpIdx1 << ", " << tpIdx2 << ", " << tpIdx3 << std::endl;
	
	bool isPhysicalOverlap = false;
	
	for (size_t i_3p = 0; i_3p < aux1.ttks.size(); ++i_3p) {
	  const reco::Track& t3p = aux1.ttks[i_3p].track();
	  
	  // Calculate deltas to see if they are physically the same track
	  double dR = reco::deltaR(t3p.eta(), t3p.phi(), opEta, opPhi);
	  double dPtRel = std::abs(t3p.pt() - opPt) / opPt;
	  
	  //	  std::cout << "    -> Daughter[" << i_3p << "] pt: " << t3p.pt() 
	  //		    << " | dR: " << dR << " | dPtRel: " << dPtRel << std::endl;
	  
	  // Check by Physics (dR < 0.001 and pT matches within 0.1%)
	  // If this is true, the track is definitely the same one.
	  if (dR < 0.001 && dPtRel < 0.001) {
	    std::cout << "    !!! PHYSICS MATCH FOUND: This 1-prong is identical to 3-prong daughter #" << i_3p << " !!!" << std::endl;
	    isPhysicalOverlap = true;
	  }
	}
	
	if (isPhysicalOverlap) {
	  std::cout << "  [RESULT] Overlap detected! Skipping kpvFitter.fit to avoid self-comparison error." << std::endl;
	  // continue; // Uncomment this if you are inside a loop to skip the fit
	}

	
        KinematicParticleVertexFitter kpvFitter;
	
	RefCountedKinematicTree ditauTree = kpvFitter.fit(kin3);
	
        pat::CompositeCandidate diTau;
        diTau.addDaughter(*tau3p, "leg1");
        diTau.addDaughter(*op.cand, "leg2");
        diTau.setP4(tau3p->p4() + op.cand->p4());
        diTau.setCharge(tau3p->charge() + op.cand->charge());

        // --- Leg 1 & 2 Info ---
        diTau.addUserFloat("leg1_pt", (float)tau3p->pt());
        diTau.addUserFloat("leg1_eta", (float)tau3p->eta());
        diTau.addUserFloat("leg1_phi", (float)tau3p->phi());
        diTau.addUserFloat("leg1_m", (float)tau3p->mass());
        diTau.addUserInt("leg1_charge", (int)tau3p->charge()); // Fix: Added charge
        for (const std::string& name : tau3p->userFloatNames()) diTau.addUserFloat("leg1_" + name, tau3p->userFloat(name));
        for (const std::string& name : tau3p->userIntNames()) diTau.addUserInt("leg1_" + name, tau3p->userInt(name));

        diTau.addUserFloat("leg2_pt", (float)op.cand->pt());
        diTau.addUserFloat("leg2_eta", (float)op.cand->eta());
        diTau.addUserFloat("leg2_phi", (float)op.cand->phi());
        diTau.addUserFloat("leg2_m", (float)op.mass);
        diTau.addUserInt("leg2_charge", (int)op.cand->charge()); // Fix: Added charge
        diTau.addUserInt("leg2_pdgId", op.pdgId);
        diTau.addUserInt("leg2_genMatchId", op.genMatchId);
        diTau.addUserInt("leg2_isFromBsTau", op.isFromBsTau);
	diTau.addUserFloat("leg2_id", op.id);
	diTau.addUserFloat("leg2_iso", op.iso);
	diTau.addUserFloat("leg2_dxy", op.dxy);
	diTau.addUserFloat("leg2_dz", op.dz);
	
        diTau.addUserInt("pairType", op.type);
        diTau.addUserFloat("m_vis", (float)diTau.mass());

        // Fix: Added DR/DZ between legs
        float dr_val = (float)reco::deltaR(tau3p->p4(), op.cand->p4());
        diTau.addUserFloat("dr_taus", dr_val);
	//        diTau.addUserFloat("dr_1p3p", dr_val); 
        float dz_val = (float)(tau3p->vz() - op.cand->vz());
        diTau.addUserFloat("dz_taus", dz_val);

	// lp, sip3d


	// --- 1. Calculate lp_doca_sv3p and sip3d_sig_1p_sv3p ---
        float lp_doca_sv3p = -1.0;
        float sip3d_sig_1p_sv3p = -1.0;
        
        GlobalPoint sv3p_pos(tau3p->userFloat("svX"), tau3p->userFloat("svY"), tau3p->userFloat("svZ"));
        reco::Vertex sv3p_vtx(reco::Vertex::Point(sv3p_pos.x(), sv3p_pos.y(), sv3p_pos.z()), 
                              AlgebraicSymMatrix33());

        // Get the TrajectoryStateOnSurface (TSOS) at the impact point
        TrajectoryStateOnSurface tsos = tt_1p.impactPointState(); // Change from stateAtBeamLine()
        
        if (tsos.isValid()) {
            AnalyticalImpactPointExtrapolator extrapolator(tt_1p.field());
            // Extrapolate the track to the 3-prong SV position
            TrajectoryStateOnSurface tsosAtSV = extrapolator.extrapolate(tsos, sv3p_pos);
            if (tsosAtSV.isValid()) {
                lp_doca_sv3p = (tsosAtSV.globalPosition() - sv3p_pos).mag();
            }
        }

        // Calculate Signed IP Significance of 1-prong wrt 3-prong SV
        std::pair<bool, Measurement1D> ip3d_sv = IPTools::signedImpactParameter3D(
            tt_1p, 
            GlobalVector(op.cand->px(), op.cand->py(), op.cand->pz()), 
            sv3p_vtx
        );
        if (ip3d_sv.first) sip3d_sig_1p_sv3p = ip3d_sv.second.significance();

	
	diTau.addUserFloat("lp_doca_sv3p", lp_doca_sv3p);
	diTau.addUserFloat("sip3d_sig_1p_sv3p", sip3d_sig_1p_sv3p);
	
        diTau.addUserFloat("iso_ratio", (float)((tau3p->pt() + op.cand->pt()) / jet.pt()));

        int is_true = 0;
        if (tau3p->userInt("isFromBsTau") && op.isFromBsTau) {
          int id1 = tau3p->userInt("genMatchId"); int id2 = op.genMatchId;
          if ((id1 == 15 && id2 == -15) || (id1 == -15 && id2 == 15)) is_true = 1;
        }
        diTau.addUserInt("is_true_signal", is_true);

        float vtxProb_bs = -1.0;
        float bsX = -999.0, bsY = -999.0, bsZ = -999.0;

        if (ditauTree->isValid()) {
          ditauTree->movePointerToTheTop();
          vtxProb_bs = (float)TMath::Prob(ditauTree->currentDecayVertex()->chiSquared(), (int)ditauTree->currentDecayVertex()->degreesOfFreedom());
          bsX = (float)ditauTree->currentDecayVertex()->position().x();
          bsY = (float)ditauTree->currentDecayVertex()->position().y();
          bsZ = (float)ditauTree->currentDecayVertex()->position().z();
        }

        // Always add these floats, even if the fit failed
        diTau.addUserFloat("vtxProb_bs", vtxProb_bs);
        diTau.addUserFloat("bsX", bsX);
        diTau.addUserFloat("bsY", bsY);
        diTau.addUserFloat("bsZ", bsZ);

        
        if (op.type == 1) outHad3P->push_back(diTau);
        else if (op.type == 2) outMu3P->push_back(diTau);
        else if (op.type == 3) outEle3P->push_back(diTau);

	//	if(debug) std::cout << "***********************"<< op.type << "1p+3p processed" << std::endl;
      }

      // for 3p+3p

      for (size_t j = i + 1; j < finalTausWithAux.size(); ++j) {
	//        if(debug) std::cout << "*********************** 3p+3p started" << std::endl;
        
        const auto* tau3p_2 = &finalTausWithAux[j].tau;
        const auto& aux2 = finalTausWithAux[j].aux;
	//      for (size_t jT = iT + 1; jT < cleanedTaus.size(); ++jT) {
	//	std::cout << "*********************** 3p+3p started" << std::endl;
	
	//        const auto* tau3p_2 = &cleanedTaus[jT];
	//        const auto& aux2 = tauAuxMap[tau3p_2];
	
	//	int idx1_2 = tau3p_2->userInt("trkIdx1");
	//	int idx2_2 = tau3p_2->userInt("trkIdx2");
	//	int idx3_2 = tau3p_2->userInt("trkIdx3");
	
	//	if(debug) std::cout << "3p: THIS IS TAU for building di-Tau!" <<  idx1_2 << " " << idx2_2 << " " << idx3_2 << std::endl;
	
	if ((tau3p->charge() + tau3p_2->charge()) != 0) continue;

	
        std::vector<reco::TransientTrack> tauTracks;
        tauTracks.push_back(ttBuilder.build(aux1.fitPart->currentState().freeTrajectoryState()));
        tauTracks.push_back(ttBuilder.build(aux2.fitPart->currentState().freeTrajectoryState()));
        TransientVertex tvBs = kvf.vertex(tauTracks);

        pat::CompositeCandidate diTau4;
        diTau4.addDaughter(*tau3p, "leg1");
        diTau4.addDaughter(*tau3p_2, "leg2");

        diTau4.addUserFloat("leg1_pt", (float)tau3p->pt());
        diTau4.addUserFloat("leg1_eta", (float)tau3p->eta());
        diTau4.addUserFloat("leg1_phi", (float)tau3p->phi());
        diTau4.addUserFloat("leg1_m", (float)tau3p->mass());
        diTau4.addUserInt("leg1_charge", (int)tau3p->charge()); // Fix: Added charge
        for (const std::string& name : tau3p->userFloatNames()) diTau4.addUserFloat("leg1_" + name, tau3p->userFloat(name));
        for (const std::string& name : tau3p->userIntNames()) diTau4.addUserInt("leg1_" + name, tau3p->userInt(name));

        diTau4.addUserFloat("leg2_pt", (float)tau3p_2->pt());
        diTau4.addUserFloat("leg2_eta", (float)tau3p_2->eta());
        diTau4.addUserFloat("leg2_phi", (float)tau3p_2->phi());
        diTau4.addUserFloat("leg2_m", (float)tau3p_2->mass());
        diTau4.addUserInt("leg2_charge", (int)tau3p_2->charge()); // Fix: Added charge
        for (const std::string& name : tau3p_2->userFloatNames()) diTau4.addUserFloat("leg2_" + name, tau3p_2->userFloat(name));
        for (const std::string& name : tau3p_2->userIntNames()) diTau4.addUserInt("leg2_" + name, tau3p_2->userInt(name));
        
        diTau4.setP4(tau3p->p4() + tau3p_2->p4());
        diTau4.setCharge(tau3p->charge() + tau3p_2->charge());
        diTau4.addUserInt("pairType", 4);
        diTau4.addUserFloat("m_vis", (float)diTau4.mass());

        // Fix: Added DR/DZ between legs
        float dr_val4 = (float)reco::deltaR(tau3p->p4(), tau3p_2->p4());
        diTau4.addUserFloat("dr_taus", dr_val4);
	//        diTau4.addUserFloat("dr_1p3p", dr_val4); 
        float dz_val4 = (float)(tau3p->vz() - tau3p_2->vz());
        diTau4.addUserFloat("dz_taus", dz_val4);

        diTau4.addUserFloat("iso_ratio", (float)((tau3p->pt() + tau3p_2->pt()) / jet.pt()));

        int is_true = 0;
        if (tau3p->userInt("isFromBsTau") && tau3p_2->userInt("isFromBsTau")) {
          int id1 = tau3p->userInt("genMatchId"); int id2 = tau3p_2->userInt("genMatchId");
          if ((id1 == 15 && id2 == -15) || (id1 == -15 && id2 == 15)) is_true = 1;
        }
        diTau4.addUserInt("is_true_signal", is_true);

	float m_exact = -1.0, x1 = -1.0, x2 = -1.0, alpha1 = -1.0, alpha2 = -1.0;

	TVector3 sv1(tau3p->userFloat("svX"), tau3p->userFloat("svY"), tau3p->userFloat("svZ"));
        TVector3 sv2(tau3p_2->userFloat("svX"), tau3p_2->userFloat("svY"), tau3p_2->userFloat("svZ"));
        
        KinRes res1 = calculate_reconstruction(tau3p->p4(), sv1, aux1.refPV);
        KinRes res2 = calculate_reconstruction(tau3p_2->p4(), sv2, aux2.refPV);

	if (res1.valid && res2.valid && res1.x > 0 && res1.x < 1 && res2.x > 0 && res2.x < 1) {
          m_exact = (float)(res1.p_full + res2.p_full).M();
          x1 = res1.x;
          x2 = res2.x;
          alpha1 = res1.alpha;
          alpha2 = res2.alpha;
        }

	// --- Kinematic Fit (1p+3p) ---
        float vtxProb_bs = -1.0;
        float bsX = -999.0, bsY = -999.0, bsZ = -999.0;

	if (tvBs.isValid()) {
          vtxProb_bs = (float)TMath::Prob(tvBs.totalChiSquared(), (int)tvBs.degreesOfFreedom());
          bsX = (float)tvBs.position().x();
          bsY = (float)tvBs.position().y();
          bsZ = (float)tvBs.position().z();
        }

	diTau4.addUserFloat("vtxProb_bs", vtxProb_bs);
        diTau4.addUserFloat("bsX", bsX);
        diTau4.addUserFloat("bsY", bsY);
        diTau4.addUserFloat("bsZ", bsZ);


	diTau4.addUserFloat("m_exact", m_exact);
        diTau4.addUserFloat("x1", x1);
        diTau4.addUserFloat("x2", x2);
        diTau4.addUserFloat("alpha1", alpha1);
        diTau4.addUserFloat("alpha2", alpha2);

	

        out3P3P->push_back(diTau4);
      }
    }
  }

  iEvent.put(std::move(outHad3P), "DiTausHad3P");
  iEvent.put(std::move(outMu3P),  "DiTausMu3P");
  iEvent.put(std::move(outEle3P), "DiTausEle3P");
  iEvent.put(std::move(out3P3P),  "DiTaus3P3P");
}

DEFINE_FWK_MODULE(TauProducer);


