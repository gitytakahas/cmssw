#include <memory>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>

// user include files
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

// for the fit
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "TMath.h"
#include "TVector3.h"

class TauProducer : public edm::stream::EDProducer<> {
public:
  explicit TauProducer(const edm::ParameterSet&);
private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const bool isMC_;
  const edm::EDGetTokenT<pat::JetCollection> jetToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pcToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;
};

TauProducer::TauProducer(const edm::ParameterSet& iConfig)
  : isMC_(iConfig.getParameter<bool>("isMC")),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
    pcToken_(consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"))),
    genToken_(isMC_ ? consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))
                    : edm::EDGetTokenT<reco::GenParticleCollection>()),
    ttBuilderToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder")))
{
  produces<pat::CompositeCandidateCollection>("SelectedTaus");    
  produces<pat::CompositeCandidateCollection>("SelectedDiTaus"); // New Collection
}

void TauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;

  auto const& Jets = iEvent.get(jetToken_);
  auto const& vertices = iEvent.get(vertexToken_);
  Handle<reco::BeamSpot> beamspotHandle; iEvent.getByToken(beamspotToken_, beamspotHandle);
  Handle<pat::PackedCandidateCollection> allPackedCandidates; iEvent.getByToken(pcToken_, allPackedCandidates);
  Handle<reco::GenParticleCollection> genParticles; if (isMC_) iEvent.getByToken(genToken_, genParticles);
  auto const& ttBuilder = iSetup.getData(ttBuilderToken_);

  if (vertices.empty()) return;
  const auto& pv = vertices.at(0); 
  
  KalmanVertexFitter kvf(true);
  AdaptiveVertexFitter avf;
  VertexDistance3D vdist;
  KinematicParticleFactoryFromTransientTrack pFactory;

  auto outputTaus = std::make_unique<pat::CompositeCandidateCollection>();
  auto outputDiTaus = std::make_unique<pat::CompositeCandidateCollection>();

  int counter_jet = -1;
  for (const auto& jet : Jets) {
    counter_jet += 1;
    float btagScore = jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb") + jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
    if (btagScore < 0.0494) continue; 

    std::vector<const pat::PackedCandidate*> tracks;
    for (unsigned int d = 0; d < jet.numberOfDaughters(); d++) {
      const auto* part = static_cast<const pat::PackedCandidate*>(jet.daughter(d));
      if (part && part->hasTrackDetails() && part->charge() != 0 && abs(part->pdgId()) == 211) {
        if (part->pt() > 0.5 && abs(part->dz()) < 1.0) tracks.push_back(part);
      }
    }
    if (tracks.size() < 3) continue;

    std::vector<pat::CompositeCandidate> tempTausInJet;
    pat::CompositeCandidate best3Prong;
    float maxVtxProb = -1.0;
    std::vector<size_t> bestTripletIdx = {0, 0, 0};
    bool foundBest = false;
    
    // For storing the intermediate objects needed for di-tau fit (best triplet only)
    std::vector<reco::TransientTrack> bestTtks;
    RefCountedKinematicParticle bestTauFitPart;
    RefCountedKinematicVertex bestTauFitVtx;
    reco::Vertex bestRefittedPV;

    for (size_t i = 0; i < tracks.size(); ++i) {
      for (size_t k = i + 1; k < tracks.size(); ++k) {
        for (size_t l = k + 1; l < tracks.size(); ++l) {
          
          const auto *p1 = tracks[i], *p2 = tracks[k], *p3 = tracks[l];
          if (abs(p1->charge() + p2->charge() + p3->charge()) != 1) continue;
          if (std::abs(p1->dz() - p2->dz()) > 0.2 || std::abs(p1->dz() - p3->dz()) > 0.2) continue;

          auto p4_3prong = p1->p4() + p2->p4() + p3->p4();
          if (p4_3prong.M() > 1.7) continue;

          std::vector<reco::TransientTrack> ttks;
          ttks.push_back(ttBuilder.build(p1->bestTrack()));
          ttks.push_back(ttBuilder.build(p2->bestTrack()));
          ttks.push_back(ttBuilder.build(p3->bestTrack()));

          TransientVertex tv = kvf.vertex(ttks);
          if (!tv.isValid()) continue;
          float vtxProb = TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
          if (vtxProb < 0.1) continue;

          // PV Refit Logic
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

          // Kinematic Fit
          float piMass = 0.13957; float piSigma = 0.0000001;
          std::vector<RefCountedKinematicParticle> kinParticles;
          for(auto const& t : ttks) kinParticles.push_back(pFactory.particle(t, piMass, 0.0, 0.0, piSigma));
          KinematicParticleVertexFitter kpvFitter;
          RefCountedKinematicTree tauTree = kpvFitter.fit(kinParticles);
          if (!tauTree->isValid()) continue;
          tauTree->movePointerToTheTop();
          RefCountedKinematicParticle tauFitPart = tauTree->currentParticle();
          RefCountedKinematicVertex tauFitVtx = tauTree->currentDecayVertex();

          // Analysis Variables
          TwoTrackMinimumDistance md; double maxDoca = -1.0; double minDoca = 99999.9;
          for (size_t m = 0; m < kinParticles.size(); ++m) {
            for (size_t n = m + 1; n < kinParticles.size(); ++n) {
              md.calculate(kinParticles[m]->currentState().freeTrajectoryState(), kinParticles[n]->currentState().freeTrajectoryState());
              if (md.distance() > maxDoca) maxDoca = md.distance();
              if (md.distance() < minDoca) minDoca = md.distance();
            }
          }
          reco::TransientTrack tauTT = ttBuilder.build(tauFitPart->currentState().freeTrajectoryState());
          std::pair<bool, Measurement1D> cur3DIP = IPTools::absoluteImpactParameter3D(tauTT, refittedPV);
          float pvips = (cur3DIP.first && cur3DIP.second.error() != 0) ? cur3DIP.second.value() / cur3DIP.second.error() : -1.0;
          TVector3 plab(tauFitPart->currentState().globalMomentum().x(), tauFitPart->currentState().globalMomentum().y(), tauFitPart->currentState().globalMomentum().z());
          TVector3 tv3diff(tauFitVtx->position().x() - refittedPV.x(), tauFitVtx->position().y() - refittedPV.y(), tauFitVtx->position().z() - refittedPV.z());
          float alpha = (plab.Mag() != 0. && tv3diff.Mag() != 0.) ? plab.Dot(tv3diff) / (plab.Mag() * tv3diff.Mag()) : -1.0;
          reco::Vertex svFit = reco::Vertex(reco::Vertex::Point(tauFitVtx->position()), tauFitVtx->error().matrix(), tauFitVtx->chiSquared(), tauFitVtx->degreesOfFreedom(), 3);
          Measurement1D dist3D = vdist.distance(refittedPV, svFit);

          if (dist3D.significance() < 3.0) continue;

          int nExtra = 0;
          for (const auto& other : tracks) {
            if (other == p1 || other == p2 || other == p3) continue;
            double dx = svFit.x() - other->vertex().x(); double dy = svFit.y() - other->vertex().y(); double dz = svFit.z() - other->vertex().z();
            if (std::sqrt(dx*dx + dy*dy + dz*dz) < 0.1) nExtra++;
          }

          // Gen Matching (Shortened for brevity but logic same as original)
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
          
          myTau.addUserFloat("vtxProb", vtxProb);
          myTau.addUserFloat("l3dSig", dist3D.significance());
          myTau.addUserFloat("flightLen", dist3D.value());
          myTau.addUserFloat("pvips", pvips);
          myTau.addUserFloat("alpha", alpha);
          myTau.addUserFloat("maxDoca", (float)maxDoca);
          myTau.addUserFloat("minDoca", (float)minDoca);
          myTau.addUserFloat("fitMass", (float)tauFitPart->currentState().mass());
          myTau.addUserFloat("fitPt", (float)tauFitPart->currentState().globalMomentum().perp());
          myTau.addUserInt("genMatchId", genMatchId);
          myTau.addUserInt("isFromBsTau", isFromBsTau ? 1 : 0);
          myTau.addUserFloat("svX", svFit.x()); myTau.addUserFloat("svY", svFit.y()); myTau.addUserFloat("svZ", svFit.z());
          myTau.addUserFloat("pvX", refittedPV.x()); myTau.addUserFloat("pvY", refittedPV.y()); myTau.addUserFloat("pvZ", refittedPV.z());
          myTau.addUserFloat("origPvX", pv.x()); myTau.addUserFloat("origPvY", pv.y()); myTau.addUserFloat("origPvZ", pv.z());
          myTau.addUserFloat("nExtra", (float)nExtra);
          myTau.addUserInt("trkIdx1", i); myTau.addUserInt("trkIdx2", k); myTau.addUserInt("trkIdx3", l);
          myTau.addUserFloat("energyFraction", p4_3prong.energy() / jet.energy());
          myTau.addUserInt("counter_jet", counter_jet);
          myTau.addUserFloat("jetPt", jet.pt()); myTau.addUserFloat("jetEta", jet.eta());
          myTau.addUserFloat("jetPhi", jet.phi()); myTau.addUserFloat("jetMass", jet.mass());
          myTau.addUserInt("jetNCharged", jet.chargedHadronMultiplicity());
          myTau.addUserInt("jetNNeutral", jet.neutralHadronMultiplicity());

          bool isTight = (jet.neutralHadronEnergyFraction() < 0.90 && jet.neutralEmEnergyFraction() < 0.90 && (jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1 && jet.chargedHadronEnergyFraction() > 0 && jet.chargedMultiplicity() > 0);
          myTau.addUserInt("jetId", isTight ? 1 : 0);
          myTau.addUserInt("puId", jet.hasUserInt("pileupJetId:fullId") ? jet.userInt("pileupJetId:fullId") : -999);
          myTau.addUserFloat("puScore", jet.hasUserFloat("pileupJetId:fullDiscriminant") ? jet.userFloat("pileupJetId:fullDiscriminant") : -999.0);
          myTau.addUserFloat("deepFlavB", jet.bDiscriminator("pfDeepFlavourJetTags:probb"));
          myTau.addUserFloat("deepFlavBB", jet.bDiscriminator("pfDeepFlavourJetTags:probbb"));
          myTau.addUserFloat("deepFlavLep", jet.bDiscriminator("pfDeepFlavourJetTags:problepb"));
          myTau.addUserFloat("jetBTag", jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb"));
          myTau.addUserInt("hadronFlavour", jet.hadronFlavour());
	  
	  myTau.addUserFloat("deltaChi2", deltaChi2);
	  myTau.addUserFloat("vweight", vweight);


          tempTausInJet.push_back(myTau);

          // Update Best Triplet for Di-Tau formation
          if (vtxProb > maxVtxProb) {
            maxVtxProb = vtxProb;
	    best3Prong = myTau;
            bestTripletIdx = {i, k, l}; bestTtks = ttks;
            bestTauFitPart = tauFitPart; bestTauFitVtx = tauFitVtx;
            bestRefittedPV = refittedPV; foundBest = true;
          }
        }
      }
    }

    // --- FORM DI-TAU from the BEST triplet only ---
// --- FORM DI-TAU from the BEST triplet only ---
    if (foundBest) {
      for (size_t m = 0; m < tracks.size(); ++m) {
        // Skip tracks already used in the best triplet
        if (m == bestTripletIdx[0] || m == bestTripletIdx[1] || m == bestTripletIdx[2]) continue;
        
        const auto* p1p = tracks[m];
        // Charge selection: Opposite Sign (OS) for Bs -> tau(3p) tau(1p)
        if ((best3Prong.charge() + p1p->charge()) != 0) continue; 

        auto p4_ditau = best3Prong.p4() + p1p->p4();
        reco::TransientTrack tt_1p = ttBuilder.build(p1p->bestTrack());
        
        // 1. 4-body Kinematic Fit (3 tracks from triplet + 1 track from 1-prong)
        float piMass = 0.13957; 
        float piSigma = 0.0000001; // Defined as variable to satisfy float& argument
        float piEmpty = 0.0;       

        std::vector<reco::TransientTrack> ttks_4b = bestTtks; 
        ttks_4b.push_back(tt_1p);

        std::vector<RefCountedKinematicParticle> kinParticles_4b;
        for(auto const& t : ttks_4b) {
            kinParticles_4b.push_back(pFactory.particle(t, piMass, piEmpty, piEmpty, piSigma));
        }

        KinematicParticleVertexFitter kpvFitter_4b;
        RefCountedKinematicTree ditauTree = kpvFitter_4b.fit(kinParticles_4b);
        
        float vtxProb_4b = -1.0; 
        float fitMass_4b = -1.0;
        if (ditauTree->isValid()) {
          ditauTree->movePointerToTheTop();
          vtxProb_4b = TMath::Prob(ditauTree->currentDecayVertex()->chiSquared(), (int)ditauTree->currentDecayVertex()->degreesOfFreedom());
          fitMass_4b = ditauTree->currentParticle()->currentState().mass();
        }

        // 2. DOCA & Signed IP relative to 3-prong SV
        TwoTrackMinimumDistance md_1p3p;
        // Dereference the pointer from freeTrajectoryState() using '*'
        md_1p3p.calculate(bestTauFitPart->currentState().freeTrajectoryState(), 
                          *(tt_1p.impactPointState().freeTrajectoryState()));

        // Create a reco::Vertex from the 3-prong Kinematic Fit result for IPTools
        reco::Vertex sv3p_reco = reco::Vertex(reco::Vertex::Point(bestTauFitVtx->position()), 
                                              bestTauFitVtx->error().matrix(), 
                                              bestTauFitVtx->chiSquared(), 
                                              bestTauFitVtx->degreesOfFreedom(), 3);

        // Calculate Signed Impact Parameter using IPTools
        // The direction vector (GlobalVector) from PV to SV defines the sign
        GlobalVector flightDir(sv3p_reco.x() - bestRefittedPV.x(), 
                               sv3p_reco.y() - bestRefittedPV.y(), 
                               sv3p_reco.z() - bestRefittedPV.z());
        
        std::pair<bool, Measurement1D> sip3d = IPTools::signedImpactParameter3D(tt_1p, flightDir, bestRefittedPV);

        // 3. Store in diTau (Inheriting all 3-prong variables via copy)
        pat::CompositeCandidate diTau = best3Prong; 
        diTau.setP4(p4_ditau);
        diTau.addDaughter(*p1p, "pion_1prong");

        // Di-Tau Kinematics
        diTau.addUserFloat("mass_ditau", (float)p4_ditau.M());
        diTau.addUserFloat("dr_1p3p", (float)reco::deltaR(best3Prong.p4(), p1p->p4()));
        diTau.addUserFloat("dz_1p3p", (float)std::abs(tracks[bestTripletIdx[0]]->dz() - p1p->dz()));
        
        // Fit and Vertexing results
        diTau.addUserFloat("vtxProb_4b", vtxProb_4b);
        diTau.addUserFloat("fitMass_4b", fitMass_4b);
        diTau.addUserFloat("lp_doca_sv3p", (float)md_1p3p.distance());
        diTau.addUserFloat("sip3d_sig_1p_sv3p", sip3d.first ? (float)sip3d.second.significance() : -999.);
        
        // Full System Isolation (Summing Pt of tracks within dR < 0.3 of Di-Tau axis)
        float iso_sumPt = 0.0;
        for (const auto& other : tracks) {
          bool isDaughter = (other == tracks[bestTripletIdx[0]] || other == tracks[bestTripletIdx[1]] || 
                             other == tracks[bestTripletIdx[2]] || other == p1p);
          if (!isDaughter && reco::deltaR(p4_ditau, other->p4()) < 0.3) {
            iso_sumPt += other->pt();
          }
        }
        diTau.addUserFloat("iso_ditau", (float)(p4_ditau.pt() / (p4_ditau.pt() + iso_sumPt)));


	int genMatchId_1p = 0;
        int isFromBsTau_1p = 0;

        // Use the same genParticles handle used for the 3-prong matching
        if (genParticles.isValid()) {
          float minDR = 0.03; 
          for (const auto& gen : *genParticles) {
            // Check if genParticle is a charged pion (or lepton if you prefer)
            if (std::abs(gen.pdgId()) != 211) continue; 
            
            float dr = reco::deltaR(gen.p4(), p1p->p4());
            if (dr < minDR) {
              minDR = dr;
              genMatchId_1p = gen.pdgId();

              // Trace back to mother to see if it comes from a Tau, then from a Bs
              const reco::Candidate* mother = gen.mother();
              bool fromTau = false;
              bool fromBs = false;
              
              while (mother) {
                if (std::abs(mother->pdgId()) == 15) fromTau = true;
                if (std::abs(mother->pdgId()) == 531) fromBs = true;
                mother = mother->mother();
              }
              if (fromTau && fromBs) isFromBsTau_1p = 1;
            }
          }
        }

        // Store 1-prong gen info in the diTau object
        diTau.addUserInt("genMatchId_1p", genMatchId_1p);
        diTau.addUserInt("isFromBsTau_1p", isFromBsTau_1p);

	
        outputDiTaus->push_back(diTau);
      }
    }


    // --- Cleaning for SelectedTaus (3-prong only) ---
    std::sort(tempTausInJet.begin(), tempTausInJet.end(), [](const pat::CompositeCandidate& a, const pat::CompositeCandidate& b) {
        return a.userFloat("vtxProb") > b.userFloat("vtxProb");
    });

    std::set<std::pair<float, float>> usedTrackIDs;
    for (const auto& tau : tempTausInJet) {
      bool hasOverlap = false;
      std::vector<std::pair<float, float>> currentTauTrackIDs;
      
      for (unsigned int d = 0; d < tau.numberOfDaughters(); ++d) {
        const auto* daug = tau.daughter(d);
        if (!daug) continue;
        std::pair<float, float> trackID = { (float)daug->pt(), (float)daug->eta() };
        if (usedTrackIDs.count(trackID)) { hasOverlap = true; break; }
        currentTauTrackIDs.push_back(trackID);
      }
      
      if (!hasOverlap) {
        outputTaus->push_back(tau);
        for (const auto& id : currentTauTrackIDs) usedTrackIDs.insert(id);
      }
    }
  } 


  iEvent.put(std::move(outputTaus), "SelectedTaus");
  iEvent.put(std::move(outputDiTaus), "SelectedDiTaus");
}

DEFINE_FWK_MODULE(TauProducer);

