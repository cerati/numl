////////////////////////////////////////////////////////////////////////
// Class:       HDF5Maker
// Plugin Type: analyzer (art v3_06_03)
// File:        HDF5Maker_module.cc
//
// Generated at Wed May  5 08:23:31 2021 by Jeremy Hewes using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

/*
1. check x correction and timing (MC, pfp, opdet/flash) => fixed MC, need to understand opdet/flash
2. check if wire positions are correct => correct based on evds (there is a 2400 ticks shift wrt evd though)
3. check truth info is filled and correct =? looks good
4. how are the MC starting positions set for showers? => seems good! (we can use end position for photons from pi0s and start position for primary electrons)
 */

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "hep_hpc/hdf5/make_ntuple.hpp"

using std::array;
using std::endl;
using std::setfill;
using std::set;
using std::setw;
using std::string;
using std::vector;

using simb::MCParticle;
using simb::MCTruth;
using sim::TrackIDE;
using recob::Hit;
using recob::OpHit;
using recob::OpFlash;
using recob::SpacePoint;
using recob::Wire;

using mf::LogInfo;

using art::ServiceHandle;
using cheat::BackTrackerService;
using cheat::ParticleInventoryService;
using detinfo::DetectorClocksService;
using detinfo::DetectorPropertiesService;

using hep_hpc::hdf5::Column;
using hep_hpc::hdf5::make_ntuple;
using hep_hpc::hdf5::make_scalar_column;
using hep_hpc::hdf5::make_column;

class HDF5Maker : public art::EDAnalyzer {
public:
  explicit HDF5Maker(fhicl::ParameterSet const& p);
  ~HDF5Maker() noexcept {}; // bare pointers are cleaned up by endSubRun

  HDF5Maker(HDF5Maker const&) = delete;
  HDF5Maker(HDF5Maker&&) = delete;
  HDF5Maker& operator=(HDF5Maker const&) = delete;
  HDF5Maker& operator=(HDF5Maker&&) = delete;

  void beginJob() override;
  void endJob() override;
  void beginSubRun(art::SubRun const& sr) override;
  void endSubRun(art::SubRun const& sr) override;
  void createFile(const art::SubRun* sr);
  void closeFile();
  void analyze(art::Event const& e) override;

private:

  string fTruthLabel;
  string fHitLabel;
  string fHitTruthLabel;
  string fSPLabel;
  string fWireLabel;
  string fPandoraLabel;
  string fOpHitLabel;
  string fOpFlashLabel;

  bool fUseMap;
  string fEventInfo;
  string fOutputName;

  float fXOffset;

  bool fFilesBySubrun;
  hep_hpc::hdf5::File fFile;  ///< Output HDF5 file

  hep_hpc::hdf5::Ntuple<Column<int, 1>     // event id (run, subrun, event)
  >* fEventNtuple; ///< event ntuple

  hep_hpc::hdf5::Ntuple<Column<int, 1>,    // event id (run, subrun, event)
                        Column<int, 1>,    // is cc
                        Column<float, 1>,  // nu energy
                        Column<float, 1>,  // lep energy
                        Column<float, 1>,  // nu dir (x, y, z)
                        Column<float, 1>,  // nu vertex position (x, y, z)
                        Column<float, 1>,  // nu vertex position, corrected (x, y, z)
			Column<int, 1>,    // nu vertex corrected wire pos
			Column<float, 1>   // nu vertex corrected wire time
  >* fEventNtupleNu; ///< event ntuple with neutrino information

  hep_hpc::hdf5::Ntuple<Column<int, 1>,    // event id (run, subrun, event)
                        Column<int, 1>,    // spacepoint id
                        Column<float, 1>,  // 3d position (x, y, z)
                        Column<int, 1>     // 2d hit (u, v, y)
  >* fSpacePointNtuple; ///< spacepoint ntuple

  hep_hpc::hdf5::Ntuple<Column<int, 1>,    // event id (run, subrun, event)
                        Column<int, 1>,    // hit id
                        Column<float, 1>,  // hit integral
                        Column<float, 1>,  // hit rms
                        Column<int, 1>,    // tpc id
                        // Column<int, 1>,    // global plane
                        // Column<float, 1>,  // global wire
                        // Column<float, 1>,  // global time
                        Column<int, 1>,    // raw plane
                        Column<float, 1>,  // raw wire
                        Column<float, 1>   // raw time
  >* fHitNtuple; ///< hit ntuple

  hep_hpc::hdf5::Ntuple<Column<int, 1>,    // event id (run, subrun, event)
                        Column<int, 1>,    // g4 id
                        Column<int, 1>,    // particle type
                        Column<int, 1>,    // parent g4 id
                        Column<float, 1>,  // momentum
                        Column<float, 1>,  // start position (x, y, z)
                        Column<float, 1>,  // end position (x, y, z)
                        Column<float, 1>,  // start position, corrected (x, y, z)
                        Column<float, 1>,  // end position, corrected (x, y, z)
			Column<int, 1>,    // start corrected wire pos
			Column<int, 1>,    // end corrected wire pos
			Column<float, 1>,  // start corrected wire time
			Column<float, 1>,  // end corrected wire time
                        Column<string, 1>, // start process
                        Column<string, 1>  // end process
  >* fParticleNtuple; ///< particle ntuple

  hep_hpc::hdf5::Ntuple<Column<int, 1>,    // event id (run, subrun, event)
                        Column<int, 1>,    // hit id
                        Column<int, 1>,    // g4 id
                        Column<float, 1>,  // deposited energy [ MeV ]
                        Column<float, 1>   // energy fraction
  >* fEnergyDepNtuple; ///< energy deposition ntuple

  hep_hpc::hdf5::Ntuple<Column<int, 1>,    // event id (run, subrun, event)
                        Column<int, 1>,    // tpc id
                        // Column<int, 1>,    // global plane
                        // Column<float, 1>,  // global wire
                        Column<int, 1>,    // raw plane
                        Column<float, 1>,  // raw wire
                        Column<float, 1>   // ADC
  >* fWireNtuple; ///< Wire ntuple

  hep_hpc::hdf5::Ntuple<Column<int, 1>,    // event id (run, subrun, event)
                        Column<int, 1>,    // hit id
			Column<int, 1>,    // hit_channel
			Column<int, 1>,    // wire pos
			Column<float, 1>,  // peaktime
			Column<float, 1>,  // width
			Column<float, 1>,  // area
			Column<float, 1>,  // amplitude
			Column<float, 1>   // pe
  >* fOpHitNtuple; ///< PMT hit ntuple

  hep_hpc::hdf5::Ntuple<Column<int, 1>,    // event id (run, subrun, event)
                        Column<int, 1>,    // flash id
			Column<int, 1>,    // wire pos
			Column<float, 1>,  // time
			Column<float, 1>,  // time width
			Column<float, 1>,  // Y center
			Column<float, 1>,  // Y width
			Column<float, 1>,  // Z center
			Column<float, 1>,  // Z width
			Column<float, 1>,  // pe in each PMT
			Column<float, 1>   // totalpe
  >* fOpFlashNtuple; ///< Flash ntuple

  hep_hpc::hdf5::Ntuple<Column<int, 1>,    // event id (run, subrun, event)
                        Column<int, 1>,    // slice id
                        Column<int, 1>,    // pdg
                        Column<float, 1>,  // neutrino score
                        Column<float, 1>,  // primary vertex position (x, y, z)
			Column<int, 1>,    // primary vertex wire pos
			Column<float, 1>   // primary vertex wire time
  >* fPandoraPrimaryNtuple; ///< Pandora primary ntuple

  hep_hpc::hdf5::Ntuple<Column<int, 1>,    // event id (run, subrun, event)
                        Column<int, 1>,    // pfp id
                        Column<int, 1>,    // pdg
                        Column<float, 1>,  // track-shower score
                        Column<float, 1>,  // pfp vertex position (x, y, z)
			Column<int, 1>,    // pfp vertex wire pos
			Column<float, 1>   // pfp vertex wire time
  >* fPandoraPfpNtuple; ///< Pandora neutrino slice ntuple

  hep_hpc::hdf5::Ntuple<Column<int, 1>,    // event id (run, subrun, event)
                        Column<int, 1>,    // hit id
                        Column<int, 1>,    // slice id
                        Column<int, 1>     // pfp id
  >* fPandoraHitNtuple; ///< Pandora hit slice and cluster ntuple


  // slice hits
  // neutrino vertex
  // pandora track hits
  // pandora shower hits

  using ProxyPfpColl_t = decltype(proxy::getCollection<std::vector<recob::PFParticle> >(
                std::declval<art::Event>(),std::declval<art::InputTag>(),
                proxy::withAssociated<larpandoraobj::PFParticleMetadata>(std::declval<art::InputTag>()),
                proxy::withAssociated<recob::Slice>(std::declval<art::InputTag>()),
                proxy::withAssociated<recob::Cluster>(std::declval<art::InputTag>()),
                proxy::withAssociated<recob::Vertex>(std::declval<art::InputTag>()) ));
  using ProxyPfpElem_t = ProxyPfpColl_t::element_proxy_t;

  // proxy to connect cluster to hit
  using ProxyClusColl_t = decltype(proxy::getCollection<std::vector<recob::Cluster>>(
      std::declval<art::Event>(), std::declval<art::InputTag>(),
      proxy::withAssociated<recob::Hit>(std::declval<art::InputTag>())));
  using ProxyClusElem_t = ProxyClusColl_t::element_proxy_t;

  void True2RecoMappingXYZ(float& t, float& x, float& y, float& z);
  void ApplySCEMappingXYZ(float& x, float& y, float& z);
  void AddDaughters(const std::map<unsigned int, unsigned int>& pfpmap, const ProxyPfpElem_t &pfp_pxy, const ProxyPfpColl_t &pfp_pxy_col, std::vector<ProxyPfpElem_t> &slice_v);
  float GetMetaData(const ProxyPfpElem_t &pfp_pxy, string metaDataName);
  int NearWire(const geo::Geometry& geo, const geo::PlaneID &ip, const float x, const float y, const float z);
};


HDF5Maker::HDF5Maker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fTruthLabel(p.get<string>("TruthLabel")),
    fHitLabel(  p.get<string>("HitLabel")),
    fHitTruthLabel(  p.get<string>("HitTruthLabel","")),
    fSPLabel(   p.get<string>("SPLabel")),
    fWireLabel(   p.get<string>("WireLabel")),
    fPandoraLabel(   p.get<string>("PandoraLabel")),
    fOpHitLabel(  p.get<string>("OpHitLabel")),
    fOpFlashLabel(  p.get<string>("OpFlashLabel")),
    fUseMap(    p.get<bool>("UseMap",false)),
    fEventInfo( p.get<string>("EventInfo")),
    fOutputName(p.get<string>("OutputName")),
    fXOffset(p.get<float>("XOffset")),
    fFilesBySubrun(p.get<bool>("FilesBySubrun",true))
{
  if (fEventInfo != "none" && fEventInfo != "nu")
    throw art::Exception(art::errors::Configuration)
      << "EventInfo must be \"none\" or \"nu\", not " << fEventInfo;
}

void HDF5Maker::analyze(art::Event const& e)
{
  cheat::BackTrackerService* bt = 0;
  if (fUseMap==false) {
    art::ServiceHandle<cheat::BackTrackerService> bt_h;
    bt = bt_h.get();
  }
  art::ServiceHandle<geo::Geometry> geo;
  auto const &detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();

  int run = e.id().run();
  int subrun = e.id().subRun();
  int event = e.id().event();

  array<int, 3> evtID { run, subrun, event };

  // Fill event table
  if (fEventInfo == "none") {
    fEventNtuple->insert( evtID.data() );
    LogInfo("HDF5Maker") << "Filling event table"
                         << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                         << ", event " << evtID[2];
  }
  if (fEventInfo == "nu") {
    // Get MC truth
    art::Handle< vector< MCTruth > > truthHandle;
    e.getByLabel(fTruthLabel, truthHandle);
    if (!truthHandle.isValid() || truthHandle->size() == 0) {
      throw art::Exception(art::errors::LogicError)
        << "Expected to find exactly one MC truth object!";
    }
    simb::MCNeutrino nutruth = truthHandle->at(0).GetNeutrino();

    array<float, 3> nuMomentum {
      (float)nutruth.Nu().Momentum().Vect().Unit().X(),
      (float)nutruth.Nu().Momentum().Vect().Unit().Y(),
      (float)nutruth.Nu().Momentum().Vect().Unit().Z()
    };

    array<float, 3> nuVtx { (float)nutruth.Nu().Vx(), (float)nutruth.Nu().Vy(), (float)nutruth.Nu().Vz() };
    
    float nuT = nutruth.Nu().T();
    array<float, 3> nuVtxCorr { (float)nutruth.Nu().Vx(), (float)nutruth.Nu().Vy(), (float)nutruth.Nu().Vz() };
    True2RecoMappingXYZ(nuT, nuVtxCorr[0], nuVtxCorr[1], nuVtxCorr[2]);

    // std::cout << "nu vtx pos=" << nuVtx[0] << " " <<  nuVtx[1] << " " << nuVtx[2] << std::endl;
    // std::cout << "nu vtx pos corr=" << nuVtxCorr[0] << " " <<  nuVtxCorr[1] << " " << nuVtxCorr[2] << std::endl;

    // double g4Ticks = detClocks->TPCG4Time2Tick(nuT) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
    // std::cout << "TPCG4Time2Tick=" << detClocks->TPCG4Time2Tick(nuT) << " GetXTicksOffset=" << detProperties->GetXTicksOffset(0, 0, 0) << " TriggerOffset=" << detProperties->TriggerOffset() << " TriggerOffsetTPC=" << detClocks->TriggerOffsetTPC() << std::endl;
    // std::cout << "nuT=" << nuT << " x->t =" << detProperties->ConvertXToTicks(nuVtxCorr[0], 0, 0, 0) << " g4Ticks=" << g4Ticks << " Offset=" << detProperties->ConvertXToTicks(fXOffset,0,0,0) << std::endl;

    vector<int> nearwires;
    for (auto p : geo->IteratePlaneIDs()) nearwires.push_back(NearWire(*geo,p,nuVtxCorr[0],nuVtxCorr[1],nuVtxCorr[2]));
    // std::cout << "nearwires=" << nearwires[0] << " " << nearwires[1] << " " << nearwires[2] << std::endl;

    fEventNtupleNu->insert( evtID.data(),
      nutruth.CCNC() == simb::kCC,
      nutruth.Nu().E(),
      nutruth.Lepton().E(),
      nuMomentum.data(),
      nuVtx.data(),
      nuVtxCorr.data(),
      nearwires.data(),
      detProperties->ConvertXToTicks(nuVtxCorr[0], 0, 0, 0)
    );
    LogInfo("HDF5Maker") << "Filling event table"
                         << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                         << ", event " << evtID[2]
                         << "\nis cc? " << (nutruth.CCNC() == simb::kCC)
                         << ", nu energy " << nutruth.Nu().E()
                         << ", lepton energy " << nutruth.Lepton().E()
                         << "\nnu momentum x " << nuMomentum[0] << ", y "
                         << nuMomentum[1] << ", z " << nuMomentum[2];
  } // if nu event info

  // Get spacepoints from the event record
  art::Handle< vector< SpacePoint > > spListHandle;
  vector< art::Ptr< SpacePoint > > splist;
  if (fSPLabel!="") {
    if (e.getByLabel(fSPLabel, spListHandle))
      art::fill_ptr_vector(splist, spListHandle);
  }

  // Get hits from the event record
  art::Handle< vector< Hit > > hitListHandle;
  vector< art::Ptr< Hit > > hitlist;
  if (fHitLabel!="") {
    if (e.getByLabel(fHitLabel, hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
  }

  // Get assocations from spacepoints to hits
  vector< vector< art::Ptr< Hit > > > sp2Hit(splist.size());
  if (splist.size()>0) {
    art::FindManyP< Hit > fmp(spListHandle, e, fSPLabel);
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    } // for spacepoint
  }

  // Fill spacepoint table
  for (size_t i = 0; i < splist.size(); ++i) {

    array<float, 3> pos {
      (float)splist[i]->XYZ()[0],
      (float)splist[i]->XYZ()[1],
      (float)splist[i]->XYZ()[2]
    };

    array<int, 3> hitID { -1, -1, -1 };
    for (size_t j = 0; j < sp2Hit[i].size(); ++j)
      hitID[sp2Hit[i][j]->View()] = sp2Hit[i][j].key();

    fSpacePointNtuple->insert(evtID.data(),
      splist[i]->ID(), pos.data(), hitID.data()
    );

    // std::cout << std::endl << "Filling spacepoint table"
    LogInfo("HDF5Maker") << "Filling spacepoint table"
                         << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                         << ", event " << evtID[2]
                         << "\nspacepoint id " << splist[i]->ID()
                         << "\nposition x " << pos[0] << ", y " << pos[1]
                         << ", z " << pos[2]
                         << "\nhit ids " << hitID[0] << ", " << hitID[1]
                         << ", " << hitID[2];

  } // for spacepoint

  std::set<int> g4id;
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> hittruth;
  if (fUseMap) {
    hittruth = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, e, fHitTruthLabel));
  }

  // Loop over hits
  for (art::Ptr< Hit > hit : hitlist) {

    // Fill hit table
    geo::WireID wireid = hit->WireID();
    size_t plane = wireid.Plane;
    size_t wire = wireid.Wire;
    double time = hit->PeakTime();
    fHitNtuple->insert(evtID.data(),
      hit.key(), hit->Integral(), hit->RMS(), wireid.TPC,
      // plane, wire, time,
      wireid.Plane, wireid.Wire, hit->PeakTime()
    );

    LogInfo("HDF5Maker") << "Filling hit table"
                         << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                         << ", event " << evtID[2]
                         << "\nhit id " << hit.key() << ", integral "
                         << hit->Integral() << ", RMS " << hit->RMS()
                         << ", TPC " << wireid.TPC
                         << "\nglobal plane " << plane << ", global wire "
                         << wire << ", global time " << time
                         << "\nlocal plane " << wireid.Plane
                         << ", local wire " << wireid.Wire
                         << ", local time " << hit->PeakTime();
                         

    // Fill energy deposit table
    if (fUseMap) {
      std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
      std::vector<anab::BackTrackerHitMatchingData const *> match_vec = hittruth->data(hit.key());
      //loop over particles
      for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
	g4id.insert(particle_vec[i_p]->TrackId());
	fEnergyDepNtuple->insert(evtID.data(),
		hit.key(), particle_vec[i_p]->TrackId(), match_vec[i_p]->energy, match_vec[i_p]->ideFraction
	);
	LogInfo("HDF5Maker") << "Filling energy deposit table"
	// std::cout << "\nFilling energy deposit table"
			     << "\nrun " << evtID[0] << ", subrun " << evtID[1]
			     << ", event " << evtID[2]
			     << "\nhit id " << hit.key() << ", g4 id "
			     << particle_vec[i_p]->TrackId() << ", energy "
			     << match_vec[i_p]->energy << ", energy fraction "
			     << match_vec[i_p]->ideFraction;
      }
    } else {
      const auto& ides = bt->HitToTrackIDEs(hit);
      for (const TrackIDE& ide : ides) {
	g4id.insert(ide.trackID);
	fEnergyDepNtuple->insert(evtID.data(),
		hit.key(), ide.trackID, ide.energy, ide.energyFrac
	);
	LogInfo("HDF5Maker") << "Filling energy deposit table"
			     << "\nrun " << evtID[0] << ", subrun " << evtID[1]
			     << ", event " << evtID[2]
			     << "\nhit id " << hit.key() << ", g4 id "
			     << ide.trackID << ", energy "
			     << ide.energy << " MeV" << ", energy fraction "
			     << ide.energyFrac;
      } // for energy deposit
    } // if using microboone map method or not
  } // for hit

  art::Handle< vector< recob::Slice > > sliceListHandle;
  vector< art::Ptr< recob::Slice > > slicelist;
  std::unique_ptr<art::FindManyP<recob::Hit> > assocSliceHit;
  if (fPandoraLabel!="") {
    if (e.getByLabel(fPandoraLabel, sliceListHandle))
      art::fill_ptr_vector(slicelist, sliceListHandle);
    assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit> >(new art::FindManyP<recob::Hit>(sliceListHandle, e, fPandoraLabel));
  }
  vector<int> hit_slice_idx(hitlist.size(),-1);
  for (art::Ptr< recob::Slice > slice : slicelist) {
    auto slicehit = assocSliceHit->at(slice.key());
    // std::cout << "Slice ID=" << is << " N slice hits=" << slicehit.size() << std::endl;
    for (unsigned int isl = 0; isl < slicehit.size(); isl++)
      {
	const auto &slhit = slicehit.at(isl);
	// this does not work if the hit collection is different from "gaushit", i.e. the one used by pandora
	hit_slice_idx[slhit.key()] = slice.key();
	// /*if (false)*/ std::cout << slhit.key() << " plane=" << slhit->View() << " wire=" << slhit->WireID().Wire << " t=" << slhit->PeakTime() << std::endl;
      }
  }
  vector<int> hit_pfp_idx(hitlist.size(),-1);
 
  if (fPandoraLabel!="") {
    // grab PFParticles in event
    art::InputTag fPandoraTag(fPandoraLabel);
    ProxyPfpColl_t const &pfp_pxy_col = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPandoraTag,
											     proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPandoraTag),
											     proxy::withAssociated<recob::Slice>(fPandoraTag),
											     proxy::withAssociated<recob::Cluster>(fPandoraTag),
											     proxy::withAssociated<recob::Vertex>(fPandoraTag)
											     );
    ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster> >(e, fPandoraTag, proxy::withAssociated<recob::Hit>(fPandoraTag));

    // Build PFP Map
    std::map<unsigned int, unsigned int> pfpmap;
    unsigned int p = 0;
    for (const auto &pfp_pxy : pfp_pxy_col) { pfpmap[pfp_pxy->Self()] = p; p++; }

    // loop through PFParticles
    for (const ProxyPfpElem_t &pfp_pxy : pfp_pxy_col)
      {

	// // get metadata for this PFP
	// const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();

	if (pfp_pxy->IsPrimary() == false) continue;

	auto slices = pfp_pxy.get<recob::Slice>();
	if (slices.size() != 1)
	  {
	    std::cout << "WRONG!!! n slices = " << slices.size() << std::endl;
	  }

	// std::cout << "slice id = " << slices[0].key() << " pdg=" << pfp_pxy->PdgCode() << std::endl;
	// std::cout << "isNeutrino = " << GetMetaData(pfp_pxy, "IsNeutrino") << std::endl;
	// std::cout << "slice topo score = " << GetMetaData(pfp_pxy, "NuScore") << std::endl;

	float nvx[3] = {-999.,-999.,-999.};
	vector<int> nearwires;
	auto nuvtx = pfp_pxy.get<recob::Vertex>();
	if (nuvtx.size() != 1)
	  {
	    if (false) std::cout << "Found primary PFP w/ != 1 associated vertices..." << std::endl;
	  }
	else {
	  // std::cout << "Vertex pos = " << nuvtx.at(0)->position() << " key=" << nuvtx.at(0).key() << std::endl;
	  nvx[0] = nuvtx.at(0)->position().X();
	  nvx[1] = nuvtx.at(0)->position().Y();
	  nvx[2] = nuvtx.at(0)->position().Z();
	  for (auto p : geo->IteratePlaneIDs()) nearwires.push_back(NearWire(*geo,p,nvx[0],nvx[1],nvx[2]));
	  // std::cout << "tick=" << detProperties->ConvertXToTicks(nvx[0], 0, 0, 0) << " u=" << nearwires[0] << " v=" << nearwires[1] << " y=" << nearwires[2] << std::endl;
	}

	fPandoraPrimaryNtuple->insert( evtID.data(),
				       slices[0].key(),
				       pfp_pxy->PdgCode(),
				       GetMetaData(pfp_pxy, "NuScore"),
				       nvx,
				       nearwires.data(),
				       detProperties->ConvertXToTicks(nvx[0], 0, 0, 0)
        );
	// LogInfo("HDF5Maker") << "Filling event table"
	//                      << "\nrun " << evtID[0] << ", subrun " << evtID[1]
	//                      << ", event " << evtID[2]
	//                      << "\nis cc? " << (nutruth.CCNC() == simb::kCC)
	//                      << ", nu energy " << nutruth.Nu().E()
	//                      << ", lepton energy " << nutruth.Lepton().E()
	//                      << "\nnu momentum x " << nuMomentum[0] << ", y "
	//                      << nuMomentum[1] << ", z " << nuMomentum[2];

	//  find neutrino candidate
	auto PDG = fabs(pfp_pxy->PdgCode());
	if ( PDG != 12 && PDG != 14 ) continue;

	auto slicehit = assocSliceHit->at(slices[0].key());
	// std::cout << "Neutrino slice ID=" << slices[0].key() << " N slice hits=" << slicehit.size() << std::endl;
 
	// collect PFParticle hierarchy originating from this neutrino candidate
	std::vector<ProxyPfpElem_t> slice_pfp_v;
	AddDaughters(pfpmap, pfp_pxy, pfp_pxy_col, slice_pfp_v);

	for (auto pfp : slice_pfp_v)
	  {
	    if (pfp->PdgCode()==12 || pfp->PdgCode()==14) continue;
	    // std::cout << "pfp pdg=" << pfp->PdgCode() << " track-shower score = " << GetMetaData(pfp, "TrackScore") << std::endl;
	    auto clus_pxy_v = pfp.get<recob::Cluster>();
	    for (auto ass_clus : clus_pxy_v)
	      {
		const auto &clus = clus_proxy[ass_clus.key()];
		auto clushit_v = clus.get<recob::Hit>();
		// std::cout << "cluster with size = " << clushit_v.size() << std::endl;
		for (unsigned int icl = 0; icl < clushit_v.size(); icl++)
		  {
		    const auto &clhit = clushit_v.at(icl);
		    hit_pfp_idx[clhit.key()] = pfp.index();
		    // if (false) std::cout << clhit.key() << std::endl;
		  }
	      }
	    auto pvtx = pfp.get<recob::Vertex>();
	    vector<int> nearwires;
	    float pvx[3] = {-999.,-999.,-999.};
	    if (pvtx.size() != 1)
	      {
		if (false) std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
	      } else {
	      // std::cout << "Vertex pos = " << pvtx.at(0)->position() << " key=" << pvtx.at(0).key() << std::endl;
	      pvx[0] = (float)pvtx.at(0)->position().X();
	      pvx[1] = (float)pvtx.at(0)->position().Y();
	      pvx[2] = (float)pvtx.at(0)->position().Z();
	      for (auto p : geo->IteratePlaneIDs()) nearwires.push_back(NearWire(*geo,p,pvx[0],pvx[1],pvx[2]));
	    }

	    fPandoraPfpNtuple->insert( evtID.data(),
				       pfp.index(),
				       pfp->PdgCode(),
				       GetMetaData(pfp, "TrackScore"),
				       pvx,
				       nearwires.data(),
				       detProperties->ConvertXToTicks(pvx[0], 0, 0, 0)
	    );
	    // LogInfo("HDF5Maker") << "Filling event table"
	    //                      << "\nrun " << evtID[0] << ", subrun " << evtID[1]
	    //                      << ", event " << evtID[2]
	    //                      << "\nis cc? " << (nutruth.CCNC() == simb::kCC)
	    //                      << ", nu energy " << nutruth.Nu().E()
	    //                      << ", lepton energy " << nutruth.Lepton().E()
	    //                      << "\nnu momentum x " << nuMomentum[0] << ", y "
	    //                      << nuMomentum[1] << ", z " << nuMomentum[2];
	  } // for all PFParticles in the slice
      } // for all PFParticles

    for (art::Ptr< Hit > hit : hitlist) {
      fPandoraHitNtuple->insert( evtID.data(),
				 hit.key(),
				 hit_slice_idx[hit.key()],
				 hit_pfp_idx[hit.key()]
      );
    }

  }

  // Get wires from the event record
  art::Handle< vector< Wire > > wireListHandle;
  vector< art::Ptr< Wire > > wirelist;
  if (fWireLabel != "") {
    if (e.getByLabel(fWireLabel, wireListHandle))
      art::fill_ptr_vector(wirelist, wireListHandle);
  }

  // Loop over wires
  for (art::Ptr< Wire > wire : wirelist) {
    auto wireid = geo->ChannelToWire(wire->Channel())[0];//fixme

    LogInfo("HDF5Maker") << "Filling wire table"
			 << "\nrun " << evtID[0] << ", subrun " << evtID[1]
			 << ", event " << evtID[2]
			 << ", TPC " << wireid.TPC
			 << "\nglobal plane " << wire->View()
			 << ", global wire " << wire->Channel()
			 << "\nlocal plane " << wireid.Plane
			 << ", local wire " << wireid.Wire;
    //for (auto w : wire->SignalROI()) if (w>0) std::cout << w << std::endl;
    //
    fWireNtuple->insert(evtID.data(),
			wireid.TPC,
			// wire->View(),
			// wire->Channel(),
			wireid.Plane,
			wireid.Wire,
			wire->Signal().data()
     );
  }

  ServiceHandle<ParticleInventoryService> pi;
  set<int> allIDs = g4id; // Copy original so we can safely modify it

  // initialize the inventory
  pi->Rebuild(e);
  pi->provider()->PrepParticleList(e);

  // Add invisible particles to hierarchy
  for (int id : g4id) { 
    const MCParticle* p = pi->TrackIdToParticle_P(abs(id));
    while (p != 0 && p->Mother() != 0) {
      allIDs.insert(abs(p->Mother()));
      p = pi->TrackIdToParticle_P(abs(p->Mother()));
    }
  }

  // Loop over true particles and fill table
  for (int id : allIDs) {
    const MCParticle* p = pi->TrackIdToParticle_P(abs(id));
    if (p==NULL) continue;
    array<float, 3> particleStart { (float)p->Vx(), (float)p->Vy(), (float)p->Vz() };
    array<float, 3> particleEnd { (float)p->EndX(), (float)p->EndY(), (float)p->EndZ() };
    //
    float startT = p->T();
    float endT = p->EndT();
    array<float, 3> particleStartCorr { (float)p->Vx(), (float)p->Vy(), (float)p->Vz() };
    array<float, 3> particleEndCorr { (float)p->EndX(), (float)p->EndY(), (float)p->EndZ() };
    True2RecoMappingXYZ(startT, particleStartCorr[0], particleStartCorr[1], particleStartCorr[2]);
    True2RecoMappingXYZ(endT, particleEndCorr[0], particleEndCorr[1], particleEndCorr[2]);

    vector<int> nearwires_start;
    for (auto ip : geo->IteratePlaneIDs()) nearwires_start.push_back(NearWire(*geo,ip,particleStartCorr[0],particleStartCorr[1],particleStartCorr[2]));
    vector<int> nearwires_end;
    for (auto ip : geo->IteratePlaneIDs()) nearwires_end.push_back(NearWire(*geo,ip,particleEndCorr[0],particleEndCorr[1],particleEndCorr[2]));

    //
    fParticleNtuple->insert(evtID.data(),
      abs(id), p->PdgCode(), p->Mother(), (float)p->P(),
      particleStart.data(), particleEnd.data(),
      particleStartCorr.data(), particleEndCorr.data(),
      nearwires_start.data(),nearwires_end.data(),
      detProperties->ConvertXToTicks(particleStartCorr[0], 0, 0, 0), detProperties->ConvertXToTicks(particleEndCorr[0], 0, 0, 0),
      p->Process(), p->EndProcess()
    );
    LogInfo("HDF5Maker") << "Filling particle table"
    // std::cout << "\nFilling particle table"
                         << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                         << ", event " << evtID[2]
                         << "\ng4 id " << abs(id) << ", pdg code "
                         << p->PdgCode() << ", parent " << p->Mother()
                         << ", momentum " << p->P()
                         << "\nparticle start x " << particleStart[0]
                         << ", y " << particleStart[1]
                         << ", z " << particleStart[2]
                         << "\nparticle end x " << particleEnd[0] << ", y "
                         << particleEnd[1] << ", z " << particleEnd[2]
                         << "\nparticle start corr x " << particleStartCorr[0]
                         << ", y " << particleStartCorr[1]
                         << ", z " << particleStartCorr[2]
	                 << ", tick " << detProperties->ConvertXToTicks(particleStartCorr[0], 0, 0, 0)
	                 << ", u " << nearwires_start[0] << ", v " << nearwires_start[1] << ", y " << nearwires_start[2]
                         << "\nparticle end corr x " << particleEndCorr[0] << ", y "
                         << particleEndCorr[1] << ", z " << particleEndCorr[2]
	                 << ", tick " << detProperties->ConvertXToTicks(particleEndCorr[0], 0, 0, 0)
	                 << ", u " << nearwires_end[0] << ", v " << nearwires_end[1] << ", y " << nearwires_end[2]
                         << "\nstart process " << p->Process()
                         << ", end process " << p->EndProcess();
  }

  // Get OpHits from the event record
  art::Handle< vector< OpHit > > opHitListHandle;
  vector< art::Ptr< OpHit > > ophitlist;
  if (fOpHitLabel != "") {
    if (e.getByLabel(fOpHitLabel, opHitListHandle))
      art::fill_ptr_vector(ophitlist, opHitListHandle);
  }
  //std::cout << "ophits size=" << ophitlist.size() << std::endl;
  // Loop over ophits
  for (art::Ptr< OpHit > ophit : ophitlist) {
    vector<int> nearwires;
    double xyz[3] = {0.,0.,0.};
    geo->OpDetGeoFromOpChannel(ophit->OpChannel()).GetCenter(xyz);
    for (auto p : geo->IteratePlaneIDs()) nearwires.push_back(NearWire(*geo,p,xyz[0],xyz[1],xyz[2]));
    // Fill ophit table
    fOpHitNtuple->insert(evtID.data(),ophit.key(),
			 nearwires.data(),
			 ophit->OpChannel(),ophit->PeakTime(),ophit->Width(),
			 ophit->Area(),ophit->Amplitude(),ophit->PE()
    );
    LogInfo("HDF5Maker") << "Filling ophit table"
			 << "\nrun " << evtID[0] << ", subrun " << evtID[1]
			 << ", event " << evtID[2]
			 << "\nhit id " << ophit.key() << ", channel "
			 << ophit->OpChannel() << ", PeakTime " << ophit->PeakTime()
			 << ", Width " << ophit->Width()
			 << "\narea " << ophit->Area() << ", amplitude "
			 << ophit->Amplitude() << ", PE " << ophit->PE();
    // std::cout << std::endl;
    // std::cout << "det center " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
    // std::cout << "nearest wires " << geo->NearestWire(xyz,0) << " " << geo->NearestWire(xyz,1) << " " << geo->NearestWire(xyz,2) << std::endl;
    // for (auto p : geo->IteratePlaneIDs()) std::cout << "planeid " <<  p << " near wire " << geo->NearestWire(xyz,p) << std::endl;
  }

  // Get OpFlashs from the event record
  art::Handle< vector< OpFlash > > opFlashListHandle;
  vector< art::Ptr< OpFlash > > opflashlist;
  if (fOpFlashLabel != "") {
    if (e.getByLabel(fOpFlashLabel, opFlashListHandle))
      art::fill_ptr_vector(opflashlist, opFlashListHandle);
  }
  // Loop over opflashs
  //std::cout << "opflashs size=" << opflashlist.size() << std::endl;
  for (art::Ptr< OpFlash > opflash : opflashlist) {
    vector<float> pes;
    for (auto pe : opflash->PEs()) pes.push_back(pe);
    vector<int> nearwires;
    double xyz[3] = {0.,opflash->YCenter(),opflash->ZCenter()};
    for (auto p : geo->IteratePlaneIDs()) nearwires.push_back(NearWire(*geo,p,xyz[0],xyz[1],xyz[2]));
    // Fill opflash table
    fOpFlashNtuple->insert(evtID.data(),
			   opflash.key(),
			   nearwires.data(),
			   opflash->Time(),opflash->TimeWidth(),
			   opflash->YCenter(),opflash->YWidth(),opflash->ZCenter(),opflash->ZWidth(),
			   pes.data(),opflash->TotalPE()
    );
    //std::cout << "PE vec size=" << opflash->PEs().size() << std::endl;
    LogInfo("HDF5Maker") << "Filling opflash table"
      // std::cout << "\nFilling opflash table"
			 << "\nrun " << evtID[0] << ", subrun " << evtID[1]
			 << ", event " << evtID[2]
			 << "\nflash id " << opflash.key() << ", Time " << opflash->Time()
			 << ", TotalPE " << opflash->TotalPE()//;
			 << "\nWireCenters size0 " << opflash->WireCenters().size()//;
			 << "\nYCenter " << opflash->YCenter()<< " ZCenter " << opflash->ZCenter()
			 << "\nYWidth " << opflash->YWidth()<< " ZWidth " << opflash->ZWidth()
			 << "\nInBeamFrame " << opflash->InBeamFrame()<< " OnBeamTime " << opflash->OnBeamTime()
			 << "\nAbsTime " << opflash->AbsTime() << " TimeWidth " << opflash->TimeWidth() << " FastToTotal " << opflash->FastToTotal();
    //std::cout << std::endl; for (auto pe : opflash->PEs()) {std::cout << pe << " ";} std::cout << std::endl;
    // for (auto p : geo->IteratePlaneIDs()) std::cout << "planeid " <<  p << " near wire " << geo->NearestWire(xyz,p) << std::endl;
    // std::cout << "\nTriggerOffset=" << detProperties->TriggerOffset() << " TriggerOffsetTPC=" << detClocks->TriggerOffsetTPC() 
    // 	      << " TriggerTime=" << detClocks->TriggerTime() << " TPCTime=" << detClocks->TPCTime()
    // 	      << " G4RefTime=" << -detClocks->G4ToElecTime(0)
    // 	      << " Time2Tick=" << detClocks->Time2Tick(opflash->Time())
    // 	      << std::endl;
  }

} // function HDF5Maker::analyze

void HDF5Maker::beginSubRun(art::SubRun const& sr) {
  if (fFilesBySubrun) createFile(&sr);
}
void HDF5Maker::beginJob() {
  if (!fFilesBySubrun) createFile(0);
}
void HDF5Maker::createFile(const art::SubRun* sr) {

  struct timeval now;
  gettimeofday(&now, NULL);

  // Open HDF5 output
  std::ostringstream fileName;
  if (sr!=0) {
  fileName << fOutputName << "_r" << setfill('0') << setw(5) << sr->run()
    << "_s" << setfill('0') << setw(5) << sr->subRun() << "_ts" << setw(6) << now.tv_usec << ".h5";
  } else {
    fileName << fOutputName << "_ts" << setfill('0') << setw(6) << now.tv_usec << ".h5";
  }

  fFile = hep_hpc::hdf5::File(fileName.str(), H5F_ACC_TRUNC);

  if (fEventInfo == "none")
    fEventNtuple = new hep_hpc::hdf5::Ntuple(
      make_ntuple({fFile, "event_table", 1000},
        make_column<int>("event_id", 3)
    ));
  if (fEventInfo == "nu")
    fEventNtupleNu = new hep_hpc::hdf5::Ntuple(
      make_ntuple({fFile, "event_table", 1000},
        make_column<int>("event_id", 3),
        make_scalar_column<int>("is_cc"),
        make_scalar_column<float>("nu_energy"),
        make_scalar_column<float>("lep_energy"),
        make_column<float>("nu_dir", 3),
        make_column<float>("nu_vtx", 3),
        make_column<float>("nu_vtx_corr", 3),
        make_column<int>("nu_vtx_wire_pos", ServiceHandle<geo::Geometry>()->Nviews()),
        make_scalar_column<float>("nu_vtx_wire_time")
    ));

  if (fSPLabel!="") {
    fSpacePointNtuple = new hep_hpc::hdf5::Ntuple(
      make_ntuple({fFile, "spacepoint_table", 1000},
	make_column<int>("event_id", 3),
	make_scalar_column<int>("spacepoint_id"),
        make_column<float>("position", 3),
        make_column<int>("hit_id", ServiceHandle<geo::Geometry>()->Nviews())
    ));
  }

  if (fHitLabel != "") {
    fHitNtuple = new hep_hpc::hdf5::Ntuple(
      make_ntuple({fFile, "hit_table", 1000},
        make_column<int>("event_id", 3),
        make_scalar_column<int>("hit_id"),
        make_scalar_column<float>("integral"),
        make_scalar_column<float>("rms"),
        make_scalar_column<int>("tpc"),
        // make_scalar_column<int>("global_plane"),
        // make_scalar_column<float>("global_wire"),
        // make_scalar_column<float>("global_time"),
        make_scalar_column<int>("local_plane"),
        make_scalar_column<float>("local_wire"),
        make_scalar_column<float>("local_time")
    ));
  }

  fParticleNtuple = new hep_hpc::hdf5::Ntuple(
    make_ntuple({fFile, "particle_table", 1000},
      make_column<int>("event_id", 3),
      make_scalar_column<int>("g4_id"),
      make_scalar_column<int>("type"),
      make_scalar_column<int>("parent_id"),
      make_scalar_column<float>("momentum"),
      make_column<float>("start_position", 3),
      make_column<float>("end_position", 3),
      make_column<float>("start_position_corr", 3),
      make_column<float>("end_position_corr", 3),
      make_column<int>("start_wire_pos", ServiceHandle<geo::Geometry>()->Nviews()),
      make_column<int>("end_wire_pos", ServiceHandle<geo::Geometry>()->Nviews()),
      make_scalar_column<float>("start_wire_time"),
      make_scalar_column<float>("end_wire_time"),
      make_scalar_column<string>("start_process"),
      make_scalar_column<string>("end_process")
  ));

  fEnergyDepNtuple = new hep_hpc::hdf5::Ntuple(
    make_ntuple({fFile, "edep_table", 1000},
      make_column<int>("event_id", 3),
      make_scalar_column<int>("hit_id"),
      make_scalar_column<int>("g4_id"),
      make_scalar_column<float>("energy"),
      make_scalar_column<float>("energy_fraction")
  ));

  if (fWireLabel != "") {
    fWireNtuple = new hep_hpc::hdf5::Ntuple(
      make_ntuple({fFile, "wire_table", 1000},
        make_column<int>("event_id", 3),
        make_scalar_column<int>("tpc"),
        // make_scalar_column<int>("global_plane"),
        // make_scalar_column<float>("global_wire"),
        make_scalar_column<int>("local_plane"),
        make_scalar_column<float>("local_wire"),
        make_column<float>("adc", ServiceHandle<DetectorPropertiesService>()->provider()->ReadOutWindowSize())
    ));
  }

  if (fOpHitLabel != "") {
    fOpHitNtuple = new hep_hpc::hdf5::Ntuple(
      make_ntuple({fFile, "ophit_table", 1000},
        make_column<int>("event_id", 3),
        make_scalar_column<int>("hit_id"),
        make_scalar_column<int>("hit_channel"),
        make_column<int>("wire_pos", ServiceHandle<geo::Geometry>()->Nviews()),
        make_scalar_column<float>("peaktime"),
        make_scalar_column<float>("width"),
        make_scalar_column<float>("area"),
        make_scalar_column<float>("amplitude"),
        make_scalar_column<float>("pe")
    ));
  }

  if (fOpFlashLabel != "") {
    fOpFlashNtuple = new hep_hpc::hdf5::Ntuple(
      make_ntuple({fFile, "opflash_table", 1000},
        make_column<int>("event_id", 3),
        make_scalar_column<int>("flash_id"),
	make_column<int>("wire_pos", ServiceHandle<geo::Geometry>()->Nviews()),
        make_scalar_column<float>("time"),
        make_scalar_column<float>("time_width"),
        make_scalar_column<float>("y_center"),
        make_scalar_column<float>("y_width"),
        make_scalar_column<float>("z_center"),
        make_scalar_column<float>("z_width"),
        make_column<float>("pe", ServiceHandle<geo::Geometry>()->NOpDets()),
        make_scalar_column<float>("totalpe")
    ));
  }

  if (fPandoraLabel != "") {
    fPandoraPrimaryNtuple = new hep_hpc::hdf5::Ntuple(
      make_ntuple({fFile, "pandoraPrimary_table", 1000},
        make_column<int>("event_id", 3),
        make_scalar_column<int>("slice_id"),
        make_scalar_column<int>("pdg"),
        make_scalar_column<float>("nu_score"),
        make_column<float>("vtx", 3),
        make_column<int>("vtx_wire_pos", ServiceHandle<geo::Geometry>()->Nviews()),
        make_scalar_column<float>("vtx_wire_time")
    ));

    fPandoraPfpNtuple = new hep_hpc::hdf5::Ntuple(
      make_ntuple({fFile, "pandoraPfp_table", 1000},
        make_column<int>("event_id", 3),
        make_scalar_column<int>("pfp_id"),
        make_scalar_column<int>("pdg"),
        make_scalar_column<float>("tkrshr_score"),
        make_column<float>("vtx", 3),
        make_column<int>("vtx_wire_pos", ServiceHandle<geo::Geometry>()->Nviews()),
        make_scalar_column<float>("vtx_wire_time")
    ));

    fPandoraHitNtuple = new hep_hpc::hdf5::Ntuple(
      make_ntuple({fFile, "pandoraHit_table", 1000},
        make_column<int>("event_id", 3),
        make_scalar_column<int>("hit_id"),
        make_scalar_column<int>("slice_id"),
        make_scalar_column<int>("pfp_id")
    ));
  }

}

void HDF5Maker::endSubRun(art::SubRun const& sr) {
  if (fFilesBySubrun) closeFile();
}
void HDF5Maker::endJob() {
  if (!fFilesBySubrun) closeFile();
}
void HDF5Maker::closeFile() {
  if (fEventInfo == "none") delete fEventNtuple;
  if (fEventInfo == "nu")   delete fEventNtupleNu;
  if (fSPLabel  != "") delete fSpacePointNtuple;
  if (fHitLabel != "") delete fHitNtuple;
  delete fParticleNtuple;
  delete fEnergyDepNtuple;
  if (fWireLabel  != "") delete fWireNtuple;
  if (fOpHitLabel != "") delete fOpHitNtuple;
  if (fOpFlashLabel != "") delete fOpFlashNtuple;
  if (fPandoraLabel != "") {
    delete fPandoraPrimaryNtuple;
    delete fPandoraPfpNtuple;
    delete fPandoraHitNtuple;
  }
  fFile.close();
}

// apply the mapping of XYZ true -> XYZ position as it would be recosntructed.
// takes into account SCE, trigger time offset, and wirecell-pandora offset.
// to be applied to truth xyz in order to compare to reconstructed variables
// e.g. used for resolution plots
void HDF5Maker::True2RecoMappingXYZ(float& t, float& x, float& y, float& z)
{
  ApplySCEMappingXYZ(x, y, z);
  auto const &detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  double g4Ticks = detClocks->TPCG4Time2Tick(t) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
  float _xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);
  x += _xtimeoffset;
  x += fXOffset;
}

// apply the mapping of XYZ true -> XYZ position after SCE-induced shift.
// to be applied to truth xyz in order to compare to reconstructed variables
// e.g. used for resolution plots
void HDF5Maker::ApplySCEMappingXYZ(float& x, float& y, float& z)
{
  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  if (SCE->EnableSimSpatialSCE() == true)
    {
      auto offset = SCE->GetPosOffsets(geo::Point_t(x, y, z));
      x -= offset.X();
      y += offset.Y();
      z += offset.Z();
    }
}

void HDF5Maker::AddDaughters(const std::map<unsigned int, unsigned int>& _pfpmap,
			     const ProxyPfpElem_t &pfp_pxy,
			     const ProxyPfpColl_t &pfp_pxy_col,
			     std::vector<ProxyPfpElem_t> &slice_v)
{

  auto daughters = pfp_pxy->Daughters();

  slice_v.push_back(pfp_pxy);

  // if (fVerbose)
  //   std::cout << "\t PFP w/ PdgCode " << pfp_pxy->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;

  for (auto const &daughterid : daughters)
  {

    if (_pfpmap.find(daughterid) == _pfpmap.end())
    {
      std::cout << "Did not find DAUGHTERID in map! error" << std::endl;
      continue;
    }

    auto pfp_pxy2 = pfp_pxy_col.begin();
    for (size_t j = 0; j < _pfpmap.at(daughterid); ++j) ++pfp_pxy2;

    AddDaughters(_pfpmap,*pfp_pxy2, pfp_pxy_col, slice_v);

  } // for all daughters

  return;
} // AddDaughters

float HDF5Maker::GetMetaData(const ProxyPfpElem_t &pfp_pxy, string metaDataName)
{
  const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();
  if (pfParticleMetadataList.size() == 0) return -999.;
  for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j)
    {
      const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
      auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
      if (!pfParticlePropertiesMap.empty())
	{
	  for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
	    {
	      if (it->first == metaDataName)
		return it->second;
	    } // for map elements
	} // if pfp metadata map not empty
    } // for list
  return -999.;
}

int HDF5Maker::NearWire(const geo::Geometry& geo, const geo::PlaneID &ip, const float x, const float y, const float z)
{
  geo::PlaneGeo const& plane = geo.Plane(ip);
  geo::WireID wireID;
  try {
    wireID = plane.NearestWireID(geo::Point_t(x,y,z));
  }
  catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wireID = plane.ClosestWireID(e.suggestedWireID());
  }
  return wireID.Wire;
}
DEFINE_ART_MODULE(HDF5Maker)
