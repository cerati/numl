////////////////////////////////////////////////////////////////////////
// Class:       PandoraNuSliceHitsProducer
// Plugin Type: producer (art v3_06_03)
// File:        PandoraNuSliceHitsProducer_module.cc
//
// Generated at Tue May 25 10:39:19 2021 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"


#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"


#include <memory>

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "lardataobj/RecoBase/SpacePoint.h"


#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


class PandoraNuSliceHitsProducer;
class PandoraNuSliceHitsProducer : public art::EDProducer {
public:
  explicit PandoraNuSliceHitsProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PandoraNuSliceHitsProducer(PandoraNuSliceHitsProducer const&) = delete;
  PandoraNuSliceHitsProducer(PandoraNuSliceHitsProducer&&) = delete;
  PandoraNuSliceHitsProducer& operator=(PandoraNuSliceHitsProducer const&) = delete;
  PandoraNuSliceHitsProducer& operator=(PandoraNuSliceHitsProducer&&) = delete;


private:

  // Declare member data here.
  std::string fOpFlashLabel;
  std::string fPfpLabel;
  std::string fSliceLabel;
  std::string fHitLabel;
  std::string fHitTruthLabel;

  std::vector<int> nearwires;

 // Required functions.
  void produce(art::Event& e) override;
  int NearWire(const geo::Geometry& geo, const geo::PlaneID &ip, const float x, const float y, const float z);
  std::vector<float> GetBarycenter(std::vector<float> hx,std::vector<float> hy,std::vector<float> hz,std::vector<float> hi); 

};


PandoraNuSliceHitsProducer::PandoraNuSliceHitsProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fOpFlashLabel( p.get<std::string>("OpFlashLabel","opflashCryoE")),
  fPfpLabel(p.get<std::string>("PfpLabel","pandoraGausCryoE")),
  fSliceLabel(p.get<std::string>("SliceLabel","pandoraGausCryoE")),
  fHitLabel(p.get<std::string>("HitLabel","cluster3DCryoE")),
  fHitTruthLabel(p.get<std::string>("HitTruthLabel","cluster3DCryoE"))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Hit>>();
   produces<std::vector<recob::SpacePoint>>();
   // produces<HitParticleAssociations>();
  produces<art::Assns<recob::Hit, recob::SpacePoint>>();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void PandoraNuSliceHitsProducer::produce(art::Event& e)
{
  art::ServiceHandle<geo::Geometry> geo;
  //Inserting the optical flashes information

  // Get OpFlashs from the event record
  art::Handle< std::vector< recob::OpFlash > > opFlashListHandle;
  std::unique_ptr<art::FindManyP<recob::OpHit> > assocFlashHit;
  std::vector< art::Ptr< recob::OpFlash > > opflashlist;
  if (fOpFlashLabel != "") {
    if (e.getByLabel(fOpFlashLabel, opFlashListHandle)) art::fill_ptr_vector(opflashlist, opFlashListHandle);
    assocFlashHit = std::unique_ptr<art::FindManyP<recob::OpHit> >(new art::FindManyP<recob::OpHit>(opFlashListHandle, e, fOpFlashLabel));
  }

 
  std::vector<std::vector<float> > Opflashbary;
 

  // Get spacepoints from the event record
  art::Handle<std::vector<recob::SpacePoint>> spListHandle;
  std::vector<art::Ptr<recob::SpacePoint>> splist;
  if (e.getByLabel(fHitLabel, spListHandle))
    art::fill_ptr_vector(splist, spListHandle);
  // Get assocations from spacepoints to hits
  art::FindManyP<recob::Hit> fmp(spListHandle, e, fHitLabel);
  std::vector<std::vector<art::Ptr<recob::Hit>>> sp2Hit(splist.size());

  for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
    sp2Hit[spIdx] = fmp.at(spIdx);
  } // for spacepoint

  

  // Implementation of required member function here.

  auto outputHits = std::make_unique<std::vector<recob::Hit> >();

  // auto outputHitPartAssns = std::make_unique<HitParticleAssociations>();
  art::PtrMaker<recob::Hit> hitPtrMaker(e);
  art::PtrMaker<recob::SpacePoint> spsPtrMaker(e);

  auto outputSpacepoints = std::make_unique<std::vector<recob::SpacePoint> >();
  auto outputSpsHitAssns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();

  art::ValidHandle<std::vector<recob::PFParticle> > inputPfp = e.getValidHandle<std::vector<recob::PFParticle> >(fPfpLabel);
  auto assocPfpSlice = std::unique_ptr<art::FindManyP<recob::Slice> >(new art::FindManyP<recob::Slice>(inputPfp, e, fPfpLabel));

  art::ValidHandle<std::vector<recob::Slice> > inputSlice = e.getValidHandle<std::vector<recob::Slice> >(fSliceLabel);
  auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit> >(new art::FindManyP<recob::Hit>(inputSlice, e, fSliceLabel));

  art::Handle< std::vector< recob::Hit > > hitListHandle;
  e.getByLabel(fHitLabel, hitListHandle);

 
  std::multimap<size_t,size_t> slc2hit;
  std::multimap<size_t,size_t> hit2slc;
  std::vector<size_t> slice_keys;

  // Loop over opflashs
  for (art::Ptr< recob::OpFlash > opflash : opflashlist) {
    std::vector<float> pes;
    for (auto pe : opflash->PEs()) pes.push_back(pe);
    std::vector<int> nearwires;
    std::vector<float> thiswire;
    double xyz[3] = {0.,opflash->YCenter(),opflash->ZCenter()};
    for (auto const& p : geo->Iterate<geo::PlaneID>()) nearwires.push_back(NearWire(*geo,p,xyz[0],xyz[1],xyz[2]));
    for (auto const& p : geo->Iterate<geo::PlaneID>()) thiswire.push_back(NearWire(*geo,p,xyz[0],xyz[1],xyz[2]));
    if(abs(opflash->Time())< 10)  Opflashbary.push_back(thiswire); // 10 microseconds interval chose (arbitrary you can modify it ot be closer to the trigger)
    thiswire.clear();
  }

  for (size_t ipfp = 0; ipfp<inputPfp->size(); ipfp++) {
    art::Ptr<recob::PFParticle> pfp(inputPfp,ipfp);
    if (pfp->IsPrimary()==false) continue;
    auto PDG = fabs(pfp->PdgCode());
    if (PDG != 12 && PDG != 14) continue;

    auto assocSlice = assocPfpSlice->at(pfp.key());
    auto sliceHits = assocSliceHit->at(assocSlice[0].key());
    slice_keys.push_back(assocSlice[0].key());
    for (size_t ihit = 0; ihit<sliceHits.size();ihit ++)
      {          
	slc2hit.insert({assocSlice[0].key(),sliceHits[ihit].key()}); 
	hit2slc.insert({sliceHits[ihit].key(),assocSlice[0].key()});
      } //loop over pfp hits
  }
  
  // Get the spacepoints
  std::multimap<int,int> slc2sps; 
  //std::cout << "SPS size=" << spListHandle->size() << std::endl;
  for (size_t i = 0; i < spListHandle->size(); ++i) {

    art::Ptr<recob::SpacePoint> spacePointPtr(spListHandle,i);
    std::vector<art::Ptr<recob::Hit>> associatedHits(fmp.at(spacePointPtr.key()));
    if (associatedHits.size() < 3)
      {
	//mf::LogDebug("SpacePointAnalysis") << "I am certain this cannot happen... but here you go, space point with " << associatedHits.size() << " hits" << std::endl;
	continue;
      }

    std::array<int, 3> hitID { -1, -1, -1 };
    for (size_t j = 0; j < associatedHits.size(); ++j)
      { 
	int plane = associatedHits[j]->WireID().Plane;
        if(associatedHits[j]->WireID().TPC==2 || associatedHits[j]->WireID().TPC==3)
          {
            if(associatedHits[j]->WireID().Plane==1) plane=2;
            else if(associatedHits[j]->WireID().Plane==2) plane=1;
          }
	//std::cout << "SLICE j=" << j << " tpc=" << associatedHits[j]->WireID().TPC << " plane=" << associatedHits[j]->WireID().Plane << " plane2=" << plane << " key=" <<  associatedHits[j].key() << std::endl;      
        hitID[plane] = associatedHits[j].key();
      }

    if(hitID[0]==-1 ||hitID[1]==-1 || hitID[2]==-1 ) continue;
    int slcID0 = hit2slc.find(hitID[0])->second;
    int slcID1 = hit2slc.find(hitID[1])->second;
    int slcID2 = hit2slc.find(hitID[2])->second;
    if(slcID0 != slcID1 || slcID2 != slcID0) continue;
    //std::cout << "SLICE id=" << slcID0 << " sp i=" << i << " hit ids=" << hitID[0] << " " << hitID[1] << " " << hitID[2] << std::endl;
    slc2sps.insert({slcID0,spacePointPtr.key()});
  }
  //std::cout << "SPS maps size=" << slc2sps.size() << std::endl;

  float mindist= 999999;
  size_t minkey= 99999;
  for (auto const& key : slice_keys)
    {
      auto p = slc2sps.equal_range(key);
      std::vector<float> slicebar0;
      std::vector<float> slicebar1;
      std::vector<float> slicebar2; 
      std::vector<float> slicebarcharge;

      for (auto& q = p.first; q != p.second; ++q)
	{
	  art::Ptr<recob::SpacePoint> spacePointPtr(spListHandle,q->second);
	  slicebar0.push_back(0);
	  slicebar1.push_back(spacePointPtr->XYZ()[1]);
	  slicebar2.push_back(spacePointPtr->XYZ()[2]);

	  std::vector<art::Ptr<recob::Hit>> associatedHits(fmp.at(spacePointPtr.key()));
	  if (associatedHits.size() < 2)
	    {
	      mf::LogDebug("SpacePointAnalysis") << "I am certain this cannot happen... but here you go, space point with " << associatedHits.size() << " hits" << std::endl;
	      continue;
	    }

	  for (size_t j = 0; j < associatedHits.size(); ++j)
	    { 
	      int plane = associatedHits[j]->WireID().Plane;
	      if(associatedHits[j]->WireID().TPC==2 || associatedHits[j]->WireID().TPC==3)
		{
		  if(associatedHits[j]->WireID().Plane==1) plane=2;
		  else if(associatedHits[j]->WireID().Plane==2) plane=1;
		}
	      if(plane==2) slicebarcharge.push_back(associatedHits[j]->Integral());
	      else continue;
	    }
	}

      std::vector<float> slcbary = GetBarycenter(slicebar0,slicebar1,slicebar2,slicebarcharge);       
      for(size_t kk=0; kk<Opflashbary.size();kk++)
	{
	  std::vector<float> localbardiff={0,9999,9999};
	  localbardiff= {abs(Opflashbary[kk][0]-slcbary[0]),abs(Opflashbary[kk][1]-slcbary[1]),abs(Opflashbary[kk][2]-slcbary[2])};
	  float auxdist= sqrt( pow(localbardiff[1],2)+ pow(localbardiff[2],2));           
	  if (auxdist<mindist)
	    { 
	      mindist=auxdist;
	      minkey=key;
	    }
	}
    }

  //std::cout << "slice minkey=" << minkey << std::endl;
  auto p = slc2hit.equal_range(minkey);  
  std::map<size_t,size_t> oldNewIdxMap;

  //std::cout << "slice hits size=" << slc2hit.count(minkey) << std::endl;
  for (auto& q = p.first; q != p.second; ++q)
    { 
      if(q->first != minkey) std::cout<< "the key is not the same, the value on this one is "<<q->first<< "while minkey is "<<minkey<<std::endl; 
      if(q->first != minkey) continue; 
      art::Ptr<recob::Hit> hit(hitListHandle,q->second);
      //std::cout << "map insert f=" << q->second << " s=" << outputHits->size() << std::endl;
      oldNewIdxMap.insert({size_t(q->second),outputHits->size()});
      outputHits->emplace_back(*hit);  
    }
  auto m = slc2sps.equal_range(minkey);    
  //std::cout << "slice sps size=" << slc2sps.count(minkey) << std::endl;
  //size_t kk = 0;
  for (auto& q = m.first; q != m.second; ++q)
    {
      art::Ptr<recob::SpacePoint> spacePointPtr(spListHandle,q->second);	 
      outputSpacepoints->emplace_back(*spacePointPtr); 

      const art::Ptr<recob::SpacePoint> sps = spsPtrMaker(outputSpacepoints->size()-1);
      std::vector<art::Ptr<recob::Hit>> associatedHits(fmp.at(spacePointPtr.key()));
      for (size_t j = 0; j < associatedHits.size(); ++j)
	{
	  const art::Ptr<recob::Hit> ahp = hitPtrMaker(oldNewIdxMap.find(associatedHits[j].key())->second);
	  outputSpsHitAssns->addSingle(ahp,sps);
	  //std::cout << "kk=" << kk << " q->first=" << q->first << " associatedHits[j].key()=" << associatedHits[j].key() << " map=" << oldNewIdxMap.find(associatedHits[j].key())->second << std::endl;
	}
      //kk++;
    }
 
  Opflashbary.clear();
 
  e.put(std::move(outputHits));
  e.put(std::move(outputSpacepoints));
  e.put(std::move(outputSpsHitAssns));

}

int PandoraNuSliceHitsProducer::NearWire(const geo::Geometry& geo, const geo::PlaneID &ip, const float x, const float y, const float z)
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

std::vector<float> PandoraNuSliceHitsProducer::GetBarycenter(std::vector<float> hx,std::vector<float> hy,std::vector<float> hz,std::vector<float> hi)
{
  float average_z_charge=0;
  float average_y_charge=0;
  float average_x_charge=0;
  float total_charge=0;
  size_t counter=9999;

  if(hz.size()<=hy.size()) counter = hz.size();
  else if(hy.size()<hz.size()) counter=hy.size();

  if(counter<=hx.size()) counter= hx.size();


  for (size_t  i = 0; i < hz.size(); i++) {
    if(isnan(hz[i])|| isnan(hx[i]) || isnan(hy[i])) continue;
    if(abs(hz[i])>9999 || abs(hy[i])>9999 || abs(hx[i])>9999 || abs(hi[i])>9999) continue;
    if(hi[i]<0) continue;
    average_z_charge+=hz[i]*hi[i];
    average_y_charge+=hy[i]*hi[i];
    average_x_charge+=hx[i]*hi[i];
    total_charge+=hi[i];

  }
  average_z_charge=average_z_charge/total_charge;
  average_y_charge=average_y_charge/total_charge;
  average_x_charge=average_x_charge/total_charge;
  std::vector<float> ThisTrackBary = {average_x_charge, average_y_charge, average_z_charge};
  return ThisTrackBary;
}

DEFINE_ART_MODULE(PandoraNuSliceHitsProducer)
