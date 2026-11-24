// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "Time.H"
#include "pointIOField.H"
#include "faceIOList.H"
#include "labelIOList.H"
#include "polyBoundaryMeshEntries.H"
#include <stdexcept>
#include <climits>
#include <algorithm>

namespace raw_resident_io {
using namespace Foam;
// Input metadata only. No native polyPatch/fvMesh or geometric calculations.
struct RawPatch {
 word name_, type_;
 label start_=0, size_=0;
 const word& name()const{return name_;}
 const word& type()const{return type_;}
 label start()const{return start_;}
 label size()const{return size_;}
 bool coupled()const{
  return type_!="patch" && type_!="wall" &&
         type_!="symmetry" && type_!="symmetryPlane";
 }
};
struct RawTopology {
 const Time* registry_=nullptr;
 pointField points_;
 faceList faces_;
 labelList owner_,neighbour_;
 List<RawPatch> patches_;
 label cells_=0;
 label nCells()const{return cells_;}
 label nFaces()const{return faces_.size();}
 label nInternalFaces()const{return neighbour_.size();}
 const pointField& points()const{return points_;}
 const faceList& faces()const{return faces_;}
 const labelList& faceOwner()const{return owner_;}
 const labelList& faceNeighbour()const{return neighbour_;}
 const List<RawPatch>& boundaryMesh()const{return patches_;}
};
inline const objectRegistry& input_registry(const RawTopology&m){return *m.registry_;}
inline bool geometric_wall(const RawPatch&p){return p.type()=="wall";}
inline RawTopology read_topology(const Time&time){
 RawTopology out;out.registry_=&time;
 const fileName directory="polyMesh";
 const word points_instance=time.findInstance(directory,"points");
 const word faces_instance=time.findInstance(directory,"faces");
 pointIOField points(IOobject("points",points_instance,directory,time,IOobject::MUST_READ,IOobject::NO_WRITE,false));
 faceCompactIOList faces(IOobject("faces",faces_instance,directory,time,IOobject::MUST_READ,IOobject::NO_WRITE,false));
 labelIOList owner(IOobject("owner",faces_instance,directory,time,IOobject::MUST_READ,IOobject::NO_WRITE,false));
 labelIOList neighbour(IOobject("neighbour",faces_instance,directory,time,IOobject::MUST_READ,IOobject::NO_WRITE,false));
 auto entries=polyBoundaryMeshEntries::readContents(IOobject("boundary",time.findInstance(directory,"boundary",IOobject::MUST_READ,faces_instance),directory,time,IOobject::MUST_READ,IOobject::NO_WRITE,false));
 out.points_=std::move(static_cast<pointField&>(points));
 out.faces_=std::move(static_cast<faceList&>(faces));
 out.owner_=std::move(static_cast<labelList&>(owner));
 out.neighbour_=std::move(static_cast<labelList&>(neighbour));
 if(out.faces_.empty()||out.owner_.size()!=out.faces_.size()||out.neighbour_.size()>out.faces_.size())
  throw std::runtime_error("raw topology extents");
 label largest=-1;
 for(label cell:out.owner_){if(cell<0||cell==INT_MAX)throw std::runtime_error("raw owner label");largest=std::max(largest,cell);}
 for(label cell:out.neighbour_){if(cell<0||cell==INT_MAX)throw std::runtime_error("raw neighbour label");largest=std::max(largest,cell);}
 out.cells_=largest+1;
 out.patches_.setSize(entries.size());
 label next=out.neighbour_.size();
 forAll(entries,p){
  const dictionary&d=entries[p].dict();auto&patch=out.patches_[p];
  patch.name_=entries[p].keyword();patch.type_=d.get<word>("type");
  patch.start_=d.get<label>("startFace");patch.size_=d.get<label>("nFaces");
  if(patch.start_!=next||patch.size_<0||patch.size_>out.faces_.size()-next||patch.coupled())
   throw std::runtime_error("raw unsupported patch or face range");
  next+=patch.size_;
 }
 if(next!=out.faces_.size())throw std::runtime_error("raw incomplete boundary ranges");
 return out;
}
}
