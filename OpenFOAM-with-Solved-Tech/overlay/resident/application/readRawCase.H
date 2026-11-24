// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "fvCFD.H"
#include "wallPolyPatch.H"
#include "XianzhiPreparation2647.hpp"
#include "readRawTopology.H"
#include <chrono>
#include <algorithm>
#include <cstdlib>
#include <cctype>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
namespace raw_resident_io {
using namespace Foam;
using namespace xk_fv_resident;
struct ParsedField {
 std::vector<double> internal;
 HostBoundary boundary;
 std::vector<int> owner_seed;
 std::vector<word> patch_type;
 bool has_spalding=false;double spalding_tolerance=.01;int spalding_iterations=10;
};
struct RawCase {
 HostResidentInput input;
 HostInitialGeometry geometry;
 double read_seconds=0,marshal_seconds=0;
};
template<class Type> struct RawField {
 dictionary metadata;
 Field<Type> internal;
 std::vector<Field<Type>> patches;
};
struct RawControls {dictionary transport,turbulence,schemes,solution,options,mrf;};
struct RawParsedCase {
 RawField<vector> U;
 RawField<scalar> p,k,omega,nut;
 RawControls controls;
 double read_seconds=0.;
};
inline std::string compact(std::string s){s.erase(std::remove_if(s.begin(),s.end(),[](unsigned char c){return std::isspace(c);}),s.end());return s;}
inline void require(bool ok,const std::string&message){if(!ok)throw std::runtime_error("raw case: "+message);}
inline void coefficient(const dictionary&dict,const word&name,double expected){if(dict.found(name))require(dict.get<scalar>(name)==expected,"unsupported coefficient "+std::string(name));}
inline const objectRegistry& input_registry(const fvMesh&mesh){return mesh;}
inline bool geometric_wall(const polyPatch&p){return isA<wallPolyPatch>(p);}
template<class Patch> inline const dictionary& patch_dictionary(const dictionary&bf,const Patch&p){const dictionary*d=bf.findDict(p.name(),keyType::REGEX);require(d!=nullptr,"missing boundary dictionary "+std::string(p.name()));return*d;}
template<class Type,class Mesh> RawField<Type> read_field(const Mesh&mesh,const Time&time,const word&name,const word&class_name,const dimensionSet&dimensions,int components){
 RawField<Type> out;
 out.metadata=IOdictionary::readContents(IOobject(name,time.timeName(),input_registry(mesh),IOobject::MUST_READ,IOobject::NO_WRITE,false),class_name);
 require(dimensionSet(out.metadata.lookup("dimensions"))==dimensions,"field dimensions "+std::string(name));
 require(!out.metadata.found("referenceLevel"),"unregistered field referenceLevel");
 out.internal=Field<Type>("internalField",out.metadata,mesh.nCells());
 out.metadata.remove("internalField");
 out.patches.resize(mesh.boundaryMesh().size());
 const dictionary&bf=out.metadata.subDict("boundaryField");
 forAll(mesh.boundaryMesh(),p){
  const auto&patch=mesh.boundaryMesh()[p];const dictionary&d=patch_dictionary(bf,patch);const word type=d.get<word>("type");
  const bool owner=type=="zeroGradient"||type=="symmetry"||type=="symmetryPlane"||type=="slip"||type=="kqRWallFunction"||type=="omegaWallFunction";
  if(owner||type=="noSlip")continue;
  if(type=="uniformFixedValue"){
   auto&stream=d.lookup("uniformValue");word function;stream>>function;
   require(function=="constant","time-dependent uniformFixedValue");Type value;stream>>value;
   require(!stream.bad()&&stream.nRemainingTokens()==0,"uniformFixedValue constant syntax");
   for(int c=0;c<components;++c)require(std::isfinite(component(value,c)),"nonfinite uniformFixedValue");
   out.patches[p].setSize(patch.size(),value);
  }else out.patches[p]=Field<Type>("value",d,patch.size(),type=="calculated"?IOobjectOption::LAZY_READ:IOobjectOption::MUST_READ);
 }
 return out;
}
template<class Type,class Mesh> ParsedField field(const Mesh&mesh,const RawField<Type>&input,const word&name,int components){
 ParsedField result;const int n=mesh.nCells(),ni=mesh.nInternalFaces(),nb=mesh.nFaces()-ni;
 result.internal.resize(size_t(n)*components);for(int i=0;i<n;++i)for(int c=0;c<components;++c)result.internal[size_t(c)*n+i]=component(input.internal[i],c);
 result.boundary.kinds.resize(nb);result.boundary.assignable.resize(nb);result.boundary.prescribed.assign(size_t(components)*nb,0.);result.owner_seed.assign(nb,0);result.patch_type.resize(mesh.boundaryMesh().size());
 const dictionary&bf=input.metadata.subDict("boundaryField");
 forAll(mesh.boundaryMesh(),p){
  const auto&patch=mesh.boundaryMesh()[p];require(!patch.coupled(),"coupled patch "+std::string(patch.name()));const dictionary&d=patch_dictionary(bf,patch);const word type=d.get<word>("type");result.patch_type[p]=type;
  const bool wall=geometric_wall(patch),omega=type=="omegaWallFunction",nut=type=="nutUSpaldingWallFunction";
  const bool zero=type=="zeroGradient"||type=="symmetry"||type=="symmetryPlane"||type=="slip"||type=="kqRWallFunction";
  const bool uniform_fixed=type=="uniformFixedValue";
  const bool fixed=type=="fixedValue"||uniform_fixed||type=="noSlip"||omega||nut||type=="calculated";
  require(zero||fixed,"unsupported field boundary "+std::string(name)+":"+std::string(type));
  require(type!="noSlip"||components==3,"scalar noSlip boundary");require(!omega||(name=="omega"&&wall),"omega wall pairing");require(!nut||(name=="nut"&&wall),"nut wall pairing");require(type!="kqRWallFunction"||(name=="k"&&wall),"k wall pairing");require(type!="calculated"||name=="nut","calculated field other than nut");
  if(type=="symmetry")require(patch.type()=="symmetry","symmetry field on incompatible patch");if(type=="symmetryPlane")require(patch.type()=="symmetryPlane","symmetryPlane field on incompatible patch");
  const bool from_owner=zero||omega;const bool assigned=!(type=="fixedValue"||uniform_fixed||type=="noSlip"||type=="slip"||omega||nut);
  const bool transform=type=="symmetry"||type=="symmetryPlane"||type=="slip";
  const auto kind=fixed?BoundaryKind::fixed_value:(transform?BoundaryKind::symmetry:BoundaryKind::zero_gradient);
  forAll(patch,j){const int at=patch.start()+j-ni;result.boundary.kinds[at]=kind;result.boundary.assignable[at]=assigned;result.owner_seed[at]=from_owner;}
  if(!from_owner&&type!="noSlip"){
   const auto&patch_values=input.patches[p];
   require(patch_values.size()==patch.size(),"parsed patch extent");
   forAll(patch,j)for(int c=0;c<components;++c)result.boundary.prescribed[size_t(c)*nb+patch.start()+j-ni]=component(patch_values[j],c);
  }
  if(omega||nut){coefficient(d,"Cmu",.09);coefficient(d,"kappa",.41);coefficient(d,"E",9.8);require(!d.found("U")||d.get<word>("U")=="U","alternate wall velocity field");}
  if(omega){coefficient(d,"beta1",.075);coefficient(d,"n",2.);require(d.getOrDefault<word>("blending","binomial")=="binomial","omega wall blending");}
  if(nut){const double tolerance=d.getOrDefault<scalar>("tolerance",.01);const int iterations=d.getOrDefault<label>("maxIter",10);require(std::isfinite(tolerance)&&tolerance>0.&&iterations>0,"Spalding controls");if(result.has_spalding)require(result.spalding_tolerance==tolerance&&result.spalding_iterations==iterations,"different per-patch Spalding controls");result.has_spalding=true;result.spalding_tolerance=tolerance;result.spalding_iterations=iterations;}
 }
 return result;
}
inline void scheme(const dictionary&dict,const word&key,const std::string&expected){require(compact(dict.lookup(key).toString())==compact(expected),"unsupported scheme "+std::string(key));}
inline void controls(const RawControls&input,HostResidentInput&h,HostInitialGeometry&geometry){
 const auto&transport=input.transport;require(transport.get<word>("transportModel")=="Newtonian","non-Newtonian transport");geometry.uniform_nu=dimensionedScalar("nu",dimensionSet(0,2,-1,0,0,0,0),transport).value();require(std::isfinite(geometry.uniform_nu)&&geometry.uniform_nu>0.,"molecular viscosity");
 const auto&turbulence=input.turbulence;require(turbulence.get<word>("simulationType")=="RAS","simulation type");const auto&ras=turbulence.subDict("RAS");require(ras.get<word>("RASModel")=="kOmegaSST"&&ras.get<Switch>("turbulence"),"RAS model");
 const dictionary&coeff=ras.optionalSubDict("kOmegaSSTCoeffs");for(auto value:{std::make_pair("alphaK1",.85),{"alphaK2",1.},{"alphaOmega1",.5},{"alphaOmega2",.856},{"gamma1",5./9.},{"gamma2",.44},{"beta1",.075},{"beta2",.0828},{"betaStar",.09},{"a1",.31},{"b1",1.},{"c1",10.},{"kInf",0.},{"omegaInf",0.}})coefficient(coeff,value.first,value.second);require(!coeff.getOrDefault<Switch>("F3",false)&&!coeff.getOrDefault<Switch>("decayControl",false),"SST switches");
 if(ras.found("kMin"))h.k_min=dimensionedScalar("kMin",dimensionSet(0,2,-2,0,0,0,0),ras).value();if(ras.found("omegaMin"))h.omega_min=dimensionedScalar("omegaMin",dimensionSet(0,0,-1,0,0,0,0),ras).value();
 const auto&schemes=input.schemes;scheme(schemes.subDict("ddtSchemes"),"default","steadyState");scheme(schemes.subDict("gradSchemes"),"default","Gauss linear");scheme(schemes.subDict("gradSchemes"),"grad(U)","cellLimited Gauss linear 1");scheme(schemes.subDict("divSchemes"),"div(phi,U)","bounded Gauss linearUpwindV grad(U)");for(auto name:{"div(phi,k)","div(phi,omega)"})scheme(schemes.subDict("divSchemes"),name,"bounded Gauss upwind");scheme(schemes.subDict("divSchemes"),"div((nuEff*dev2(T(grad(U)))))","Gauss linear");scheme(schemes.subDict("laplacianSchemes"),"default","Gauss linear corrected");scheme(schemes.subDict("interpolationSchemes"),"default","linear");scheme(schemes.subDict("snGradSchemes"),"default","corrected");const auto&wd=schemes.subDict("wallDist");require(wd.get<word>("method")=="meshWave"&&wd.getOrDefault<Switch>("correctWalls",true)&&!wd.getOrDefault<Switch>("nRequired",false),"wall-distance controls");
 const auto&solution=input.solution;const auto&simple=solution.subDict("SIMPLE");require(simple.getOrDefault<label>("nNonOrthogonalCorrectors",0)==0&&simple.getOrDefault<Switch>("consistent",false)&&simple.getOrDefault<Switch>("momentumPredictor",true),"SIMPLE controls");
 const auto&solvers=solution.subDict("solvers");for(auto method:{std::make_pair("U",15),{"k",8},{"omega",8}}){const auto&s=solvers.subDict(method.first);word name=s.get<word>("solver");if(name=="BenchLinear")name=s.get<word>("nativeSolver");require(name=="smoothSolver"&&s.get<word>("smoother")=="GaussSeidel"&&s.getOrDefault<label>("nSweeps",1)>0&&method.second%s.getOrDefault<label>("nSweeps",1)==0&&s.getOrDefault<label>("minIter",0)<=method.second&&s.get<label>("maxIter")==method.second&&s.get<scalar>("tolerance")==0.&&s.get<scalar>("relTol")==0.,"fixed-sweep solver "+std::string(method.first));}
 const auto&relax=solution.subDict("relaxationFactors");const auto&equations=relax.subDict("equations");coefficient(equations,"U",.9);coefficient(equations,"k",.6);coefficient(equations,"omega",.6);if(auto*fields=relax.findDict("fields"))coefficient(*fields,"p",1.);
 const auto&options=input.options;forAllConstIter(dictionary,options,it){if(!it().isDict())continue;const auto&option=it().dict();if(!option.getOrDefault<Switch>("active",true))continue;require(option.get<word>("type")=="velocityDampingConstraint","fvOption");const auto&c=option.optionalSubDict("velocityDampingConstraintCoeffs");require(c.get<word>("selectionMode")=="all"&&!c.found("UNames")&&c.getOrDefault<word>("U","U")=="U"&&h.damping_maximum==0.,"damping selection");h.damping_maximum=c.get<scalar>("UMax");h.damping_coefficient=c.getOrDefault<scalar>("C",1.);}
 const auto&mrf=input.mrf;forAllConstIter(dictionary,mrf,it)if(it().isDict())require(!it().dict().getOrDefault<Switch>("active",true),"active MRF");
}
template<class Mesh> inline RawParsedCase read_input(const Mesh&mesh,const Time&time){
 const auto begin=std::chrono::steady_clock::now();RawParsedCase input;
 input.U=read_field<vector>(mesh,time,"U","volVectorField",dimensionSet(0,1,-1,0,0,0,0),3);
 input.p=read_field<scalar>(mesh,time,"p","volScalarField",dimensionSet(0,2,-2,0,0,0,0),1);
 input.k=read_field<scalar>(mesh,time,"k","volScalarField",dimensionSet(0,2,-2,0,0,0,0),1);
 input.omega=read_field<scalar>(mesh,time,"omega","volScalarField",dimensionSet(0,0,-1,0,0,0,0),1);
 input.nut=read_field<scalar>(mesh,time,"nut","volScalarField",dimensionSet(0,2,-1,0,0,0,0),1);
 input.controls.transport=IOdictionary(IOobject("transportProperties",time.constant(),input_registry(mesh),IOobject::MUST_READ,IOobject::NO_WRITE,false));
 input.controls.turbulence=IOdictionary(IOobject("turbulenceProperties",time.constant(),input_registry(mesh),IOobject::MUST_READ,IOobject::NO_WRITE,false));
 input.controls.schemes=IOdictionary(IOobject("fvSchemes",time.system(),input_registry(mesh),IOobject::MUST_READ,IOobject::NO_WRITE,false));
 input.controls.solution=IOdictionary(IOobject("fvSolution",time.system(),input_registry(mesh),IOobject::MUST_READ,IOobject::NO_WRITE,false));
 input.controls.options=IOdictionary(IOobject("fvOptions",time.system(),input_registry(mesh),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE,false));
 input.controls.mrf=IOdictionary(IOobject("MRFProperties",time.constant(),input_registry(mesh),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE,false));
 input.read_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();return input;
}
template<class Mesh> inline RawCase marshal(const Mesh&mesh,const RawParsedCase&input){
 const auto begin=std::chrono::steady_clock::now();RawCase out;auto&h=out.input;auto&g=out.geometry;
 auto U=field<vector>(mesh,input.U,"U",3);auto p=field<scalar>(mesh,input.p,"p",1);
 auto k=field<scalar>(mesh,input.k,"k",1);auto w=field<scalar>(mesh,input.omega,"omega",1);auto nut=field<scalar>(mesh,input.nut,"nut",1);
 h.velocity=std::move(U.internal);h.pressure=std::move(p.internal);h.k=std::move(k.internal);h.omega=std::move(w.internal);h.nut=std::move(nut.internal);h.U=std::move(U.boundary);h.p=std::move(p.boundary);h.K=std::move(k.boundary);h.W=std::move(w.boundary);h.Nut=std::move(nut.boundary);g.k_owner_seed=std::move(k.owner_seed);g.omega_owner_seed=std::move(w.owner_seed);g.nut_owner_seed=std::move(nut.owner_seed);h.spalding_tolerance=nut.spalding_tolerance;h.spalding_iterations=nut.spalding_iterations;
 h.cells=mesh.nCells();h.internal_faces=mesh.nInternalFaces();h.owner.assign(mesh.faceOwner().begin(),mesh.faceOwner().end());h.neighbour.assign(mesh.faceNeighbour().begin(),mesh.faceNeighbour().end());const int nf=mesh.nFaces(),nb=nf-h.internal_faces;h.use_velocity_predictor_boundary.resize(nb);for(int b=0;b<nb;++b)h.use_velocity_predictor_boundary[b]=!h.U.assignable[b];
 h.face_offsets.resize(size_t(nf)+1);for(int f=0;f<nf;++f){for(int v:mesh.faces()[f])h.face_vertices.push_back(v);require(h.face_vertices.size()<=INT32_MAX,"polygon incidence exceeds int32");h.face_offsets[f+1]=int(h.face_vertices.size());}h.points.reserve(size_t(mesh.points().size())*3);for(auto p:mesh.points())for(int a=0;a<3;++a)h.points.push_back(p[a]);
 forAll(mesh.boundaryMesh(),p){const auto&patch=mesh.boundaryMesh()[p];require(patch.type()=="patch"||patch.type()=="wall"||patch.type()=="symmetry"||patch.type()=="symmetryPlane","patch geometry override");const bool wall=geometric_wall(patch),omega=w.patch_type[p]=="omegaWallFunction",spalding=nut.patch_type[p]=="nutUSpaldingWallFunction";require(omega==spalding&&(!omega||wall),"unpaired turbulence wall functions");forAll(patch,j){const int f=patch.start()+j;if(wall)g.geometric_wall_faces.push_back(f);if(omega)h.wall_faces.push_back(f-h.internal_faces);}}
 controls(input.controls,h,g);out.read_seconds=input.read_seconds;out.marshal_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
 const char*gs=std::getenv("RESIDENT_GS_CONFIG"),*amg=std::getenv("RESIDENT_PRESSURE_CONFIG");require(gs&&amg,"solver JSON paths missing");h.one_sweep_config=gs;h.pressure_config=amg;return out;
}
template<class Mesh> inline RawCase read(const Mesh&mesh,const Time&time){return marshal(mesh,read_input(mesh,time));}
}
