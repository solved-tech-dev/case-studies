// SPDX-License-Identifier: AGPL-3.0-only
// Copyright 2026 Solved Tech
#include "model.h"
#include "grid.h"
#include "pangulu_fwi_adapter.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Z = SolverComplex;
using Real = std::vector<double>;
using Wave = std::vector<Z>;

constexpr int nx = FWI_EX+1, ny = FWI_EY+1, nz = FWI_EZ+1;
constexpr int nodes = nx*ny*nz, ndof = 3*nodes;
constexpr int nnz = ndof+42*(FWI_EX-1)*(FWI_EY-1)*(FWI_EZ-1);
constexpr int iterations = 12;
constexpr double h = FWI_SPACING, ref_h = 85, pi = 3.141592653589793;
constexpr std::array<double, 2> frequencies{3., 5.};
int node(int x, int y, int z) { return (z*ny+y)*nx+x; }
void require(bool ok, const char* msg) { if (!ok) throw std::runtime_error(msg); }

int mapped_node(int x, int y, int z) {
    auto nearest=[](int i){return int(std::lround(i*ref_h/h));};
    const int i=nearest(x),j=nearest(y),k=nearest(z);
    require(i>0 && i<nx-1 && j>0 && j<ny-1 && k>0 && k<nz-1,
            "mapped acquisition node is outside the model interior");
    return node(i,j,k);
}

struct Receiver { int point; bool held_out; int ref_x,ref_y; };
struct Survey {
    std::vector<int> shots;
    std::vector<Receiver> receivers;
    Survey() {
        shots.push_back(mapped_node(49,7,2));  // Map the 85 m central vertical point force.
        for (int y : {3, 5, 7, 9, 11})
            for (int x = 5; x <= 95; x += 2)
                receivers.push_back({mapped_node(x,y,3), ((x-5)/2+2*y)%5 == 0, x, y});
        require(shots.size() == kRhs, "wrong number of source columns");
    }
    bool used(int, Receiver r) const {
        const int dx=49-r.ref_x, dy=7-r.ref_y;
        return dx*dx+dy*dy >= 25;  // retain the 85 m trace mask for the central shot
    }
    Wave forcing() const {
        Wave rhs(size_t(ndof)*kRhs);
        for (int s=0; s<kRhs; ++s) rhs[size_t(s)*ndof+3*shots[s]+2]=1.;
        return rhs;
    }
    Wave sample(const Wave& wave) const {
        Wave values;
        values.reserve(kRhs*receivers.size()*3);
        for (int s=0;s<kRhs;++s)
            for (auto r:receivers)
                for (int c=0;c<3;++c)
                    values.push_back(wave[size_t(s)*ndof+3*r.point+c]);
        return values;
    }
};

template<class T> void save(const std::string& path, const std::vector<T>& values) {
    std::ofstream f(path,std::ios::binary);
    f.write(reinterpret_cast<const char*>(values.data()), values.size()*sizeof(T));
    require(bool(f),"write output");
}
template<class T> std::vector<T> read(const std::string& path, size_t count) {
    std::ifstream f(path,std::ios::binary);
    std::vector<T> data(count);
    f.read(reinterpret_cast<char*>(data.data()), data.size()*sizeof(T));
    require(bool(f) && f.peek()==EOF,"observation file length");
    return data;
}
struct System {
    Matrix a{};
    explicit System(double freq, const Real& vp) {
        require(build_matrix(freq,vp.data(),&a)==ndof && a.nnz==nnz,"unexpected matrix dimensions");
    }
    ~System() { release_matrix(&a); }
    System(const System&) = delete;
};

Wave solve(const Matrix& a, Wave& rhs, bool adjoint, double frequency, const std::string& stage) {
    std::cout<<"FWI_SYSTEM stage="<<stage<<" frequency="<<frequency
             <<" orientation="<<(adjoint?"adjoint":"forward")<<" nrhs="<<kRhs<<std::endl;
    Real imag;
    if (adjoint) {
        imag.assign(a.ai,a.ai+a.nnz);
        for (double& v:imag) v=-v;
    }
    const int* r=adjoint?a.c:a.r;
    const int* c=adjoint?a.r:a.c;
    const double* ai=adjoint?imag.data():a.ai;
    Wave out(rhs.size());
    int status=case_factor(ndof,a.nnz,r,c,a.ar,ai);
    if (!status) status=case_solve_many(rhs.data(),out.data());
    const int cleanup=case_cleanup();
    require(!status && !cleanup,"native factorization, solve or cleanup failed");
    std::array<double,kRhs> errors{};
    #pragma omp parallel for num_threads(kRhs)
    for(int s=0;s<kRhs;++s) {
        const Z* x=out.data()+size_t(s)*ndof;
        const Z* b=rhs.data()+size_t(s)*ndof;
        Wave residual(b,b+ndof);
        for(int64_t i=0;i<a.nnz;++i) residual[r[i]]-=Z(a.ar[i],ai[i])*x[c[i]];
        double u=0,v=0;
        for(int i=0;i<ndof;++i){u+=std::norm(residual[i]);v+=std::norm(b[i]);}
        errors[s]=std::sqrt(u/(v>0?v:1.));
    }
    for(double error:errors) require(std::isfinite(error)&&error<=kTolerance,"an original COO RHS residual failed");
    const double worst=*std::max_element(errors.begin(),errors.end());
    std::cout<<"FWI_RESIDUAL stage="<<stage<<" frequency="<<frequency
             <<" orientation="<<(adjoint?"adjoint":"forward")<<" columns="<<kRhs<<" max_relative="<<worst<<std::endl;
    return out;
}

Real smooth(Real v, double sx, double sy, double sz) {
    // Widths are given in 85 m cells to preserve physical smoothing scales.
    const std::array<double,3> sigma{sx*ref_h/h,sy*ref_h/h,sz*ref_h/h};
    const std::array<int,3> size{nx,ny,nz}, stride{1,nx,nx*ny};
    Real tmp(nodes);
    for(int axis=0;axis<3;++axis){
        const int radius=int(std::ceil(2.7*sigma[axis]));
        Real kernel(2*radius+1);
        for(int d=-radius;d<=radius;++d) kernel[d+radius]=std::exp(-.5*std::pow(d/sigma[axis],2));
        const double sum=std::accumulate(kernel.begin(),kernel.end(),0.);
        for(double& w:kernel)w/=sum;
        #pragma omp parallel for num_threads(12)
        for(int q=0;q<nodes;++q) {
            const int c=(q/stride[axis])%size[axis];
            double val=0;
            for(int d=-radius;d<=radius;++d){
                const int other=std::clamp(c+d,0,size[axis]-1);
                val+=kernel[d+radius]*v[q+(other-c)*stride[axis]];
            }
            tmp[q]=val;
        }
        v.swap(tmp);
    }
    return v;
}
bool free_node(int q) {
    const int x=q%nx,y=(q/nx)%ny,z=q/(nx*ny);
    return x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1;
}
Real vp_to_q(const Real& vp) {
    Real q(nodes);for(int i=0;i<nodes;++i)q[i]=std::pow(vp[i]/1000,2);return q;
}
Real q_to_vp(const Real& q) {
    Real vp(nodes);for(int i=0;i<nodes;++i)vp[i]=1000*std::sqrt(q[i]);return vp;
}
double dot(const Real& a,const Real& b) { return std::inner_product(a.begin(),a.end(),b.begin(),0.); }

void accumulate_gradient(Real& g, const Real& q, const Matrix& a, const Wave& u,
                         const Wave& lambda, double freq) {
    // A is affine in q=(Vp/1000)^2 for fixed density and Vs=0.67*Vp.
    // Each axial term contributes to an off-diagonal and the opposite diagonal.
    const double scale=1e6/std::pow(h*2*pi*freq,2), shear=.67*.67;
    std::vector<Real> shot(kRhs, Real(nodes));
    #pragma omp parallel for num_threads(kRhs)
    for(int s=0;s<kRhs;++s){
        const Z* x=u.data()+size_t(s)*ndof;
        const Z* l=lambda.data()+size_t(s)*ndof;
        auto& accum=shot[s];
        for(int64_t k=0;k<a.nnz;++k){
            const int r=a.r[k],c=a.c[k],i=r/3,j=c/3;
            if(r==c) continue;
            const Z left=std::conj(l[r]);
            if(r%3==c%3){
                const int stride=std::abs(i-j);
                const int axis=stride==1?0:stride==nx?1:2;
                const double contribution=-.5*scale*(r%3==axis?1.:shear)*std::real(left*(x[r]-x[c]));
                accum[i]+=contribution;
                accum[j]+=contribution;
            } else {
                accum[i]-=std::real(left*x[c])*(a.ar[k]/q[i]);
            }
        }
    }
    for(const auto& v:shot)
        for(int i=0;i<nodes;++i) if(free_node(i))g[i]+=v[i];
}

struct Objective {
    const Survey& survey;
    std::array<Wave,frequencies.size()> observed;
    std::array<double,frequencies.size()> train_norm{}, test_norm{};
    Objective(const Survey& s,const std::string& dir):survey(s) {
        for(size_t f=0;f<frequencies.size();++f){
            observed[f]=read<Z>(dir+"/data-"+std::to_string(f)+".bin",kRhs*s.receivers.size()*3);
            size_t p=0;
            for(int shot=0;shot<kRhs;++shot)
                for(auto r:s.receivers) {
                    if(s.used(shot,r))
                        (r.held_out?test_norm[f]:train_norm[f])+=
                            std::norm(observed[f][p])+std::norm(observed[f][p+1])+std::norm(observed[f][p+2]);
                    p+=3;
                }
            require(std::isfinite(train_norm[f]) && train_norm[f]>0 &&
                    std::isfinite(test_norm[f]) && test_norm[f]>0,"observation energy");
        }
    }
    struct Evaluation {
        double loss=0, held_out=0;
        Real g;
        std::array<Wave,frequencies.size()> predicted;
    };
    Evaluation evaluate(const Real& q, bool gradient, const std::string& stage) const {
        Evaluation result;
        const Real vp=q_to_vp(q);
        if(gradient)result.g.assign(nodes,0.);
        for(size_t f=0;f<frequencies.size();++f){
            System system(frequencies[f],vp);
            auto rhs=survey.forcing();
            const Wave field=solve(system.a,rhs,false,frequencies[f],stage);
            auto& samples=result.predicted[f];
            samples=survey.sample(field);
            if(gradient) std::fill(rhs.begin(),rhs.end(),Z{});
            size_t p=0;
            const double a=1./(frequencies.size()*train_norm[f]);
            const double b=1./(frequencies.size()*test_norm[f]);
            for(int s=0;s<kRhs;++s)
                for(auto r:survey.receivers)
                    for(int c=0;c<3;++c,++p){
                        if(!survey.used(s,r))continue;
                        const Z mismatch=samples[p]-observed[f][p];
                        if(r.held_out) result.held_out+=.5*b*std::norm(mismatch);
                        else {
                            result.loss+=.5*a*std::norm(mismatch);
                            if(gradient)rhs[size_t(s)*ndof+3*r.point+c]=a*mismatch;
                        }
                    }
            if(gradient){
                auto adj=solve(system.a,rhs,true,frequencies[f],stage);
                accumulate_gradient(result.g,q,system.a,field,adj,frequencies[f]);
            }
        }
        require(std::isfinite(result.loss)&&std::isfinite(result.held_out),"nonfinite objective");
        return result;
    }
};

Real gradient_gate(const Objective& objective,const Real& q,const Real& grad) {
    Real direction(nodes);
    for(int i=0;i<nodes;++i){
        const double x=(i%nx)*h/ref_h,y=((i/nx)%ny)*h/ref_h,z=(i/(nx*ny))*h/ref_h;
        if(free_node(i)) direction[i]=.8*std::sin(.17*x+.21*z)*std::cos(.3*y)+.2*std::sin(.51*z+.31*y);
    }
    const double derivative=dot(direction,grad);
    require(std::isfinite(derivative)&&std::abs(derivative)>1e-9,"degenerate gradient check direction");
    double last=1.;
    for(double eps:{.03,.01,.003}){
        Real plus=q, minus=q;
        for(int i=0;i<nodes;++i){plus[i]+=eps*direction[i];minus[i]-=eps*direction[i];}
        const double fd=(objective.evaluate(plus,false,"gradient_gate").loss-
                         objective.evaluate(minus,false,"gradient_gate").loss)/(2*eps);
        last=std::abs(fd-derivative)/std::max(std::abs(fd),std::abs(derivative));
        std::cout<<"FWI_GRAD_CHECK epsilon="<<eps<<" analytic="<<derivative
                 <<" finite_difference="<<fd<<" relative_error="<<last<<std::endl;
    }
    require(last<1e-3,"directional gradient check did not pass");
    return direction;
}

struct Pair { Real step, delta; double inverse; };
Real direction_for(const Real& g,const std::vector<Pair>& memory) {
    Real d=g;
    Real alpha(memory.size());
    for(size_t n=memory.size();n>0;--n) {
        auto& p=memory[n-1];alpha[n-1]=p.inverse*dot(p.step,d);
        for(int i=0;i<nodes;++i)d[i]-=alpha[n-1]*p.delta[i];
    }
    d=smooth(d,1.4,1.05,1.15);
    for(int i=0;i<nodes;++i)if(!free_node(i))d[i]=0.;
    if(!memory.empty()){
        const auto& p=memory.back();
        const double scale=1./(p.inverse*dot(p.delta,p.delta));
        for(double& v:d)v*=scale;
    }
    for(size_t n=0;n<memory.size();++n){
        auto& p=memory[n];const double beta=p.inverse*dot(p.delta,d);
        for(int i=0;i<nodes;++i)d[i]+=p.step[i]*(alpha[n]-beta);
    }
    for(double& v:d)v=-v;
    return d;
}
double model_error(const Real& q,const Real& truth) {
    const auto vp=q_to_vp(q);
    double numerator=0,denominator=0;
    for(int i=0;i<nodes;++i) if(free_node(i)){
        numerator+=std::pow(vp[i]-truth[i],2);denominator+=truth[i]*truth[i];
    }
    return std::sqrt(numerator/denominator);
}
void emit_iteration(int n,const Objective::Evaluation& e,const Real& q,const Real& truth,
                    double step,int trials) {
    std::cout<<"FWI_ITER step="<<n<<" objective="<<e.loss<<" waveform_relative="<<std::sqrt(2*e.loss)
             <<" heldout_relative="<<std::sqrt(2*e.held_out)<<" vp_relative="<<model_error(q,truth)
             <<" step_length="<<step<<" trials="<<trials<<std::endl;
}
int main(int argc,char** argv) try {
    require(argc==4,"usage: executable generate|invert MODEL OBS_DIR");
    std::cout<<std::setprecision(15);
    Real truth(nodes);
    require(!load_truth(argv[2],truth.data()),"could not load true Marmousi Vp");
    const Survey survey;
    const std::string dir=argv[3];
    std::cout<<"FWI_CONFIG grid="<<nx<<"x"<<ny<<"x"<<nz<<" spacing_m="<<h
             <<" n="<<ndof<<" nnz="<<nnz<<" shots="<<survey.shots.size()
             <<" receivers="<<survey.receivers.size()
             <<" heldout="<<std::count_if(survey.receivers.begin(),survey.receivers.end(),[](auto r){return r.held_out;})
             <<" frequencies_hz=3,5 parameter=Vp Vs_over_Vp=0.67 rho=2000 reuse=off"<<std::endl;
    if(std::string(argv[1])=="generate") {
        save(dir+"/vp-true.bin",truth);
        Real initial=smooth(truth,5.,2.5,3.5);
        save(dir+"/vp-initial.bin",initial);
        for(size_t f=0;f<frequencies.size();++f){
            System system(frequencies[f],truth);
            auto rhs=survey.forcing();
            auto u=solve(system.a,rhs,false,frequencies[f],"synthetic_data");
            save(dir+"/data-"+std::to_string(f)+".bin",survey.sample(u));
        }
        std::cout<<"FWI_OBSERVATIONS_COMPLETE"<<std::endl;
        return 0;
    }
    require(std::string(argv[1])=="invert","invalid mode");
    require(truth==read<double>(dir+"/vp-true.bin",nodes),"truth mismatch between branches");
    const Real initial=read<double>(dir+"/vp-initial.bin",nodes);
    Real q=vp_to_q(initial), last_q,last_g;
    const Objective objective(survey,dir);
    {
        System system(5.,truth);
        auto rhs=survey.forcing();
        const Wave reference=objective.observed[1];
        const Wave actual=survey.sample(solve(system.a,rhs,false,5.,"matched_forward"));
        double error=0., energy=0.;
        for(size_t i=0;i<actual.size();++i){
            error+=std::norm(actual[i]-reference[i]);
            energy+=std::norm(reference[i]);
        }
        const double relative=std::sqrt(error/energy);
        std::cout<<"FWI_MATCHED_FORWARD relative_receiver_l2="<<relative<<std::endl;
        require(std::isfinite(relative)&&relative<=kTolerance,"matched forward wavefield check failed");
    }
    auto current=objective.evaluate(q,true,"initial_gradient");
    gradient_gate(objective,q,current.g);
    emit_iteration(0,current,q,truth,0,0);
    std::vector<Pair> memory;
    int accepted=0;
    for(int iteration=1;iteration<=iterations;++iteration){
        if(iteration>1){
            current=objective.evaluate(q,true,"iteration_gradient");
            Pair pair{Real(nodes),Real(nodes),0};
            for(int i=0;i<nodes;++i){pair.step[i]=q[i]-last_q[i];pair.delta[i]=current.g[i]-last_g[i];}
            const double cross=dot(pair.step,pair.delta);
            if(cross>1e-10*std::sqrt(dot(pair.step,pair.step)*dot(pair.delta,pair.delta))){
                pair.inverse=1./cross;memory.push_back(std::move(pair));
                if(memory.size()>5) memory.erase(memory.begin());
            }
        }
        auto d=direction_for(current.g,memory);
        if(dot(d,current.g)>=0.){memory.clear();d=direction_for(current.g,memory);}
        const double max_change=*std::max_element(d.begin(),d.end(),[](double a,double b){return std::abs(a)<std::abs(b);});
        require(std::isfinite(max_change) && std::abs(max_change)>0,"invalid descent direction");
        const double limit=.85/std::abs(max_change);
        double step=memory.empty()?limit:std::min(1.,limit);
        bool took_step=false;
        for(int trial=1;trial<=7;++trial){
            Real candidate=q, change(nodes);
            for(int i=0;i<nodes;++i){
                if(free_node(i))candidate[i]=std::clamp(q[i]+step*d[i],1.44,25.);
                change[i]=candidate[i]-q[i];
            }
            const double slope=dot(change,current.g);
            require(slope<0,"bounded direction is not descent");
            auto value=objective.evaluate(candidate,false,"line_search");
            std::cout<<"FWI_TRIAL step="<<iteration<<" trial="<<trial<<" loss="<<value.loss
                     <<" step_length="<<step<<std::endl;
            if(value.loss<=current.loss+1e-4*slope){
                last_q=std::move(q);last_g=std::move(current.g);q=std::move(candidate);
                emit_iteration(iteration,value,q,truth,step,trial);
                accepted=iteration;
                for(size_t f=0;f<frequencies.size();++f)
                    save("predicted-"+std::to_string(f)+".bin",value.predicted[f]);
                save("vp-final.bin",q_to_vp(q));
                took_step=true;
                break;
            }
            step*=.5;
        }
        if(!took_step){std::cout<<"FWI_STOP line_search_stalled iteration="<<iteration<<std::endl;break;}
    }
    require(accepted==iterations,"the inversion stopped before 12 model updates");
    std::cout<<"FWI_COMPLETE accepted_iterations="<<accepted<<" requested_iterations="<<iterations<<std::endl;
    return 0;
} catch(const std::exception& e) {
    case_cleanup();
    std::cerr<<"FWI_ERROR "<<e.what()<<std::endl;
    return 1;
}
