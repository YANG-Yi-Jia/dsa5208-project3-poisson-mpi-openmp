#include "options.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char** argv) {
    try {
        Options o;
        if(!parse_options(argc,argv,o)) return 0;
        poisson::Solver s(o.dim,o.n);
        int completed=0;
        bool converged=false;
        double reduction_sec=0;
        auto start=std::chrono::steady_clock::now();
        int limit=o.mode=="benchmark"?o.iterations:o.max_iters;
        for(int t=0;t<limit;++t) {
            s.step(); completed=t+1;
            if(o.mode=="validation" && (completed%o.check_every==0 || completed==limit)) {
                auto check=std::chrono::steady_clock::now();
                converged=s.metrics().relative_residual<o.tol;
                reduction_sec+=std::chrono::duration<double>(std::chrono::steady_clock::now()-check).count();
                if(converged) break;
            }
        }
        double runtime=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        auto m=s.metrics();
        if(!std::isfinite(m.relative_residual) || !std::isfinite(m.l2_rms_error)) throw std::runtime_error("nonfinite solution metrics");
        std::ifstream existing(o.output,std::ios::binary|std::ios::ate);
        bool header=!existing || existing.tellg()==0;
        existing.close();
        std::ofstream out(o.output,std::ios::app);
        if(!out) throw std::runtime_error("cannot open output: "+o.output);
        const std::string columns="version,dim,n,mpi_ranks,omp_threads,total_cores,iterations,mode,runtime_sec,compute_sec,halo_sec,reduction_sec,relative_residual,l2_rms_error,linf_error,boundary_max,converged";
        if(header) out<<columns<<'\n';
        auto record=[&](std::ostream& stream) {
            stream<<std::setprecision(17)<<"serial,"<<o.dim<<','<<o.n<<",1,1,1,"<<completed<<','<<o.mode<<','<<runtime<<','<<runtime-reduction_sec<<",0,"<<reduction_sec<<','<<m.relative_residual<<','<<m.l2_rms_error<<','<<m.linf_error<<','<<m.boundary_max<<','<<(o.mode=="validation"?(converged?"true":"false"):"NA")<<'\n';
        };
        record(out); out.close(); if(!out) throw std::runtime_error("failed writing output");
        if(!o.field.empty()) {
            std::ofstream field(o.field);
            if(!field) throw std::runtime_error("cannot open field output");
            field<<std::setprecision(17);
            s.visit([&](std::size_t p,int,int,int) { field<<s.old[p]<<'\n'; });
            field.close(); if(!field) throw std::runtime_error("failed writing field");
        }
        std::cout<<columns<<'\n'; record(std::cout);
        return o.mode=="validation" && !converged ? 2 : 0;
    } catch(const std::exception& e) { std::cerr<<"error: "<<e.what()<<'\n'; return 1; }
}
