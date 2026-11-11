#include "options.hpp"
#include "mpi_solver.hpp"
#include <fstream>
#include <iomanip>

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    int rank=0,ranks=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try {
        Options o;
        if(!parse_options(argc,argv,o)) { MPI_Finalize(); return 0; }
        poisson::MPISolver s(o.dim,o.n);
        int completed=0; bool converged=false;
        double compute=0,halo=0,diagnostics=0;
        MPI_Barrier(MPI_COMM_WORLD);
        double start=MPI_Wtime();
        int limit=o.mode=="benchmark"?o.iterations:o.max_iters;
        for(int t=0;t<limit;++t) {
            double tick=MPI_Wtime(); s.exchange(); halo+=MPI_Wtime()-tick;
            tick=MPI_Wtime(); s.step(); compute+=MPI_Wtime()-tick; completed=t+1;
            if(o.mode=="validation" && (completed%o.check_every==0 || completed==limit)) {
                tick=MPI_Wtime(); s.exchange(); halo+=MPI_Wtime()-tick;
                tick=MPI_Wtime(); converged=s.metrics().relative_residual<o.tol; diagnostics+=MPI_Wtime()-tick;
                if(converged) break;
            }
        }
        double local_times[4]={MPI_Wtime()-start,compute,halo,diagnostics}, times[4];
        MPI_Reduce(local_times,times,4,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        // Final diagnostics deliberately lie outside the iteration timer.
        s.exchange(); auto m=s.metrics();
        if(!std::isfinite(m.relative_residual) || !std::isfinite(m.l2_rms_error)) throw std::runtime_error("nonfinite metrics");
        if(!o.field.empty()) {
            std::ofstream field;
            if(rank==0) { field.open(o.field); if(!field) throw std::runtime_error("cannot open field output"); field<<std::setprecision(17); }
            std::vector<double> buffer(s.layer);
            for(int owner=0;owner<ranks;++owner) {
                int slabs=o.n/ranks+(owner<o.n%ranks);
                for(int i=1;i<=slabs;++i) {
                    if(rank==owner && rank!=0) MPI_Send(s.old.data()+i*s.layer,static_cast<int>(s.layer),MPI_DOUBLE,0,20,MPI_COMM_WORLD);
                    if(rank==0) {
                        const double* values;
                        if(owner==0) values=s.old.data()+i*s.layer;
                        else { MPI_Recv(buffer.data(),static_cast<int>(s.layer),MPI_DOUBLE,owner,20,MPI_COMM_WORLD,MPI_STATUS_IGNORE); values=buffer.data(); }
                        for(int j=1;j<=o.n;++j) {
                            if(o.dim==2) field<<values[j]<<'\n';
                            else for(int k=1;k<=o.n;++k) field<<values[j*s.side+k]<<'\n';
                        }
                    }
                }
            }
            if(rank==0) { field.close(); if(!field) throw std::runtime_error("failed writing field"); }
        }
        if(rank==0) {
            std::ifstream existing(o.output,std::ios::binary|std::ios::ate); bool header=!existing || existing.tellg()==0; existing.close();
            std::ofstream out(o.output,std::ios::app); if(!out) throw std::runtime_error("cannot open CSV output");
            const std::string columns="version,dim,n,mpi_ranks,omp_threads,total_cores,iterations,mode,runtime_sec,compute_sec,halo_sec,reduction_sec,relative_residual,l2_rms_error,linf_error,boundary_max,converged";
            auto record=[&](std::ostream& stream) {
                stream<<std::setprecision(17)<<"mpi,"<<o.dim<<','<<o.n<<','<<ranks<<",1,"<<ranks<<','<<completed<<','<<o.mode<<','<<times[0]<<','<<times[1]<<','<<times[2]<<','<<times[3]<<','<<m.relative_residual<<','<<m.l2_rms_error<<','<<m.linf_error<<','<<m.boundary_max<<','<<(o.mode=="validation"?(converged?"true":"false"):"NA")<<'\n';
            };
            if(header) out<<columns<<'\n';
            record(out); out.close(); if(!out) throw std::runtime_error("failed writing CSV");
            std::cout<<columns<<'\n'; record(std::cout);
        }
        int code=o.mode=="validation" && !converged?2:0;
        MPI_Finalize(); return code;
    } catch(const std::exception& e) {
        std::cerr<<"rank "<<rank<<" error: "<<e.what()<<'\n';
        MPI_Abort(MPI_COMM_WORLD,1); return 1;
    }
}
