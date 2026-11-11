#pragma once
#include "solver.hpp"
#include <mpi.h>
#include <omp.h>

namespace poisson {
class HybridSolver {
public:
    int dim,n,rank,ranks,local_n,start,left,right;
    std::size_t side,layer;
    double h;
    std::vector<double> old,next,rhs;
    HybridSolver(int dimension,int interior,MPI_Comm communicator=MPI_COMM_WORLD)
        : dim(dimension),n(interior),comm(communicator) {
        MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
        if((dim!=2 && dim!=3) || n<ranks) throw std::invalid_argument("dim must be 2 or 3; n must be >= MPI ranks");
        local_n=n/ranks+(rank<n%ranks);
        start=1+rank*(n/ranks)+std::min(rank,n%ranks);
        left=rank==0?MPI_PROC_NULL:rank-1; right=rank==ranks-1?MPI_PROC_NULL:rank+1;
        side=static_cast<std::size_t>(n)+2; layer=side;
        if(dim==3) {
            if(side>std::numeric_limits<std::size_t>::max()/side) throw std::length_error("layer overflow");
            layer*=side;
        }
        if(layer>static_cast<std::size_t>(std::numeric_limits<int>::max())) throw std::length_error("halo exceeds MPI int count");
        if(layer>std::vector<double>().max_size()/(static_cast<std::size_t>(local_n)+2)) throw std::length_error("local allocation overflow");
        std::size_t count=layer*(static_cast<std::size_t>(local_n)+2);
        old.assign(count,0); next.assign(count,0); rhs.assign(count,0);
        h=1.0/(static_cast<double>(n)+1);
        visit([&](std::size_t p,int i,int j,int k) { rhs[p]=dim*Solver::pi()*Solver::pi()*exact(i,j,k); });
    }
    std::size_t index(int i,int j,int k=0) const { return static_cast<std::size_t>(i)*layer+(dim==2?j:static_cast<std::size_t>(j)*side+k); }
    double exact(int i,int j,int k) const {
        double v=std::sin(Solver::pi()*(start+i-1)*h)*std::sin(Solver::pi()*j*h);
        return dim==3?v*std::sin(Solver::pi()*k*h):v;
    }
    template<class F> void visit(F f) const {
        for(int i=1;i<=local_n;++i) for(int j=1;j<=n;++j) {
            if(dim==2) f(index(i,j),i,j,0);
            else for(int k=1;k<=n;++k) f(index(i,j,k),i,j,k);
        }
    }
    void exchange() {
        MPI_Sendrecv(old.data()+layer,static_cast<int>(layer),MPI_DOUBLE,left,10,
                     old.data()+(local_n+1)*layer,static_cast<int>(layer),MPI_DOUBLE,right,10,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(old.data()+local_n*layer,static_cast<int>(layer),MPI_DOUBLE,right,11,
                     old.data(),static_cast<int>(layer),MPI_DOUBLE,left,11,comm,MPI_STATUS_IGNORE);
    }
    double neighbors(std::size_t p) const {
        if(dim==2) return old[p-layer]+old[p+layer]+old[p-1]+old[p+1];
        return old[p-layer]+old[p+layer]+old[p-side]+old[p+side]+old[p-1]+old[p+1];
    }
    void step() {
        // Independent Jacobi updates; MPI communication stays on the main thread.
        #pragma omp parallel for schedule(static)
        for (int i=1; i<=local_n; ++i) {
            for (int j=1; j<=n; ++j) {
                if (dim==2) {
                    const std::size_t p=index(i,j);
                    next[p]=(neighbors(p)+h*h*rhs[p])/(2.0*dim);
                } else {
                    for (int k=1; k<=n; ++k) {
                        const std::size_t p=index(i,j,k);
                        next[p]=(neighbors(p)+h*h*rhs[p])/(2.0*dim);
                    }
                }
            }
        }
        old.swap(next);
    }
    // Caller must refresh current halos before this operation.
    Metrics metrics() const {
        long double sums[3]={0,0,0},global_sums[3];
        double maxima[2]={0,0},global_maxima[2];
        visit([&](std::size_t p,int i,int j,int k) {
            double r=rhs[p]-(2.0*dim*old[p]-neighbors(p))/(h*h),e=old[p]-exact(i,j,k);
            sums[0]+=static_cast<long double>(r)*r; sums[1]+=static_cast<long double>(rhs[p])*rhs[p];
            sums[2]+=static_cast<long double>(e)*e; maxima[0]=std::max(maxima[0],std::abs(e));
        });
        for(int i=0;i<=local_n+1;++i) for(std::size_t q=0;q<layer;++q) {
            bool physical=(i==0 && rank==0)||(i==local_n+1 && rank==ranks-1);
            physical=physical || (dim==2 ? q==0 || q==side-1 : q/side==0 || q/side==side-1 || q%side==0 || q%side==side-1);
            if(physical) maxima[1]=std::max({maxima[1],std::abs(old[i*layer+q]),std::abs(next[i*layer+q])});
        }
        MPI_Allreduce(sums,global_sums,3,MPI_LONG_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(maxima,global_maxima,2,MPI_DOUBLE,MPI_MAX,comm);
        long double points=1; for(int q=0;q<dim;++q) points*=n;
        return {std::sqrt(static_cast<double>(global_sums[0]/global_sums[1])),std::sqrt(static_cast<double>(global_sums[2]/points)),global_maxima[0],global_maxima[1]};
    }
private:
    MPI_Comm comm;
};
}
