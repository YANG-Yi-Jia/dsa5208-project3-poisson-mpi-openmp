#include "mpi_solver.hpp"
#include <iostream>
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    double maximum=0;
    for(int dim:{2,3}) for(int n:{ranks,2*ranks+1}) {
        poisson::Solver serial(dim,n);
        poisson::MPISolver parallel(dim,n);
        auto value=[](int i,int j,int k) { return 0.1*i*i+3*j+7*k*k+0.25*i*j; };
        serial.visit([&](std::size_t p,int i,int j,int k){serial.old[p]=value(i,j,k);});
        parallel.visit([&](std::size_t p,int i,int j,int k){parallel.old[p]=value(parallel.start+i-1,j,k);});
        for(int iteration=0;iteration<5;++iteration) {
            serial.step(); parallel.exchange(); parallel.step();
            parallel.visit([&](std::size_t p,int i,int j,int k){
                maximum=std::max(maximum,std::abs(parallel.old[p]-serial.old[serial.index(parallel.start+i-1,j,k)]));
            });
        }
        parallel.exchange(); auto m=parallel.metrics();
        maximum=std::max(maximum,m.boundary_max);
    }
    double global=0;MPI_Allreduce(&maximum,&global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(rank==0) std::cout<<"asymmetric field max difference: "<<global<<'\n';
    MPI_Finalize();return global<1e-12?0:1;
}
