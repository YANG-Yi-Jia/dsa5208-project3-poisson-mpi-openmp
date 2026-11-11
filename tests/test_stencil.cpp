#include "solver.hpp"
#include <cassert>
int main() {
    for(int dim : {2,3}) {
        poisson::Solver s(dim,4);
        s.visit([&](std::size_t p,int i,int j,int k) { s.old[p]=i*i+3*j+7*k*k+0.1*i*j; });
        auto before=s.old;
        s.step();
        s.visit([&](std::size_t p,int i,int j,int k) {
            double sum=before[s.index(i-1,j,k)]+before[s.index(i+1,j,k)]
                      +before[s.index(i,j-1,k)]+before[s.index(i,j+1,k)];
            if(dim==3) sum+=before[s.index(i,j,k-1)]+before[s.index(i,j,k+1)];
            double expected=(sum+s.h*s.h*s.rhs[p])/(2*dim);
            assert(std::abs(s.old[p]-expected)<1e-12);
        });
        assert(s.metrics().boundary_max==0);
    }
}
