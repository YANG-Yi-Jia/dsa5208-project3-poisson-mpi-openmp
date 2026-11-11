#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace poisson {
struct Metrics { double relative_residual, l2_rms_error, linf_error, boundary_max; };
class Solver {
public:
    int dim, n;
    std::size_t side, plane, count;
    double h;
    std::vector<double> old, next, rhs;
    Solver(int dimension, int interior) : dim(dimension), n(interior), side(0), plane(0), count(0), h(0) {
        if ((dim != 2 && dim != 3) || n < 1) throw std::invalid_argument("dim must be 2 or 3; n must be positive");
        side = static_cast<std::size_t>(n) + 2;
        count = 1;
        for (int q = 0; q < dim; ++q) {
            if (count > std::vector<double>().max_size() / side) throw std::length_error("grid size overflow");
            count *= side;
        }
        plane = side * side;
        h = 1.0 / (static_cast<double>(n) + 1);
        old.assign(count, 0); next.assign(count, 0); rhs.assign(count, 0);
        visit([&](std::size_t p, int i, int j, int k) { rhs[p] = dim * pi() * pi() * exact(i,j,k); });
    }
    static double pi() { return std::acos(-1.0); }
    std::size_t index(int i, int j, int k = 0) const {
        return dim == 2 ? static_cast<std::size_t>(i)*side+j : (static_cast<std::size_t>(i)*side+j)*side+k;
    }
    double exact(int i, int j, int k) const {
        double v = std::sin(pi()*i*h)*std::sin(pi()*j*h);
        return dim == 3 ? v*std::sin(pi()*k*h) : v;
    }
    template<class F> void visit(F f) const {
        for (int i=1; i<=n; ++i) for (int j=1; j<=n; ++j) {
            if (dim==2) f(index(i,j),i,j,0);
            else for (int k=1; k<=n; ++k) f(index(i,j,k),i,j,k);
        }
    }
    double neighbors(std::size_t p) const {
        if (dim==2) return old[p-side]+old[p+side]+old[p-1]+old[p+1];
        return old[p-plane]+old[p+plane]+old[p-side]+old[p+side]+old[p-1]+old[p+1];
    }
    void step() {
        const double h2=h*h, denom=2.0*dim;
        visit([&](std::size_t p,int,int,int) { next[p]=(neighbors(p)+h2*rhs[p])/denom; });
        old.swap(next);
    }
    Metrics metrics() const {
        long double residual2=0, rhs2=0, error2=0;
        double max_error=0, boundary=0;
        visit([&](std::size_t p,int i,int j,int k) {
            double r=rhs[p]-(2.0*dim*old[p]-neighbors(p))/(h*h);
            double e=old[p]-exact(i,j,k);
            residual2+=static_cast<long double>(r)*r; rhs2+=static_cast<long double>(rhs[p])*rhs[p];
            error2+=static_cast<long double>(e)*e; max_error=std::max(max_error,std::abs(e));
        });
        for (std::size_t p=0; p<count; ++p) {
            std::size_t a=p;
            bool edge=false;
            for (int q=0; q<dim; ++q) { auto c=a%side; a/=side; edge=edge || c==0 || c==side-1; }
            if(edge) boundary=std::max({boundary,std::abs(old[p]),std::abs(next[p])});
        }
        long double points=1; for(int q=0;q<dim;++q) points*=n;
        return {std::sqrt(static_cast<double>(residual2/rhs2)),std::sqrt(static_cast<double>(error2/points)),max_error,boundary};
    }
};
}
