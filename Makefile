CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude
.PHONY: all test clean
all: poisson_serial
poisson_serial: src/serial.cpp include/solver.hpp include/options.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/serial.cpp -o $@
test: poisson_serial
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) tests/test_stencil.cpp -o /tmp/poisson_test_stencil
	/tmp/poisson_test_stencil
	python3 tests/test_serial.py
clean:
	rm -f poisson_serial poisson_mpi poisson_hybrid test_mpi_stencil

MPICXX ?= mpicxx
.PHONY: mpi test-mpi
mpi: poisson_mpi
poisson_mpi: src/mpi.cpp include/mpi_solver.hpp include/options.hpp include/solver.hpp
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -DOMPI_SKIP_MPICXX src/mpi.cpp -o $@
test-mpi: poisson_serial poisson_mpi
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -DOMPI_SKIP_MPICXX tests/test_mpi_stencil.cpp -o test_mpi_stencil
	python3 tests/test_mpi.py

.PHONY: hybrid test-hybrid
hybrid: poisson_hybrid
poisson_hybrid: src/hybrid.cpp include/hybrid_solver.hpp include/options.hpp include/solver.hpp
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -fopenmp -DOMPI_SKIP_MPICXX src/hybrid.cpp -o $@
test-hybrid: poisson_serial poisson_hybrid
	python3 tests/test_hybrid.py
