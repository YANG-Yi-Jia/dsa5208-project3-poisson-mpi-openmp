#pragma once

#include "solver.hpp"
#include <filesystem>
#include <iostream>
#include <string>

struct Options {
    int dim = 2, n = 64, iterations = 1000, max_iters = 100000, check_every = 100;
    double tol = 1e-8;
    std::string mode = "benchmark", output = "results.csv", field = "";
};

inline int integer(const std::string& s) {
    std::size_t used = 0;
    long long value = std::stoll(s, &used);
    if (used != s.size() || value < 1 || value > std::numeric_limits<int>::max()) {
        throw std::invalid_argument("expected positive integer: " + s);
    }
    return static_cast<int>(value);
}

inline bool parse_options(int argc, char** argv, Options& o) {
    for (int a = 1; a < argc; ++a) {
        std::string key = argv[a];

        if (key == "--help") {
            std::cout << "--dim 2|3 --n N --mode benchmark|validation "
                         "--iterations K --tol EPS --max-iters K --check-every K "
                         "--output PATH --field PATH\n";
            return false;
        }

        if (a + 1 >= argc) {
            throw std::invalid_argument("missing value for " + key);
        }

        std::string val = argv[++a];
        if (key == "--dim") o.dim = integer(val);
        else if (key == "--n") o.n = integer(val);
        else if (key == "--iterations") o.iterations = integer(val);
        else if (key == "--max-iters") o.max_iters = integer(val);
        else if (key == "--check-every") o.check_every = integer(val);
        else if (key == "--mode") o.mode = val;
        else if (key == "--output") o.output = val;
        else if (key == "--field") o.field = val;
        else if (key == "--tol") {
            std::size_t used = 0;
            o.tol = std::stod(val, &used);
            if (used != val.size() || !std::isfinite(o.tol) || o.tol <= 0) {
                throw std::invalid_argument("invalid tolerance");
            }
        } else {
            throw std::invalid_argument("unknown option: " + key);
        }
    }

    if (o.mode != "benchmark" && o.mode != "validation") {
        throw std::invalid_argument("invalid mode");
    }

    if (!o.field.empty()) {
        namespace fs = std::filesystem;
        const auto output = fs::weakly_canonical(fs::absolute(o.output));
        const auto field = fs::weakly_canonical(fs::absolute(o.field));

        if (output == field ||
            (fs::exists(output) && fs::exists(field) && fs::equivalent(output, field))) {
            throw std::invalid_argument("--output and --field must refer to different files");
        }
    }

    return true;
}
