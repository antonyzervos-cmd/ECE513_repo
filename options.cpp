// options.cpp
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>

struct RunOptions {
    bool use_spd = false;
    bool use_custom = false;
    bool do_dc_sweep = false;

    std::string sweep_source;
    double sweep_start = 0;
    double sweep_end = 0;
    double sweep_step = 0;

    std::vector<std::string> plot_nodes;
};

static inline std::string to_lower_str(const std::string& s) {
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(), ::tolower);
    return t;
}

RunOptions parse_options_from_file(const std::string& filename) {
    RunOptions opts;
    std::ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "parse_options_from_file: cannot open " << filename << "\n";
        return opts;
    }

    std::string line;
    while (std::getline(ifs, line)) {
        std::string trimmed = line;
        trimmed.erase(0, trimmed.find_first_not_of(" \t"));

        if (trimmed.empty()) 
            continue;
        if (trimmed[0] == '*') 
            continue;

        std::string low = to_lower_str(trimmed);

        if (low.rfind(".options", 0) == 0) {
            std::istringstream iss(low);
            std::string dot, tok;
            iss >> dot;
            while (iss >> tok) {
                if (tok == "spd") opts.use_spd = true;
                else if (tok == "custom") opts.use_custom = true;
            }
        }
        else if (low.rfind(".dc", 0) == 0) {
            std::istringstream iss(trimmed);
            std::string dot, var;
            double s, e, step;
            if (iss >> dot >> var >> s >> e >> step) {
                opts.do_dc_sweep = true;
                opts.sweep_source = to_lower_str(var);
                opts.sweep_start = s;
                opts.sweep_end = e;
                opts.sweep_step = step;
            }
        }
        else if (low.rfind(".plot", 0) == 0 || low.rfind(".print", 0) == 0) {
            std::istringstream iss(trimmed);
            std::string cmd, tok;
            iss >> cmd; // skip .plot
            while (iss >> tok) {
                if (!tok.empty() && tok.back() == ',') tok.pop_back();
                opts.plot_nodes.push_back(tok);
            }
        }
    }

    return opts;
}
