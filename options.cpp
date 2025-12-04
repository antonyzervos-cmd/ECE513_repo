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

    bool use_iter = false;     // iterative
    double itol = 1e-3;

    std::string sweep_source; // onoma phghs sarwshs
    double sweep_start = 0;
    double sweep_end = 0;
    double sweep_step = 0;

    std::vector<std::string> plot_nodes;  // nodes to plot 
};

static inline std::string to_lower_str(const std::string& s) {
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(), ::tolower);
    return t;
}

std::vector<RunOptions> parse_all_options(const std::string& filename)
{
    std::ifstream file(filename);
    std::vector<RunOptions> all;

    RunOptions base; 
    base.do_dc_sweep = false;
    // default values
    base.use_iter = false;
    base.use_spd = false;
    base.itol = 1e-3;

    std::string line;
    while (std::getline(file, line))
    {
        std::string low = line;
        to_lowercase(low);
        if (low.rfind(".options", 0) == 0) {
            std::istringstream iss(low);
            std::string token;
            iss >> token; // skip .options

            while (iss >> token) {
                if (token == "iter") {
                    base.use_iter = true;
                }
                else if (token == "spd") {
                    base.use_spd = true;
                }
                else if (token.find("itol=") == 0) {
                    // Parse itol value e.g., itol=1e-4
                    std::string val_str = token.substr(5); // skip "itol="
                    try {
                        base.itol = std::stod(val_str);
                    } catch (...) {
                        std::cerr << "Warning: Invalid ITOL format, using default 1e-3\n";
                    }
                }
            }
            continue;
        }
        // -----------------------
        // .print / .plot
        // -----------------------
        if (low.rfind(".print", 0) == 0 || low.rfind(".plot", 0) == 0)
        {
            std::istringstream iss(low);
            std::string cmd;
            iss >> cmd;         // .print
            std::string token;

            while (iss >> token)
                base.plot_nodes.push_back(token);

            continue;
        }

        // -----------------------
        // .dc
        // -----------------------
        if (low.rfind(".dc", 0) == 0)
        {
            RunOptions dc = base;  

            std::istringstream iss(low);
            std::string cmd;
            iss >> cmd;  // .dc

            iss >> dc.sweep_source;
            iss >> dc.sweep_start;
            iss >> dc.sweep_end;
            iss >> dc.sweep_step;

            dc.do_dc_sweep = true;

            all.push_back(dc);
            continue;
        }
    }

    if (all.empty())
        all.push_back(base);

    return all;
}