// options.cpp
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>


struct CommandDC {
    std::string source_name; // onoma phghs sarwshs
    double start = 0;
    double end = 0;
    double step = 0;
    std::vector<std::string> plot_nodes;  // nodes to plot 
};

struct RunOptions {
    bool use_spd = false;
    bool use_custom = false;
    bool do_dc_sweep = false;

    bool use_sparse = false;
    bool use_iter = false;  // iterative
    double itol = 1e-3;

    std::vector<CommandDC> dc_commands; 
};

static inline std::string to_lower_str(const std::string& s) {
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(), ::tolower);
    return t;
}


extern void to_lowercase(std::string& str); 

RunOptions parse_options_from_file(const std::string& filename) {
    std::ifstream file(filename);
    RunOptions opts;

    opts.use_iter = false;
    opts.use_spd = false;
    opts.use_sparse = false;
    opts.itol = 1e-3;
    opts.do_dc_sweep = false;

    std::string line;
    while (std::getline(file, line))
    {
        std::string low = line;
        to_lowercase(low);

        if (low.rfind(".options", 0) == 0) {
            std::istringstream iss(low);
            std::string token;
            iss >> token; // .options

            while (iss >> token) {
                if (token == "iter") 
                    opts.use_iter = true;
                else if (token == "spd") {
                    opts.use_spd = true;
                }
                else if (token == "sparse") 
                    opts.use_sparse = true;
                else if (token.find("itol=") == 0)
                    opts.itol = std::stod(token.substr(5));
            }
        }
        else if (low.rfind(".print", 0) == 0 || low.rfind(".plot", 0) == 0)
        {
            std::istringstream iss(low);
            std::string token;
            iss >> token; // .print
            if (!opts.dc_commands.empty()) {
                while (iss >> token)
                    opts.dc_commands.back().plot_nodes.push_back(token);
            }
        }
        else if (low.rfind(".dc", 0) == 0) {
            // NEW DC
            CommandDC new_cmd;

            std::istringstream iss(low);
            std::string token;
            iss >> token;          // .dc
            iss >> new_cmd.source_name; // sto new struct
            iss >> new_cmd.start;
            iss >> new_cmd.end;
            iss >> new_cmd.step;

            opts.dc_commands.push_back(new_cmd); // sto vector

            opts.do_dc_sweep = true;
        }
    }

    return opts;
}