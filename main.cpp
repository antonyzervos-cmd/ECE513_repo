// main.cpp
// g++ main.cpp -I /usr/include/eigen3 -lgsl -lgslcblas -O2 -std=c++17 -o main

#include <iostream>
#include "parser.cpp"
#include "circuit_equations.cpp"
#include "options.cpp"
#include "solver.cpp"
#include "dc_sweep.cpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "No args!\n";
        return 1;
    }

    std::string filename = argv[1];

    auto [head, num_nodes, node_map] = parse_netlist(filename);

    std::vector<RunOptions> all_opts = parse_all_options(filename);

    print_the_list(head);

    auto [A, rhs, vsrc_index_map, isrc_index_map] = build_MNA_DC(head, num_nodes);

    SolveResult sol = solve_system(head, num_nodes, all_opts[0], A, rhs);
    write_dc_op("dc_op.txt", head, sol);

    for (const RunOptions& opts : all_opts)
    {
        if (!opts.do_dc_sweep)
            continue;

        run_dc_sweep( head,num_nodes, opts,A,rhs, vsrc_index_map, isrc_index_map, node_map);
    }

    // free memory
    while (head) {
        element* t = head;
        head = head->next;
        delete t;
    }

    return 0;
}
