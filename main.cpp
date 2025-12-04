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

    // Parse netlist
    auto [head, num_nodes, node_map] = parse_netlist(filename);

    // Parse .options, .dc, .plot
    RunOptions opts = parse_options_from_file(filename);

    // Print circuit
    print_the_list(head);

    // Build MNA equations
    auto [A, rhs, vsrc_index_map, isrc_index_map] = build_MNA_DC(head, num_nodes);

    // Solve DC operating point
    SolveResult sol = solve_system(head, num_nodes, opts, A, rhs);

    write_dc_op("dc_op.txt", head, sol);

    // DC sweep 
    if (opts.do_dc_sweep) {
        VectorXd rhs_backup = rhs; 
        for (const auto& cmd : opts.dc_commands) {
            rhs = rhs_backup;
            run_dc_sweep(head,num_nodes,cmd,opts,A,rhs,vsrc_index_map, isrc_index_map, node_map);
        }
    }

    // Free memory
    while (head) {
        element* t = head;
        head = head->next;
        delete t;
    }

    return 0;
}