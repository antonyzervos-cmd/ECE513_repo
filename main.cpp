//g++ main.cpp -I /usr/include/eigen3 -o main
#include "parser.cpp"
#include "circuit_equations.cpp"
#include <iostream>
#include <utility>

int main() {
    auto [head, num_nodes] = parse_netlist("new 1.cir");

    if (head == nullptr) {
        std::cerr << "Error while parsing netlist.\n";
        return 1;
    }

    print_the_list(head);

    auto mna = build_MNA_DC(head, num_nodes);
    const auto& A = mna.first;
    const auto& rhs = mna.second;

    
    // FREE MEMORY
    element* temp;
    while (head != nullptr) {
        temp = head;
        head = head->next;
        delete temp;
    }
    return 0;
}