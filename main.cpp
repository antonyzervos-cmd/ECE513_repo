//g++ main.cpp -I /usr/include/eigen3 -o main
#include "parser.cpp"
#include "circuit_equations.cpp"

int main() {
    auto [head, num_nodes] = parse_netlist("part1_simple.cir");

    if (head == nullptr) {
        std::cerr << "Error while parsing netlist.\n";
        return 1;
    }

    print_the_list(head);

    build_MNA_DC(head, num_nodes);


    // FREE MEMORY
    element* temp;
    while (head != nullptr) {
        temp = head;
        head = head->next;
        delete temp;
    }
    return 0;
}