#include "Header.h"
#include <iostream>
#include <ctime>

int main() {
    Graph gr(R"(D:\Nikita\Programms\JetBrains\CLionProjects\algorithms\input.txt)");
    std::vector<int> new_path, best_path;
    srand(time(0));
    gr.Greedy(rand() % gr.GetVertexes());
    gr.LocalSearch();

    std::vector<int> path = gr.getPathList();
    for (size_t i = 0; i < path.size(); i++) {
        std::cout << path[i] << " ";
    }
    std::cout << '\n' << gr.GetPathSize();
    return 0;
}