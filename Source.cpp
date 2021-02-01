#include "Header.h"
#include <ctime>
#include <cmath>
#include <string>

double dist(int x1, int y1, int x2, int y2) {
    double dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    return dist;
}

Graph::Graph(const std::string& filename) {
    int ID, x, y;
    std::ifstream fin(filename);

    fin >> this->numOfVertices;

    this->matrix.resize(this->numOfVertices);
    for (auto & i : matrix) {
        i.resize(this->numOfVertices);
    }
    this->visited.resize(this->numOfVertices);
    this->path.resize(this->numOfVertices);

    this->edges = (this->numOfVertices * (this->numOfVertices - 1)) / 2;

    for (size_t i = 0; i < this->numOfVertices; i++) {
        fin >> ID;
        fin >> x >> y;
        vertices.emplace_back(ID, std::make_pair(x, y));
    }

    fin.close();

    for (size_t i = 0; i < this->numOfVertices; i++) {
        for (size_t j = 0; j < this->numOfVertices; j++) {
            if (i != j)
                matrix[i][j] = dist(vertices[i].second.first, vertices[i].second.second,
                                    vertices[j].second.first, vertices[j].second.second);
        }
    }
}

Graph::~Graph() = default;

int min(std::vector<double> edges, std::vector<bool> visited) {
    double min = 999999.0;
    int min_ind = -1;
    for (size_t i = 0; i < edges.size(); i++) {
        if (edges[i] < min && edges[i] != 0 && visited[i] == 0) {
            min = edges[i];
            min_ind = i;
        }
    }
    return min_ind;
}

bool Graph::AllVisited() {
    for (int i = 0; i < this->visited.size(); i++) {
        if (!visited[i]) return false;
    }
    return true;
}

void Graph::Greedy(int start) {
    int i = start;
    int newInd;

    this->pathSize = 0;
    this->path.erase(this->path.begin(), this->path.end());

    this->visited[start] = 1;
    this->path.push_back(this->vertices[i].first);

    while (min(this->matrix[i], this->visited) != -1) {
        newInd = min(this->matrix[i], this->visited);
        this->path.push_back(this->vertices[newInd].first);
        this->pathSize += this->matrix[i][newInd];
        i = newInd;
        this->visited[i] = 1;
    }

    if (this->AllVisited()) {
        this->pathSize += this->matrix[i][start];
    }

    for (size_t j = 0; j < visited.size(); j++) {
        this->visited[j] = 0;
    }
}

std::vector<std::string> twoOpt(std::vector<std::string> route, int i, int k) {
    int start(i), end(k);
    std::vector<std::string> new_route(route.size());
    for (size_t count = 0; count < i; count++) {
        new_route[count] = route[count];
    }
    for (size_t count = i; count < k + 1; count++) {
        new_route[start++] = route[end--];
    }
    for (size_t count = k + 1; count < route.size(); count++) {
        new_route[count] = route[count];
    }
    return new_route;
}

void Graph::switch_edges(int i, int j) {
    int tmp = this->path[j];
    this->path[j] = this->path[i + 1];
    this->path[i + 1] = tmp;
    i += 2;
    j--;
    while (i < j) {
        tmp = this->path[j];
        this->path[j] = this->path[i];
        this->path[i] = tmp;
        i++;
        j--;
    }
}

double Graph::CalcGain(int a, int b, int c, int d) {
    double currCost =
            matrix[this->path[a] - 1][this->path[b] - 1] +
            matrix[this->path[c] - 1][this->path[d % this->path.size()] - 1];
    double newCost =
            matrix[this->path[a] - 1][this->path[c] - 1] +
            matrix[this->path[b] - 1][this->path[d % this->path.size()] - 1];
    return currCost - newCost;
}

void Graph::LocalSearch() {
    this->bestSize = this->pathSize;
    _2OPT();
    while (this->bestSize != this->pathSize) {
        if (this->bestSize > this->pathSize) {
            this->bestSize = this->pathSize;
            this->bestPath = this->path;
        }
        _2OPT();
    }
}

double Graph::GetPathSize() const {
    return this->pathSize;
}

std::vector<int> Graph::getPathList() {
    return this->path;
}

int Graph::GetVertexes() const {
    return this->numOfVertices;
}

void Graph::_2OPT() {
    struct best_move {
        int i, j;
        double gain;
    };
    bool locallyOptimal = false;
    int i, j;
    int x1, x2, y1, y2;
    int counter2Limit;
    double gainExpected;
    best_move move{0};

    while (!locallyOptimal) {
        locallyOptimal = true;

        for (int counter_1 = 0; counter_1 < this->numOfVertices - 2; counter_1++) {
            i = counter_1;
            x1 = this->path[i];
            x2 = this->path[(i + 1) % this->numOfVertices];

            if (i == 0)
                counter2Limit = this->numOfVertices - 2;
            else
                counter2Limit = this->numOfVertices - 1;

            for (int counter_2 = (i + 2); counter_2 <= counter2Limit; counter_2++) {
                j = counter_2;
                y1 = this->path[j];
                y2 = this->path[(j + 1) % this->numOfVertices];

                gainExpected = CalcGain(i, i + 1, j, j + 1);
                if (gainExpected > move.gain) {
                    move = {i, j, gainExpected};
                    locallyOptimal = false;
                }
            }
        }
        if (!locallyOptimal) {
            switch_edges(move.i, move.j);
            this->pathSize -= move.gain;
        }
    }
}
