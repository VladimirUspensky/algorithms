#ifndef HEADER_H
#define HEADER_H

#include <vector>
#include <fstream>

double dist(int x1, int y1, int x2, int y2);
int min(std::vector<double>, std::vector<bool>);

class Graph {
public:
    explicit Graph(const std::string &);

    ~Graph();

    void Greedy(int start);

    void LocalSearch();

    bool AllVisited();

    double CalcGain(int a, int b, int c, int d);

    void _2OPT();

    double GetPathSize() const;

    int GetVertexes() const;

    std::vector<int> getPathList();

    void switch_edges(int i, int j);

private:
    double pathSize = 0;
    double bestSize = 0;
    int numOfVertices = 0;
    int edges = 0;
    std::vector<bool> visited;
    std::vector<std::vector<double>> matrix;
    std::vector<std::pair<int, std::pair<int, int>>> vertices;
    std::vector<int> path;
    std::vector<int> bestPath;
};

#endif
