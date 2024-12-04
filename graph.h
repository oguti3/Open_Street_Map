#pragma once

#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

using namespace std;

/// @brief Simple directed graph using an adjacency list.
/// @tparam VertexT vertex type
/// @tparam WeightT edge weight type
template <typename VertexT, typename WeightT>
class graph {
   private:
        map<VertexT, map<VertexT, WeightT>> vertices; // every vertex has a map of their adjacent vertex and weight of edge between them
        size_t numEdges;
        
   public:
    /// Default constructor
    graph() {
        numEdges = 0;
    }

    /// @brief Add the vertex `v` to the graph, must run in at most O(log |V|).
    /// @param v
    /// @return true if successfully added; false if it existed already
    bool addVertex(VertexT v) {
        if((vertices.find(v)) != vertices.end()) {
            return false;
        }
        vertices[v];
        return true;
    }

    /// @brief Add or overwrite directed edge in the graph, must run in at most O(log |V|).
    /// @param from starting vertex
    /// @param to ending vertex
    /// @param weight edge weight / label
    /// @return true if successfully added or overwritten;
    ///         false if either vertices isn't in graph
    bool addEdge(VertexT from, VertexT to, WeightT weight) {
        if((vertices.find(from) == vertices.end()) || (vertices.find(to) == vertices.end())) {
            return false;
        }
        if(vertices[from].find(to) == vertices[from].end()) { // new directed edge
            ++numEdges;
        }
        vertices[from][to] = weight;
        return true;
    }

    /// @brief Maybe get the weight associated with a given edge, must run in at most O(log |V|).
    /// @param from starting vertex
    /// @param to ending vertex
    /// @param weight output parameter
    /// @return true if the edge exists, and `weight` is set;
    ///         false if the edge does not exist
    bool getWeight(VertexT from, VertexT to, WeightT& weight) const {   
        if((vertices.find(from) == vertices.end()) || (vertices.at(from).find(to) == vertices.at(from).end())) {
            return false;
        }
        weight = vertices.at(from).at(to);
        return true;
    }

    /// @brief Get the out-neighbors of `v`. Must run in at most O(|V|).
    /// @param v
    /// @return vertices that v has an edge to
    set<VertexT> neighbors(VertexT v) const {
        set<VertexT> S;
        for(const auto& outNeighbor : vertices.at(v)) {
            S.emplace(outNeighbor.first);
        }
        return S;
    }

    /// @brief Return a vector containing all vertices in the graph
    vector<VertexT> getVertices() const {
        vector<VertexT> vertices;
        for(const auto& v : this->vertices) {
            vertices.push_back(v.first);
        }
        return vertices;
    }

    /// @brief Get the number of vertices in the graph. Runs in O(1).
    size_t NumVertices() const {
        return vertices.size();
    }

    /// @brief Get the number of directed edges in the graph. Runs in at most O(|V|).
    size_t NumEdges() const {
        return numEdges;
    }
};
