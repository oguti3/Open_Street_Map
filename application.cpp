#include "application.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip> /*setprecision*/
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "dist.h"
#include "graph.h"
#include "osm.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

double INF = numeric_limits<double>::max();

class prioritize {
   public:
    bool operator()(const pair<long long, double>& p1,
                    const pair<long long, double>& p2) const {
        return p1.second > p2.second;
    }
};

/// @brief Builds a graph using given nodes, footways, and buildings.
/// @param Nodes, Footways, Buildings
/// @return A graph with all the vertexes and edges given the nodes, footways, and buildings.
graph<long long, double> buildGraph(
    const map<long long, Coordinates>& Nodes,
    const vector<FootwayInfo>& Footways,
    const vector<BuildingInfo>& Buildings) {
    graph<long long, double> G;  
    for(const auto& node : Nodes) {
        G.addVertex(node.first);
    }
    for(const auto& building : Buildings) {
        G.addVertex(building.Coords.ID); // add current building to graph
        for(const auto& node : Nodes) {
            double dist = distBetween2Points(building.Coords.Lat, building.Coords.Lon, node.second.Lat, node.second.Lon);
            if((dist <= 0.041) && node.second.OnFootway) { // dist between current building and node is within 0.041 and node is on footway
                G.addEdge(building.Coords.ID, node.first, dist); // add directed edge from current building to current node
                G.addEdge(node.first, building.Coords.ID, dist); // add directed edge from current node to current building
            }
        }
    }
    for(const auto& footway : Footways) {
        for(size_t i = 0; i < (footway.Nodes.size()-1); ++i) {
            double dist = distBetween2Points(
                            Nodes.at(footway.Nodes.at(i)).Lat, 
                            Nodes.at(footway.Nodes.at(i)).Lon, 
                            Nodes.at(footway.Nodes.at(i+1)).Lat, 
                            Nodes.at(footway.Nodes.at(i+1)).Lon); 
            G.addEdge(footway.Nodes.at(i), footway.Nodes.at(i+1), dist); // add directed edge from footways current node to it's next node
            G.addEdge(footway.Nodes.at(i+1), footway.Nodes.at(i), dist); // add directed edge from footways next node to it's current node
        }
    }
    return G;
}

vector<long long> dijkstra(
    const graph<long long, double>& G,
    long long start,
    long long target,
    const set<long long>& ignoreNodes) {
    vector<long long> path;

    map<long long, pair<long long, double>> vertices; // key vertex maps to a pair of predV(first) and dist(second)
    priority_queue<pair<long long, double>, 
               vector<pair<long long, double>>, 
               prioritize>
    worklist;
    pair<long long, double> currV;
    double distToAdd;
    for(const auto& v : G.getVertices()) { // sets vertices to default with 0 predV and INF dist
        vertices[v] = make_pair(0, INF);
    }
    vertices.at(start) = make_pair(0, 0); // initializes start's predV and dist to 0
    worklist.emplace(make_pair(start, 0)); // adds start and current dist as 0
    while(!worklist.empty()) {
        currV = worklist.top(); // saves top of worklist before popping
        worklist.pop();
        for(const auto& neighbor : G.neighbors(currV.first)) { // loops through currV's neighbors
            if(((neighbor == target) || (ignoreNodes.find(neighbor) == ignoreNodes.end())) && // target is the neighbor or not in ignoreNodes
                (G.getWeight(currV.first, neighbor, distToAdd) && // an edge exists between the currV and neighbor
                ((currV.second + distToAdd) < vertices.at(neighbor).second))) { // currV's current dist with the added dist is less than neighbor's current dist
                vertices.at(neighbor) = make_pair(currV.first, currV.second + distToAdd); // neighbor's predV is currV and dist is currV's current dist with added dist
                worklist.emplace(make_pair(neighbor, vertices.at(neighbor).second)); // add neighbor and it's distance to worklist
            }
        }  
    }
    if((start == target) || (vertices.at(target).first != 0)) { // there is a path from start to target
        long long curr = target;
        while(curr != start) { // loops through predecessors until reaching start
            path.push_back(curr);
            curr = vertices.at(curr).first;
        }
        path.push_back(curr); // pushes start to the end of path
        for(size_t i = 0; i < (path.size()/2); ++i) { // reverses the path to go from start to target
            long long temp = path.at(i);
            path.at(i) = path.at(path.size()-i-1);
            path.at(path.size()-i-1) = temp;
        }
    }
    return path;
}

double pathLength(const graph<long long, double>& G, const vector<long long>& path) {
    double length = 0.0;
    double weight;
    for (size_t i = 0; i + 1 < path.size(); i++) {
        bool res = G.getWeight(path.at(i), path.at(i + 1), weight);
        assert(res);
        length += weight;
    }
    return length;
}

void outputPath(const vector<long long>& path) {
    for (size_t i = 0; i < path.size(); i++) {
        cout << path.at(i);
        if (i != path.size() - 1) {
            cout << "->";
        }
    }
    cout << endl;
}

void application(
    const vector<BuildingInfo>& Buildings,
    const graph<long long, double>& G) {
    string person1Building, person2Building;

    set<long long> buildingNodes;
    for (const auto& building : Buildings) {
        buildingNodes.insert(building.Coords.ID);
    }

    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);

    while (person1Building != "#") {
        cout << "Enter person 2's building (partial name or abbreviation)> ";
        getline(cin, person2Building);

        //
        // find the building coordinates
        //
        bool foundP1 = false;
        bool foundP2 = false;
        Coordinates P1Coords, P2Coords;
        string P1Name, P2Name;

        for (const BuildingInfo& building : Buildings) {
            if (building.Abbrev == person1Building) {
                foundP1 = true;
                P1Name = building.Fullname;
                P1Coords = building.Coords;
            }
            if (building.Abbrev == person2Building) {
                foundP2 = true;
                P2Name = building.Fullname;
                P2Coords = building.Coords;
            }
        }

        for (const BuildingInfo& building : Buildings) {
            if (!foundP1 &&
                building.Fullname.find(person1Building) != string::npos) {
                foundP1 = true;
                P1Name = building.Fullname;
                P1Coords = building.Coords;
            }
            if (!foundP2 && building.Fullname.find(person2Building) != string::npos) {
                foundP2 = true;
                P2Name = building.Fullname;
                P2Coords = building.Coords;
            }
        }

        if (!foundP1) {
            cout << "Person 1's building not found" << endl;
        } else if (!foundP2) {
            cout << "Person 2's building not found" << endl;
        } else {
            cout << endl;
            cout << "Person 1's point:" << endl;
            cout << " " << P1Name << endl;
            cout << " (" << P1Coords.Lat << ", " << P1Coords.Lon << ")" << endl;
            cout << "Person 2's point:" << endl;
            cout << " " << P2Name << endl;
            cout << " (" << P2Coords.Lat << ", " << P2Coords.Lon << ")" << endl;

            string destName;
            Coordinates destCoords;

            Coordinates centerCoords = centerBetween2Points(
                P1Coords.Lat, P1Coords.Lon, P2Coords.Lat, P2Coords.Lon);

            double minDestDist = numeric_limits<double>::max();

            for (const BuildingInfo& building : Buildings) {
                double dist = distBetween2Points(
                    centerCoords.Lat, centerCoords.Lon,
                    building.Coords.Lat, building.Coords.Lon);
                if (dist < minDestDist) {
                    minDestDist = dist;
                    destCoords = building.Coords;
                    destName = building.Fullname;
                }
            }

            cout << "Destination Building:" << endl;
            cout << " " << destName << endl;
            cout << " (" << destCoords.Lat << ", " << destCoords.Lon << ")" << endl;

            vector<long long> P1Path = dijkstra(G, P1Coords.ID, destCoords.ID, buildingNodes);
            vector<long long> P2Path = dijkstra(G, P2Coords.ID, destCoords.ID, buildingNodes);

            // This should NEVER happen with how the graph is built
            if (P1Path.empty() || P2Path.empty()) {
                cout << endl;
                cout << "At least one person was unable to reach the destination building. Is an edge missing?" << endl;
                cout << endl;
            } else {
                cout << endl;
                cout << "Person 1's distance to dest: " << pathLength(G, P1Path);
                cout << " miles" << endl;
                cout << "Path: ";
                outputPath(P1Path);
                cout << endl;
                cout << "Person 2's distance to dest: " << pathLength(G, P2Path);
                cout << " miles" << endl;
                cout << "Path: ";
                outputPath(P2Path);
            }
        }

        //
        // another navigation?
        //
        cout << endl;
        cout << "Enter person 1's building (partial name or abbreviation), or #> ";
        getline(cin, person1Building);
    }
}
