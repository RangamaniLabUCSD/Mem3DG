
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <iostream>
#include <bitset>

using HalfedgeMesh = geometrycentral::surface::HalfedgeMesh;
using PolygonSoupMesh = geometrycentral::PolygonSoupMesh;
using Vector3 = geometrycentral::Vector3;

template<typename T>
::std::ostream& operator<<(::std::ostream& output, const std::vector<T>& v) {
    output << "[";
    for (T item : v) {
        output << item << ", ";
    }
    output << "]";
    return output;
}


// strip unused vertices from face-vertex lists
void stripUnusedVertices(std::vector<Vector3>& positions, std::vector<std::vector<size_t>>& faceIndices) {
    size_t nVert = positions.size();


    // Find any unused vertices
    std::vector<size_t> vertexDegreeCount(nVert, 0);
    size_t nUsedVerts = 0;
    for (auto& f : faceIndices) {
        for (auto& i : f) {
            // Make sure we can safely index positions
            GC_SAFETY_ASSERT(i < positions.size(),
                "face index list has a vertex index which is greater than the number of vertices");
            vertexDegreeCount[i]++;
            if (vertexDegreeCount[i] == 1) {
                nUsedVerts++;
            }
        }
    }


    // Early exit if dense
    if (nUsedVerts == nVert) {
        return;
    }


    // Else: strip unused vertices and re-index faces
    size_t nNewVertices = 0;
    std::vector<size_t> oldToNewVertexInd(nVert);
    for (size_t iV = 0; iV < nVert; iV++) {
        if (vertexDegreeCount[iV] > 0) {
            oldToNewVertexInd[iV] = nNewVertices;
            positions[nNewVertices] = positions[iV];
            nNewVertices++;
        }
    }
    positions.resize(nNewVertices);
    for (auto& f : faceIndices) {
        for (auto& i : f) {
            i = oldToNewVertexInd[i];
        }
    }
}



int main() {
    //std::unique_ptr<HalfedgeMesh> mesh = new HalfedgeMesh();
    // Construct an icosahedron
    // Construct a std::vector of Vector3
    std::vector<Vector3> coords;
    std::vector<std::vector<std::size_t>> polygons;
    // [[1,2,3], [2,3,4],[..] ]
    coords.push_back(Vector3{ 1,0,0 });
    coords.push_back(Vector3{ -1,0,0 });
    coords.push_back(Vector3{ 0,1,0 });
    coords.push_back(Vector3{ 0,0,1 });
    for (Vector3 v : coords)
        std::cout << v << std::endl;
    polygons.push_back(std::vector<std::size_t>{0, 1, 2});
    polygons.push_back(std::vector<std::size_t>{0, 1, 3});
    polygons.push_back(std::vector<std::size_t>{0, 2, 3});
    polygons.push_back(std::vector<std::size_t>{1, 2, 3});
    for (std::vector<std::size_t> f : polygons)
        std::cout << f << std::endl;
    // std::cout << static_cast<std::size_t>(3.14159) << std::endl;
    // int i = std::pow(2,31)-1;
    // std::cout << i << ": " << std::bitset<32>(i) << std::endl;
    // std::cout << i+1 << ": " << std::bitset<32>(i+1) << std::endl;
    // std::size_t t = std::pow(2,31)-1;
    // std::cout << t << ": " << std::bitset<32>(t) << std::endl;
    // std::cout << t+1 << ": " << std::bitset<32>(t+1) << std::endl;
    PolygonSoupMesh soup(polygons, coords);
    soup.mergeIdenticalVertices();
    // stripUnusedVertices(soup.vertexCoordinates, soup.polygons);
    std::unique_ptr<HalfedgeMesh> mesh;
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> vpg;
    std::tie(mesh, vpg) = geometrycentral::surface::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);
    return 0;
}