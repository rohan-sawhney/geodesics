#ifndef VERTEX_H
#define VERTEX_H

#include "Types.h"

class Vertex {
public:
    // outgoing halfedge
    HalfEdgeIter he;
    
    // location in 3d
    Eigen::Vector3d position;
    
    // geodesic distance from source point
    double phi;
    
    // id between 0 and |V|-1
    int index;
           
    // checks if vertex is contained in any edge or face
    bool isIsolated() const;
    
    // returns area of barycentric dual cell associated with the vertex
    double dualArea() const;
    
    bool onBoundary() const;
};

#endif
