#ifndef MESH_H
#define MESH_H

#include "Types.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "HalfEdge.h"
#include <Eigen/SparseCore>

class Mesh {
public:
    // default constructor
    Mesh();
    
    // copy constructor
    Mesh(const Mesh& mesh);
        
    // read mesh from file
    bool read(const std::string& fileName);
    
    // write mesh to file
    bool write(const std::string& fileName) const;
    
    // computes geodesic distances from source point
    void computeGeodesics(const int vIdx);
        
    // member variables
    std::vector<HalfEdge> halfEdges;
    std::vector<Vertex> vertices;
    std::vector<Eigen::Vector3d> uvs;
    std::vector<Eigen::Vector3d> normals;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::vector<HalfEdgeIter> boundaries;

private:
    // builds Laplacian operator
    void buildLaplacian(Eigen::SparseMatrix<double>& L) const;
    
    // builds area matrix
    void buildAreaMatrix(Eigen::SparseMatrix<double>& A) const;
    
    // computes time step
    double computeTimeStep() const;
    
    // computes gradient for faces
    void computeFaceGradients(Eigen::MatrixXd& gradients, const Eigen::VectorXd& u) const;
    
    // computes integrated gradient for vertices
    void computeIntegratedDivergence(Eigen::VectorXd& integratedDivs,
                                     const Eigen::MatrixXd& gradients) const;
    
    // center mesh about origin and rescale to unit radius
    void normalize();
};

#endif