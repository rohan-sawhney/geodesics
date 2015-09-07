#include "Mesh.h"
#include "MeshIO.h"
#include <Eigen/SparseCholesky>

Mesh::Mesh()
{
    
}

Mesh::Mesh(const Mesh& mesh)
{
    *this = mesh;
}

bool Mesh::read(const std::string& fileName)
{
    std::ifstream in(fileName.c_str());

    if (!in.is_open()) {
        std::cerr << "Error: Could not open file for reading" << std::endl;
        return false;
    }
    
    bool readSuccessful = false;
    if ((readSuccessful = MeshIO::read(in, *this))) {
        normalize();
    }
    
    return readSuccessful;
}

bool Mesh::write(const std::string& fileName) const
{
    std::ofstream out(fileName.c_str());
    
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file for writing" << std::endl;
        return false;
    }
    
    MeshIO::write(out, *this);
    
    return false;
}

void Mesh::buildLaplacian(Eigen::SparseMatrix<double>& L) const
{
    std::vector<Eigen::Triplet<double>> LTriplet;
    
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        
        HalfEdgeCIter he = v->he;
        double sumCoefficients = 0.0;
        do {
            // (cotA + cotB) / 2A
            double coefficient = 0.5 * (he->cotan() + he->flip->cotan());
            sumCoefficients += coefficient;
            
            LTriplet.push_back(Eigen::Triplet<double>(v->index, he->flip->vertex->index, coefficient));
            
            he = he->flip->next;
        } while (he != v->he);
        
        LTriplet.push_back(Eigen::Triplet<double>(v->index, v->index, -sumCoefficients));
    }
    
    L.setFromTriplets(LTriplet.begin(), LTriplet.end());
}

void Mesh::buildAreaMatrix(Eigen::SparseMatrix<double>& A) const
{
    std::vector<Eigen::Triplet<double>> ATriplet;
    
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        ATriplet.push_back(Eigen::Triplet<double>(v->index, v->index, v->dualArea()));
    }
    
    A.setFromTriplets(ATriplet.begin(), ATriplet.end());
}

double Mesh::computeTimeStep() const
{
    // t = avg edge length ^ 2
    double avgLength = 0.0;
    for (EdgeCIter e = edges.begin(); e != edges.end(); e++) {
        avgLength += e->length();
    }
    avgLength /= (double)edges.size();
    
    return avgLength * avgLength;
}

void Mesh::computeFaceGradients(Eigen::MatrixXd& gradients, const Eigen::VectorXd& u) const
{
    for (FaceCIter f = faces.begin(); f != faces.end(); f++) {
        
        Eigen::Vector3d gradient;
        gradient.setZero();
        Eigen::Vector3d normal = f->normal();
        
        HalfEdgeCIter he = f->he;
        do {
            double ui = u(he->next->next->vertex->index);
            Eigen::Vector3d ei = he->next->vertex->position - he->vertex->position;
            
            gradient += ui * normal.cross(ei);
            
            he = he->next;
        } while (he != f->he);
        
        gradient /= (2.0 * f->area());
        gradient.normalize();
        
        gradients.row(f->index) = -gradient;
    }
}

void Mesh::computeIntegratedDivergence(Eigen::VectorXd& integratedDivs,
                                       const Eigen::MatrixXd& gradients) const
{
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        double integratedDiv = 0.0;
        Eigen::Vector3d p = v->position;
        
        HalfEdgeCIter he = v->he;
        do {
            Eigen::Vector3d gradient = gradients.row(he->face->index);
            
            Eigen::Vector3d p1 = he->next->vertex->position;
            Eigen::Vector3d p2 = he->next->next->vertex->position;
            
            Eigen::Vector3d e1 = p1 - p;
            Eigen::Vector3d e2 = p2 - p;
            Eigen::Vector3d ei = p2 - p1;
            
            double theta1 = acos((-e2).dot(-ei) / (e2.norm() * ei.norm()));
            double cot1 = 1.0 / tan(theta1);
            
            double theta2 = acos((-e1).dot(ei) / (e1.norm() * ei.norm()));
            double cot2 = 1.0 / tan(theta2);
            
            integratedDiv += e1.dot(gradient) * cot1 + e2.dot(gradient) * cot2;
            
            he = he->flip->next;
            
        } while (he != v->he);
        
        integratedDivs[v->index] = 0.5 * integratedDiv;
    }
}

void Mesh::computeGeodesics(const int vIdx)
{
    int v = (int)vertices.size();
    
    // build Laplacian
    Eigen::SparseMatrix<double> L(v, v);
    buildLaplacian(L);
    
    // build diagonal area matrix
    Eigen::SparseMatrix<double> A(v, v);
    buildAreaMatrix(A);
    
    // set random point on mesh to 1
    Eigen::VectorXd u(v);
    u.setZero();
    u(vIdx) = 1.0;
    
    // compute time step
    double t = computeTimeStep();
    
    // 1. compute heat flow for time t
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> heatSolver(A - t*L);
    u = heatSolver.solve(u);
    
    // 2. evaluate face gradients
    Eigen::MatrixXd gradients((int)faces.size(), 3);
    computeFaceGradients(gradients, u);
    
    // 3. solve poisson equation
    Eigen::VectorXd integratedDivs(v);
    computeIntegratedDivergence(integratedDivs, gradients);
    
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> poissonSolver(L);
    Eigen::VectorXd phis = poissonSolver.solve(integratedDivs);
    
    // compute max and min phis
    double minPhi = INFINITY;
    double maxPhi = -INFINITY;
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        if (minPhi > phis(v->index)) minPhi = phis(v->index);
        if (maxPhi < phis(v->index)) maxPhi = phis(v->index);
    }
    double range = maxPhi - minPhi;
    
    // set phis
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->phi = 1 - ((phis(v->index) - minPhi) / range);
    }
}

void Mesh::normalize()
{
    // compute center of mass
    Eigen::Vector3d cm = Eigen::Vector3d::Zero();
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        cm += v->position;
    }
    cm /= (double)vertices.size();
    
    // translate to origin
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position -= cm;
    }
    
    // determine radius
    double rMax = 0;
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        rMax = std::max(rMax, v->position.norm());
    }
    
    // rescale to unit sphere
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position /= rMax;
    }
}
