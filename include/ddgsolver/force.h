#pragma once

#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

using HalfedgeMesh = geometrycentral::surface::HalfedgeMesh;

Eigen::Matrix<double, Eigen::Dynamic, 3> bending_force(std::unique_ptr<HalfedgeMesh>& mesh,
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry>& vpg, double Kb, double H0);

/*
Eigen::Matrix<double, Eigen::Dynamic, 1> stretching_force(std::unique_ptr<HalfedgeMesh> mesh,
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> vpg);

Eigen::Matrix<double, Eigen::Dynamic, 1> pressure_force(std::unique_ptr<HalfedgeMesh> mesh,
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> vpg);
*/
