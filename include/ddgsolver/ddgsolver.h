#pragma once
#include <geometrycentral/surface/surface_mesh.h>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

int viewer(std::string fileName);

int genIcosphere(size_t nSub, std::string path, double R);

int driver(std::string inputMesh, std::string refMesh, double Kb, double H0,
						double Kse, double Ksl, double Ksg,
						double Kv, double Vt, double gamma,
						double kt, size_t ptInd, double extF,
						double conc, double h, double T,
						double eps, double closeZone, double increment, 
						double tSave, std::string outputDir);