import flags
import numpy
import scipy.sparse
from _typeshed import Incomplete
from typing import Callable, overload

__short_version__: str
__version__: str

class Adsorption:
    epsilon: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class Aggregation:
    chi: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class Bending:
    D: float
    H0c: float
    Kb: float
    Kbc: float
    Kd: float
    Kdc: float
    Kg: float
    Kgc: float
    alpha: float
    dA0: float
    relation: str
    def __init__(self) -> None:
        """__init__(self: pymem3dg._core.Bending) -> None"""

class Boundary:
    proteinBoundaryCondition: str
    shapeBoundaryCondition: str
    def __init__(self) -> None:
        """__init__(self: pymem3dg._core.Boundary) -> None"""

class ConjugateGradient:
    c1: float
    constraintTolerance: float
    ifAdaptiveStep: bool
    ifJustGeometryPly: bool
    ifOutputMeshFile: bool
    ifOutputTrajFile: bool
    ifPrintToConsole: bool
    isAugementedLagrangian: bool
    isBacktrack: bool
    outputDirectory: str
    restartPeriod: int
    rho: float
    trajFileName: str
    def __init__(
        self,
        system,
        characteristicTimeStep: float,
        totalTime: float,
        savePeriod: float,
        tolerance: float,
        outputDirectory: str,
        frame: int = ...,
    ) -> None:
        """__init__(self: pymem3dg._core.ConjugateGradient, system: mem3dg::solver::System, characteristicTimeStep: float, totalTime: float, savePeriod: float, tolerance: float, outputDirectory: str, frame: int = 0) -> None


        Conjugate Gradient optimizer constructor

        """
    def integrate(self) -> bool:
        """integrate(self: pymem3dg._core.ConjugateGradient) -> bool


        integrate

        """
    def march(self) -> None:
        """march(self: pymem3dg._core.ConjugateGradient) -> None


        stepping forward

        """
    def saveData(
        self, ifOutputTrajFile: bool, ifOutputMeshFile: bool, ifPrintToConsole: bool
    ) -> None:
        """saveData(self: pymem3dg._core.ConjugateGradient, ifOutputTrajFile: bool, ifOutputMeshFile: bool, ifPrintToConsole: bool) -> None


        save data to output directory

        """
    def status(self) -> None:
        """status(self: pymem3dg._core.ConjugateGradient) -> None


        status computation and thresholding

        """
    def step(self, n: int) -> None:
        """step(self: pymem3dg._core.ConjugateGradient, n: int) -> None


        step for n iterations

        """
    @property
    def characteristicTimeStep(self) -> float: ...
    @property
    def savePeriod(self) -> float: ...
    @property
    def tolerance(self) -> float: ...
    @property
    def totalTime(self) -> float: ...

class DPD:
    gamma: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class Dirichlet:
    eta: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class Energy:
    def __init__(self) -> None:
        """__init__(self: pymem3dg._core.Energy) -> None"""
    @property
    def adsorptionEnergy(self) -> float: ...
    @property
    def aggregationEnergy(self) -> float: ...
    @property
    def areaDifferenceEnergy(self) -> float: ...
    @property
    def deviatoricCurvatureEnergy(self) -> float: ...
    @property
    def dirichletEnergy(self) -> float: ...
    @property
    def edgeSpringEnergy(self) -> float: ...
    @property
    def entropyEnergy(self) -> float: ...
    @property
    def externalWork(self) -> float: ...
    @property
    def faceSpringEnergy(self) -> float: ...
    @property
    def kineticEnergy(self) -> float: ...
    @property
    def lcrSpringEnergy(self) -> float: ...
    @property
    def potentialEnergy(self) -> float: ...
    @property
    def pressureEnergy(self) -> float: ...
    @property
    def proteinInteriorPenalty(self) -> float: ...
    @property
    def spontaneousCurvatureEnergy(self) -> float: ...
    @property
    def surfaceEnergy(self) -> float: ...
    @property
    def totalEnergy(self) -> float: ...

class Entropy:
    xi: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class Euler:
    c1: float
    ifAdaptiveStep: bool
    ifJustGeometryPly: bool
    ifOutputMeshFile: bool
    ifOutputTrajFile: bool
    ifPrintToConsole: bool
    isBacktrack: bool
    outputDirectory: str
    rho: float
    trajFileName: str
    @overload
    def __init__(
        self,
        system,
        characteristicTimeStep: float,
        totalTime: float,
        savePeriod: float,
        tolerance: float,
        outputDirectory: str,
        frame: int = ...,
    ) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pymem3dg._core.Euler, system: mem3dg::solver::System, characteristicTimeStep: float, totalTime: float, savePeriod: float, tolerance: float, outputDirectory: str, frame: int = 0) -> None


                Euler integrator (steepest descent) constructor


        2. __init__(self: pymem3dg._core.Euler, system: mem3dg::solver::System, characteristicTimeStep: float, tolerance: float, outputDirectory: str) -> None


                Euler integrator (steepest descent) constructor

        """
    @overload
    def __init__(
        self,
        system,
        characteristicTimeStep: float,
        tolerance: float,
        outputDirectory: str,
    ) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pymem3dg._core.Euler, system: mem3dg::solver::System, characteristicTimeStep: float, totalTime: float, savePeriod: float, tolerance: float, outputDirectory: str, frame: int = 0) -> None


                Euler integrator (steepest descent) constructor


        2. __init__(self: pymem3dg._core.Euler, system: mem3dg::solver::System, characteristicTimeStep: float, tolerance: float, outputDirectory: str) -> None


                Euler integrator (steepest descent) constructor

        """
    def closeMutableNetcdfFile(self) -> None:
        """closeMutableNetcdfFile(self: pymem3dg._core.Euler) -> None


        close netcdf file

        """
    def createMutableNetcdfFile(self, isContinue: bool) -> None:
        """createMutableNetcdfFile(self: pymem3dg._core.Euler, isContinue: bool) -> None


        create netcdf file

        """
    def integrate(self) -> bool:
        """integrate(self: pymem3dg._core.Euler) -> bool


        integrate

        """
    def march(self) -> None:
        """march(self: pymem3dg._core.Euler) -> None


        stepping forward

        """
    def saveData(
        self, ifOutputTrajFile: bool, ifOutputMeshFile: bool, ifPrintToConsole: bool
    ) -> None:
        """saveData(self: pymem3dg._core.Euler, ifOutputTrajFile: bool, ifOutputMeshFile: bool, ifPrintToConsole: bool) -> None


        save data to output directory

        """
    def saveMutableNetcdfData(self) -> None:
        """saveMutableNetcdfData(self: pymem3dg._core.Euler) -> None


        write to netcdf file

        """
    def status(self) -> None:
        """status(self: pymem3dg._core.Euler) -> None


        status computation and thresholding

        """
    def step(self, n: int) -> None:
        """step(self: pymem3dg._core.Euler, n: int) -> None


        step for n iterations

        """
    @property
    def EXIT(self) -> bool: ...
    @property
    def characteristicTimeStep(self) -> float: ...
    @property
    def savePeriod(self) -> float: ...
    @property
    def timeStep(self) -> float: ...
    @property
    def tolerance(self) -> float: ...
    @property
    def totalTime(self) -> float: ...

class External:
    form: Callable[
        [
            numpy.ndarray[numpy.float64[m, 3]],
            numpy.ndarray[numpy.float64[m, 1]],
            float,
            numpy.ndarray[numpy.float64[m, 1]],
        ],
        numpy.ndarray[numpy.float64[m, 3]],
    ]
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class Forces:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def getAdsorptionForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getAdsorptionForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the adsorption force

        """
    def getAdsorptionPotential(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getAdsorptionPotential(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 1]]


        get the adsorption potential

        """
    def getAggregationForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getAggregationForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the aggregation force

        """
    def getAggregationPotential(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getAggregationPotential(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 1]]


        get the aggregation potential

        """
    def getAreaDifferenceForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getAreaDifferenceForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the area difference force of the system

        """
    def getCapillaryForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getCapillaryForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the tension-induced capillary Force

        """
    def getChemicalPotential(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getChemicalPotential(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 1]]


        get the chemical Potential

        """
    def getDeviatoricCurvatureForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getDeviatoricCurvatureForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the deviatoric curvature force of the system

        """
    def getDeviatoricCurvaturePotential(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getDeviatoricCurvaturePotential(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 1]]


        get the deviatoric curvature potential

        """
    def getDirichletPotential(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getDirichletPotential(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 1]]


        get the dirichlet Potential

        """
    def getEntropyForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getEntropyForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the entropy force

        """
    def getEntropyPotential(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getEntropyPotential(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 1]]


        get the entropy Potential

        """
    def getExternalForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getExternalForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the externally-applied Force

        """
    def getInteriorPenaltyPotential(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getInteriorPenaltyPotential(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 1]]


        get the interior point potential

        """
    def getLineCapillaryForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getLineCapillaryForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the interfacial line tension

        """
    def getMechanicalForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getMechanicalForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the the total mechanical force

        """
    def getOsmoticForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getOsmoticForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the osmotic force

        """
    def getOsmoticPressure(self) -> float:
        """getOsmoticPressure(self: pymem3dg._core.Forces) -> float


        get the osmotic pressure

        """
    def getSpontaneousCurvatureForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getSpontaneousCurvatureForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the spontaneous curvature force of the system

        """
    def getSpontaneousCurvatureForceVec_areaGrad(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getSpontaneousCurvatureForceVec_areaGrad(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the area gradient component of the spontaneous curvature force of the system

        """
    def getSpontaneousCurvatureForceVec_gaussVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getSpontaneousCurvatureForceVec_gaussVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the the gaussian curvature vector component of the spontaneous curvature force of the system

        """
    def getSpontaneousCurvatureForceVec_schlafliVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getSpontaneousCurvatureForceVec_schlafliVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the Schlaflic (smoothing) component of the spontaneous curvature force of the system

        """
    def getSpontaneousCurvaturePotential(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getSpontaneousCurvaturePotential(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 1]]


        get the spontaneous curvature potential

        """
    def getSpringForceVec(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getSpringForceVec(self: pymem3dg._core.Forces) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the spring force

        """
    def getSurfaceTension(self) -> float:
        """getSurfaceTension(self: pymem3dg._core.Forces) -> float


        get the Surface tension

        """

class Geometry:
    @overload
    def __init__(
        self, inputMesh: str, notableVertex: numpy.ndarray[bool[m, 1]]
    ) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pymem3dg._core.Geometry, inputMesh: str, notableVertex: numpy.ndarray[bool[m, 1]]) -> None


                Geometry constructor with .ply files.


        2. __init__(self: pymem3dg._core.Geometry, inputMesh: str) -> None


                  Geometry constructor with .ply files.


        3. __init__(self: pymem3dg._core.Geometry, faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]], notableVertex: numpy.ndarray[bool[m, 1]]) -> None


                Geometry constructor with Matrices


        4. __init__(self: pymem3dg._core.Geometry, faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]]) -> None


                Geometry constructor with Matrices


        5. __init__(self: pymem3dg._core.Geometry, trajFile: str, startingFrame: int) -> None


                Geometry constructor with NetCDF trajectory file

        """
    @overload
    def __init__(self, inputMesh: str) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pymem3dg._core.Geometry, inputMesh: str, notableVertex: numpy.ndarray[bool[m, 1]]) -> None


                Geometry constructor with .ply files.


        2. __init__(self: pymem3dg._core.Geometry, inputMesh: str) -> None


                  Geometry constructor with .ply files.


        3. __init__(self: pymem3dg._core.Geometry, faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]], notableVertex: numpy.ndarray[bool[m, 1]]) -> None


                Geometry constructor with Matrices


        4. __init__(self: pymem3dg._core.Geometry, faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]]) -> None


                Geometry constructor with Matrices


        5. __init__(self: pymem3dg._core.Geometry, trajFile: str, startingFrame: int) -> None


                Geometry constructor with NetCDF trajectory file

        """
    @overload
    def __init__(
        self,
        faceMatrix: numpy.ndarray[numpy.uint64[m, 3]],
        vertexMatrix: numpy.ndarray[numpy.float64[m, 3]],
        notableVertex: numpy.ndarray[bool[m, 1]],
    ) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pymem3dg._core.Geometry, inputMesh: str, notableVertex: numpy.ndarray[bool[m, 1]]) -> None


                Geometry constructor with .ply files.


        2. __init__(self: pymem3dg._core.Geometry, inputMesh: str) -> None


                  Geometry constructor with .ply files.


        3. __init__(self: pymem3dg._core.Geometry, faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]], notableVertex: numpy.ndarray[bool[m, 1]]) -> None


                Geometry constructor with Matrices


        4. __init__(self: pymem3dg._core.Geometry, faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]]) -> None


                Geometry constructor with Matrices


        5. __init__(self: pymem3dg._core.Geometry, trajFile: str, startingFrame: int) -> None


                Geometry constructor with NetCDF trajectory file

        """
    @overload
    def __init__(
        self,
        faceMatrix: numpy.ndarray[numpy.uint64[m, 3]],
        vertexMatrix: numpy.ndarray[numpy.float64[m, 3]],
    ) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pymem3dg._core.Geometry, inputMesh: str, notableVertex: numpy.ndarray[bool[m, 1]]) -> None


                Geometry constructor with .ply files.


        2. __init__(self: pymem3dg._core.Geometry, inputMesh: str) -> None


                  Geometry constructor with .ply files.


        3. __init__(self: pymem3dg._core.Geometry, faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]], notableVertex: numpy.ndarray[bool[m, 1]]) -> None


                Geometry constructor with Matrices


        4. __init__(self: pymem3dg._core.Geometry, faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]]) -> None


                Geometry constructor with Matrices


        5. __init__(self: pymem3dg._core.Geometry, trajFile: str, startingFrame: int) -> None


                Geometry constructor with NetCDF trajectory file

        """
    @overload
    def __init__(self, trajFile: str, startingFrame: int) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pymem3dg._core.Geometry, inputMesh: str, notableVertex: numpy.ndarray[bool[m, 1]]) -> None


                Geometry constructor with .ply files.


        2. __init__(self: pymem3dg._core.Geometry, inputMesh: str) -> None


                  Geometry constructor with .ply files.


        3. __init__(self: pymem3dg._core.Geometry, faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]], notableVertex: numpy.ndarray[bool[m, 1]]) -> None


                Geometry constructor with Matrices


        4. __init__(self: pymem3dg._core.Geometry, faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]]) -> None


                Geometry constructor with Matrices


        5. __init__(self: pymem3dg._core.Geometry, trajFile: str, startingFrame: int) -> None


                Geometry constructor with NetCDF trajectory file

        """
    @overload
    def computeGeodesicDistance(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """computeGeodesicDistance(*args, **kwargs)
        Overloaded function.

        1. computeGeodesicDistance(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 1]]


                  compute the geodesic distance centered around Center cached in System

                  return: NDarray[double]


        2. computeGeodesicDistance(self: pymem3dg._core.Geometry, points: list[int]) -> numpy.ndarray[numpy.float64[m, 1]]


                    compute the geodesic distance from a set of vertex indices

                    return: NDarray[double]

        """
    @overload
    def computeGeodesicDistance(
        self, points: list[int]
    ) -> numpy.ndarray[numpy.float64[m, 1]]:
        """computeGeodesicDistance(*args, **kwargs)
        Overloaded function.

        1. computeGeodesicDistance(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 1]]


                  compute the geodesic distance centered around Center cached in System

                  return: NDarray[double]


        2. computeGeodesicDistance(self: pymem3dg._core.Geometry, points: list[int]) -> numpy.ndarray[numpy.float64[m, 1]]


                    compute the geodesic distance from a set of vertex indices

                    return: NDarray[double]

        """
    def getCotanLaplacian(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """getCotanLaplacian(self: pymem3dg._core.Geometry) -> scipy.sparse.csc_matrix[numpy.float64]


        get the Cotan Laplacian matrix of the mesh

        """
    def getEdgeAdjacencyMatrix(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """getEdgeAdjacencyMatrix(self: pymem3dg._core.Geometry) -> scipy.sparse.csc_matrix[numpy.float64]


        get the signed F-E edge adjacency matrix, equivalent of d1 operator

        """
    def getEdgeLengths(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getEdgeLengths(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 1]]


        get edge lengths

        """
    def getFaceAreas(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getFaceAreas(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 1]]


        get face area

        """
    def getFaceMatrix(self) -> numpy.ndarray[numpy.uint64[m, n]]:
        """getFaceMatrix(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.uint64[m, n]]


        get the face matrix

        """
    def getGeodesicDistance(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getGeodesicDistance(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 1]]


        get the geodesic distance from notable vertices of the mesh

        """
    def getLumpedMassMatrix(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """getLumpedMassMatrix(self: pymem3dg._core.Geometry) -> scipy.sparse.csc_matrix[numpy.float64]


        get the lumped mass matrix of the mesh

        """
    def getNotableVertex(self) -> numpy.ndarray[bool[m, 1]]:
        """getNotableVertex(self: pymem3dg._core.Geometry) -> numpy.ndarray[bool[m, 1]]


        get vertex data to track the notable vertex, which may or may not be a single vertex

        """
    def getPolyscopeEdgeOrientations(self) -> numpy.ndarray[numpy.int8[m, 1]]:
        """getPolyscopeEdgeOrientations(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.int8[m, 1]]


        get polyscope edge orientation

        """
    def getPolyscopePermutations(self, *args, **kwargs):
        """getPolyscopePermutations(self: pymem3dg._core.Geometry) -> Annotated[list[tuple[list[int], int]], FixedSize(5)]


        get polyscope permutation

        """
    def getSurfaceArea(self) -> float:
        """getSurfaceArea(self: pymem3dg._core.Geometry) -> float


        get the surface area of the mesh

        """
    def getVertexAdjacencyMatrix(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """getVertexAdjacencyMatrix(self: pymem3dg._core.Geometry) -> scipy.sparse.csc_matrix[numpy.float64]


        get the signed E-V vertex adjacency matrix, equivalent of d0 operator

        """
    def getVertexDualAreas(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getVertexDualAreas(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 1]]


        get vertex dual area

        """
    def getVertexGaussianCurvatureVectors(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getVertexGaussianCurvatureVectors(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the integrated vector Gaussian Curvature

        """
    def getVertexGaussianCurvatures(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getVertexGaussianCurvatures(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 1]]


        get the integrated scalar Gaussian Curvature

        """
    def getVertexMatrix(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getVertexMatrix(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the vertex matrix

        """
    def getVertexMeanCurvatureVectors(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getVertexMeanCurvatureVectors(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the integrated vector Mean Curvature

        """
    def getVertexMeanCurvatures(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getVertexMeanCurvatures(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 1]]


        get the integrated scalar mean curvature

        """
    def getVertexNormals(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getVertexNormals(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get angle-weighted normal on vertices

        """
    def getVertexSchlafliLaplacianMeanCurvatureVectors(
        self, arg0: numpy.ndarray[numpy.float64[m, 1]]
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getVertexSchlafliLaplacianMeanCurvatureVectors(self: pymem3dg._core.Geometry, arg0: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the vertex Schlafli based Laplacian of mean curvature Vectors

        """
    def getVertexVolumeVariationVectors(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getVertexVolumeVariationVectors(self: pymem3dg._core.Geometry) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]


        get the integrated vector volume variation (dual area)

        """
    def getVolume(self) -> float:
        """getVolume(self: pymem3dg._core.Geometry) -> float


        get the enclosed volume of the mesh

        """
    def isBoundary(self, index: int) -> bool:
        """isBoundary(self: pymem3dg._core.Geometry, index: int) -> bool


        check if the vertex is on the boundary

        return: bool

        """
    def nEdges(self) -> int:
        """nEdges(self: pymem3dg._core.Geometry) -> int"""
    def nFaces(self) -> int:
        """nFaces(self: pymem3dg._core.Geometry) -> int"""
    def nVertices(self) -> int:
        """nVertices(self: pymem3dg._core.Geometry) -> int"""
    def setInputVertexPositions(self, arg0: numpy.ndarray[numpy.float64[m, 3]]) -> None:
        """setInputVertexPositions(self: pymem3dg._core.Geometry, arg0: numpy.ndarray[numpy.float64[m, 3]]) -> None


        set the vertex matrix

        """

class MeshMutator:
    collapseFlat: bool
    collapseFlatScaleFactor: float
    collapseSkinny: bool
    collapseSmall: bool
    curvTol: float
    flipNonDelaunay: bool
    flipNonDelaunayRequireFlat: bool
    isShiftVertex: bool
    isSmoothenMesh: bool
    maximumEdgeLength: float
    maximumFaceArea: float
    minimumEdgeLength: float
    minimumFaceArea: float
    mutateMeshPeriod: int
    splitCurved: bool
    splitCurvedScaleFactor: float
    splitFat: bool
    splitLarge: bool
    splitLong: bool
    splitSharp: bool
    splitSkinnyDelaunay: bool
    def __init__(self) -> None:
        """__init__(self: pymem3dg._core.MeshMutator) -> None


        Constructor of mesh mutator object

        """
    @property
    def isChangeTopology(self) -> bool: ...
    @property
    def isCollapseEdge(self) -> bool: ...
    @property
    def isFlipEdge(self) -> bool: ...
    @property
    def isSplitEdge(self) -> bool: ...

class MeshProcessor:
    meshMutator: MeshMutator
    def __init__(self) -> None:
        """__init__(self: pymem3dg._core.MeshProcessor) -> None


        MeshProcessor constructor

        """
    @property
    def isMeshMutate(self) -> bool: ...

class Osmotic:
    form: Callable[[float], tuple[float, float]]
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class Parameters:
    adsorption: Adsorption
    aggregation: Aggregation
    bending: Bending
    boundary: Boundary
    damping: float
    dirichlet: Dirichlet
    dpd: DPD
    entropy: Entropy
    external: External
    osmotic: Osmotic
    point: Point
    protein: Protein
    proteinMobility: float
    selfAvoidance: SelfAvoidance
    spring: spring
    temperature: float
    tension: Tension
    variation: Variation
    def __init__(self) -> None:
        """__init__(self: pymem3dg._core.Parameters) -> None"""

class Point:
    prescribeNotableVertex: Callable[[Geometry], numpy.ndarray[bool[m, 1]]]
    updateGeodesicsPeriod: float
    updateNotableVertexPeriod: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class Protein:
    prescribeProteinDensityDistribution: Callable[
        [
            Geometry,
            float,
            numpy.ndarray[numpy.float64[m, 1]],
            numpy.ndarray[numpy.float64[m, 1]],
        ],
        numpy.ndarray[numpy.float64[m, 1]],
    ]
    proteinInteriorPenalty: float
    updateProteinDensityDistributionPeriod: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class SelfAvoidance:
    d: float
    mu: float
    n: int
    p: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

class System:
    meshProcessor: MeshProcessor
    parameters: Incomplete
    time: float
    @overload
    def __init__(
        self,
        geometry: Geometry,
        proteinDensity: numpy.ndarray[numpy.float64[m, 1]],
        velocity: numpy.ndarray[numpy.float64[m, 3]],
        parameters,
        time: float = ...,
    ) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pymem3dg._core.System, geometry: pymem3dg._core.Geometry, proteinDensity: numpy.ndarray[numpy.float64[m, 1]], velocity: numpy.ndarray[numpy.float64[m, 3]], parameters: mem3dg::solver::Parameters, time: float = 0) -> None


                System constructor with Matrices


        2. __init__(self: pymem3dg._core.System, geometry: pymem3dg._core.Geometry, parameters: mem3dg::solver::Parameters, time: float = 0) -> None


                System constructor with Matrices


        3. __init__(self: pymem3dg._core.System, geometry: pymem3dg._core.Geometry, trajFile: str, startingFrame: int, parameters: mem3dg::solver::Parameters) -> None


                System constructor with NetCDF trajectory file

        """
    @overload
    def __init__(self, geometry: Geometry, parameters, time: float = ...) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pymem3dg._core.System, geometry: pymem3dg._core.Geometry, proteinDensity: numpy.ndarray[numpy.float64[m, 1]], velocity: numpy.ndarray[numpy.float64[m, 3]], parameters: mem3dg::solver::Parameters, time: float = 0) -> None


                System constructor with Matrices


        2. __init__(self: pymem3dg._core.System, geometry: pymem3dg._core.Geometry, parameters: mem3dg::solver::Parameters, time: float = 0) -> None


                System constructor with Matrices


        3. __init__(self: pymem3dg._core.System, geometry: pymem3dg._core.Geometry, trajFile: str, startingFrame: int, parameters: mem3dg::solver::Parameters) -> None


                System constructor with NetCDF trajectory file

        """
    @overload
    def __init__(
        self, geometry: Geometry, trajFile: str, startingFrame: int, parameters
    ) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pymem3dg._core.System, geometry: pymem3dg._core.Geometry, proteinDensity: numpy.ndarray[numpy.float64[m, 1]], velocity: numpy.ndarray[numpy.float64[m, 3]], parameters: mem3dg::solver::Parameters, time: float = 0) -> None


                System constructor with Matrices


        2. __init__(self: pymem3dg._core.System, geometry: pymem3dg._core.Geometry, parameters: mem3dg::solver::Parameters, time: float = 0) -> None


                System constructor with Matrices


        3. __init__(self: pymem3dg._core.System, geometry: pymem3dg._core.Geometry, trajFile: str, startingFrame: int, parameters: mem3dg::solver::Parameters) -> None


                System constructor with NetCDF trajectory file

        """
    def addNonconservativeForcing(self, timeStep: float = ...) -> None:
        """addNonconservativeForcing(self: pymem3dg._core.System, timeStep: float = 0) -> None


        Compute and append all non-conservative forces, update mechanicalForce(Vec) and mechErrorNorm

        """
    def computeConservativeForcing(self) -> None:
        """computeConservativeForcing(self: pymem3dg._core.System) -> None


        Compute and update all conservative forces, update mechanicalForce(Vec) with conservativeForce(Vec)

        """
    def computeInPlaneFluxForm(
        self, chemicalPotential: numpy.ndarray[numpy.float64[m, 1]]
    ) -> numpy.ndarray[numpy.float64[m, 1]]:
        """computeInPlaneFluxForm(self: pymem3dg._core.System, chemicalPotential: numpy.ndarray[numpy.float64[m, 1]]) -> numpy.ndarray[numpy.float64[m, 1]]


        Compute in plane flux form from chemical potential

        Args:
            chemicalPotential (npt.NDArray[np.float64]): components (or sum) of chemical potential

        Returns
            np.NDArray[np.float64]: in plane flux form on edges

        """
    def computeIntegratedPower(self, dt: float) -> float:
        """computeIntegratedPower(self: pymem3dg._core.System, dt: float) -> float


        Intermediate function to integrate the power

        """
    def computeTotalEnergy(self) -> float:
        """computeTotalEnergy(self: pymem3dg._core.System) -> float


        compute the total energy, where total energy = kinetic energy + potential energy - external work

        """
    def getEnergy(self, *args, **kwargs):
        """getEnergy(self: pymem3dg._core.System) -> mem3dg::solver::Energy


        Get the energy

        Args:
            System (:py:class:`System`): System of interest

        Returns:
            :py:class:`float`: Energy

        """
    def getForces(self) -> Forces:
        """getForces(self: pymem3dg._core.System) -> pymem3dg._core.Forces


        Get the pointwise forces

        Args:
            System (:py:class:`System`): System of interest

        Returns:
            :py:class:`pymem3dg.Forces`: Forces

        """
    def getGeometry(self) -> Geometry:
        """getGeometry(self: pymem3dg._core.System) -> pymem3dg._core.Geometry


        Get the geometry

        Args:
            System (:py:class:`System`): System of interest

        Returns:
            :py:class:`pymem3dg.Geometry`: Geometry

        """
    def getProteinDensity(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getProteinDensity(self: pymem3dg._core.System) -> numpy.ndarray[numpy.float64[m, 1]]

        get the protein Density
        """
    def getProteinRateOfChange(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getProteinRateOfChange(self: pymem3dg._core.System) -> numpy.ndarray[numpy.float64[m, 1]]

        get the protein rate of change
        """
    def getSpontaneousCurvature(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """getSpontaneousCurvature(self: pymem3dg._core.System) -> numpy.ndarray[numpy.float64[m, 1]]

        get the pointwise spontaneous curvature
        """
    def getVelocity(
        self,
    ) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]:
        """getVelocity(self: pymem3dg._core.System) -> numpy.ndarray[numpy.float64[m, 3], flags.writeable, flags.c_contiguous]

        get the vertex velocity matrix
        """
    def initialize(self, ifMutateMesh: bool = ...) -> None:
        """initialize(self: pymem3dg._core.System, ifMutateMesh: bool = False) -> None


        initialize the system

        """
    def mutateMesh(self, nMutation: int = ...) -> None:
        """mutateMesh(self: pymem3dg._core.System, nMutation: int = 1) -> None


        mutate the mesh

        """
    def prescribeExternalForce(self) -> numpy.ndarray[numpy.float64[m, 3]]:
        """prescribeExternalForce(self: pymem3dg._core.System) -> numpy.ndarray[numpy.float64[m, 3]]


        prescribe the External Force

        """
    def saveRichData(self, pathToSave: str, isJustGeometry: bool = ...) -> None:
        """saveRichData(self: pymem3dg._core.System, pathToSave: str, isJustGeometry: bool = False) -> None


        save snapshot data to directory

        """
    def smoothenMesh(
        self, initStep: float, target: float, maxIteration: int
    ) -> numpy.ndarray[bool[m, 1]]:
        """smoothenMesh(self: pymem3dg._core.System, initStep: float, target: float, maxIteration: int) -> numpy.ndarray[bool[m, 1]]


        smoothen the mesh using spontaneous curvature force

        """
    def testConservativeForcing(self, timeStep: float) -> bool:
        """testConservativeForcing(self: pymem3dg._core.System, timeStep: float) -> bool


        test conservative force computation by validating energy decrease

        """
    def updateConfigurations(self) -> None:
        """updateConfigurations(self: pymem3dg._core.System) -> None


        update the system configuration due to changes in state variables (e.g vertex positions or protein density)

        """

class Tension:
    form: Callable[[float], tuple[float, float]]
    def __init__(self) -> None:
        """__init__(self: pymem3dg._core.Tension) -> None"""

class Variation:
    geodesicMask: float
    isProteinConservation: bool
    isProteinVariation: bool
    isShapeVariation: bool
    updateMaskPeriod: float
    def __init__(self) -> None:
        """__init__(self: pymem3dg._core.Variation) -> None"""

class VelocityVerlet:
    c1: float
    ifAdaptiveStep: bool
    ifOutputMeshFile: bool
    ifOutputTrajFile: bool
    ifPrintToConsole: bool
    isBacktrack: bool
    isCapEnergy: bool
    outputDirectory: str
    rho: float
    trajFileName: str
    def __init__(
        self,
        system,
        characteristicTimeStep: float,
        totalTime: float,
        savePeriod: float,
        tolerance: float,
        outputDirectory: str,
        frame: int = ...,
    ) -> None:
        """__init__(self: pymem3dg._core.VelocityVerlet, system: mem3dg::solver::System, characteristicTimeStep: float, totalTime: float, savePeriod: float, tolerance: float, outputDirectory: str, frame: int = 0) -> None


        Velocity Verlet integrator constructor

        """
    def integrate(self) -> bool:
        """integrate(self: pymem3dg._core.VelocityVerlet) -> bool


        integrate

        """
    def march(self) -> None:
        """march(self: pymem3dg._core.VelocityVerlet) -> None


        stepping forward

        """
    def saveData(
        self, ifOutputTrajFile: bool, ifOutputMeshFile: bool, ifPrintToConsole: bool
    ) -> None:
        """saveData(self: pymem3dg._core.VelocityVerlet, ifOutputTrajFile: bool, ifOutputMeshFile: bool, ifPrintToConsole: bool) -> None


        save data to output directory

        """
    def status(self) -> None:
        """status(self: pymem3dg._core.VelocityVerlet) -> None


        status computation and thresholding

        """
    def step(self, n: int) -> None:
        """step(self: pymem3dg._core.VelocityVerlet, n: int) -> None


        step for n iterations

        """

class spring:
    Kse: float
    Ksl: float
    Kst: float
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

def getCylinder(
    radius: float,
    radialSubdivision: int,
    axialSubdivision: int,
    frequency: float = ...,
    amplitude: float = ...,
) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]:
    """getCylinder(radius: float, radialSubdivision: int, axialSubdivision: int, frequency: float = 1, amplitude: float = 0) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]

    get face and vertex matrix of a non-capped cylinder
    """

def getData(
    plyName: str, elementName: str, propertyName: str
) -> numpy.ndarray[numpy.float64[m, 1]]:
    """getData(plyName: str, elementName: str, propertyName: str) -> numpy.ndarray[numpy.float64[m, 1]]

    read richData from .ply file

    """

def getDataElementName(plyName: str) -> list[str]:
    """getDataElementName(plyName: str) -> list[str]

    retrieve all richData element name from .ply file. Namely the list of the places where data live in, such as vertex, edge or face.

    """

def getDataPropertyName(plyName: str, elementName: str) -> list[str]:
    """getDataPropertyName(plyName: str, elementName: str) -> list[str]

    retrieve all richData property name from .ply file. Namely the list of the data where data on the particular element, such as vertex, edge or face.

    """

def getDiamond(
    dihedral: float,
) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]:
    """getDiamond(dihedral: float) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]

    get face and vertex matrix of diamond
    """

def getFaceAndVertexMatrix(
    plyName: str,
) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]:
    """getFaceAndVertexMatrix(plyName: str) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]


    read face matrix and vertex matrix from .ply file

    """

def getFaceSurfacePointClosestToEmbeddedCoordinate(*args, **kwargs):
    """getFaceSurfacePointClosestToEmbeddedCoordinate(faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]], embeddedCoordinate: Annotated[list[float], FixedSize(3)], filter: numpy.ndarray[bool[m, 1]], accountedCoordinate: Annotated[list[bool], FixedSize(3)] = [True, True, True]) -> tuple[int, Annotated[list[float], FixedSize(3)]]

    find the face surface point closest to a embedded coordinate
    """

def getHexagon(
    radius: float, subdivision: int = ...
) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]:
    """getHexagon(radius: float, subdivision: int = 0) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]

    get face and vertex matrix of Hexagon
    """

def getIcosphere(
    radius: float, subdivision: int = ...
) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]:
    """getIcosphere(radius: float, subdivision: int = 0) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]

    get face and vertex matrix of icosphere
    """

def getTetrahedron() -> (
    tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]
):
    """getTetrahedron() -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]

    get face and vertex matrix of tetrahedron
    """

def getVertexClosestToEmbeddedCoordinate(
    vertexMatrix: numpy.ndarray[numpy.float64[m, 3]],
    embeddedCoordinate,
    filter: numpy.ndarray[bool[m, 1]],
    accountedCoordinate=...,
) -> int:
    """getVertexClosestToEmbeddedCoordinate(vertexMatrix: numpy.ndarray[numpy.float64[m, 3]], embeddedCoordinate: Annotated[list[float], FixedSize(3)], filter: numpy.ndarray[bool[m, 1]], accountedCoordinate: Annotated[list[bool], FixedSize(3)] = [True, True, True]) -> int

    find the index of vertex closest to a embedded coordinate
    """

def getVertexFurthestFromBoundary(
    faceMatrix: numpy.ndarray[numpy.uint64[m, 3]],
    vertexMatrix: numpy.ndarray[numpy.float64[m, 3]],
) -> int:
    """getVertexFurthestFromBoundary(faceMatrix: numpy.ndarray[numpy.uint64[m, 3]], vertexMatrix: numpy.ndarray[numpy.float64[m, 3]]) -> int

    find the vertex furthest away from the boundaries
    """

def linearSubdivide(
    face: numpy.ndarray[numpy.uint64[m, 3]],
    vertex: numpy.ndarray[numpy.float64[m, 3]],
    nSub: int,
) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]:
    """linearSubdivide(face: numpy.ndarray[numpy.uint64[m, 3]], vertex: numpy.ndarray[numpy.float64[m, 3]], nSub: int) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]


    subdivide the mesh with linear interpolation

    """

def loopSubdivide(
    face: numpy.ndarray[numpy.uint64[m, 3]],
    vertex: numpy.ndarray[numpy.float64[m, 3]],
    nSub: int,
) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]:
    """loopSubdivide(face: numpy.ndarray[numpy.uint64[m, 3]], vertex: numpy.ndarray[numpy.float64[m, 3]], nSub: int) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]


    subdivide the mesh in loop scheme

    """

def processSoup(
    meshName: str,
) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]:
    """processSoup(meshName: str) -> tuple[numpy.ndarray[numpy.uint64[m, 3]], numpy.ndarray[numpy.float64[m, 3]]]

    process soup data in .ply and return face and vertex matrices.

    """

def startProfiler(filename: str) -> None:
    """startProfiler(filename: str) -> None

    start profiler
    """

def stopProfiler() -> None:
    """stopProfiler() -> None

    stop profiler
    """
