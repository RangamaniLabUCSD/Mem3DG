import numpy as np
import matplotlib.pyplot as plt

import pymem3dg as dg
import pymem3dg.read as dg_read
import pymem3dg.util as dg_util

import polyscope as ps
import polyscope.imgui as psim


def matplotlibStyle(small=6, medium=8, large=10):
    """ Formatting style of matplotlib """
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.family'] = "sans-serif"
    # mpl.rcParams.update({'font.size': 8})
    SMALL_SIZE = small
    MEDIUM_SIZE = medium
    BIGGER_SIZE = large
    plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    plt.rc('pdf', fonttype=42)


def getColorMap(fig, range, orient="horizontal", num_points=100):
    x = np.linspace(0, 1, num=num_points)
    y = np.linspace(0, 1, num=num_points)
    fig.gca().set_visible(False)
    Temp = np.linspace(range[0], range[1], num=num_points)
    plt.scatter(x, y, c=Temp)
    plt.colorbar(orientation=orient)
    plt.tight_layout()
    return fig


def animate(trajNc, **kwargs):
    buildSystem = kwargs is not None
    hasParameters = "parameters" in kwargs
    if hasParameters:
        parameters = kwargs.get("parameters")

    ps.init()
    ps.set_transparency_mode("pretty")

    prevFrame = 0
    currFrame = 0
    isPlay = False
    transparency = 0.4
    maxFrame = dg_read.sizeOf(trajNc) - 1

    if buildSystem:
        system = dg.System(trajNc, currFrame)
        if hasParameters:
            system = dg.System(trajNc, currFrame, parameters, True)

    def show(trajNc):
        nonlocal currFrame, system
        face, vertex = dg_read.readMeshByNc(trajNc, currFrame)
        proteinDensity = dg_read.readMeshDataByNc(
            trajNc, currFrame, "Trajectory", "proteindensity", 1)
        limit = (np.min(proteinDensity), np.max(proteinDensity))
        mesh = ps.register_surface_mesh(
            "mesh", vertex, face, transparency=transparency)
        mesh.add_scalar_quantity(
            "proteinDensity", proteinDensity, vminmax=limit)

        if buildSystem:
            if hasParameters:
                system = dg.System(trajNc, currFrame, parameters, True)
                system.computePhysicalForcing()
                if (kwargs.get('bendingForce')):
                    bendingForce = system.forces.getBendingForce()
                    mesh.add_vector_quantity("bendingForce", bendingForce)
                    mesh.add_scalar_quantity("<bendingForce,n>", dg_util.rowwiseDotProduct(
                        bendingForce, system.getVertexNormal()))
                if (kwargs.get('externalForce')):
                    externalForce = system.forces.getExternalForce()
                    mesh.add_vector_quantity("externalForce", externalForce)
                    mesh.add_scalar_quantity("<externalForce,n>", dg_util.rowwiseDotProduct(
                        externalForce, system.getVertexNormal()))
                if (kwargs.get('mechanicalForce')):
                    mechanicalForce = system.forces.getMechanicalForce()
                    mesh.add_vector_quantity("mechanicalForce", mechanicalForce)
                    mesh.add_scalar_quantity("<mechanicalForce,n>", dg_util.rowwiseDotProduct(
                        mechanicalForce, system.getVertexNormal()))
                if (kwargs.get('capillaryForce')):
                    capillaryForce = system.forces.getCapillaryForce()
                    mesh.add_vector_quantity("capillaryForce", capillaryForce)
                    mesh.add_scalar_quantity("<capillaryForce,n>", dg_util.rowwiseDotProduct(
                        capillaryForce, system.getVertexNormal()))
                if (kwargs.get('lineCapillaryForce')):
                    lineCapillaryForce = system.forces.getLineCapillaryForce()
                    mesh.add_vector_quantity("lineCapillaryForce", lineCapillaryForce)
                    mesh.add_scalar_quantity("<lineCapillaryForce,n>", dg_util.rowwiseDotProduct(
                        lineCapillaryForce, system.getVertexNormal()))
                if (kwargs.get('osmoticForce')):
                    osmoticForce = system.forces.getOsmoticForce()
                    mesh.add_vector_quantity("osmoticForce", osmoticForce)
                    mesh.add_scalar_quantity("<osmoticForce,n>", dg_util.rowwiseDotProduct(
                        osmoticForce, system.getVertexNormal()))
            else:
                system = dg.System(trajNc, currFrame)
                
            if (kwargs.get('meanCurvature')):
                meanCurvature = system.getVertexMeanCurvature()
                mesh.add_scalar_quantity("meanCurvature", meanCurvature)
            if (kwargs.get('gaussianCurvature')):
                gaussianCurvature = system.getVertexGaussianCurvature()
                mesh.add_scalar_quantity("gaussianCurvature", gaussianCurvature)

    show(trajNc)
    
    def callback():
        psim.PushItemWidth(100)
        nonlocal prevFrame, currFrame, isPlay, transparency
        changed, currFrame = psim.SliderInt("frame", currFrame, 0, maxFrame)
        changed, transparency = psim.SliderFloat(
            "transparency", transparency, 0, 1)
        if psim.Button("Play/Pause"):
            isPlay = ~isPlay
        if (prevFrame != currFrame | changed):
            show(trajNc)
        if isPlay:
            currFrame = currFrame + 1
            if currFrame >= maxFrame:
                currFrame = 0

    ps.set_up_dir("z_up")
    ps.set_user_callback(callback)
    ps.show()
