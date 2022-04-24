import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import pymem3dg as dg
import pymem3dg.read as dg_read
import pymem3dg.util as dg_util

import polyscope as ps
import polyscope.imgui as psim


def matplotlibStyle(s=6, m=8, l=10):
    """ Formatting style of matplotlib """
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.family'] = "sans-serif"
    # mpl.rcParams.update({'font.size': 8})
    SMALL_SIZE = s
    MEDIUM_SIZE = m
    BIGGER_SIZE = l
    plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    plt.rc('pdf', fonttype=42)


def getColorMap(fig, range, orientation="horizontal", label="",  num_points=100):
    x = np.linspace(0, 1, num=num_points)
    y = np.linspace(0, 1, num=num_points)
    fig.gca().set_visible(False)
    Temp = np.linspace(range[0], range[1], num=num_points)
    plt.scatter(x, y, c=Temp)
    cbar = plt.colorbar(orientation=orientation)
    cbar.set_label(label)
    plt.tight_layout()
    return fig


def overlayColorMap(pngFile, clim, label, orientation="horizontal"):
    img = mpimg.imread(pngFile)

    fig, ax = plt.subplots(1)
    matplotlibStyle()
    inchPerPixel = 1 / plt.rcParams['figure.dpi']
    fig.set_size_inches(492 * inchPerPixel, 276 * inchPerPixel)
    ax.axis("off")

    plt.imshow(img)
    plt.clim(clim[0], clim[1])
    cbar = plt.colorbar(fraction=0.046, pad=0.04, orientation=orientation)
    # cbar.formatter.set_powerlimits((0, 0))
    cbar.set_label(label)
    # plt.title("$t=$"+"{time: .1f}".format(time = tSave * frame))

    plt.tight_layout()
    return fig


def buildSystemAndAddQuantities(psMesh, face, vertex, proteinDensity, velocity, **kwargs):
    hasParameters = "parameters" in kwargs
    if hasParameters:
        parameters = kwargs.get("parameters")
    if hasParameters:
        system = dg.System(face, vertex, proteinDensity, velocity, parameters)
        system.initialize(nMutation=0, ifMute=True)
        system.computePhysicalForcing()
        if (kwargs.get('mechanicalForce')):
            mechanicalForce = system.getForces().getMechanicalForceVec()
            psMesh.add_vector_quantity(
                "mechanicalForce", mechanicalForce)
            psMesh.add_scalar_quantity("<mechanicalForce,n>", dg_util.rowwiseDotProduct(
                mechanicalForce, system.getVertexNormals()))
        if (kwargs.get('bendingForce')):
            bendingForce = system.getForces().getBendingForceVec()
            psMesh.add_vector_quantity("bendingForce", bendingForce)
            psMesh.add_scalar_quantity("<bendingForce,n>", dg_util.rowwiseDotProduct(
                bendingForce, system.getVertexNormals()))
        if (kwargs.get('externalForce')):
            externalForce = system.getForces().getExternalForceVec()
            psMesh.add_vector_quantity("externalForce", externalForce)
            psMesh.add_scalar_quantity("<externalForce,n>", dg_util.rowwiseDotProduct(
                externalForce, system.getVertexNormals()))
        if (kwargs.get('capillaryForce')):
            capillaryForce = system.getForces().getCapillaryForceVec()
            psMesh.add_vector_quantity("capillaryForce", capillaryForce)
            psMesh.add_scalar_quantity("<capillaryForce,n>", dg_util.rowwiseDotProduct(
                capillaryForce, system.getVertexNormals()))
        if (kwargs.get('lineCapillaryForce')):
            lineCapillaryForce = system.getForces().getLineCapillaryForceVec()
            psMesh.add_vector_quantity(
                "lineCapillaryForce", lineCapillaryForce)
            psMesh.add_scalar_quantity("<lineCapillaryForce,n>", dg_util.rowwiseDotProduct(
                lineCapillaryForce, system.getVertexNormals()))
        if (kwargs.get('osmoticForce')):
            osmoticForce = system.getForces().getOsmoticForceVec()
            psMesh.add_vector_quantity("osmoticForce", osmoticForce)
            psMesh.add_scalar_quantity("<osmoticForce,n>", dg_util.rowwiseDotProduct(
                osmoticForce, system.getVertexNormals()))
        if (kwargs.get('adsorptionForce')):
            adsorptionForce = system.getForces().getAdsorptionForceVec()
            psMesh.add_vector_quantity(
                "adsorptionForce", adsorptionForce)
            psMesh.add_scalar_quantity("<adsorptionForce,n>", dg_util.rowwiseDotProduct(
                adsorptionForce, system.getVertexNormals()))
        if (kwargs.get('aggregationForce')):
            aggregationForce = system.getForces().getAggregationForceVec()
            psMesh.add_vector_quantity(
                "aggregationForce", aggregationForce)
            psMesh.add_scalar_quantity("<aggregationForce,n>", dg_util.rowwiseDotProduct(
                aggregationForce, system.getVertexNormals()))

        if (kwargs.get('chemicalPotential')):
            chemicalPotential = system.getForces().getChemicalPotential()
            psMesh.add_scalar_quantity(
                "chemicalPotential", chemicalPotential)
        if (kwargs.get('bendingPotential')):
            bendingPotential = system.getForces().getBendingPotential()
            psMesh.add_scalar_quantity(
                "bendingPotential", bendingPotential)
        if (kwargs.get('aggregationPotential')):
            aggregationPotential = system.getForces().getAggregationPotential()
            psMesh.add_scalar_quantity(
                "aggregationPotential", aggregationPotential)
        if (kwargs.get('diffusionPotential')):
            diffusionPotential = system.getForces().getDiffusionPotential()
            psMesh.add_scalar_quantity(
                "diffusionPotential", diffusionPotential)
        if (kwargs.get('adsorptionPotential')):
            adsorptionPotential = system.getForces().getAdsorptionPotential()
            psMesh.add_scalar_quantity(
                "adsorptionPotential", adsorptionPotential)
    else:
        system = dg.System(face, vertex, proteinDensity, velocity)
        system.initialize(nMutation=0, ifMute=True)

    if (kwargs.get('meanCurvature')):
        meanCurvature = system.getVertexMeanCurvatures()/system.getVertexDualAreas()
        psMesh.add_scalar_quantity("meanCurvature", meanCurvature)
    if (kwargs.get('gaussianCurvature')):
        gaussianCurvature = system.getVertexGaussianCurvatures()/system.getVertexDualAreas()
        psMesh.add_scalar_quantity(
            "gaussianCurvature", gaussianCurvature)


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
    transparency = 1
    maxFrame = dg_read.sizeOf(trajNc) - 1

    if buildSystem:
        system = dg.System(trajNc, currFrame)
        system.initialize(nMutation=0, ifMute=True)
        if hasParameters:
            system = dg.System(trajNc, currFrame, parameters)
            system.initialize(nMutation=0, ifMute=True)

    def show(trajNc):
        nonlocal currFrame, system
        face, vertex = dg_read.readMeshByNc(trajNc, currFrame)
        proteinDensity = dg_read.readMeshDataByNc(
            trajNc, currFrame, "Trajectory", "proteindensity", 1)
        velocity = dg_read.readMeshDataByNc(
            trajNc, currFrame, "Trajectory", "velocities", 3)
        mesh = ps.register_surface_mesh(
            "mesh", vertex, face, transparency=transparency)
        limit = (np.min(proteinDensity), np.max(proteinDensity))
        mesh.add_scalar_quantity(
            "proteinDensity", proteinDensity, vminmax=limit)
        mesh.add_vector_quantity(
            "velocity", velocity)
        if buildSystem:
            buildSystemAndAddQuantities(
                mesh, face, vertex, proteinDensity, velocity, **kwargs)

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
            changed = False
            prevFrame = currFrame
        if isPlay:
            currFrame = currFrame + 1
            if currFrame >= maxFrame:
                currFrame = 0

    show(trajNc)
    ps.set_up_dir("z_up")
    ps.set_user_callback(callback)
    ps.show()
