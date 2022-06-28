import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns

import netCDF4 as nc

import pymem3dg as dg
import pymem3dg.read as dg_read
import pymem3dg.util as dg_util

import polyscope as ps
import polyscope.imgui as psim


def matplotlibStyle(s=6, m=8, l=10):
    plt.rcParams["font.sans-serif"] = "Arial"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["lines.linewidth"] = 2
    plt.rcParams["savefig.dpi"] = 600
    # mpl.rcParams.update({'font.size': 8})
    plt.rc("font", size=l)  # controls default text sizes
    plt.rc("axes", titlesize=l)  # fontsize of the axes title
    plt.rc("axes", labelsize=m)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=m)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=m)  # fontsize of the tick labels
    plt.rc("legend", fontsize=s, frameon=False)  # legend fontsize
    plt.rc("figure", titlesize=l)  # fontsize of the figure title
    plt.rc("pdf", fonttype=42)


def polyscopeStyle(setLengthScale=False):
    ps.set_transparency_mode("pretty")
    ps.set_up_dir("z_up")
    ps.set_ground_plane_mode("none")
    ps.set_autocenter_structures(True)
    ps.set_autoscale_structures(False)
    ps.set_view_projection_mode("orthographic")
    ps.set_give_focus_on_show(True)

    if setLengthScale:
        ps.set_automatically_compute_scene_extents(False)
        ps.set_length_scale(1.0)
        low = np.array((-1, -1.0, -1.0))
        high = np.array((1.0, 1.0, 1.0))
        ps.set_bounding_box(low, high)


def getColorMap(fig, range, orientation="horizontal", label="", num_points=100):
    x = np.linspace(0, 1, num=num_points)
    y = np.linspace(0, 1, num=num_points)
    fig.gca().set_visible(False)
    Temp = np.linspace(range[0], range[1], num=num_points)
    plt.scatter(x, y, c=Temp)
    cbar = plt.colorbar(orientation=orientation)
    cbar.set_label(label)
    plt.tight_layout()
    return fig


def plotScreenshot(ax, pngFile, clim, label, orientation="horizontal"):
    img = mpimg.imread(pngFile)

    ax.axis("off")

    # inchPerPixel = 1 / plt.rcParams["figure.dpi"]
    # fig.set_size_inches(492 * inchPerPixel, 276 * inchPerPixel)

    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)

    plt.imshow(img)
    plt.clim(clim[0], clim[1])
    cbar = plt.colorbar(fraction=0.046, pad=0.04, orientation=orientation)
    # cbar.formatter.set_powerlimits((0, 0))
    cbar.set_label(label)
    # plt.title("$t=$"+"{time: .1f}".format(time = tSave * frame))

    plt.tight_layout()


def overlayColormap():
    frameList = np.arange(0, 100, 1)
    for frame in frameList:
        matplotlibStyle()
        fig, ax = plt.subplots(1, 1)
        fig.set_size_inches(4, 3)

        png = "video/frame{}.png".format(dg_util.zeroPadding(frame))
        figureName = "colormapped/frame{}.png".format(frame)

        plotScreenshot(ax, png, [0, 1], "$\phi$", orientation="vertical")
        fig.savefig(figureName, transparent=True, dpi=1500)


def visualizePly(plyFile, **kwargs):
    face, vertex = dg.readMesh(plyFile)
    # print(dg.readData(ply))
    # print(dg.readData(ply, 'vertex'))

    ps.init()
    polyscopeStyle()

    ps_mesh = ps.register_surface_mesh("mesh", vertex, face, smooth_shade=True)
    if kwargs.get("meanCurvature"):
        meanCurvature = dg.readData(plyFile, "vertex", "meanCurvature")
        ps_mesh.add_scalar_quantity("meanCurvature", meanCurvature, enabled=True)
    if kwargs.get("bendingForce"):
        bendingForce = dg.readData(plyFile, "vertex", "bendingForce")
        ps_mesh.add_scalar_quantity("bendingForce", bendingForce, enabled=True)
    # ps.show()


def plotProteinDensity(ax, trajFile, parameters, **kwargs):
    if "frames" in kwargs:
        ncFrames = kwargs.get("frames")
    else:
        ncFrames = range(dg_read.sizeOf(trajFile))
    frameNum = np.size(ncFrames)
    proteinDensity = np.array([])
    proteinDensity_time = np.array([])
    area_weight = np.array([])
    totalProtein = np.zeros(frameNum)
    time = np.zeros(frameNum)

    for plotFrame in range(frameNum):
        ncFrame = ncFrames[plotFrame]
        with nc.Dataset(trajFile) as ds:
            proteinDensity_now = np.array(
                ds.groups["Trajectory"].variables["proteindensity"][ncFrame]
            )
            proteinDensity = np.append(proteinDensity, proteinDensity_now)
            time_now = (
                np.ones(np.size(proteinDensity_now))
                * ds.groups["Trajectory"].variables["time"][ncFrame]
            )
            proteinDensity_time = np.append(proteinDensity_time, time_now)
            time[plotFrame] = ds.groups["Trajectory"].variables["time"][ncFrame]
        system = dg.System(trajFile, ncFrame, parameters)
        system.initialize(nMutation=0, ifMute=True)
        area_weight_now = system.getVertexDualAreas()
        area_weight = np.append(area_weight, area_weight_now)
        totalProtein[plotFrame] = np.sum(area_weight_now * proteinDensity_now)

    ax.set_xlabel("$t$")
    ax.set_ylabel("$\phi$", color="tab:blue")
    cmap = sns.cubehelix_palette(start=0, light=1, as_cmap=True)
    # sns.kdeplot(
    #     x=proteinDensity_time,
    #     y=proteinDensity,
    #     n_levels=5,
    #     cmap=cmap,
    #     shade=True,
    #     ax=ax,
    #     weight=area_weight,
    #     # cbar=True,
    # )
    sns.histplot(
        x=proteinDensity_time,
        y=proteinDensity,
        bins=100,
        pthresh=0.01,
        pmax=0.99,
        weights=area_weight,
        # cbar=True,
        # color='black',
        ax=ax,
        # kde=True,
    )
    ax.tick_params(axis="y", labelcolor="tab:blue")
    ax.set_ylim([0, 1])

    ax2 = ax.twinx()
    ax2.set_ylim([0, 2])
    ax2.set_ylabel(
        "$\int \phi ~dA$", color="grey"
    )  # we already handled the x-label with ax1
    ax2.plot(time, totalProtein, ":", color="grey")
    ax2.tick_params(axis="y", labelcolor="grey")

    # sns.jointplot(proteinDensity_time, proteinDensity, kind="kde", ax=secAx)


def plotChemicalPotentials(ax, trajFile, parameters, logScale=False, **kwargs):
    if "frames" in kwargs:
        ncFrames = kwargs.get("frames")
    else:
        ncFrames = range(dg_read.sizeOf(trajFile))
    frameNum = np.size(ncFrames)
    time = np.zeros(frameNum)

    bendingPotential = np.zeros(frameNum)
    deviatoricPotential = np.zeros(frameNum)
    dirichletPotential = np.zeros(frameNum)
    adsorptionPotential = np.zeros(frameNum)
    aggregationPotential = np.zeros(frameNum)
    entropyPotential = np.zeros(frameNum)
    chemicalPotential = np.zeros(frameNum)

    for plotFrame in range(frameNum):
        ncFrame = ncFrames[plotFrame]
        system = dg.System(trajFile, ncFrame, parameters)
        system.initialize(nMutation=0, ifMute=True)
        time[plotFrame] = system.time
        system.computePhysicalForcing()
        bendingPotential[plotFrame] = (
            np.linalg.norm(system.getForces().getBendingPotential()) ** 2
        )
        deviatoricPotential[plotFrame] = (
            np.linalg.norm(system.getForces().getDeviatoricPotential()) ** 2
        )
        dirichletPotential[plotFrame] = (
            np.linalg.norm(system.getForces().getDirichletPotential()) ** 2
        )
        adsorptionPotential[plotFrame] = (
            np.linalg.norm(system.getForces().getAdsorptionPotential()) ** 2
        )
        aggregationPotential[plotFrame] = (
            np.linalg.norm(system.getForces().getAggregationPotential()) ** 2
        )
        entropyPotential[plotFrame] = (
            np.linalg.norm(system.getForces().getEntropyPotential()) ** 2
        )
        chemicalPotential[plotFrame] = (
            np.linalg.norm(system.getForces().getChemicalPotential()) ** 2
        )
    ms = 4  # marker size
    markerEvery = round(frameNum / 10)
    ax.set_xlabel("$t$")
    ax.set_ylabel("$\|\mu\|^2 [k_BT]$")

    def actualPlotting(ax_local):
        ax_local.plot(
            time,
            dg_util.femtoJToKBT(chemicalPotential, 296),
            "-x",
            label="$\mu$",
            markersize=ms,
            markevery=markerEvery,
        )
        if kwargs.get("bendingPotential"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(bendingPotential, 296),
                ":x",
                label="$\mu_b$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("adsorptionPotential"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(adsorptionPotential, 296),
                ":x",
                label="$\mu_a$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("dirichletPotential"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(dirichletPotential, 296),
                ":x",
                label="$\mu_d$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("deviatoricPotential"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(deviatoricPotential, 296),
                ":x",
                label="$\mu_{dev}$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("aggregationPotential"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(aggregationPotential, 296),
                ":x",
                label="$\mu_{agg}$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("entropyPotential"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(entropyPotential, 296),
                ":x",
                label="$\mu_{entr}$",
                markersize=ms,
                markevery=markerEvery,
            )

    actualPlotting(ax)
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    if logScale:
        ax.set_yscale("log")
    # axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
    # actualPlotting(axins)
    # axins.set_ylim(0, 10)
    # ax.indicate_inset_zoom(axins, edgecolor="grey")
    # axins.set_xticklabels([])


def plotMechanicalForces(ax, trajFile, parameters, logScale=False, **kwargs):
    if "frames" in kwargs:
        ncFrames = kwargs.get("frames")
    else:
        ncFrames = range(dg_read.sizeOf(trajFile))
    frameNum = np.size(ncFrames)
    time = np.zeros(frameNum)

    bendingForce = np.zeros(frameNum)
    capillaryForce = np.zeros(frameNum)
    externalForce = np.zeros(frameNum)
    osmoticForce = np.zeros(frameNum)
    adsorptionForce = np.zeros(frameNum)
    aggregationForce = np.zeros(frameNum)
    entropyForce = np.zeros(frameNum)
    springForce = np.zeros(frameNum)
    lineCapillaryForce = np.zeros(frameNum)
    mechanicalForce = np.zeros(frameNum)

    for plotFrame in range(frameNum):
        ncFrame = ncFrames[plotFrame]
        system = dg.System(trajFile, ncFrame, parameters)
        system.initialize(nMutation=0, ifMute=True)
        time[plotFrame] = system.time
        system.computePhysicalForcing()
        bendingForce[plotFrame] = (
            np.linalg.norm(system.getForces().getBendingForceVec()) ** 2
        )
        capillaryForce[plotFrame] = (
            np.linalg.norm(system.getForces().getCapillaryForceVec()) ** 2
        )
        osmoticForce[plotFrame] = (
            np.linalg.norm(system.getForces().getOsmoticForceVec()) ** 2
        )
        adsorptionForce[plotFrame] = (
            np.linalg.norm(system.getForces().getAdsorptionForceVec()) ** 2
        )
        aggregationForce[plotFrame] = (
            np.linalg.norm(system.getForces().getAggregationForceVec()) ** 2
        )
        externalForce[plotFrame] = (
            np.linalg.norm(system.getForces().getExternalForceVec()) ** 2
        )
        entropyForce[plotFrame] = (
            np.linalg.norm(system.getForces().getEntropyForceVec()) ** 2
        )
        springForce[plotFrame] = (
            np.linalg.norm(system.getForces().getSpringForceVec()) ** 2
        )
        lineCapillaryForce[plotFrame] = (
            np.linalg.norm(system.getForces().getLineCapillaryForceVec()) ** 2
        )
        mechanicalForce[plotFrame] = (
            np.linalg.norm(system.getForces().getMechanicalForceVec()) ** 2
        )

    ms = 4  # marker size
    markerEvery = round(frameNum / 10)
    ax.set_xlabel("$t$")
    ax.set_ylabel("$\|f\|^2 [k_BT]$")

    def actualPlotting(ax_local):
        ax_local.plot(
            time,
            dg_util.femtoJToKBT(mechanicalForce, 296),
            "-o",
            label="$f$",
            markersize=ms,
            markevery=markerEvery,
        )
        if kwargs.get("bendingForce"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(bendingForce, 296),
                ":o",
                label="$f_b$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("capillaryForce"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(capillaryForce, 296),
                ":o",
                label="$f_s$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("osmoticForce"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(osmoticForce, 296),
                ":o",
                label="$f_p$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("adsorptionForce"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(adsorptionForce, 296),
                ":o",
                label="$f_a$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("lineCapillaryForce"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(lineCapillaryForce, 296),
                ":o",
                label="$f_d$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("aggregationForce"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(aggregationForce, 296),
                ":o",
                label="$f_{agg}$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("externalForce"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(externalForce, 296),
                ":o",
                label="$f_{ext}$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("entropyForce"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(entropyForce, 296),
                ":o",
                label="$f_{entr}$",
                markersize=ms,
                markevery=markerEvery,
            )
        if kwargs.get("springForce"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(springForce, 296),
                ":o",
                label="$f_{spr}$",
                markersize=ms,
                markevery=markerEvery,
            )

    actualPlotting(ax)
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    if logScale:
        ax.set_yscale("log")
    # axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
    # actualPlotting(axins)
    # axins.set_ylim(0, 50)
    # ax.indicate_inset_zoom(axins, edgecolor="grey")
    # axins.set_xticklabels([])


def plotEnergy(ax, trajFile, parameters, zeroing=False, logScale=False, **kwargs):
    if "frames" in kwargs:
        ncFrames = kwargs.get("frames")
    else:
        ncFrames = range(dg_read.sizeOf(trajFile))
    frameNum = np.size(ncFrames)
    time = np.zeros(frameNum)

    kineticEnergy = np.zeros(frameNum)
    potentialEnergy = np.zeros(frameNum)
    externalWork = np.zeros(frameNum)
    totalEnergy = np.zeros(frameNum)
    bendingEnergy = np.zeros(frameNum)
    deviatoricEnergy = np.zeros(frameNum)
    surfaceEnergy = np.zeros(frameNum)
    pressureEnergy = np.zeros(frameNum)
    adsorptionEnergy = np.zeros(frameNum)
    aggregationEnergy = np.zeros(frameNum)
    entropyEnergy = np.zeros(frameNum)
    edgeSpringEnergy = np.zeros(frameNum)
    faceSpringEnergy = np.zeros(frameNum)
    lcrSpringEnergy = np.zeros(frameNum)
    dirichletEnergy = np.zeros(frameNum)

    for plotFrame in range(frameNum):
        ncFrame = ncFrames[plotFrame]
        system = dg.System(trajFile, ncFrame, parameters)
        system.initialize(nMutation=0, ifMute=True)
        time[plotFrame] = system.time
        system.computeTotalEnergy()
        energy = system.getEnergy()
        kineticEnergy[plotFrame] = energy.kineticEnergy
        potentialEnergy[plotFrame] = energy.potentialEnergy
        bendingEnergy[plotFrame] = energy.bendingEnergy
        deviatoricEnergy[plotFrame] = energy.deviatoricEnergy
        surfaceEnergy[plotFrame] = energy.surfaceEnergy
        pressureEnergy[plotFrame] = energy.pressureEnergy
        adsorptionEnergy[plotFrame] = energy.adsorptionEnergy
        aggregationEnergy[plotFrame] = energy.aggregationEnergy
        entropyEnergy[plotFrame] = energy.entropyEnergy
        edgeSpringEnergy[plotFrame] = energy.edgeSpringEnergy
        lcrSpringEnergy[plotFrame] = energy.lcrSpringEnergy
        faceSpringEnergy[plotFrame] = energy.faceSpringEnergy
        dirichletEnergy[plotFrame] = energy.dirichletEnergy
        if plotFrame != 0:
            externalWork[plotFrame] = externalWork[
                plotFrame - 1
            ] + system.computeIntegratedPower(time[plotFrame] - time[plotFrame - 1])
    if zeroing:
        kineticEnergy = kineticEnergy - kineticEnergy[-1]
        potentialEnergy = potentialEnergy - potentialEnergy[-1]
        bendingEnergy = bendingEnergy - bendingEnergy[-1]
        deviatoricEnergy = deviatoricEnergy - deviatoricEnergy[-1]
        surfaceEnergy = surfaceEnergy - surfaceEnergy[-1]
        pressureEnergy = pressureEnergy - pressureEnergy[-1]
        adsorptionEnergy = adsorptionEnergy - adsorptionEnergy[-1]
        aggregationEnergy = aggregationEnergy - aggregationEnergy[-1]
        entropyEnergy = entropyEnergy - entropyEnergy[-1]
        edgeSpringEnergy = edgeSpringEnergy - edgeSpringEnergy[-1]
        lcrSpringEnergy = lcrSpringEnergy - lcrSpringEnergy[-1]
        faceSpringEnergy = faceSpringEnergy - faceSpringEnergy[-1]
        dirichletEnergy = dirichletEnergy - dirichletEnergy[-1]

    totalEnergy = potentialEnergy + kineticEnergy - externalWork

    ax.set_xlabel("$t$")
    ax.set_ylabel("$E[k_BT]$")

    def actualPlotting(ax_local):
        if kwargs.get("kineticEnergy"):
            ax_local.plot(time, dg_util.femtoJToKBT(kineticEnergy, 296), label="$K$")
        if kwargs.get("potentialEnergy"):
            ax_local.plot(time, dg_util.femtoJToKBT(potentialEnergy, 296), label="$E$")
        if kwargs.get("externalWork"):
            ax_local.plot(time, dg_util.femtoJToKBT(externalWork, 296), label="$W$")
        if kwargs.get("totalEnergy"):
            ax_local.plot(time, dg_util.femtoJToKBT(totalEnergy, 296), label="$W$")
        if kwargs.get("bendingEnergy"):
            ax_local.plot(
                time, dg_util.femtoJToKBT(bendingEnergy, 296), ":", label="$E_b$"
            )
        if kwargs.get("surfaceEnergy"):
            ax_local.plot(
                time, dg_util.femtoJToKBT(surfaceEnergy, 296), ":", label="$E_s$"
            )
        if kwargs.get("pressureEnergy"):
            ax_local.plot(
                time, dg_util.femtoJToKBT(pressureEnergy, 296), ":", label="$E_p$"
            )
        if kwargs.get("adsorptionEnergy"):
            ax_local.plot(
                time, dg_util.femtoJToKBT(adsorptionEnergy, 296), ":", label="$E_a$"
            )
        if kwargs.get("dirichletEnergy"):
            ax_local.plot(
                time, dg_util.femtoJToKBT(dirichletEnergy, 296), ":", label="$E_{d}$"
            )
        if kwargs.get("deviatoricEnergy"):
            ax_local.plot(
                time, dg_util.femtoJToKBT(deviatoricEnergy, 296), ":", label="$E_{dev}$"
            )
        if kwargs.get("aggregationEnergy"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(aggregationEnergy, 296),
                ":",
                label="$E_{agg}$",
            )
        if kwargs.get("entropyEnergy"):
            ax_local.plot(
                time, dg_util.femtoJToKBT(entropyEnergy, 296), ":", label="$E_{entr}$"
            )
        if kwargs.get("edgeSpringEnergy"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(edgeSpringEnergy, 296),
                ":",
                label="$E_{edge}$",
            )
        if kwargs.get("faceSpringEnergy"):
            ax_local.plot(
                time,
                dg_util.femtoJToKBT(faceSpringEnergy, 296),
                ":",
                label="$E_{face}$",
            )
        if kwargs.get("lcrSpringEnergy"):
            ax_local.plot(
                time, dg_util.femtoJToKBT(lcrSpringEnergy, 296), ":", label="$E_{lcr}$"
            )

    actualPlotting(ax)
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    if logScale:
        ax.set_yscale("log")
    # axins = ax.inset_axes([0.5, 0.5, 0.43, 0.6])
    # actualPlotting(axins)
    # axins.set_xlim(0, 0.1)
    # ax.indicate_inset_zoom(axins, edgecolor="grey")
    # axins.set_xticklabels([])


def interlaceImages(fileList, videoName):
    """example: fileList = [frame1.png, frame2.png]"""
    from moviepy.editor import ImageSequenceClip
    clip = ImageSequenceClip(fileList, fps=20)
    # clip.write_videofile("video.mp4", fps=24)
    clip.write_gif(videoName + ".gif")

def buildSystemAndAddQuantities(psMesh, trajNc, frame, isPointwiseValue=False, **kwargs):
    hasParameters = "parameters" in kwargs

    # Build system
    if hasParameters:
        parameters = kwargs.get("parameters")
        system = dg.System(trajNc, frame, parameters)
        system.initialize(nMutation=0, ifMute=True)
        system.computePhysicalForcing()
    else:
        system = dg.System(trajNc, frame)
        system.initialize(nMutation=0, ifMute=True)
    polyscopePermutations = system.getPolyscopePermutations()
    setPolyscopePermutations(psMesh, polyscopePermutations)

    # Add Quantities
    vertexDualAreas = system.getVertexDualAreas()
    proteinDensity = system.getProteinDensity()
    if isPointwiseValue:
        psMesh.add_scalar_quantity(
                "proteinDensity", proteinDensity, vminmax=(-1, 1), cmap = 'coolwarm', enabled=True
            )
    else:
        proteinDensity = proteinDensity * vertexDualAreas
        psMesh.add_scalar_quantity(
                "proteinDensity", proteinDensity, enabled=True, cmap = 'viridis',
            )
    if kwargs.get("meanCurvature"):
        meanCurvature = system.getVertexMeanCurvatures()
        if isPointwiseValue:
            meanCurvature = meanCurvature / vertexDualAreas
        psMesh.add_scalar_quantity("meanCurvature", meanCurvature)
    if kwargs.get("gaussianCurvature"):
        gaussianCurvature = system.getVertexGaussianCurvatures()
        if isPointwiseValue:
            gaussianCurvature = gaussianCurvature / vertexDualAreas
        psMesh.add_scalar_quantity("gaussianCurvature", gaussianCurvature)
    if hasParameters:
        if kwargs.get("mechanicalForce"):
            mechanicalForce = system.getForces().getMechanicalForceVec()
            if isPointwiseValue:
                mechanicalForce = dg_util.rowwiseScaling(
                    1 / vertexDualAreas, mechanicalForce
                )
            psMesh.add_vector_quantity("mechanicalForce", mechanicalForce)
            psMesh.add_scalar_quantity(
                "<mechanicalForce,n>",
                dg_util.rowwiseDotProduct(mechanicalForce, system.getVertexNormals()),
            )
        if kwargs.get("bendingForce"):
            bendingForce = system.getForces().getBendingForceVec()
            if isPointwiseValue:
                bendingForce = dg_util.rowwiseScaling(1 / vertexDualAreas, bendingForce)
            psMesh.add_vector_quantity("bendingForce", bendingForce)
            psMesh.add_scalar_quantity(
                "<bendingForce,n>",
                dg_util.rowwiseDotProduct(bendingForce, system.getVertexNormals()),
            )
        if kwargs.get("externalForce"):
            externalForce = system.getForces().getExternalForceVec()
            if isPointwiseValue:
                externalForce = dg_util.rowwiseScaling(
                    1 / vertexDualAreas, externalForce
                )
            psMesh.add_vector_quantity("externalForce", externalForce)
            psMesh.add_scalar_quantity(
                "<externalForce,n>",
                dg_util.rowwiseDotProduct(externalForce, system.getVertexNormals()),
            )
        if kwargs.get("capillaryForce"):
            capillaryForce = system.getForces().getCapillaryForceVec()
            if isPointwiseValue:
                capillaryForce = dg_util.rowwiseScaling(
                    1 / vertexDualAreas, capillaryForce
                )
            psMesh.add_vector_quantity("capillaryForce", capillaryForce)
            psMesh.add_scalar_quantity(
                "<capillaryForce,n>",
                dg_util.rowwiseDotProduct(capillaryForce, system.getVertexNormals()),
            )
        if kwargs.get("lineCapillaryForce"):
            lineCapillaryForce = system.getForces().getLineCapillaryForceVec()
            if isPointwiseValue:
                lineCapillaryForce = dg_util.rowwiseScaling(
                    1 / vertexDualAreas, lineCapillaryForce
                )
            psMesh.add_vector_quantity("lineCapillaryForce", lineCapillaryForce)
            psMesh.add_scalar_quantity(
                "<lineCapillaryForce,n>",
                dg_util.rowwiseDotProduct(
                    lineCapillaryForce, system.getVertexNormals()
                ),
            )
        if kwargs.get("osmoticForce"):
            osmoticForce = system.getForces().getOsmoticForceVec()
            if isPointwiseValue:
                osmoticForce = dg_util.rowwiseScaling(1 / vertexDualAreas, osmoticForce)
            psMesh.add_vector_quantity("osmoticForce", osmoticForce)
            psMesh.add_scalar_quantity(
                "<osmoticForce,n>",
                dg_util.rowwiseDotProduct(osmoticForce, system.getVertexNormals()),
            )
        if kwargs.get("adsorptionForce"):
            adsorptionForce = system.getForces().getAdsorptionForceVec()
            if isPointwiseValue:
                adsorptionForce = dg_util.rowwiseScaling(
                    1 / vertexDualAreas, adsorptionForce
                )
            psMesh.add_vector_quantity("adsorptionForce", adsorptionForce)
            psMesh.add_scalar_quantity(
                "<adsorptionForce,n>",
                dg_util.rowwiseDotProduct(adsorptionForce, system.getVertexNormals()),
            )
        if kwargs.get("aggregationForce"):
            aggregationForce = system.getForces().getAggregationForceVec()
            if isPointwiseValue:
                aggregationForce = dg_util.rowwiseScaling(
                    1 / vertexDualAreas, aggregationForce
                )
            psMesh.add_vector_quantity("aggregationForce", aggregationForce)
            psMesh.add_scalar_quantity(
                "<aggregationForce,n>",
                dg_util.rowwiseDotProduct(aggregationForce, system.getVertexNormals()),
            )
        if kwargs.get("entropyForce"):
            entropyForce = system.getForces().getEntropyForceVec()
            if isPointwiseValue:
                entropyForce = dg_util.rowwiseScaling(1 / vertexDualAreas, entropyForce)
            psMesh.add_vector_quantity("entropyForce", entropyForce)
            psMesh.add_scalar_quantity(
                "<entropyForce,n>",
                dg_util.rowwiseDotProduct(entropyForce, system.getVertexNormals()),
            )
        if kwargs.get("springForce"):
            springForce = system.getForces().getSpringForceVec()
            if isPointwiseValue:
                springForce = dg_util.rowwiseScaling(1 / vertexDualAreas, springForce)
            psMesh.add_vector_quantity("springForce", springForce)
            psMesh.add_scalar_quantity(
                "<springForce,n>",
                dg_util.rowwiseDotProduct(springForce, system.getVertexNormals()),
            )
        if kwargs.get("chemicalPotential"):
            chemicalPotential = system.getForces().getChemicalPotential()
            psMesh.add_scalar_quantity("chemicalPotential", chemicalPotential)
        if kwargs.get("bendingPotential"):
            bendingPotential = system.getForces().getBendingPotential()
            psMesh.add_scalar_quantity("bendingPotential", bendingPotential)
        if kwargs.get("aggregationPotential"):
            aggregationPotential = system.getForces().getAggregationPotential()
            psMesh.add_scalar_quantity("aggregationPotential", aggregationPotential)
        if kwargs.get("dirichletPotential"):
            dirichletPotential = system.getForces().getDirichletPotential()
            psMesh.add_scalar_quantity("dirichletPotential", dirichletPotential)
        if kwargs.get("adsorptionPotential"):
            adsorptionPotential = system.getForces().getAdsorptionPotential()
            psMesh.add_scalar_quantity("adsorptionPotential", adsorptionPotential)
        if kwargs.get("entropyPotential"):
            entropyPotential = system.getForces().getEntropyPotential()
            psMesh.add_scalar_quantity("entropyPotential", entropyPotential)
        if kwargs.get("inPlaneFluxForm"):
            inPlaneFluxForm = system.getForces().getInPlaneFluxForm()
            psMesh.add_one_form_vector_quantity(
                "inPlaneFluxForm",
                inPlaneFluxForm,
                system.getPolyscopeEdgeOrientations(),
            )


def setPolyscopePermutations(psmesh, polyscopePermutations):
    psmesh.set_all_permutations(
        vertex_perm=np.array(polyscopePermutations[0][0]),
        vertex_perm_size=polyscopePermutations[0][1],
        face_perm=np.array(polyscopePermutations[1][0]),
        face_perm_size=polyscopePermutations[1][1],
        edge_perm=np.array(polyscopePermutations[2][0]),
        edge_perm_size=polyscopePermutations[2][1],
        corner_perm=np.array(polyscopePermutations[4][0]),
        corner_perm_size=polyscopePermutations[4][1],
        halfedge_perm=np.array(polyscopePermutations[3][0]),
        halfedge_perm_size=polyscopePermutations[3][1],
    )


def animate(trajNc, **kwargs):
    buildSystem = kwargs is not None
    hasParameters = "parameters" in kwargs
    if hasParameters:
        parameters = kwargs.get("parameters")

    if "frames" in kwargs:
        frames = kwargs.get("frames")
    else:
        frames = range(dg_read.sizeOf(trajNc))

    maxFrameInd = np.size(frames) - 1
    prevFrameInd = 0
    currFrameInd = 0
    time = dg_read.readMeshDataByNc(
        trajNc, frames[currFrameInd], "Trajectory", "time", 1
    )
    isPointwiseValue = False
    recordingDir = ""

    ps.init()
    polyscopeStyle()

    isPlay = False
    transparency = 1
    isRecord = False

    if buildSystem:
        if hasParameters:
            system = dg.System(trajNc, frames[currFrameInd], parameters)
            system.initialize(nMutation=0, ifMute=True)
        else:
            system = dg.System(trajNc, frames[currFrameInd])
            system.initialize(nMutation=0, ifMute=True)

    def show(trajNc):
        nonlocal currFrameInd, system, time, isPointwiseValue
        time = dg_read.readMeshDataByNc(
            trajNc, frames[currFrameInd], "Trajectory", "time", 1
        )
        face, vertex = dg_read.readMeshByNc(trajNc, frames[currFrameInd])
        proteinDensity = dg_read.readMeshDataByNc(
            trajNc, frames[currFrameInd], "Trajectory", "proteindensity", 1
        )
        velocity = dg_read.readMeshDataByNc(
            trajNc, frames[currFrameInd], "Trajectory", "velocities", 3
        )
        psmesh = ps.register_surface_mesh(
            "mesh", vertex, face, transparency=transparency, smooth_shade=True
        )
        psmesh.add_scalar_quantity(
            "proteinDensity", proteinDensity, vminmax=(-1, 1), cmap = 'coolwarm', enabled=True
        )
        psmesh.add_vector_quantity("velocity", velocity)
        if buildSystem:
            buildSystemAndAddQuantities(
                psmesh,
                trajNc,
                frames[currFrameInd],
                isPointwiseValue=isPointwiseValue,
                **kwargs
            )

    def callback():
        psim.PushItemWidth(100)
        nonlocal prevFrameInd, currFrameInd, time, isPlay, isRecord, transparency, isPointwiseValue, recordingDir
    
        if psim.Button("Play/Pause"):
            isPlay = ~isPlay
        psim.SameLine()
        changed, currFrameInd = psim.SliderInt("", currFrameInd, 0, maxFrameInd)
        psim.SameLine()
        psim.TextUnformatted(
            "Frames: {:d}, time: ".format(frames[currFrameInd])
            + "{0:.3g}".format(time)
        )
        
        if psim.Button("Record/Stop"):
            isRecord = ~isRecord
            isPlay = ~isPlay
        psim.SameLine() 
        changed, recordingDir = psim.InputText("Directory", recordingDir)
        if isRecord:
            if ~isPlay:
                isPlay = True
            defaultName = "frame{}.png".format(currFrameInd)
            ps.screenshot(recordingDir + defaultName)
        if isPlay:
            currFrameInd = currFrameInd + 1
            if currFrameInd >= maxFrameInd:
                currFrameInd = maxFrameInd
                isPlay = ~isPlay
                
        changed, transparency = psim.SliderFloat("Transparency", transparency, 0, 1)
        psim.SameLine()
        changed, isPointwiseValue = psim.Checkbox(
            "Pointwise (vs. integrated)", isPointwiseValue
        )
        
        if (prevFrameInd != currFrameInd) or changed:
            show(trajNc)
            changed = False
            prevFrameInd = currFrameInd
            
    show(trajNc)
    ps.set_user_callback(callback)
    # ps.show()
    