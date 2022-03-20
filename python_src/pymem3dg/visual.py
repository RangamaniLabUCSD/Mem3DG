import numpy as np
import matplotlib.pyplot as plt

import pymem3dg as dg
import pymem3dg.read as dg_read

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


def animate():
    trajNc = "../../results/temp/traj.nc"
    ps.init()
    prevFrame = 0
    currFrame = 0
    isPlay = False
    transparency = 0.4
    maxFrame = dg_read.sizeOf(trajNc) - 1
    face, vertex = dg_read.readMeshByNc(trajNc, currFrame)
    mesh = ps.register_surface_mesh("mesh", vertex, face)
    ps.set_transparency_mode("pretty")

    def callback():
        psim.PushItemWidth(100)
        nonlocal prevFrame, currFrame, isPlay, face, vertex, mesh, transparency
        changed, currFrame = psim.SliderInt("frame", currFrame, 0, maxFrame)
        changed, transparency = psim.SliderFloat(
            "transparency", transparency, 0, 1)
        if psim.Button("Play/Pause"):
            isPlay = ~isPlay
        if (prevFrame != currFrame | changed):
            face, vertex = dg_read.readMeshByNc(trajNc, currFrame)
            mesh = ps.register_surface_mesh(
                "mesh", vertex, face, transparency=transparency
            )
        if isPlay:
            currFrame = currFrame + 1
            if currFrame >= maxFrame:
                currFrame = 0

    ps.set_up_dir("z_up")
    ps.set_user_callback(callback)
    ps.show()
