# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Copyright (c) 2020:
#     Laboratory for Computational Cellular Mechanobiology
#     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
#     Christopher T. Lee (ctlee@ucsd.edu)
#     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
#     Padmini Rangamani (prangamani@eng.ucsd.edu)
#

"""PyMem3DG: Membrane Dynamics in 3D using Discrete Differential Geometry

PyMem3DG performs membrane simulations.
"""

import os
import sys
import subprocess
import re


def git_version():
    """Get the version from git describe

    Returns:
        string: version string or None if invalid
    """
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ["SYSTEMROOT", "PATH", "HOME"]:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env["LANGUAGE"] = "C"
        env["LANG"] = "C"
        env["LC_ALL"] = "C"
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, env=env)
        return out

    try:
        out = _minimal_ext_cmd(["git", "describe", "--tags", "--dirty", "--always"])
        GIT_REVISION = out.strip().decode("ascii")
    except (subprocess.SubprocessError, OSError):
        GIT_REVISION = None

    return GIT_REVISION

def standardize_version(version_string):
    """Standardize the version string

    Args:
        version_string (str): input version string
    Returns:
        str: standardized version
    """
    VERSION_PATTERN = r"""
        v?
        (?:
            (?:(?P<epoch>[0-9]+)!)?                           # epoch
            (?P<release>[0-9]+(?:\.[0-9]+)*)                  # release segment
            (?P<pre>                                          # pre-release
                [-_\.]?
                (?P<pre_l>(alpha|beta|a|b|c|rc))
                (?P<pre_n>[0-9]+)?
            )?
            (?P<dev>                                          # dev release
                [-_\.]?
                (?P<dev_l>dev)
                [-_\.]?
                (?P<dev_n>[0-9]+)?
            )?
            (?P<meta>
                [-_\.]?
                (?P<commits_since>[0-9]+)?
                [-_\.]?
                (?P<sha>[a-z0-9]*)?
                [-_\.]?
                (?P<dirty>dirty)?
            )?
        )
    """

    _regex = re.compile(
        r"^\s*" + VERSION_PATTERN + r"\s*$",
        re.VERBOSE | re.IGNORECASE,
    )

    match = _regex.match(version_string)
    if match:
        if match.group("release"):
            version = match.group("release")
            if match.group("pre_l"):
                version += match.group("pre_l")
                if match.group("pre_n"):
                    version += match.group("pre_n")
                else:
                    version += "0"
            if match.group("dev"):
                version += "dev"
                if match.group("dev_n"):
                    version += match.group("dev_n")
                else:
                    version += "0"
    else:
        version = "0.0.0"
    return version

version = git_version()
if version is None:
    # Git describe failed... read version from file
    with open("VERSION", "r") as f:
        version = f.readline()
    version = standardize_version(version)
else:
    version = standardize_version(version)


cmake_args=['-DBUILD_PYMEM3DG=ON', '-DSUITESPARSE=OFF']

if('CONDA_PREFIX' in os.environ):
    print("Setting library search path (CMAKE_PREFIX_PATH): %s"%(os.environ['CONDA_PREFIX']))
    cmake_args.append('-DCMAKE_PREFIX_PATH=%s'%(os.environ['CONDA_PREFIX']))

DOCLINES = __doc__.split("\n")

CLASSIFIERS = """\
Development Status :: 3 - Alpha
Environment :: Console
Intended Audience :: Science/Research
License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)
Natural Language :: English
Operating System :: OS Independent
Programming Language :: C++
Programming Language :: Python :: 3 :: Only
Programming Language :: Python :: Implementation :: CPython
Topic :: Scientific/Engineering :: Chemistry
Topic :: Scientific/Engineering :: Mathematics
Topic :: Scientific/Engineering :: Physics
Topic :: Scientific/Engineering :: Visualization
"""

# If building on readthedocs.io build with unix makefiles
_on_rtd = os.environ.get("READTHEDOCS", None) == "True"
if _on_rtd:
    sys.argv.extend(["-G", "Unix Makefiles"])
    cmake_args.append("-DBUILD_MEM3DG_DOCS=ON")

try:
    from skbuild import setup
except ImportError:
    print('\nERROR: scikit-build is required to build from source.', file=sys.stderr)
    print('Please run:', file=sys.stderr)
    print('', file=sys.stderr)
    print('  python -m pip install scikit-build\n', file=sys.stderr)
    print('  -- or --\n', file=sys.stderr)
    print('  conda install scikit-build', file=sys.stderr)
    sys.exit(1)


from setuptools import find_packages

setup(
    name="pymem3dg",
    version=version,
    maintainer="Cuncheng Zhu and Christopher T. Lee",
    maintainer_email="cuzhu@ucsd.edu, ctlee@ucsd.edu",
    author="The Mem3DG Team",
    author_email="cuzhu@ucsd.edu, ctlee@ucsd.edu",
    url="https://github.com/RangamaniLabUCSD/Mem3DG",
    packages=find_packages(where="python_src"),
    package_dir={"": "python_src"},
    cmake_install_dir="python_src/pymem3dg",
    include_package_data=True,
    extras_require={"test": ["pytest"]},

    description=DOCLINES[0],
    long_description=open("README.md", encoding="utf8").read(),
    long_description_content_type="text/markdown",
    platforms=["Windows", "Linux", "Mac OS-X", "Unix"],
    classifiers=[c for c in CLASSIFIERS.split("\n") if c],
    keywords="meshing membrane mechanics",
    cmake_args=cmake_args,
    zip_safe=False,
)
