import sys
import json
import argparse

# parse the config file 
parser = argparse.ArgumentParser()
parser.add_argument("config", help = "configuration file (.json) used for the simulation", type = str)
parser.add_argument("meshFile", help = "Mesh file (.ply, .obj, .etc.) used for visualization", type = str)
args = parser.parse_args()

# read the config file
with open(args.config) as f:
    data = json.load(f)
geo = data["parameters"]["geometry"]
prop = data["parameters"]["properties"]
var = data["parameters"]["variables"]
inte = data["parameters"]["integration"]
dep = data["parameters"]["dependencies"]

# find dependency
sys.path.append(dep["pymem3dg"])
import pymem3dg

# run viewer 
# pyddg.viewer(args.meshFile)
pymem3dg.view_animation(args.meshFile)
