import sys
import json
import os
import argparse

# parse the config file 
parser = argparse.ArgumentParser()
parser.add_argument("config", help = "configuration file (.json) used for the simulation", type = str)
args = parser.parse_args()

# read the config file
with open(args.config) as f:
    data = json.load(f)
geo = data["parameters"]["geometry"]
prop = data["parameters"]["properties"]
var = data["parameters"]["variables"]
inte = data["parameters"]["integration"]
dep = data["parameters"]["dependencies"]
opt = data["parameters"]["options"]

# find dependency
sys.path.append(dep["pymem3dg"])
import pymem3dg

# make I/O directories
cwd = os.getcwd()
IDir = os.path.join(cwd, os.path.split(geo["inputMesh"])[0])
ODir = os.path.join(cwd, geo["outputDir"])
if not os.path.exists(IDir):
    os.mkdir(IDir)
if not os.path.exists(ODir):
    os.mkdir(ODir)

# create starting mesh
if geo["generateGeometry"] == True:
    pymem3dg.genIcosphere( nSub = geo["nSub"] 
                        , path = geo["refMesh"], R = geo["R"])
# elif geo["inputMesh"] == "UVsphere.ply":
#     pymem3dg.genUVsphere( nSub = geo["nSub"] 
#                     , path = geo["inputMesh"])

# run simulation
pymem3dg.driver(inputMesh = geo["inputMesh"],
            outputDir = geo["outputDir"],
            refMesh = geo["refMesh"],
            radius = geo["radiusOfIntegration"],

            isTuftedLaplacian = opt["isTuftedLaplacian"],
            mollifyFactor = opt["mollifyFactor"],
            isVertexShift = opt["isVertexShift"],

            H0 = var["H0*R"] / geo["R"],
            Vt = var["Vt"],
            ptInd = var["ptInd"],  
            Kf = var["Kf"],
            conc = var["conc"],
            height = var["height"],
            
            Kb  = prop["Kb"],
            Kse = prop["Kse"],     
            Ksl = prop["Ksl"],		
            Ksg = prop["Ksg"],		
            Kv  = prop["Kv"], 	
            kt = prop["kt"], 
            gamma = prop["gamma"], 
            
            h = inte["h"],
            T = inte["T"],
            eps = inte["eps"],
            closeZone = inte["closeZone"],
            increment = inte["increment"],
            tSave = inte["tSave"],
            tMollify = inte["tMollify"])
