import sys
import json
import os

with open("config.json") as f:
    data = json.load(f)

geo = data["parameters"]["geometry"]
prop = data["parameters"]["properties"]
var = data["parameters"]["variables"]
inte = data["parameters"]["integration"]
dep = data["parameters"]["dependencies"]

sys.path.append(dep["pyddg"])
import pyddg

cwd = os.getcwd()
IDir = os.path.join(cwd, os.path.split(geo["inputMesh"])[0])
ODir = os.path.join(cwd, geo["outputDir"])
if not os.path.exists(IDir):
    os.mkdir(IDir)
if not os.path.exists(ODir):
    os.mkdir(ODir)

pyddg.genIcosphere( nSub = geo["nSub"] 
                    , path = geo["inputMesh"])

# pyddg.genUVsphere( nSub = geo["nSub"] 
#                    , path = geo["inputMesh"])

pyddg.driver(inputMesh = geo["inputMesh"],
            outputDir = geo["outputDir"],

            H0 = var["H0"],
            Vt = var["Vt"],
            ptInd = var["ptInd"],  
            extF = var["extF"],
            conc = var["conc"],
            
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
            tSave = inte["tSave"])