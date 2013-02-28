import arcpy as ARCPY
import arcpy.management as DM
import os as OS
import sys as SYS
import subprocess as SUB
import arcpyWithR as RARC

#### Parameter Dictionaries ####
clusterDict = {"KMEANS_HARTIGAN": "kmeansHartigan", "CLARA": "clara", 
               "B_CLUST": "bclust", "M_CLUST": "Mclust", 
               "KCCA_KMEANS": "kccaKmeans", 
               "CMEANS": "cmeans"} 

def PointClusters():
    #### Get User Provided Inputs ####
    inputFC = ARCPY.GetParameterAsText(0) 
    outputFC = ARCPY.GetParameterAsText(1) 
    numClusters = ARCPY.GetParameterAsText(2)
    clusterMethod = ARCPY.GetParameterAsText(3)
    clusterMethodStr = clusterDict[clusterMethod]
    varNames = ARCPY.GetParameterAsText(4)
    varNames = [ str(i) for i in varNames.split(";") ]
    varNames = ";".join(varNames)
    useLocation = ARCPY.GetParameterAsText(5)
    if useLocation == 'true':
        useLocation = "1"
    else:
        useLocation = "0"

    #### Create R Command ####
    pyScript = SYS.argv[0]
    toolDir = OS.path.dirname(pyScript)
    rScript = OS.path.join(toolDir, "PointClusters.r")
    ARCPY.SetProgressor("default", "Executing R Script...")
    args = ["R", "--slave", "--vanilla", "--args",
            inputFC, outputFC, numClusters, clusterMethodStr,
            varNames, useLocation] 

    #### Uncomment Next Two Lines to Print/Create Command Line Args ####
    #cmd = RARC.createRCommand(args, rScript)
    #ARCPY.AddWarning(cmd)

    #### Execute Command ####
    scriptSource = open(rScript, 'rb')
    rCommand = SUB.Popen(args, 
                         stdin = scriptSource,
                         stdout = SUB.PIPE, 
                         stderr = SUB.PIPE,
                         shell=True)

    #### Print Result ####
    resString, errString = rCommand.communicate()

    #### Push Output to Message Window ####
    if errString and "Calculations Complete..." not in resString:
        ARCPY.AddError(errString)
    else:
        resOutString = RARC.printRMessages(resString)
        ARCPY.AddMessage(resOutString)

        #### Project the Data ####
        DM.DefineProjection(outputFC, inputFC)

        #### Render the Results ####
        params = ARCPY.gp.GetParameterInfo()
        renderFile = OS.path.join(toolDir, "RenderClusters.lyr")
        params[1].Symbology = renderFile

if __name__ == '__main__':

    test = PointClusters() 
