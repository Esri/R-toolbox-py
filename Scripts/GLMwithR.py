import arcpy as ARCPY
import arcpy.management as DM
import SSDataObject as SSDO
import SSUtilities as UTILS
import arcpyWithR as RARC
import numpy as NUM
import ErrorUtils as ERROR
import subprocess as SUB
import os as OS
import sys as SYS
import locale as LOCALE
LOCALE.setlocale(LOCALE.LC_ALL, '')

def setupGLM():
    #### Get User Provided Inputs ####
    inputFC = ARCPY.GetParameterAsText(0)    
    outputFC = ARCPY.GetParameterAsText(1) 
    depVarName = str(ARCPY.GetParameterAsText(2))
    indVarNames = ARCPY.GetParameterAsText(3) 
    indVarNames = [ str(i) for i in indVarNames.split(";") ]
    indVarNames = ";".join(indVarNames)

    coefTableIn = ARCPY.GetParameterAsText(4) 
    coefTable, dbf = UTILS.returnTableName(coefTableIn)

    diagTableIn = ARCPY.GetParameterAsText(5) 
    diagTable, dbf = UTILS.returnTableName(diagTableIn)

    #### Create R Command ####
    pyScript = SYS.argv[0]
    toolDir = OS.path.dirname(pyScript)
    rScript = OS.path.join(toolDir, "GLMwithR.r")
    ARCPY.SetProgressor("default", "Executing R Script...")
    args = [RARC.findRExecutable(), "--slave", "--vanilla", "--args",
            inputFC, outputFC, depVarName, indVarNames, 
            coefTable, diagTable]

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

        #### Create SSDO ####
        ssdo = SSDO.SSDataObject(outputFC)

        #### Display Symbology ####
        params = ARCPY.gp.GetParameterInfo()
        try:
            renderType = UTILS.renderType[ssdo.shapeType.upper()]
            if renderType == 0:
                renderLayerFile = "StdResidPoints.lyr"
            elif renderType == 1:
                renderLayerFile = "StdResidPolylines.lyr"
            else:
                renderLayerFile = "StdResidPolygons.lyr"
            fullRLF = OS.path.join(ARCPY.GetInstallInfo()['InstallDir'], 
                                   "ArcToolbox", "Templates", "Layers",
                                   renderLayerFile)
            params[1].Symbology = fullRLF 
        except:
            ARCPY.AddIDMessage("WARNING", 973)

        #### Add to TOC ####
        ARCPY.SetParameterAsText(4, coefTable)


        #### Add to TOC ####
        ARCPY.SetParameterAsText(5, diagTable)

if __name__ == '__main__':
    test = setupGLM() 
