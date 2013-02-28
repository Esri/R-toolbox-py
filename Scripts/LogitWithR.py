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

def setupLogit():
    #### Get User Provided Inputs ####
    inputFC = ARCPY.GetParameterAsText(0)    
    outputFC = ARCPY.GetParameterAsText(1) 
    depVarName = str(ARCPY.GetParameterAsText(2))
    indVarNames = ARCPY.GetParameterAsText(3) 
    indVarNames = [ str(i) for i in indVarNames.split(";") ]
    indVarNames = ";".join(indVarNames)
    usePenalty = ARCPY.GetParameterAsText(4) 
    if usePenalty == 'true':
        usePenalty = "1"
    else:
        usePenalty = "0"

    coefTableIn = ARCPY.GetParameterAsText(5) 
    coefTable, dbf = UTILS.returnTableName(coefTableIn)

    diagTableIn = ARCPY.GetParameterAsText(6) 
    diagTable, dbf = UTILS.returnTableName(diagTableIn)

    #### Create R Command ####
    pyScript = SYS.argv[0]
    toolDir = OS.path.dirname(pyScript)
    rScript = OS.path.join(toolDir, "logitWithR.r")
    ARCPY.SetProgressor("default", "Executing R Script...")
    args = ["R", "--slave", "--vanilla", "--args",
            inputFC, outputFC, depVarName, indVarNames, 
            usePenalty, coefTable, diagTable]

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

        #### Print Coef Output Table ####
        try:
            rows = ARCPY.SearchCursor(coefTable)
        except:
            ARCPY.AddIDMessage("ERROR", 204)
            raise ERROR.ScriptError()

        labels = ["Variable", "Coef", "StdError", "Wald", "Prob"]
        header = "Logistic Regression Coefficient Table"
        res = [ labels ]
        for row in rows:
            rowRes = []
            for i, val in enumerate(labels):
                if i == 0:
                    rowRes.append(row.getValue(val))
                else:
                    rowRes.append(LOCALE.format("%0.6f", row.getValue(val)))
            res.append(rowRes)
        del rows

        coefTextTab = UTILS.outputTextTable(res, header = header)
        ARCPY.AddMessage("\n")
        ARCPY.AddMessage(coefTextTab)

        #### Add to TOC ####
        ARCPY.SetParameterAsText(5, coefTable)

        #### Print Diag Table (In Two Parts) ####
        try:
            rows = ARCPY.SearchCursor(diagTable)
        except:
            ARCPY.AddIDMessage("ERROR", 204)
            raise ERROR.ScriptError()

        labels = ["Diag_Name", "Diag_Value"]
        header = "Logistic Regression Diagnostic Table"
        resLab1 = []
        resVal1 = []
        resLab2 = []
        resVal2 = []
        c = 0
        for row in rows:
            for i, val in enumerate(labels):
                if i == 0:
                    cellVal = row.getValue(val)
                    if c <= 6:
                        resLab1.append(cellVal)
                    else:
                        resLab2.append(cellVal)
                else:
                    cellVal = LOCALE.format("%0.6f", row.getValue(val))
                    if c <= 6:
                        resVal1.append(cellVal)
                    else:
                        resVal2.append(cellVal)
            c += 1
        del rows

        diagTextTab1 = UTILS.outputTextTable([resLab1, resVal1], header = header)
        ARCPY.AddMessage("\n")
        ARCPY.AddMessage(diagTextTab1)
        ARCPY.AddMessage("\n")
        diagTextTab2 = UTILS.outputTextTable([resLab2, resVal2])
        ARCPY.AddMessage(diagTextTab2)
        ARCPY.AddMessage("\n")

        #### Add to TOC ####
        ARCPY.SetParameterAsText(6, diagTable)

if __name__ == '__main__':
    test = setupLogit() 
