
def printRMessages(resString):
    initLine = '[1] "Loading Libraries...."'
    outInitLine = "Loading Libraries...."
    resString = resString.replace(initLine, outInitLine)
    resString = resString.replace("\n[1] ", "\n")
    resString = resString.replace('"\n', "\n")
    return resString.replace('\n"', "\n")

def createRCommand(args, rScript):
    argString = args + ["<", rScript]
    return " ".join(argString)

