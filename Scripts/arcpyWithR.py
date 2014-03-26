import _winreg

def findRExecutable():
    handle = None
    try:
        handle =_winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE,
                                r"SOFTWARE\R-Core\R")
        i = 0
        while True:
            key, value, _ = _winreg.EnumValue(handle, i)
            if key.lower() == "installpath":
                r_executable_path = os.path.join(value, "bin", "R.exe")
                return r_executable_path
            i += 1
    except:
        return "R"
    finally:
        if handle is not None:
            _winreg.CloseKey(handle)
    # Fall back to default (hope it's on the %PATH%)
    return "R"

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

