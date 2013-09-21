#!/usr/bin/env python
"""Build and prettify fortran source code"""

import sys,re,os,commands,time,shutil
import prettify

defaultDirectives={"normalize-use":0,"upcase-keywords":0,"clean":0,
            "replace":0,"synopsis":0,"prettify-cvs":0}

def build(srcRoot,buildType="",buildDir=None,clean=None,logFilePath=None):

    if buildDir is None:
        buildDir= os.path.join(srcRoot,"build-"+time.strftime("%y%m%d-%H:%M")+"-"+buildType)

    if not logFilePath:
        logFilePath=os.path.join(buildDir,"build-"+time.strftime("%y%m%d-%H:%M")+"-"+buildType+".log")

    if clean and os.access(buildDir,os.W_OK):
        commands.getoutput('rm -rf "'+buildDir+'"')

    if not os.path.exists(buildDir):
        os.mkdir(buildDir)

    os.chdir(buildDir)

    if os.access(os.path.expanduser("~/cmake/bin/cmake"),os.X_OK):
        makeCmd=os.path.expanduser("~/cmake/bin/cmake")
    else:
        makeCmd="cmake"
    pipe=os.popen("{ { FC=ifort "+makeCmd+' -DCMAKE_Fortran_FLAGS="-mp1 -fp-model precise -fp-model source -fimf-arch-consistency=true" '+os.path.abspath(srcRoot)+"; } 2>&1 ; } >>"+logFilePath)
    if (pipe.close()):
        logFile=open(logFilePath,'a')
        logFile.write("\n+++ ERROR, cmake "+buildType+" FAILED! +++\n")
        logFile.close()
        return None

    if os.access("/usr/bin/gnumake",os.X_OK): makeCmd="/usr/bin/gnumake"
    else: makeCmd="gmake"
    pipe=os.popen("{ { "+makeCmd+"; } 2>&1 ; } >>"+logFilePath)

    if (pipe.close()):
        logFile=open(logFilePath,'a')
        logFile.write("\n+++ ERROR, build "+buildType+" FAILED! +++\n")
        logFile.close()
        return None
    else:
        logFile=open(logFilePath,'a')
        logFile.write("\n+++ build "+buildType+" SUCESSFULLY! +++\n")
        logFile.close()
        return 1

def prettify(srcRoot,buildType="",logDirPath=None,mainLog=sys.stdout,
                 directives=defaultDirectives):

    if not logDirPath:
        logDirPath=os.path.join(srcRoot,"prettify-"+time.strftime("%y%m%d-%H:%M"))

    mainLog.write("====== prettifying files ======\n")
    mainLog.flush()

    logFile=open(os.path.join(logDirPath,"prettify.log"),'w')
    buildDir= os.path.join(srcRoot,"test-"+buildType)
    outDir=os.path.join(logDirPath,"prettifyDir")

    if os.access(outDir,os.W_OK): commands.getoutput('rm -rf "'+outDir+'"')
    os.mkdir(outDir)

    filesToPret2=[] # glob.glob(os.path.join(srcRoot,"src","*.F90"))
    filesToPret=templateInstances
    baseNames=map(os.path.basename,filesToPret)
    for fileToPret in filesToPret2:
        if not os.path.basename(fileToPret) in baseNames:
            filesToPret.append(fileToPret)

    # if requested adds cvs modified files to files to prettify
    if directives["prettify-cvs"]:
        mainLog.write("+ adding cvs modified files to the files to prettify\n")
        mainLog.flush()
        os.chdir(os.path.join(srcRoot,"src"))
        filesC=commands.getoutput("cvs -n update")
        shouldUpdate=0
        filesToPret2=[]
        conflicts=0
        conflictsDir=os.path.join(outDir,"conflicts")
        fileCRe=re.compile(r"([ACRMUP]?) ([a-zA-Z_\.\-]+\.F)$")

        for line in filesC.splitlines():
            m=fileCRe.match(line)
            if m:
                if m.groups()[0] in ["A","M","C"]:
                    filesToPret2.append(os.path.join(srcRoot,"src",m.groups()[1]))
                if m.groups()[0]=="C":
                    conflicts=conflicts+1
                    if not os.path.isdir(conflictsDir): os.mkdir(conflictsDir)
                    shutil.copyfile(m.groups()[1],
                                    os.path.join(conflictsDir,m.groups()[1]))
                    logFile.write(" copied "+m.groups()[1]+" to "+
                                  conflictsDir+"\n")
                if m.groups()[0] in ["U","P","C"]:
                    shouldUpdate=1

        baseNames=map(os.path.basename,filesToPret)
        for fileToPret in filesToPret2:
            if not os.path.basename(fileToPret) in baseNames:
                filesToPret.append(fileToPret)

        logFile.write("cvs modified file to prettify:\n"+str(filesToPret2)+"\n")
        if shouldUpdate:
            mainLog.write("++ WARNING cvs update is needed\n")
            mainLog.write("++ there were "+str(conflicts)+" conflicts\n")
            if conflicts:
                mainLog.write("++ consider restoring the original files from\n"+
                              "++   "+conflictsDir+"\n")
            mainLog.flush()

    # start prettyfication
    logFile.write("Files to prettify:\n"+str(filesToPret)+"\n")
    logFile.write("\nStarting prettification\n")
    logFile.flush()
    mainLog.write("+ starting prettify process\n")
    mainLog.flush()

    if not directives["synopsis"]: buildDir=None
    errors=0
    for fileP in filesToPret:
        try:
            logFile.write("\n+ processing file '"+os.path.basename(fileP)+"'\n")
            infile=open(fileP,"r")
            outfile=open(os.path.join(outDir,os.path.basename(fileP)),"w")
            prettify.prettifyFile(infile,outfile,
                                  normalize_use=directives["normalize-use"],
                                  upcase_keywords=directives["upcase-keywords"],
                                  interfaces_dir=buildDir,
                                  replace=directives["replace"],
                                  logFile=logFile)
            infile.close()
            outfile.close()
        except:
            logFile.write("\n** ERROR prettifying the file "+fileP+"\n")
            import traceback
            logFile.write('-'*60+"\n")
            traceback.print_exc(file=logFile)
            logFile.write('-'*60+"\n")
            logFile.flush()
            errors=errors+1
            if fileP in templateInstances:
                shutil.copyfile(fileP,os.path.join(outDir,os.path.basename(fileP)))
                logFile.write("+ used unprettified template instance for file "+
                              os.path.basename(fileP)+"\n")
    os.mkdir(os.path.join(outDir,"orig"))
    os.chdir(outDir)
    for fileP in os.listdir(outDir):
        try:
            if os.path.isfile(fileP) and fileP[-4:]!=".err":
                logFile.write("moving of file "+fileP)
                origF=os.path.join(srcRoot,"src",os.path.basename(fileP))
                if os.path.isfile(origF):
                    if commands.getoutput('diff "'+origF+'" "'+fileP+'"'):
                        os.rename(origF,os.path.join(outDir,"orig",
                                                     os.path.basename(fileP)))
                        os.rename(fileP,origF)
                    else:
                        logFile.write(" NOT")
                else:
                    os.rename(fileP,origF)
                logFile.write(" done\n")
        except:
            logFile.write("\n** ERROR moving the file "+
                          os.path.basename(fileP)+"\n")
            import traceback
            logFile.write('-'*60+"\n")
            traceback.print_exc(file=logFile)
            logFile.write('-'*60+"\n")
            logFile.flush()
            errors=errors+1

    logFile.write(" ***** prettification finished *****\n")
    logFile.write("+ there were "+str(errors)+" errors\n")
    logFile.close()
    mainLog.write("+ there were "+str(errors)+" errors while prettifying\n")
    mainLog.write("  prettification logFile in '%s'\n"%os.path.basename(logFile.name))
    mainLog.flush()
    # os.rename returns before the operation is complete, so wait a little 
    # (use command.getoutput intread of os.rename?)
    time.sleep(3.0)

if __name__ == '__main__':
    if len(sys.argv)>2:
      print "Usage: ", os.path.basename(sys.argv[0])," [build type]"
    else:
        if len(sys.argv)==2: buildType=sys.argv[1]
        else: buildType=''

    SrcRoot=os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),".."))
    BuildRoot="build-"+time.strftime("%y%m%d-%H:%M")
    print "BuildRoot:",os.path.join(SrcRoot,BuildRoot)
    if not os.path.exists(os.path.join(SrcRoot,BuildRoot)):
        os.mkdir(os.path.join(SrcRoot,BuildRoot))

    build(srcRoot=SrcRoot,buildType=buildType,buildDir=BuildDir)
