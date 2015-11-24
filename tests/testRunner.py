#!/usr/bin/env python
"""Prepare fortran source code for cvs check-in: building, testing, prettifying"""

import sys,re,os,commands,time,unittest
import buildUtils
import testUtils

# read directives
directives={"clean":0,"test":1,"synopsis":0,"upcase-keywords":0
            ,"prettify-cvs":0,"replace":0,"normalize-use":0}

directiveRe=re.compile(r"--(no-)?(clean|test|synopsis|upcase-keywords|prettify-cvs|replace|normalize-use)$")

descStr=("Usage: "+os.path.basename(sys.argv[0])+"""
  [--[no-]clean] [--[no-]test] [--help]
  [--[no-]synopsis] [--[no-]upcase-keywords] [--[no-]prettify-cvs]
  [--[no-]replace] [--[no-]normalize-use]

  Prepares for checkin the source.
  defaults="""+str(directives)
  )

if "--help" in sys.argv[1:]:
    print descStr
    sys.exit(0)

for directive in sys.argv[1:]:
    m=directiveRe.match(directive)
    if m:
        directives[m.groups()[1]]=not m.groups()[0]
    else:
        print " ** ERROR **\nUnknown argument",directive
        print descStr
        sys.exit(-1)

Root=os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),".."))

logDirPath=os.path.join(Root,"test-"+time.strftime("%y%m%d"))
if not os.path.exists(logDirPath):
    os.mkdir(logDirPath)

buildDir=os.path.join(logDirPath,"build")

mainLog=open(os.path.join(logDirPath,"main.log"),'w')
mainLog.write(" ******* "+time.strftime("%y%m%d-%H:%M")+
              " prepare check-in BEGAN *******\n")
mainLog.flush()

print "Root directory: "+Root
print "Test directory: "+logDirPath
print "Main log: "+mainLog.name

# prettify
# if directives["prettify-cvs"]:
#     buildUtils.prettify(srcRoot=Root,buildType="",logDirPath=logDirPath
#                     ,mainLog=mainLog,directives=directives)

# build
mainLog.write("====== clean="+str(directives["clean"])+" compile topmon ======\n")
mainLog.write("  compilation logFile in 'build.log'\n")
mainLog.flush()

if not buildUtils.build(srcRoot=Root,buildType="",buildDir=buildDir,clean=directives["clean"],logFilePath=os.path.join(logDirPath,"build.log")):
    mainLog.write("+++ ERROR, build FAILED! +++\n")
else:
    mainLog.write("+++ build SUCESSFULL! +++\n")
mainLog.flush()

# test
if directives["test"]:
    mainLog.write("====== tests ======\n")
    suite=testUtils.simpleTests(srcRoot=Root,exeRoot=os.path.join(buildDir,"src"),testsRoot=logDirPath,command="topmon")
    testRunner=unittest.TextTestRunner(stream=mainLog,verbosity=2)
    retCode=testRunner.run(suite).wasSuccessful()
    if retCode:
        exitCode=0
    else:
        exitCode=-1

mainLog.write(" ******* "+time.strftime("%y%m%d-%H:%M")+" prepare check-in FINISHED *******\n")
mainLog.close()
sys.exit(exitCode)

