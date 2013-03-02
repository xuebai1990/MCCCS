#! /usr/bin/env python
# Perform test runs

import sys,re,os,commands,time,shutil,unittest,glob
import diffEpsilon

defaultSrcRoot=os.path.expanduser("~/MCCCS/topmon")
defaultTestsRoot=os.path.join(defaultSrcRoot,"test-"+time.strftime("%y%m%d-%H:%M"))
defaultExeRoot=os.path.join(defaultTestsRoot,"build","src")
defaultCommand="topmon"

class simpleRunTest(unittest.TestCase):
    def __init__(self,testName,srcRoot=defaultSrcRoot,exeRoot=defaultExeRoot,testsRoot=defaultTestsRoot,command=defaultCommand):
        unittest.TestCase.__init__(self)
        self.testName=testName
        self.srcRoot=srcRoot
        self.testRoot=os.path.join(testsRoot,testName+"-"+os.path.basename(command))
        self.command=os.path.join(exeRoot,command)
        self.outputFile="run1a.dat"

    def setUp(self):
        if not os.access(os.path.abspath(os.path.join(self.testRoot,"..")),os.W_OK):
            os.mkdir(os.path.abspath(os.path.join(self.testRoot,"..")))
        #if not os.path.exists(self.testRoot):
        #    os.mkdir(self.testRoot)
        commands.getoutput('cp -rL "'+os.path.join(self.srcRoot,"tests","SHORT",self.testName)+'" '+self.testRoot)

    def tearDown(self):pass

    def runTest(self):
        os.chdir(self.testRoot)
        pipe=os.popen("{ { "+self.command+"; } 2>&1 ; } >>"+self.outputFile)

        if (pipe.close()):
            self.fail('error, the command returned an error')
        else:
            file1=open(os.path.join(self.testRoot,self.outputFile))
            if os.access(os.path.join(self.testRoot,"spec.out"),os.R_OK):
                diffs=open(os.path.join(self.testRoot,os.path.basename(self.command)+".spec.diff"),"w")
                file2=open(os.path.join(self.testRoot,"spec.out"))
                if diffEpsilon.compareFiles(file1,file2,logFile=diffs)!=0:
                    arch_spec_test="failed"
                else:
                    arch_spec_test="succeded"
            else:
                arch_spec_test=0

            if os.access(os.path.join(self.testRoot,"generic.out"),os.R_OK):
                diffs=open(os.path.join(self.testRoot,os.path.basename(self.command)+".gen.diff"),"w")
                file2=open(os.path.join(self.testRoot,"generic.out"))
                diffVal=diffEpsilon.compareFiles(file1,file2,logFile=diffs)
            elif not arch_spec_test:
                self.fail('no generic output for this test')
            else:
                diffVal=0
                print 'no generic output for this test'

            if diffVal>1.0e-9:
                if arch_spec_test=="succeded":
                    print "** generic output to be updated? Too different from validated output"
                    print self.id(),"diff=",str(diffVal),'see',diffs.name
                elif arch_spec_test:
                    self.fail('output was different than expected for this platform (and too different from generic output), see '+diffs.name)
                else:
                    self.fail('output was too different from generic output ('+str(diffVal)+'), see '+diffs.name)
            else:
                if arch_spec_test and arch_spec_test!="succeded":
                    self.fail('output was different than expected for this platform (but generic output within error", see '+diffs.name)
                print self.id(),"generic_diff=",diffVal
    
    def id(self):
        return "%s.%s-%s" % (self.__class__, self.testName,os.path.basename(self.command))

    def shortDescription(self):
        return "interactive run of the test %s with %s" % (self.testName,os.path.basename(self.command))

def simpleTests(srcRoot=defaultSrcRoot,exeRoot=defaultExeRoot,testsRoot=defaultTestsRoot,command=defaultCommand):
    suite=unittest.TestSuite()
    testDirectories=glob.glob(os.path.join(srcRoot,"tests","SHORT","TEST*"))
    for testDirectory in testDirectories:
        suite.addTest(simpleRunTest(testName=os.path.basename(testDirectory),srcRoot=srcRoot,exeRoot=exeRoot,testsRoot=testsRoot,command=command))
    return suite

simpleTestsDefault=lambda:simpleTests()

for test in simpleTestsDefault()._tests:
    globals()[test.testName.replace("-","_")]=lambda :test

if __name__=="__main__":
    defaultSrcRoot=os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),".."))
    defaultTestsRoot="test-"+time.strftime("%y%m%d-%H:%M")
    print "testsRoot:",os.path.join(defaultSrcRoot,defaultTestsRoot)
    if not os.path.exists(os.path.join(defaultSrcRoot,defaultTestsRoot)):
        os.mkdir(os.path.join(defaultSrcRoot,defaultTestsRoot))

    unittest.TestProgram(defaultTest="simpleTests")
