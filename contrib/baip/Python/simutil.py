import os,re
import fileutil as fu

class ProcSim:
    def __init__(self,sim,runDir):
        self.runDir=runDir
        self.sim=sim
        self.run=os.walk(sim)
        self.runIter=iter(self.run)
        self.count=0
    def __del__(self):
        pass
    def __iter__(self):
        return self
    def next(self):
        for self.curRunTuple in self.runIter:
            self.curRun=self.curRunTuple[0]
            self.runMatch=self.runDir.search(self.curRun)
            if self.runMatch is not None:
                self.count+=1
                return self.curRun
        self.curRun=''
        raise StopIteration
    def getCurRun(self):
        return self.curRun
    def getCount(self):
        return self.count
    def getMatchedGroups(self):
        return self.runMatch.groups()

class ProcDir:
    def __init__(self,simRegex,runRegex=r'RUN\d+$'):
        self.simDir=re.compile(simRegex)
        self.runDir=re.compile(runRegex)
        self.root=os.getcwd()
        self.sim=os.walk(self.root)
        self.simIter=iter(self.sim)
    def __del__(self):
        pass
    def __iter__(self):
        return self
    def next(self):
        for self.curSimTuple in self.simIter:
            self.curSim=self.curSimTuple[0]
            self.simMatch=self.simDir.search(self.curSim)
            if self.simMatch is not None:
                return self.curSim
        self.curSim=''
        raise StopIteration
    def getCurSim(self):
        return self.curSim
    def getRuns(self):
        return ProcSim(self.curSim,self.runDir)
    def getMatchedGroups(self):
        return self.simMatch.groups()

def modPBS(pathname,prefix=None,phaselist=None,proc=None,start=None):
    if prefix is None:
        prefix=pathname
    fcfg_in=open(os.path.join(prefix,'sub.script'))
    cfg_in=fcfg_in.readlines()
    fcfg_in.close()
    fcfg_out=open(os.path.join(pathname,'sub.script'),'w')
    for line in cfg_in:
        if re.match(r'#PBS -N ',line) is not None:
            line='#PBS -N %s\n' % os.path.basename(pathname)
        elif re.match(r'DATA=',line) is not None:
            line='DATA=%s\n' % pathname
        elif phaselist is not None and re.match(r'for FCUR in ',line) is not None:
            line='for FCUR in %s;do\n' % phaselist
        elif proc is not None and '$PROG' in line:
            if start is not None: #Koronis-specific
                processorList='dplace -c%d-%d -s1 ' % (start,start+proc-1)
            else:
                processorList=''
            ip=line.find('$PROG')
            im=line.find('mpirun')
            if im<0:
                im=ip
            if proc>1:
                line=line[0:im]+('mpirun -np %d '%proc)+processorList+line[ip:]
            else:
                line=line[0:im]+line[ip:]
        fcfg_out.write(line)
    fcfg_out.close()
    os.chmod(os.path.join(pathname,'sub.script'),0700);

class qSubmitter:
    def __init__(self,jobsPerGroup=-1,jobsPerNode=8,walltime=24*3600,coresPerNode=8,rootDir=None,iBatchBase=100,koronis=False):
        self.jobsPerGroup=jobsPerGroup
        self.jobsPerNode=jobsPerNode
        self.hours=int(walltime/3600)
        self.minutes=int((walltime-self.hours*3600)/60)
        self.seconds=int(walltime-self.hours*3600-self.minutes*60)
        self.coresPerNode=coresPerNode
        if rootDir is None:
            self.rootDir=os.getcwd()
        else:
            self.rootDir=rootDir
        self.iBatchBase=iBatchBase
        self.koronis=koronis
        self.iJob=0
        self.qBatch=[]
        self.iBatch=self.iBatchBase
        self.qBatchList=[]

    def writeBatchFile(self):
        import math
        self.iBatch+=1
        cwd=os.getcwd()
        os.chdir(self.rootDir)
        fname='qbatch-%d.script'%self.iBatch
        self.qBatchList.append(fname)
        fqbatch=open(fname,'w')
        if self.koronis:
            fqbatch.write('#!/bin/bash -l\n#PBS -l select=%d:ncpus=%d,place=excl:group=socket,walltime=%d:%d:%d\n#PBS -N %s-%d\n#PBS -e stderr-%d\n#PBS -o stdout-%d\n#PBS -q uv1000\n\nDATA=$PBS_O_WORKDIR\ncd $DATA\n\nNODELIST=$(cat $PBS_NODEFILE | uniq)\necho $NODELIST > tmp-$PBS_JOBID\nread -a NODE < tmp-$PBS_JOBID\nrm tmp-$PBS_JOBID\n\ninode=-1\nijob=0\n\n' % (math.ceil(float(self.iJob)/self.jobsPerNode),self.coresPerNode,self.hours,self.minutes,self.seconds,os.path.basename(self.rootDir),self.iBatch,self.iBatch,self.iBatch))
        else:
            fqbatch.write('#!/bin/bash -l\n#PBS -l nodes=%d:ppn=%d,walltime=%d:%d:%d\n#PBS -N %s-%d\n#PBS -e stderr-%d\n#PBS -o stdout-%d\n###PBS -q devel\n\nDATA=$PBS_O_WORKDIR\ncd $DATA\n\nNODELIST=$(cat $PBS_NODEFILE | uniq)\necho $NODELIST > tmp-$PBS_JOBID\nread -a NODE < tmp-$PBS_JOBID\nrm tmp-$PBS_JOBID\n\ninode=-1\nijob=0\n\n' % (math.ceil(float(self.iJob)/self.jobsPerNode),self.coresPerNode,self.hours,self.minutes,self.seconds,os.path.basename(self.rootDir),self.iBatch,self.iBatch,self.iBatch))
        fqbatch.writelines(self.qBatch)
        fqbatch.write('wait\n\nexit 0\n')
        fqbatch.close()
        self.qBatch=[]
        self.iJob=0
        os.chdir(cwd)

    def addJob(self,runString):
        if self.jobsPerGroup>0:
            self.iJob+=1
            if self.koronis:
                self.qBatch.append('bash %s &\n\n' % (runString))
            else:
                self.qBatch.append('[ $((ijob++ %% %d)) -eq 0 ] && ((inode++))\nssh ${NODE[inode]} %s &\n\n' % (self.jobsPerNode,runString))
            if self.iJob % self.jobsPerGroup == 0:
                self.writeBatchFile()
        else:
            self.iBatch+=1
            self.qBatchList.append(runString)

    def submit(self):
        cwd=os.getcwd()
        os.chdir(self.rootDir)
        if self.iJob>0:
            self.writeBatchFile()
        submit=raw_input('Submit all %d jobs? yes/[no]: ' % (self.iBatch-self.iBatchBase))
        if submit=='yes':
            for fname in self.qBatchList:
                if self.jobsPerGroup<=0:
                    os.chdir(os.path.dirname(fname))
                    fname=os.path.basename(fname)
                os.system('qsub %s' %fname)
        os.chdir(cwd)


def batchMove(suffix,last,nbox=3,suffixPrev=None,prev=None):
    import shutil
    try:
        if suffixPrev is not None and prev is not None:
            shutil.move('config.'+suffixPrev+str(prev),'config.'+suffix+str(last))
            shutil.move('run.'+suffixPrev+str(prev),'run.'+suffix+str(last))
            shutil.move('movie.'+suffixPrev+str(prev),'movie.'+suffix+str(last))
            shutil.move('fort12.'+suffixPrev+str(prev),'fort12.'+suffix+str(last))
            for j in range(nbox):
                shutil.move('box'+str(j+1)+'movie.'+suffixPrev+str(prev),'box'+str(j+1)+'movie.'+suffix+str(last))
                shutil.move('box'+str(j+1)+'config.'+suffixPrev+str(prev),'box'+str(j+1)+'config.'+suffix+str(last))
        elif suffixPrev is None and prev is None:
            shutil.copy('fort.77','config.'+suffix+str(last))
            shutil.move('run1a.dat','run.'+suffix+str(last))
            shutil.move('movie1a.dat','movie.'+suffix+str(last))
            shutil.move('fort.12','fort12.'+suffix+str(last))
            for j in range(nbox):
                shutil.move('box'+str(j+1)+'movie1a.xyz','box'+str(j+1)+'movie.'+suffix+str(last))
        else:
            return -1
    except IOError,X:
        print(type(X),X.args)


def prodToEquil(prodSuffix,prodLast,equilSuffix,equilString,nbox=3):
    """Rename all production run output files to equilibration run output files"""
    lastfile=re.compile(equilString)
    equilLast=0
    for cur_file in os.listdir('.'):
        filematch=lastfile.match(cur_file)
        if filematch is not None and int(filematch.group(1)) > equilLast:
            equilLast=int(filematch.group(1))
    for i in range(prodLast):
        batchMove(suffix=equilSuffix,last=equilLast+i+1,nbox=nbox,suffixPrev=prodSuffix,prev=i+1)
    prodLast+=equilLast
    prodSuffix=equilSuffix
    print "%s%d->%d" % (equilSuffix,equilLast,prodLast)
    return (prodSuffix,prodLast)

class procRun:
    def __init__(self,inputname,eps=1e-7):
        self.volTooSmall=False
        ibox=0

        for line in open(inputname):
######################################################################
# Accumulating block averages
######################################################################
            if '------------ box:' in line:
                ibox = int((line.split())[-1])
            elif ibox>0:
                fields=line.split()
                try:
                    self.energyTrajectory[ibox].append(float(fields[1]))
                    self.densityTrajectory[ibox].append(float(fields[2]))
                    self.pressureTrajectory[ibox].append(float(fields[3]))
                    self.surfaceTensionTrajectory[ibox].append(float(fields[4]))
                    self.molarFractionTrajectory[ibox].append(float(fields[5]))
                except ValueError:
                    continue
######################################################################
            elif 'number of cycles:' in line:
                self.ncycles = int((line.split())[-1])
                #print 'ncycles: ', self.ncycles
            elif 'number of boxes' in line:
                self.nbox = int((line.split())[-1])
                self.chempot=map(lambda x:[0] , range(self.nbox+1))
                self.nmol=map(lambda x:[0] , range(self.nbox+1))
                self.numdens=map(lambda x:[0] , range(self.nbox+1))
                self.molfrac=map(lambda x:[0] , range(self.nbox+1))
                self.eta=[[0] for i in range(self.nbox+1)]
                self.press=[0]
                self.molvol=[0]
                self.Hvap=[0]
                self.dens=[0]
                self.energy=[0]
                self.energyTrajectory=map(lambda x:[] , range(self.nbox+1))
                self.densityTrajectory=map(lambda x:[] , range(self.nbox+1))
                self.pressureTrajectory=map(lambda x:[] , range(self.nbox+1))
                self.surfaceTensionTrajectory=map(lambda x:[] , range(self.nbox+1))
                self.molarFractionTrajectory=map(lambda x:[] , range(self.nbox+1))
                #print 'nbox: ', self.nbox,'chempot=',self.chempot
            elif '  temperature:' in line:
                fields=line.split()
                if len(fields) == 3:
                    self.temp = float(fields[1])
                    #print 'temp=',self.temp
            elif 'number of molecule types' in line:
                self.nmolty = int((line.split())[-1])
                #print 'nmolty: ',self.nmolty
            elif 'pmvol:' in line:
                self.pmvol = float((line.split())[-1])
                #print 'pmvol: ',self.pmvol
            elif 'pmcb:' in line:
                self.pmcb = float((line.split())[-1])
                if abs(self.pmcb)<eps:
                    self.pmcb=self.pmswap
                #print 'pmcb: ',self.pmcb
######################################################################
# Details for swap & identity switch moves
######################################################################
            elif 'pmswat:' in line:
                self.pmswat = float((line.split())[-1])
                if abs(self.pmswat)<eps:
                    self.pmswat=self.pmvol
                self.pmsatc=[0]
                self.pmswtcb=[[]]
                #print 'pmswat: ',self.pmswat
            elif 'probability of each swatch pair' in line:
                self.pmsatc.append(float((line.split())[-1]))
                if abs(self.pmsatc[-1])<eps:
                    self.pmsatc[-1]=self.pmsatc[-2]
            elif 'number of swatch box pairs:' in line:
                self.pmswtcb.append([0])
            elif 'probability of the swatch box pair:' in line:
                self.pmswtcb[-1].append(float(line.split()[-1]))
                if abs(self.pmswtcb[-1][-1])<eps:
                    self.pmswtcb[-1][-1]=self.pmswtcb[-1][-2]
            elif 'pmswap:' in line:
                self.pmswap = float((line.split())[-1])
                if abs(self.pmswap)<eps:
                    self.pmswap=self.pmswat
                self.pmswmt=[0]
                self.pmswapb=[[]]
                #print 'pmswap: ',self.pmswap
            elif 'swap probability for molecule type' in line:
                self.pmswmt.append(float((line.split())[-1]))
                if abs(self.pmswmt[-1])<eps:
                    self.pmswmt[-1]=self.pmswmt[-2]
            elif 'number of swap box pairs for molecule type' in line:
                self.pmswapb.append([0])
            elif 'pmswapb' in line:
                self.pmswapb[-1].append(float((line.split())[-1]))
                if abs(self.pmswapb[-1][-1])<eps:
                    self.pmswapb[-1][-1]=self.pmswapb[-1][-2]
######################################################################

######################################################################
# Counting acceptance statistics for swap & identity swap moves
######################################################################
            elif 'Molecule swap' in line:
                ltag = 1
                self.ntswap=0.0
                self.nswap=[0]
            elif 'Molecule swatch' in line:
                ltag = 2
                self.ntswat=0.0
                self.nswat=[0]
            elif 'accepted =' in line and 'between box' in line:
                if ltag == 1:
                    try:
                        self.nswap.append(float((line.split())[-1].lstrip('=')))
                    except ValueError:
                        self.nswap.append(0.0)
                    self.ntswap+=self.nswap[-1]
                    #print 'swap: ',len(self.nswap)-1,', accepted: ',self.nswap[-1]
                elif ltag == 2:
                    self.nswat.append(float((line.split())[-1].lstrip('=')))
                    self.ntswat+=self.nswat[-1]
                    #print 'swatch: ',len(self.nswat)-1,', accepted: ',self.nswat[-1]
######################################################################
            elif 'chemical potential  itype' in line:
                fields=line.split()
                for i in range(1,self.nbox+1):
                    self.chempot[int(fields[5])].append(float(fields[-3]))
                #print 'moltyp: ',fields[4],', chem potential: ',self.chempot[1][int(fields[4])],self.chempot[2][int(fields[4])],self.chempot[3][int(fields[4])]
            elif 'no. of chains' in line:
                fields=line.split()
                for i in range(1,self.nbox+1):
                    self.nmol[i].append(float(fields[6+i]))
                    #print 'nmol[',i,'][',len(box[i])-1,']=',self.nmol[ibox][-1]
            elif 'number density      itype' in line:
                fields=line.split()
                try:
                    self.numdens[int(fields[5])].append(float(fields[-3]))
                except ValueError:
                    self.numdens[int(fields[5])].append(0.0)
            elif 'mole fraction       itype' in line:
                fields=line.split()
                try:
                    self.molfrac[int(fields[5])].append(float(fields[-3]))
                except ValueError:
                    self.molfrac[int(fields[5])].append(0.0)
            elif 'energy offset for box' in line:
                fields=line.split()
                self.eta[int(fields[4].rstrip(':'))].append(float(fields[5]))
            elif 'A move was attempted' in line:
                self.volTooSmall=True
            elif 'pressure         box' in line:
                try:
                    self.press.append(float(line.split()[-3]))
                except ValueError:
                    self.press.append(0.0)
                #print 'press[',len(self.press)-1,']=',self.press[-1]
            elif 'molar volume' in line:
                fields=line.split()
                for i in fields[-3:]:
                    try:
                        self.molvol.append(float(i))
                    except ValueError:
                        self.molvol.append(0.0)
            elif 'H_vap      [kJ/mol]' in line:
                try:
                    self.Hvap.append(float(line.split()[-3]))
                except ValueError:
                    self.Hvap.append(0.0)
                #print 'H_vap=',self.Hvap[-1]
            elif 'specific density box' in line:
                try:
                    self.dens.append(float(line.split()[-3]))
                except ValueError:
                    self.dens.append(0.0)
            elif 'Total energy' in line and 'kJ' in line:
                fields=line.split()
                for i in range(self.nbox):
                    try:
                        self.energy.append(float(fields[i-2*self.nbox])/120.27) # K to kJ
                    except ValueError:
                        self.energy.append(0.0)
            elif re.match(r'molecule\s+\d+:',line) is not None:
                fields=line.split()
                imol=int(fields[1][:-1])
                for i in range(1,self.nbox+1):
                   self.eta[i][imol]=float(fields[i-self.nbox-1])
                #print 'eta ',len(self.eta)-1,': ',self.eta[-1]

        if self.volTooSmall:
            print "volume too small!"

        self.pmcb-=self.pmswap
        self.pmswap-=self.pmswat;
        self.pmswat-=self.pmvol;
        for i in range(1,len(self.pmsatc)):
            self.pmsatc[-i]-=self.pmsatc[-i-1];
            #print 'prob. of swatch pair ',len(self.pmsatc)-i,': ',self.pmsatc[-i]
        for i in range(1,len(self.pmswtcb)):
            for j in range(1,len(self.pmswtcb[i])):
                self.pmswtcb[i][-j]-=self.pmswtcb[i][-j-1];
            #print 'prob. of swatch for molecule ',i,': ',self.pmswtcb[i]
        for i in range(1,len(self.pmswmt)):
            self.pmswmt[-i]-=self.pmswmt[-i-1];
            #print 'prob. of swap ',len(self.pmswmt)-i,': ',self.pmswmt[-i]
        for i in range(1,len(self.pmswapb)):
            for j in range(1,len(self.pmswapb[i])):
                self.pmswapb[i][-j]-=self.pmswapb[i][-j-1];
            #print 'prob. of swap for molecule ',i,': ',self.pmswapb[i]

    def __del__(self):
        pass

class procMov:
    def __init__(self,inputname):
        self.fin = fu.ReadFile(inputname)
        self.nframe=self.fin.readNum()
        self.nchain=self.fin.readNum()
        self.nmolty=self.fin.readNum()
        self.nbox=self.fin.readNum()
        self.nhere = self.fin.readNum()
        #print self.nframe,self.nchain,self.nmolty,self.nbox,self.nhere
        self.rcut=[0] * (self.nbox+1)
        for ibox in range(1,self.nbox+1):
            self.rcut[ibox]=self.fin.readNum()
            #print self.rcut[ibox]
        self.temphere = [0] *(self.nhere+1)
        for i in range(1,self.nhere+1):
            self.temphere[i]=self.fin.readNum()
            #print self.temphere[i]
        self.nunit = [0] *(self.nmolty+1)
        self.nvib=[0]*(self.nmolty+1)
        self.ijvib=[0]*(self.nmolty+1)
        self.ntor=[0]*(self.nmolty+1)
        self.itor2=[0]*(self.nmolty+1)
        self.itor3=[0]*(self.nmolty+1)
        self.itor4=[0]*(self.nmolty+1)
        for imolty in range(1,self.nmolty+1):
            self.nunit[imolty]=self.fin.readNum()
            #print self.nunit[imolty]
            self.nvib[imolty]=[0]*(self.nunit[imolty]+1)
            self.ijvib[imolty]=[0]*(self.nunit[imolty]+1)
            self.ntor[imolty]=[0]*(self.nunit[imolty]+1)
            self.itor2[imolty]=[0]*(self.nunit[imolty]+1)
            self.itor3[imolty]=[0]*(self.nunit[imolty]+1)
            self.itor4[imolty]=[0]*(self.nunit[imolty]+1)
            for iunit in range(1,self.nunit[imolty]+1):
                self.nvib[imolty][iunit]=self.fin.readNum()
                #print self.nvib[imolty][iunit]
                self.ijvib[imolty][iunit]=[0]*(self.nvib[imolty][iunit]+1)
                for i in range(1,self.nvib[imolty][iunit]+1):
                    self.ijvib[imolty][iunit][i]=self.fin.readNum()
                    #print self.ijvib[imolty][iunit][i]
            for iunit in range(1,self.nunit[imolty]+1):
                self.ntor[imolty][iunit]=self.fin.readNum()
                #print self.ntor[imolty][iunit]
                self.itor2[imolty][iunit]=[0]*(self.ntor[imolty][iunit]+1)
                self.itor3[imolty][iunit]=[0]*(self.ntor[imolty][iunit]+1)
                self.itor4[imolty][iunit]=[0]*(self.ntor[imolty][iunit]+1)
                for i in range(1,self.ntor[imolty][iunit]+1):
                    self.itor2[imolty][iunit][i]=self.fin.readNum()
                    self.itor3[imolty][iunit][i]=self.fin.readNum()
                    self.itor4[imolty][iunit][i]=self.fin.readNum()
                    #print self.itor2[imolty][iunit][i],self.itor3[imolty][iunit][i],self.itor4[imolty][iunit][i]
        self.moltyp=[0]*(self.nchain+1)
        self.ncmt=[0]*(self.nbox+1);
        for ibox in range(1,self.nbox+1):
            self.ncmt[ibox]=[0]*(self.nmolty+1)
        self.boxlx=[0]*(self.nbox+1);
        self.boxly=[0]*(self.nbox+1);
        self.boxlz=[0]*(self.nbox+1);
        self.molbox=[0]*(self.nchain+1);
        self.xcom=[0]*(self.nchain+1);
        self.ycom=[0]*(self.nchain+1);
        self.zcom=[0]*(self.nchain+1);
        self.xbead=[0]*(self.nchain+1);
        self.ybead=[0]*(self.nchain+1);
        self.zbead=[0]*(self.nchain+1);
        self.qbead=[0]*(self.nchain+1);
        self.beadtype=[0]*(self.nchain+1);
        self.iframe=0

    def readFrame(self):
        self.iframe+=1
        if self.iframe>self.nframe:
            raise StopIteration
        self.ncycles=self.fin.readNum()
        #print self.ncycles
        for ibox in range(1,self.nbox+1):
            for imolty in range(1,self.nmolty+1):
                self.ncmt[ibox][imolty] = self.fin.readNum()
                #print self.ncmt[ibox][imolty]
            self.boxlx[ibox]=self.fin.readNum()
            self.boxly[ibox]=self.fin.readNum()
            self.boxlz[ibox]=self.fin.readNum()
            #print self.boxlx[ibox],self.boxly[ibox],self.boxlz[ibox]
        for i in range(1,self.nchain+1):
            jchain=self.fin.readNum()
            self.moltyp[jchain]=self.fin.readNum()
            self.nunit[self.moltyp[jchain]]=self.fin.readNum()
            self.molbox[jchain]=self.fin.readNum()
            self.xcom[jchain]=self.fin.readNum()
            self.ycom[jchain]=self.fin.readNum()
            self.zcom[jchain]=self.fin.readNum()
            self.xbead[jchain]=[0]*(self.nunit[self.moltyp[jchain]]+1)
            self.ybead[jchain]=[0]*(self.nunit[self.moltyp[jchain]]+1)
            self.zbead[jchain]=[0]*(self.nunit[self.moltyp[jchain]]+1)
            self.qbead[jchain]=[0]*(self.nunit[self.moltyp[jchain]]+1)
            self.beadtype[jchain]=[0]*(self.nunit[self.moltyp[jchain]]+1)
            #print jchain,self.moltyp[jchain],self.nunit[self.moltyp[jchain]],self.molbox[jchain],self.xcom[jchain],self.ycom[jchain],self.zcom[jchain]
            for iunit in range(1,self.nunit[self.moltyp[jchain]]+1):
                self.xbead[jchain][iunit]=self.fin.readNum()
                self.ybead[jchain][iunit]=self.fin.readNum()
                self.zbead[jchain][iunit]=self.fin.readNum()
                self.qbead[jchain][iunit]=self.fin.readNum()
                self.beadtype[jchain][iunit]=self.fin.readNum()
                #print self.xbead[jchain][iunit],self.ybead[jchain][iunit],self.zbead[jchain][iunit],self.qbead[jchain][iunit],self.beadtype[jchain][iunit]
        return self.iframe

class procTraj:
    def __init__(self,outputFileName,trajFileName,last,first=None,logFile=None,excludeMolty=[],numDensity=False,massFraction=False,nblock=1,nbox=None):
        import numpy as np
        if nbox is None:
            newformat=True
        else:
            newformat=False
        excludeMolty=[i-1 for i in excludeMolty]
        if first is None and logFile is None:
            print 'Both first and logFile are not given!'
            raise SystemExit
        if logFile is not None:
            try:
                flast=file(logFile)
                st=int(flast.readline())+1
                flast.close()
            except IOError:
                st=1
        if first is not None:
            if logFile is not None:
                if first!=st:
                    print 'logFile and first inconsistent. Start from st = ',st
            else:
                st=first
        if st==1:
            readflag='w'
        else:
            readflag='a'
        fou=file(outputFileName,readflag)
        self.nPressure=0
        self.pressure=0
        for i in range(st,last+1):
            fname=trajFileName+str(i)
            print "\tprocessing ", fname
            try:
                fin = fu.ReadFile(fname)
                nstep = fin.readNum()
                if newformat:
                    iratp = fin.readNum()
                    nbox = fin.readNum()
                nmolty = fin.readNum()
                press=np.zeros((nbox),dtype=np.float64)
                molfrac=np.zeros((nbox,nmolty),dtype=np.float64)
                ncmt=np.zeros((nmolty),dtype=np.float64)
                M=np.zeros((nmolty),dtype=np.float64)
                for imolty in range(nmolty):
                    M[imolty]=fin.readNum()
                    if not massFraction:
                        M[imolty]=1.0
                eof=False
                for istep in range(1,nstep+1):
                    try:
                        for ibox in range(nbox):
                            boxlx=fin.readNum()
                            boxly=fin.readNum()
                            boxlz=fin.readNum()
                            v=fin.readNum()
                            if newformat:
                                press[ibox]=fin.readNum()
                            ntbox=0
                            for imolty in range(nmolty):
                                ncmt[imolty]=M[imolty]*fin.readNum()
                                if not imolty in excludeMolty:
                                    ntbox+=ncmt[imolty]
                            if numDensity:
                                ntbox=boxlx*boxly*boxlz
                            molfrac[ibox,:]+=ncmt/ntbox
                    except EOFError:
                        print "EOF Error"
                        eof=True
                    finally:
                        if newformat and istep % iratp == 0:
                            self.nPressure+=1
                            if numDensity:
                                vapbox=np.sum(molfrac,axis=1).argmin(0)
                            else:
                                vapbox=nbox-1
                            self.pressure=float(self.nPressure-1)/float(self.nPressure)*self.pressure+press[vapbox]/float(self.nPressure)
                        if eof or istep % nblock == 0 or istep == nstep: #reaching file end prematurely or normally, or time to output
                            if eof or istep%nblock != 0: #reaching file end prematurely, or normally (but the number of steps in this file is not multiple of nblock so the current block is not complete)
                                print "nstep = ",nstep,"; istep = ",istep,"; nblock = ",nblock
                            if eof:
                                nblock = (istep-1) % nblock
                                if nblock == 0: #this is the first step in a new block; nothing needs to be output
                                    break
                                #otherwise try to properly average the data that have been read successfully
                                for iibox in range(ibox):
                                    molfrac[iibox,:]/=np.float64(nblock+1)
                                for iibox in range(ibox,nbox):
                                    molfrac[iibox,:]/=np.float64(nblock)
                            elif istep % nblock != 0: #this is when istep==nstep
                                nblock = istep % nblock #reset the number of steps in a block to be the actual number of steps in the last block
                            if not eof:
                                molfrac/=np.float64(nblock)
                            for iibox in range(nbox):
                                for imolty in range(nmolty):
                                    if not imolty in excludeMolty:
                                        fou.write("%.5G " % (molfrac[iibox,imolty]))
                            fou.write("\n")
                            if eof:
                                break
                            molfrac=np.zeros((nbox,nmolty),dtype=np.float64)
            except IOError:
                print '\tCannot open ',fname
                continue
        flast=file(logFile,'w')
        flast.write('%d' % last)
        flast.close()
        fou.close()

def avgFrac(baseName,runs,nEntry,fRunAvg,fAvg,outputPrefix='',useLastNSteps=-1,shuffle=None):
    import numpy as np
    fracAvg=np.zeros((nEntry),dtype=np.float64)
    nAvg=0
    for irun in runs:
        fname=baseName+str(irun)+".dat"
        if not os.path.exists(fname):
            print fname,' not exists'
            continue
        print "\tprocessing ", fname
        if useLastNSteps>0:
            nstepMin=fu.fileLen(fname)-useLastNSteps
            if nstepMin<0:
                nstepMin=0
        else:
            nstepMin=0
        frac=np.zeros((nEntry),dtype=np.float64)
        fracNew=np.zeros((nEntry),dtype=np.float64)
        fin = fu.ReadFile(fname)
        nstep = 0
        while True:
            try:
                for i in range(nEntry):
                    fracNew[i]=fin.readNum()
                nstep+=1
                if nstep > nstepMin:
                    if shuffle is not None:
                        shuffle(fracNew)
                    frac+=fracNew
            except EOFError:
                frac/=(nstep-nstepMin)
                break

        fRunAvg.write(outputPrefix+" "+str(irun))
        for i in range(nEntry):
            fRunAvg.write(' %.5G' % frac[i])
        fRunAvg.write('\n')
        fracAvg+=frac
        nAvg+=1
    fracAvg/=nAvg
    fAvg.write(outputPrefix)
    for i in range(nEntry):
        fAvg.write(' %.5G' % fracAvg[i])
    fAvg.write('\n')

def incHist(frac,hist,nbin,sep):
    for i in range(len(nbin)):
        if (frac<sep[i+1]):
            hist[sum(nbin[0:i])+(frac-sep[i])*nbin[i]/(sep[i+1]-sep[i])]+=1.0
            break

def histogram(inputname,nEntry,element,sep=[0.0,0.5,1.0],nbin=[100,100],nstepMin=[0]):
    import numpy as np
    hist=map(lambda x:np.zeros((sum(nbin)),dtype=np.float64),nstepMin)
    fracNew=np.zeros((nEntry),dtype=np.float64)
    fin = fu.ReadFile(inputname)
    nstep = 0
    while True:
        try:
            for i in range(nEntry):
                fracNew[i]=fin.readNum()
            nstep+=1
            for i in range(len(nstepMin)):
                if nstep > nstepMin[i]:
                    for x in element:
                        incHist(fracNew[x],hist[i],nbin,sep)
        except EOFError:
            for i in range(len(nstepMin)):
                hist[i]/=(nstep-nstepMin[i])
            break
    return hist

def locatePeak(histFile,pos,val,peakLocator):
    import subprocess
    lines=subprocess.check_output([peakLocator,histFile]).splitlines()
    if len(lines)!=len(pos):
        if len(lines)!=2:
            print 'Found ',len(lines),' peaks in ',histFile
        for i in range(len(pos),len(lines)):
            pos.append(0.0)
            val.append(0.0)
    for i in range(len(lines)):
        fields=[float(x) for x in lines[i].split()]
        pos[i]=fields[0]
        val[i]=fields[1]
    return (pos,val)

def outputHistogramPeaks(fi,hist,sep,nbin,peakLocator=None):
    pos=[]
    val=[]
    last=0
    cur=0
    for i in range(len(nbin)):
        pos.append(sep[i])
        val.append(0.0)
        while cur-last<nbin[i]:
            conc=sep[i]+(cur-last+0.5)*(sep[i+1]-sep[i])/nbin[i]
            fi.write('%.5f %.12G\n' % (conc,hist[cur]))
            if (hist[cur]>val[i]):
                val[i]=hist[cur]
                pos[i]=conc
            cur+=1
        last=cur
    fi.close()
    if peakLocator is not None:
        (pos,val)=locatePeak(fi.name,pos,val,peakLocator)
    return (pos)

def getTime(inputname):
    try:
        i=2
        while True:
            fields=fu.readLastN(inputname,i).split('user')
            if len(fields)!=2:
                i+=1
            else:
                break
        return float(fields[0])+float((fields[1].split('system'))[0])
    except IOError:
        return -1

def getSuffix(inputname):
    try:
        i=1
        while True:
            fields=fu.readLastN(inputname,i).split(':')
            if len(fields)!=2 or not fields[1].startswith(' '):
                i+=1
            else:
                break
        line=fields[0]
        first=-1
        try:
            while line[first].isdigit():
                first-=1
        except IndexError:
            first-=1
        first+=1
        try:
            cur=int(line[first:])
            if i!=3:
                last=cur-1
            else:
                last=cur
        except ValueError:
            raise
        suf=line[:first]
    except IOError:
        suf='melt'
        last=0
        cur=last
    return (suf,last,cur)

def procCFG(cfglast,cfgprev):
    fin=open(cfglast)
    #with open(cfglast) as fin:
    nlast=int(fin.readline().split()[0])
    fin.close()
    fin=open(cfgprev)
    #with open(cfgprev) as fin:
    nprev=int(fin.readline().split()[0])
    fin.close()
    return nlast-nprev

class avgResult:
    def __init__(self,nbox,nmolty,nmoltyeta,nswapbp,nswatbp,npairbox,triclinic=[]):
        self.length=[]
        self.ncycles=[]
        self.pmswat=[]
        self.pmswap=[]
        self.pmcb=[]
        self.pmtra=[]
        self.eta=[[[] for i in range(nmoltyeta)] for ibox in range(nbox)]
        self.pmswmt=[[] for i in range(nswapbp)]
        self.pmswapb=[[[] for ibox in range(npairbox)] for i in range(nswapbp)]
        self.pmsatc=[[] for i in range(nswatbp)]
        self.pmswtcb=[[[] for ibox in range(npairbox)] for i in range(nswatbp)]
        self.rmtra=[[[[] for i in range(3)] for imolty in range(nmolty)] for ibox in range(nbox)]
        self.rmrot=[[[[] for i in range(3)] for imolty in range(nmolty)] for ibox in range(nbox)]
        self.rmvol=[[] for ibox in range(nbox)]
        self.rmhmat=[[] for i in range(9)]
        self.triclinic=triclinic
        self.lineNo={}
    def __str__(self):
        prtstr='length: '+str(self.length)+'\n\n' \
            + 'ncycles: '+str(self.ncycles)+'\n\n' \
            + 'pmswat: '+str(self.pmswat)+'\n\n' \
            + 'pmswap: '+str(self.pmswap)+'\n\n' \
            + 'pmcb: '+str(self.pmcb)+'\n\n' \
            + 'pmtra: '+str(self.pmtra)+'\n\n' \
            + 'eta: '+str(self.eta)+'\n\n' \
            + 'pmswmt: '+str(self.pmswmt)+'\n\n' \
            + 'pmswapb: '+str(self.pmswapb)+'\n\n' \
            + 'pmsatc: '+str(self.pmsatc)+'\n\n' \
            + 'pmswtcb: '+str(self.pmswtcb)+'\n\n' \
            + 'rmtra: '+str(self.rmtra)+'\n\n' \
            + 'rmrot: '+str(self.rmrot)+'\n\n' \
            + 'rmvol: '+str(self.rmvol)+'\n\n' \
            + 'rmhmat: '+str(self.rmhmat)+'\n\n' \
            + 'lineNo: '+str(self.lineNo)+'\n\n'
        return prtstr
    def average(self,fname=None):
        doAvg=False
        if fname is not None:
            import pickle
            try:
                fbak=open(fname,'rb')
                oldAvg=pickle.load(fbak)
                fbak.close()
                doAvg=True
            except IOError:
                doAvg=False

        if len(self.length)>0:
            self.length=sum(self.length)/len(self.length)
            if doAvg:
                this_length=self.length
                self.length=this_length+oldAvg.length
        if len(self.ncycles)>0:
            self.ncycles=min(self.ncycles)
        if len(self.pmswat)>0:
            self.pmswat=sum(self.pmswat)/len(self.pmswat)
            if doAvg:
                self.pmswat=float(this_length)/self.length*self.pmswat+float(oldAvg.length)/self.length*oldAvg.pmswat
        if len(self.pmswap)>0:
            self.pmswap=sum(self.pmswap)/len(self.pmswap)
            if doAvg:
                self.pmswap=float(this_length)/self.length*self.pmswap+float(oldAvg.length)/self.length*oldAvg.pmswap
        if len(self.pmcb)>0:
            self.pmcb=sum(self.pmcb)/len(self.pmcb)
            if doAvg:
                self.pmcb=float(this_length)/self.length*self.pmcb+float(oldAvg.length)/self.length*oldAvg.pmcb
        if len(self.pmtra)>0:
            self.pmtra=sum(self.pmtra)/len(self.pmtra)
            if doAvg:
                self.pmtra=float(this_length)/self.length*self.pmtra+float(oldAvg.length)/self.length*oldAvg.pmtra
        for ibox in range(len(self.eta)):
            for i in range(len(self.eta[ibox])):
                if len(self.eta[ibox][i])>0:
                    self.eta[ibox][i]=sum(self.eta[ibox][i])/len(self.eta[ibox][i])
                    if doAvg:
                        self.eta[ibox][i]=float(this_length)/self.length*self.eta[ibox][i]+float(oldAvg.length)/self.length*oldAvg.eta[ibox][i]
        for i in range(len(self.pmswmt)):
            if len(self.pmswmt[i])>0:
                self.pmswmt[i]=sum(self.pmswmt[i])/len(self.pmswmt[i])
                if doAvg:
                    self.pmswmt[i]=float(this_length)/self.length*self.pmswmt[i]+float(oldAvg.length)/self.length*oldAvg.pmswmt[i]
            for ibox in range(len(self.pmswapb[i])):
                if len(self.pmswapb[i][ibox])>0:
                    self.pmswapb[i][ibox]=sum(self.pmswapb[i][ibox])/len(self.pmswapb[i][ibox])
                    if doAvg:
                        self.pmswapb[i][ibox]=float(this_length)/self.length*self.pmswapb[i][ibox]+float(oldAvg.length)/self.length*oldAvg.pmswapb[i][ibox]
        for i in range(len(self.pmsatc)):
            if len(self.pmsatc[i])>0:
                self.pmsatc[i]=sum(self.pmsatc[i])/len(self.pmsatc[i])
                if doAvg:
                    self.pmsatc[i]=float(this_length)/self.length*self.pmsatc[i]+float(oldAvg.length)/self.length*oldAvg.pmsatc[i]
            for ibox in range(len(self.pmswtcb[i])):
                if len(self.pmswtcb[i][ibox])>0:
                    self.pmswtcb[i][ibox]=sum(self.pmswtcb[i][ibox])/len(self.pmswtcb[i][ibox])
                    if doAvg:
                        self.pmswtcb[i][ibox]=float(this_length)/self.length*self.pmswtcb[i][ibox]+float(oldAvg.length)/self.length*oldAvg.pmswtcb[i][ibox]
        for ibox in range(len(self.rmtra)):
            for imolty in range(len(self.rmtra[ibox])):
                for i in range(3):
                    if len(self.rmtra[ibox][imolty][i])>0:
                        self.rmtra[ibox][imolty][i]=sum(self.rmtra[ibox][imolty][i])/len(self.rmtra[ibox][imolty][i])
                        if doAvg:
                            self.rmtra[ibox][imolty][i]=float(this_length)/self.length*self.rmtra[ibox][imolty][i]+float(oldAvg.length)/self.length*oldAvg.rmtra[ibox][imolty][i]
                    if len(self.rmrot[ibox][imolty][i])>0:
                        self.rmrot[ibox][imolty][i]=sum(self.rmrot[ibox][imolty][i])/len(self.rmrot[ibox][imolty][i])
                        if doAvg:
                            self.rmrot[ibox][imolty][i]=float(this_length)/self.length*self.rmrot[ibox][imolty][i]+float(oldAvg.length)/self.length*oldAvg.rmrot[ibox][imolty][i]
        for ibox in range(len(self.rmvol)):
            if len(self.rmvol[ibox])>0:
                self.rmvol[ibox]=sum(self.rmvol[ibox])/len(self.rmvol[ibox])
                if doAvg:
                    self.rmvol[ibox]=float(this_length)/self.length*self.rmvol[ibox]+float(oldAvg.length)/self.length*oldAvg.rmvol[ibox]
        for i in range(9):
            if len(self.rmhmat[i])>0:
                self.rmhmat[i]=sum(self.rmhmat[i])/len(self.rmhmat[i])
                if doAvg:
                    self.rmhmat[i]=float(this_length)/self.length*self.rmhmat[i]+float(oldAvg.length)/self.length*oldAvg.rmhmat[i]

        if fname is not None:
            fbak=open(fname,'wb')
            pickle.dump(self,fbak,pickle.HIGHEST_PROTOCOL)
            fbak.close()

    def readMaximumDisplacement(self,cfg):
        fcfg=open(cfg)
        fcfg.readline()
        fcfg.readline()
        nbox=len(self.rmtra)
        nmolty=len(self.rmtra[0])
        for ibox in range(nbox):
            for imolty in range(nmolty):
                fieldsTra=fcfg.readline().split()
                fieldsRot=fcfg.readline().split()
                for i in range(3):
                    self.rmtra[ibox][imolty][i].append(float(fieldsTra[i]))
                    self.rmrot[ibox][imolty][i].append(float(fieldsRot[i]))
        line=0
        for ibox in range(nbox):
            fields=fcfg.readline().split()
            j=0
            for imolty in range(nmolty):
                if len(fields)-j<=0:
                    fields=fcfg.readline().split()
                    line+=1
                    j=0
                #rmflcq[ibox][imolty].append(float(fields[j]))
                j+=1
        self.lineNo['rmvol']=2+nbox*(nmolty*2+1)+line
        line=0
        fields=fcfg.readline().split()
        j=0
        for ibox in range(nbox):
            if len(fields)-j<=0:
                fields=fcfg.readline().split()
                line+=1
                j=0
            self.rmvol[ibox].append(float(fields[j]))
            j+=1
        self.lineNo['rmhmat']=self.lineNo['rmvol']+1+line
        if len(self.triclinic)==1:
            fields=fcfg.readline().split()
            j=0
            for i in range(9):
                if len(fields)-j<=0:
                    fields=fcfg.readline().split()
                    j=0
                self.rmhmat[i].append(float(fields[j]))
                j+=1
        elif len(self.triclinic)!=0:
            raise SystemExit
        fcfg.close()
        return

    def writeMaximumDisplacement(self,cfg,skipBox):
        fcfg=open(cfg)
        fdat=fcfg.readlines()
        fcfg.close()
        iLine=1
        nbox=len(self.rmtra)
        nmolty=len(self.rmtra[0])
        for ibox in range(nbox):
            if ibox in skipBox:
                iLine+=2*nmolty
                continue
            for imolty in range(nmolty):
                iLine+=1
                fdat[iLine]=str(self.rmtra[ibox][imolty][0])+'  '+str(self.rmtra[ibox][imolty][1])+'  '+str(self.rmtra[ibox][imolty][2])+'\n'
                iLine+=1
                fdat[iLine]=str(self.rmrot[ibox][imolty][0])+'  '+str(self.rmrot[ibox][imolty][1])+'  '+str(self.rmrot[ibox][imolty][2])+'\n'
        if len(self.triclinic)==1:
            fdat[self.lineNo['rmhmat']]='  '.join(str(self.rmhmat[i]) for i in range(3))+'\n'
            fdat[self.lineNo['rmhmat']+1]='  '.join(str(self.rmhmat[i]) for i in range(3,6))+'\n'
            fdat[self.lineNo['rmhmat']+2]='  '.join(str(self.rmhmat[i]) for i in range(6,9))+'\n'
        fdat[self.lineNo['rmvol']:self.lineNo['rmhmat']]='  '.join(str(self.rmvol[ibox]) for ibox in range(nbox))+'\n'
        fcfg=open(cfg,'w')
        fcfg.writelines(fdat)
        fcfg.close()
        return

def arrangeLattice(nmol):
    x=[int(nmol**(1.0/3))]*3
    i=0
    while (x[0]*x[1]*x[2]<nmol):
        x[i]+=1
        i+=1
    return x

def estimateRuntime(timeTarget,ncycles,timeElapsed,minIncCycles=500):
    return minIncCycles*int(timeTarget*ncycles/timeElapsed/minIncCycles)

def estimateMolSwapProb(Range,nswap,pmswmt):
    p=[]
    for imol in Range:
        prod=1
        for i in Range:
            if i!=imol:
                prod*=(nswap[i*6+1]+nswap[i*6+2]+nswap[i*6+3]+nswap[i*6+4]+nswap[i*6+5]+nswap[i*6+6]+1)/pmswmt[i+1]
        p.append(prod)
    total=sum(p)
    for i in range(len(p)):
        p[i]=p[i]/total
        if p[i]==0:
            p[i]+=0.05
        elif p[i]==1:
            p[i]-=0.05*(len(p)-1)
    q=[sum(p[:i+1]) for i in range(len(p))]
    return q

def estimateSwapBoxProb(nswap,pmswapb,boxpair=range(3)):
    p=[]
    for iboxpair in boxpair:
        prod=1
        for i in boxpair:
            if i!=iboxpair:
                prod*=(nswap[i*2+1]+nswap[i*2+2]+1)/pmswapb[i+1]
        p.append(prod)
    total=sum(p)
    for i in range(len(p)):
        p[i]=p[i]/total
        if p[i]==1:
            p[i]-=0.05*(len(p)-1)
    q=[sum(p[:i+1]) for i in range(len(p))]
    return q

def estimateSwatchPairProb(Range,nswat,pmsatc):
    p=[]
    for ipair in Range:
        prod=1
        for i in Range:
            if i!=ipair:
                prod*=(nswat[6*i+1]+nswat[6*i+2]+nswat[6*i+3]+nswat[6*i+4]+nswat[6*i+5]+nswat[6*i+6]+1)/pmsatc[i+1]
        p.append(prod)
    total=sum(p)
    for i in range(len(p)):
        p[i]=p[i]/total
        if p[i]==0:
            p[i]=0.05
        elif p[i]==1:
            p[i]-=0.05*(len(p)-1)
    q=[sum(p[:i+1]) for i in range(len(p))]
    return q

def estimateSwatchBoxProb(nswat,pmswtcb,boxpair=range(3)):
    p=[]
    for iboxpair in boxpair:
        prod=1
        for i in boxpair:
            if i!=iboxpair:
                prod*=(nswat[i*2+1]+nswat[i*2+2]+1)/pmswtcb[i+1]
        p.append(prod)
    total=sum(p)
    for i in range(len(p)):
        p[i]=p[i]/total
        if p[i]==1:
            p[i]-=0.05*(len(p)-1)
    q=[sum(p[:i+1]) for i in range(len(p))]
    return q
