#!/usr/bin/env python

import os,shutil,random
import simutil as su

def modFort4(flag,suf=None,last=None,info=None):
    avgParam=info.avgParam
    lineNo=info.lineNo
    lineNoSwat=info.lineNoSwat
    temp=info.temp
    press=info.press
    nmol=info.nmol
    boxlen=info.boxlen
    rcut=[i*0.4 for i in boxlen]
    rcut[0]=14.0

    for suffix in ['melt','cool','volume','equil1','prod']:
#    for suffix in ['equil1','prod']:
        fcfg_in=open("fort.4."+suffix)
        fdat=fcfg_in.readlines()
        fcfg_in.close()

        #>>> Generate a random number to use as the random number seed
        mseed=random.random()*1e6
        fdat[1]="     seed=%d\n" % mseed

        #>>> Set the box lengths, cut-offs, temperatures, and pressures for the two simulation boxes
        if suffix == 'melt':
            fdat[lineNo['ncmt']-2]='%.5e %.5e %.5e %.2f 0.23d0 0.0d0 3  F  T  F  F  %.2f  %.5e\n' % (boxlen[0],boxlen[0],boxlen[0],rcut[0],5000,press)
            fdat[lineNo['ncmt']+4]='%.5e %.5e %.5e %.2f 0.23d0 0.0d0 3  F  T  F  F  %.2f  %.5e\n' % (boxlen[1],boxlen[1],boxlen[1],rcut[1],5000,press)
        else:
            fdat[lineNo['ncmt']-2]='%.5e %.5e %.5e %.2f 0.23d0 0.0d0 3  F  T  F  F  %.2f  %.5e\n' % (boxlen[0],boxlen[0],boxlen[0],rcut[0],temp,press)
            fdat[lineNo['ncmt']+4]='%.5e %.5e %.5e %.2f 0.23d0 0.0d0 3  F  T  F  F  %.2f  %.5e\n' % (boxlen[1],boxlen[1],boxlen[1],rcut[1],temp,press)

        #>>> Set total number of molecules, and number of molecules of each type in each box
        fdat[4]='     nchain = %d\n' % sum(nmol)
        fdat[lineNo['ncmt']]='%d        0\n' % nmol[0]
        (x,y,z)=su.arrangeLattice(nmol[0])
        fdat[lineNo['ncmt']+2]='%d %d %d 0 0  0.0  0.000 F 0.0\n' % (x,y,z)
        fdat[lineNo['ncmt']+6]='%d        0\n' % nmol[1]
        (x,y,z)=su.arrangeLattice(nmol[1])
        fdat[lineNo['ncmt']+8]='%d %d %d 0 0  0.0  0.000 F 0.0\n' % (x,y,z)

        #>>> Volume move probabilities
        if suffix=='volume' or suffix=='equil1' or suffix=='prod':
            fdat[lineNo['pmvol']]="     pmvol=  %.4e\n" % (2.0/sum(nmol))

        #>>> Probabilities for CBMC (if needed) and translation moves
        fdat[lineNo['pmcb']]="     pmcb=  0.34d0\n"
        fdat[lineNo['pmtra']]="     pmtra=  0.67d0\n"
        if suffix=='equil1' or suffix=='prod':
            #>>> Set swap move probabilities and adjust the rest
            fdat[lineNo['pmswap']]="     pmswap=  0.34d0\n"
            fdat[lineNo['pmcb']]="     pmcb=  0.56d0\n"
            fdat[lineNo['pmtra']]="     pmtra=  0.78d0\n"

        if flag==1:
            #>>> If the simulation was killed due to overtime,
            #   estimate the number of cycles successfully completed
            #   from the difference between fort.77 and save-config.1
            ncyclesn=su.procCFG("save-config.1","fort.77")
            #>>> Estimate number of cycles for the next run so that it can finish in 22 hours
            #   ncyclesn cycles were done in 24 hours (when the job was killed)
            ncyclesn=su.estimateRuntime(timeTarget=22,ncycles=ncyclesn,timeElapsed=24)
        elif flag==2:
            #>>> If the simulation has completed and has generated an
            #   output file successfully, then process the output
            result=su.procRun("run."+suf+str(last))
            #   and read the time it used from the log file
            timeElapsed=su.getTime("log")
            #   and estimate the number of cycles for the next run
            ncyclesn=su.estimateRuntime(timeTarget=22,ncycles=result.ncycles,timeElapsed=timeElapsed/3600)

        if flag==1 or flag==2:
            avgParam.ncycles.append(ncyclesn)
            fdat[5]="     nstep              = %d\n" % ncyclesn
            ncyclesn/=10
            fdat[17]="     iprint       = %d\n" % ncyclesn
            fdat[19]="     iblock       = %d\n" % ncyclesn
            fdat[18]="     imv          = %d\n" % ncyclesn

        if flag==2 and suf=='equil': #>>> optimize simulation parameters, but only when it is equilibration
            if result.pmswat!=0:
                #>>> If swatch moves are used, try to get one accepted per cycle
                result.pmswat=result.pmswat*result.ncycles/result.ntswat+result.pmvol
            else:
                result.pmswat=result.pmvol

            #>>> Adjust swap move probabilities so that one is accepted per cycle
            result.pmswap=result.pmswap*result.ncycles/result.ntswap+result.pmswat

            if result.pmswap>0.6: #>>> Always reserve at least 40% probabilities for the remaining types of moves
                print "%.5G,%.5G" % (result.pmswat,result.pmswap)
                result.pmswat=result.pmvol+(result.pmswat-result.pmvol)*(0.6-result.pmvol)/(result.pmswap-result.pmvol)
                result.pmswap=0.6

            #>>> Evenly divide between CBMC, translation, and rotation
            if result.pmcb!=0: #>>> Only if CBMC is used
                result.pmcb=(1+2*result.pmswap)/3
            else:
                result.pmcb=result.pmswap
            pmtra=0.5+result.pmcb/2

            avgParam.pmswat.append(result.pmswat)
            avgParam.pmswap.append(result.pmswap)
            avgParam.pmcb.append(result.pmcb)
            avgParam.pmtra.append(pmtra)

            #fdat[lineNo['pmswat']]="     pmswat= %.5e\n" % result.pmswat
            fdat[lineNo['pmswap']]="     pmswap= %.5e\n" % result.pmswap
            fdat[lineNo['pmcb']]="     pmcb   = %.5e\n" % result.pmcb
            fdat[lineNo['pmtra']]="     pmtra = %.5e\n" % pmtra

            #>>> Adjust the probabilities for 2 swatch pairs
            # p=su.estimateSwatchPairProb(range(2),result.nswat,result.pmsatc)
            # for i in range(len(p)):
            #     avgParam.pmsatc[i].append(p[i])
            # fdat[lineNo['pmsatc']]="     pmsatc= %.5e 1.0d0\n" % (p[0])

            #>>> Adjust the probabilities for doing swatching moves between the different box pairs
            #   The following assumes that there are three box pairs, 1<->2, 1<->3, 2<->3
            # for ipair in range(2):
            #     p=su.estimateSwatchBoxProb(result.nswat[ipair*6:ipair*6+7],result.pmswtcb[ipair+1],boxpair=[0,1,2])
            #     for i in range(3):
            #         avgParam.pmswtcb[ipair][i].append(p[i])
            #     fdat[lineNoSwat[ipair]]="   3     %.5e %.5e 1.0d0\n" % (p[0],p[1]) # swatch box pair prob

            #>>> In this hypothetical example, all molecules are also swapped directly.
            #   It similarly assumes that there are three box pairs, 1<->2, 1<->3, 2<->3
            #   but that we are not swapping between boxes 1 and 2
            #   Also try to achieve that the accepted moves are proportional to the number of molecules of that type
            # for imol in range(result.nmolty):
            #     p=su.estimateSwapBoxProb(result.nswap[imol*6:imol*6+7],result.pmswapb[imol+1],boxpair=[1,2])
            #     avgParam.pmswapb[imol][1].append(p[0])
            #     fdat[lineNo['pmswapb']+imol*6]="   3     0.0d0 %.5e 1.0d0\n" % p[0] # swap box pair prob
            #     result.nswap[imol*6:imol*6+7]=[float(x)/nmol[imol] for x in result.nswap[imol*6:imol*6+7]]

            #>>> Adjust the probabilities for swapping each molecule type
            # p=su.estimateMolSwapProb(range(result.nmolty),result.nswap,result.pmswmt) # swap water,furan-0.25
            # for i in range(len(p)):
            #     avgParam.pmswmt[i].append(p[i])
            # fdat[lineNo['pmswap']+1]="     pmswmt= %.5e 1.0d0\n" % (p[0])

        fcfg_out=open("fort.4."+suffix,'w')
        fcfg_out.writelines(fdat)
        fcfg_out.close()

def modFort4Prod(info=None):
    avgParam=info.avgParam
    lineNo=info.lineNo
    lineNoSwat=info.lineNoSwat

    for suffix in ['equil1','prod']:
        fcfg_in=open('fort.4.'+suffix)
        fdat=fcfg_in.readlines()
        fcfg_in.close()

        fdat[5]="     nstep              = %d\n" % avgParam.ncycles
        ncyclesn=avgParam.ncycles/10
        fdat[17]="     iprint       = %d\n" % ncyclesn
        fdat[19]="     iblock       = %d\n" % ncyclesn
        fdat[18]="     imv          = %d\n" % ncyclesn

        #fdat[lineNo['pmswat']]="     pmswat= %.5e\n" % avgParam.pmswat
        fdat[lineNo['pmswap']]="     pmswap= %.5e\n" % avgParam.pmswap
        fdat[lineNo['pmcb']]="     pmcb   = %.5e\n" % avgParam.pmcb
        fdat[lineNo['pmtra']]="     pmtra = %.5e\n" % avgParam.pmtra

        #fdat[lineNo['pmsatc']]="     pmsatc= %.5e 1.0d0\n" % (avgParam.pmsatc[0])
        #for ipair in range(2):
        #    fdat[lineNoSwat[ipair]]="   3     %.5e %.5e 1.0d0\n" % (avgParam.pmswtcb[ipair][0],avgParam.pmswtcb[ipair][1]) # swatch box pair prob

        #fdat[lineNo['pmswap']+1]="     pmswmt= %.5e 1.0d0\n" % (avgParam.pmswmt[0])
        #for imol in range(2):
        #    fdat[lineNo['pmswapb']+imol*6]="   3     0.0d0 %.5e 1.0d0\n" % avgParam.pmswapb[imol][0] # swap box pair prob

        fcfg_out=open('fort.4.'+suffix,'w')
        fcfg_out.writelines(fdat)
        fcfg_out.close()

class getLineNo:
    def __init__(self,curSim):
        #>>> Line numbers for the input file are assigned according to characteristic strings in the directory path; NOTE that Python arrays start with index 0
        if 'methanol' in curSim:
            self.lineNo={'ncmt':99,'pmvol':35,'pmswat':44,'pmsatc':46,'pmswap':51,'pmswapb':205,'pmcb':58,'pmtra':89}
            #self.lineNoSwat=[133,155]
            self.avgParam=su.avgResult(nbox=2,nmolty=1,nmoltyeta=0,nswapbp=0,nswatbp=0,npairbox=1)

#>>> If you want to group N sub.script together as one big job, set jobsPerGroup=N
#   Depending on how many cores one sub.script uses, set jobsPerNode=8 for serial jobs
#   or, for example, jobsPerNode=2 for parallel jobs that use 4 cores each
qbatch=su.qSubmitter(jobsPerGroup=-1,jobsPerNode=8,walltime=24*3600,coresPerNode=8)
#>>> The directory is expected to be organized as VLE/methanol/300K/RUNxx
sim=su.ProcDir(r'/VLE/.*/([\.\d]+)K$',r'RUN(\d+)$')
for curSim in sim:
    info=getLineNo(curSim)
    info.temp=float(sim.simMatch.group(1))
    info.press=0.101325 # in MPa
    info.nmol=[270,30]
    info.boxlen=[35.0,(info.nmol[1]*76.0944*10/6.022/0.0061403)**(1.0/3)] #g/mL
    print ("temp -> %f, press -> %f, boxlen -> " % (info.temp,info.press)),info.boxlen

    print "#####New Directory#####"
    run=sim.getRuns()
    for pathname in run:
        os.chdir(pathname)
        #>>> Suppose a shared topmon.inp file exists in the VLE/ directory
        os.system('ln -sf ../../../topmon.inp topmon.inp')
        #os.system('ln -sf ../../../TIP4P.xyz input_struc.xyz')

        (suffix,last,cur)=su.getSuffix('log')
        print "%s -> %s%d -> %d" % (pathname,suffix,last,cur)

        if os.path.exists('final-config'):
        #if False:
            print "volume too small!"
            last=cur
            shutil.move('final-config','fort.77')
            su.batchMove(suffix,last,nbox=2)
        elif os.path.exists('run1a.dat'):
            print "not finished!"
            modFort4(1,suffix,last,info)
            if suffix=='equil':
                last=cur
                shutil.move('save-config.1','fort.77')
                su.batchMove(suffix,last,nbox=2)
        else:
            #>>> To set up the input files intially
            #modFort4(0,suffix,last,info)
            #   Or to optimize simulation parameters for subsequent equilibration runs
            modFort4(2,suffix,last,info)

        if suffix=='equil':
            shutil.copy('fort.4.equil1','fort.4'); #equil
        elif suffix=='prod':
            shutil.copy('fort.4.prod','fort.4'); #prod

        #>>> A template sub.script file is expected to be found in /path/to/VLE,
        #   For the very first run, uncomment the next line
        #su.modPBS(pathname,'/path/to/VLE','melt cool volume equil1',1)
        #   For subsequent runs, use the following
        su.modPBS(pathname,'/path/to/VLE',suffix+str(last+1),1)

        qbatch.addJob(pathname+'/sub.script')

    #>>> The following part averages simulation parameters and maximum displacements for all the independent simulations
    #   Usually you would want to uncomment them after you have started equilibration runs
    ### START HERE ###
    #     info.avgParam.readMaximumDisplacement('fort.77')

    # print "#####Averaging#####"
    # try:
    #     info.avgParam.average()
    #     run=sim.getRuns()
    #     for pathname in run:
    #         os.chdir(pathname)
    #         (suffix,last,cur)=su.getSuffix('log')
    #         print "%s -> %s%d -> %d" % (pathname,suffix,last,cur)
    #         modFort4Prod(info)
    #         if suffix=='equil':
    #             shutil.copy('fort.4.equil1','fort.4'); #equil
    #         elif suffix=='prod':
    #             shutil.copy('fort.4.prod','fort.4'); #prod
    #         info.avgParam.writeMaximumDisplacement('fort.77',skipBox=[])
    # except (TypeError, ValueError) as e:
    #     print 'skip averaging..'
    #     continue
    # finally:
    #     print "###################\n"
    ### END HERE ###

qbatch.submit()
