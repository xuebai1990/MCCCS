#!/usr/bin/env python

import sys,getopt,os,shutil

def distance(a,b):
    r1=[float(i) for i in a.split()]
    r2=[float(i) for i in b.split()]
    return (sum((p-q)*(p-q) for p,q in zip(r1,r2)))**0.5

def extract(finname,founame,beadList,scanList):
    '''Extract energies from Gaussian scan output files:
    /path/to/exGaussian.py --qmenergy -i gaussian.log -o qmenergy.txt [-x index.txt/--nbead N] -s scanList.txt
    REMEMBER: Gaussian uses the protein convention for dihedral angles (trans is 180 deg)'''

    EAU2K=315775.043382 #a.u. to K

    path=os.path.dirname(founame)
    fou=open(founame,'w')
    fin=open(finname)
    for line in fin:
        # Currently only one scan per run is supported
        if 'The following ModRedundant input section has been read:' in line:
            line=next(fin)
            while 'GradGrad' not in line:
                fields=line.split()
                if len(fields)>=8 and fields[0]=='D' and fields[5]=='S':
                    # expect the ModRedundant section to contain something like:
                    # D       1       3       4       5 S  36 5.0000
                    nStep=int(fields[6])
                    step=float(fields[7])
                    print 'Torsion scan: '
                elif len(fields)>=7 and fields[0]=='A' and fields[4]=='S':
                    # A       1       3       4 S  10 1.0000
                    nStep=int(fields[5])
                    step=float(fields[6])
                    print 'Angle scan: '
                line=next(fin)
        elif 'Initial Parameters' in line:
            for i in range(4):
                next(fin)
            line=next(fin)
            while '-------------' not in line:
                if 'Scan' in line:
                    fields=line.split()
                    idx=fields.index('Scan')
                    val0=float(fields[idx-1])
                    if fields[1] not in scanList:
                        scanList.append(fields[1])
                    print nStep,' steps of size ',step,', starting at ',val0
                    iStep=0
                line=next(fin)
        elif 'Input orientation' in line:
            for i in range(4):
                next(fin)
            mol=[]
            line=next(fin)
            while '-------------' not in line:
                fields=line.split()
                mol.append(' '.join(fields[3:])+'\n')
                line=next(fin)
        elif 'SCF Done:' in line:
            QMEnergy=float(line.split()[4])*EAU2K
        elif 'Optimization completed' in line or 'Optimization stopped' in line:
            val=val0+iStep*step
            iStep+=1
            fou.write('%.2f ' % (val))

            for i in range(7):
                next(fin)
            line=next(fin)
            while '-------------' not in line:
                fields=line.split()
                if fields[1] in scanList:
                    fou.write('%s ' % fields[3])
                line=next(fin)

            #>>> Output any other distances not monitored by Gaussian <<<
            #fou.write('%e %e ' % (distance(mol[0],mol[5]),distance(mol[1],mol[4]))) #EG
            #fou.write('%e %e ' % (distance(mol[0],mol[6]),distance(mol[1],mol[4]))) #PG and up

            fou.write('%.16f\n' % (QMEnergy))

            fname=os.path.join(path,str(iStep)+'.xyz')
            fcfg=open(fname,'w')
            fcfg.write('%3d out of %3d steps of size %.2f, starting at %.2f, current value is %.2f\n' %(iStep,nStep,step,val0,val))
            if len(beadList)==0:
                fcfg.writelines(mol)
            else:
                for i in beadList:
                    fcfg.write(mol[i])
            fcfg.close()
    fin.close()
    fou.close()
    print iStep,' nenergies extracted.'

def subtract(finname,founame,fexename):
    '''Subtract energies from MM runs using topmon:
    /path/to/exGaussian.py --mmenergy -i qmenergy.txt -o mmenergy.txt -e /path/to/topmon'''

    energies=[]
    lineNo=0
    for line in open(finname):
        lineNo+=1
        shutil.copy(str(lineNo)+'.xyz','input_struc.xyz')
        os.system(fexename)
        for line2 in open('run1a.dat'):
            if 'total energy' in line2:
                MMEnergy=float(line2.split()[-1])
                break
        diffEnergy=float(line.split()[-1])-MMEnergy
        if lineNo==1:
            minEnergy=diffEnergy
        else:
            minEnergy=min(minEnergy,diffEnergy)
        energies.append([line.replace('\n',' '),MMEnergy,diffEnergy])
    fou=open(founame,'w')
    for item in energies:
        fou.write('%s%.16f %.16f\n' %(item[0],item[1],item[2]-minEnergy))
    fou.close()

if __name__=='__main__':
    try:
        opts,args=getopt.getopt(sys.argv[1:],"i:o:x:s:e:",["qmenergy","mmenergy","input=","output=","index=","nbead=","scan=","exec="])
    except getopt.GetoptError,err:
        print str(err)
        sys.exit(-1)

    nbead=None
    fidxname=None
    fscanname=None
    for o,a in opts:
        if o == "--qmenergy":
            mode=1
        elif o == "--mmenergy":
            mode=2
        elif o in ("-i","--input"):
            finname=a
        elif o in ("-o","--output"):
            founame=a
        elif o in ("-x","--index"):
            fidxname=a
        elif o in ("--nbead"):
            nbead=int(a)
        elif o in ("-s","--scan"):
            fscanname=a
        elif o in ("-e","--exec"):
            fexename=a
        else:
            assert False,"ignored option"

    if mode is None or finname is None or founame is None or (mode==2 and fexename is None):
        assert False,"input error"

    if mode==1:
        beadList=[]
        if fidxname is not None:
            for line in open(fidxname):
                for item in line.split():
                    beadList.append(int(item)-1)
        elif nbead is not None:
            beadList=[i for i in range(nbead)]
        print 'beadList = ',beadList
        scanList=[]
        if fscanname is not None:
            for line in open(fscanname):
                for item in line.split():
                    scanList.append(item)
        print 'scanList = ',scanList
        extract(finname,founame,beadList,scanList)
    elif mode==2:
        subtract(finname,founame,fexename)
