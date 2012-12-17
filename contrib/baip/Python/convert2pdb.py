def checkMolecule(molSpec,printMask):
    failed=False
    if len(molSpec)!=len(printMask):
        print "Number of molecules don't match in molSpec and printMask"
        failed=True
    else:
        for i,j in zip(molSpec,printMask):
            if len(i)!=len(j):
                print "Number of atoms don't match: molSpec ",i,"; printMask ",j
                failed=True
    return failed

def xyz2pdb(argv,molSpec,printMask):
    failed=checkMolecule(molSpec,printMask)

    if len(argv) != 3:
        print "Usage: %s /path/to/xyz /path/to/pdb\n" % argv[0]
        failed=True

    if failed:
        quit()

    fou=file(argv[2],"w")

    try:
        for lineStr in open(argv[1]):
            line=lineStr.split()
            if len(line)==0:
                pass
            elif len(line)==1:
                natom=0
                nmol=0
                atoms=[]
                coords=[]
            elif len(line)==4:
                atoms.append(line[0])
                coords.append([float(line[1]),float(line[2]),float(line[3])])
                for i,j in enumerate(molSpec):
                    if j==atoms:
                        printed=False
                        for a,m,c in zip(atoms,printMask[i],coords):
                            if m!=0:
                                printed=True
                                natom+=1
                                fou.write("ATOM  %5d  %-4s%-5d   %-5d%8.3f%8.3f%8.3f  1.00  0.00\n" % (natom,a,i+1,nmol+1,c[0],c[1],c[2]))
                        if printed:
                            nmol+=1
                        atoms=[]
                        coords=[]
            else:
                print "something goes wrong!\n"
                quit()
    finally:
        fou.write("END\n")
        fou.close()
        print nmol," molecules, ",natom," atoms written."
        if len(atoms)==0:
            print "Finished."
        else:
            print "Not done! ",atoms

def cfg2pdb(argv,molSpec,printMask,nbox=3):
    import numpy as np
    import fileutil as fu

    failed=checkMolecule(molSpec,printMask)

    if len(argv) != 2:
        print "Usage: %s /path/to/cfg\n" % argv[0]
        failed=True

    if failed:
        quit()

    fou=[]
    for ibox in range(nbox):
        fname=argv[1]+'.box'+str(ibox+1)+'.pdb'
        print fname
        fou.append(file(fname,'w'))

    nmolty=len(molSpec)

    try:
        fin=fu.ReadFile(argv[1])
        ncycle=fin.readline()
        armtra=fin.readline()
        for ibox in range(nbox):
            for imol in range(nmolty):
                rmtra=fin.readline()
                rmrot=fin.readline()
        for ibox in range(nbox):
            for imol in range(nmolty):
                rmflcq=fin.readNum()
        for ibox in range(nbox):
            rmvol=fin.readNum()
        for ibox in range(nbox):
            boxl=fin.readline()
        nchain=fin.readNum()
        nmolty_tmp=fin.readNum()
        print "nchain = ",nchain,"; nmolty = ",nmolty_tmp
        if nmolty_tmp!=nmolty:
            print "Number of molecule types don't match"
            quit()
        nunit=np.zeros((nmolty),dtype=np.int)
        moltyp=np.zeros((nchain),dtype=np.int)
        nboxi=np.zeros((nchain),dtype=np.int)
        for imol in range(nmolty):
            nunit[imol]=fin.readNum()
            if nunit[imol]!=len(molSpec[imol]):
                print "Number of atoms don't match: molSpec ",i,"; printMask ",j,"; nunit = ",nunit[imol]
                quit()
        for i in range(nchain):
            moltyp[i]=fin.readNum() - 1
        for i in range(nchain):
            nboxi[i]=fin.readNum() - 1
        natom=np.zeros((nbox),dtype=np.int)
        nmol=np.zeros((nbox),dtype=np.int)
        for i in range(nchain):
            printed=False
            for j in range(nunit[moltyp[i]]):
                rxu=fin.readNum()
                ryu=fin.readNum()
                rzu=fin.readNum()
                qqu=fin.readNum()
                if printMask[moltyp[i]][j]!=0:
                    printed=True
                    natom[nboxi[i]]+=1
                    fou[nboxi[i]].write("ATOM  %5d  %-4s%-5d   %-5d%8.3f%8.3f%8.3f  1.00  0.00\n" % (natom[nboxi[i]],molSpec[moltyp[i]][j],moltyp[i]+1,nmol[nboxi[i]]+1,rxu,ryu,rzu))
            if printed:
                nmol[nboxi[i]]+=1
    except IndexError:
        print i,j,nboxi[i]+1,moltyp[i]+1,nunit[moltyp[i]]
    except EOFError:
        print "Unexpected end of file\n"
    finally:
        for ibox in range(nbox):
            fou[ibox].write("END\n")
            fou[ibox].close()
            print "Box ",ibox+1,": ",nmol[ibox]," molecules, ",natom[ibox]," atoms written."
        print "Finished."
