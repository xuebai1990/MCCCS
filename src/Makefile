#########################
### Compiler definitions
#########################

### Intel Compiler for x86 architecture
#COMPFLAGS      =  -w90  -check all -traceback 
#COMPFLAGS     =     -w90 -w95 -cm -check all -traceback#-check_bounds# -g
#COMPFLAGS     =  -O2 -w90 -w95 -cm -p#-check_bounds# -g
## only for debugging
#COMPFLAGS     = -w90 -w95 -cm -CB -CU -g

#COMPFLAGS     =  -O2 -w90 -w95 -cm -p
COMPFLAGS     =  -O3 -no-prec-div -ipo #-p  #-ipo
LD            = ifort
LDS           = ifort -c $(COMPFLAGS)
LDFLAGS       = -static $(COMPFLAGS)

### NAG Compiler for OSX
#COMPFLAGS     = -O1 -gline -f77 -mismatch_all
## only for debugging
#COMPFLAGS     = -O2 
#LD            = f95
#LDS           = f95 -c $(COMPFLAGS)
#LDFLAGS       = $(COMPFLAGS)

### SGI Altix: Intel Compiler for ia64 architecture
#COMPFLAGS     = -O2 -tpp2 -w90 -w95 -cm #-check all -traceback
# only for debugging
#COMPFLAGS     =  -check all -traceback -w90 -w95 -cm -CB -CU -g
#LD            = ifort
#LDS           = ifort -c $(COMPFLAGS)
#LDFLAGS       = -static $(COMPFLAGS)

### GNU Compiler
#COMPFLAGS     =  -ftrace=full -pg -fbounds-check
#COMPFLAGS     = -O3 -fomit-frame-pointer -finline-functions -ffast-math 
#LD	      = g95
#LDS 	      = g95 -c $(COMPFLAGS)
#LDFLAGS	      = $(COMPFLAGS)

### Portland Compiler for x86 architecture
#COMPFLAGS     =  -O2 -tp k8-64 -Mvect=sse#-check_bounds# -g
#COMPFLAGS     =  -tp k8-64 -Mbounds -g -Mchkptr#-check_bounds# -g
#LD	      = /opt/pgi/linux86-64/6.1/bin/pgf95
#LDS 	      = /opt/pgi/linux86-64/6.1/bin/pgf95 -c $(COMPFLAGS)
#LDFLAGS	      = $(COMPFLAGS)

### SGI Compiler for Origin
#COMPFLAGS     = -O2 #-march=opteron # -mips4 -r12000 #-check_bounds# -g
#LD            = f77 
#LDS           = f77 -c $(COMPFLAGS)
#LDFLAGS       = $(COMPFLAGS)

### IBM AIX Compiler for SP
#COMPFLAGS     = -O3 -qstrict  -bmaxdata:2000000000 #-check_bounds# -g
#LD            = xlf
#LDS           = xlf -c $(COMPFLAGS)
#LDFLAGS	      = $(COMPFLAGS)

#############################
### End compiler definitions 
#############################

LIBS	      =

LINTLIBS      =

MAKEFILE      = Makefile

OBJS	      = Atom_traxyz.o \
		Atom_energy.o \
		analysis.o \
		atomtype.o \
		anes.o \
		boltz.o \
		bondlength.o \
		calcsolpar.o \
		calctor.o \
		charge.o \
		close.o \
		cone.o \
		coneangle.o \
		config.o \
		corp.o \
		coru.o \
		coruz.o \
		ctrmas.o \
		dipole.o \
		dump.o \
		energy.o \
		erfunc.o \
		expand.o \
		exsix.o \
                exgrph.o \
		exmuir.o \
		explct.o \
		exsami.o \
		exzeo.o \
		exzeof.o \
		flucq.o \
		forcefield.o \
		geometry.o \
		inclus.o \
		initia.o \
		linkcell.o \
		ljmuir.o \
		ljpsur.o \
		ljsami.o \
		matops.o \
                mbessel.o \
		mmff.o \
		mimage.o \
		molsetup.o \
		monola.o \
		monper.o \
                ninesix.o \
		nrtab.o \
		place.o \
		pgrid.o \
		pressure.o \
		prvolume.o \
		random.o \
		readdat.o \
		recip.o \
		recip_atom.o \
		recippress.o \
		recipsum.o \
		rigfix.o \
		rigrot.o \
		rosenbluth.o \
		rotxyz.o \
		rwztab.o \
		safecbmc.o \
		safeschedule.o \
		setnn.o \
		setpbc.o \
		schedule.o \
		slitpore.o \
		sphere.o \
		suijtab.o \
		sumup.o \
		susami.o \
		suvibe.o \
		suzeo.o \
		swatch.o \
		swap.o \
		trim_swap.o \
		topmon.o \
		traxyz.o \
		update.o \
		update1.o \
		updnn.o \
		virial.o\
		volume.o \
		vtorso.o \
		writepdb.o \
		zeocoord.o \
		ztab.o \
		ztest.o \
		Intra_energy.o  \
                deriv.o \
                ipswsetup.o \
                eemove.o \
                eesetup.o \
		ee_index_swap.o \
                ee_energy.o \
		ee_recip.o \
		v_elect_field.o \
		torsion.o \
		garofalini.o \
		garo3body.o

PROGRAM       = topmon4 

SHELL	      = /bin/sh

SRCS	     =  Atom_traxyz.f \
		Atom_energy.f \
		analysis.f \
		atomtype.f \
		anes.f \
		angle.f \
		boltz.f \
		bondlength.f \
		calcsolpar.f \
		calctor.f \
		charge.f \
		close.f \
		cone.f \
		coneangle.f \
		config.f \
		corp.f \
		coru.f \
		coruz.f \
		ctrmas.f \
		dipole.f \
		dump.f \
		energy.f \
		erfunc.f \
		expand.f \
                exgrph.f \
		exsix.f \
		exmuir.f \
		exsami.f \
		explct.f \
		exzeo.f \
		exzeof.f \
		flucq.f \
		forcefield.f \
		geometry.f \
		initia.f \
		inclus.f \
		linkcell.f \
		ljmuir.f \
		ljpsur.f \
		ljsami.f \
		matops.f \
                mbessel.f \
		mmff.f \
		mimage.f \
		molsetup.f \
		monola.f \
		monper.f \
                ninesix.f \
		nrtab.f \
		pgrid.f \
		place.f \
		pressure.f \
		prvolume.f \
		random.f \
		rigfix.f \
		rigrot.f \
		readdat.f \
		recip.f \
		recip_atom.f \
		recippress.f \
		recipsum.f \
		rosenbluth.f \
		rotxyz.f \
		rwztab.f \
		safecbmc.f \
		safeschedule.f \
		setnn.f \
		setpbc.f \
		schedule.f \
		slitpore.f \
		sphere.f \
		suijtab.f \
		sumup.f \
		susami.f \
		suvibe.f \
		suzeo.f \
		swatch.f \
		swap.f \
		trim_swap.f \
		topmon.f \
		traxyz.f \
		update.f \
		update1.f \
		updnn.f \
		virial.f\
		volume.f \
		vtorso.f \
		writepdb.f \
		zeocoord.f \
		ztab.f \
		ztest.f \
		Intra_energy.f \
                deriv.f \
                ipswsetup.f \
                eemove.f \
                eesetup.f \
                ee_index_swap.f \
                ee_energy.f \
                ee_recip.f \
		v_elect_field.f \
		torsion.f \
		garofalini.f \
		garo3body.f

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS) $(MAKEFILE)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "done"

###
Atom_traxyz.o: coord.inc coord2.inc ensemble.inc neigh2.inc system.inc \
        control.inc inputdata.inc bnbsma.inc fepsi.inc Atom_traxyz.f
	$(LDS) Atom_traxyz.f
Atom_energy.o: control.inc coord.inc system.inc neigh.inc poten.inc coord2.inc \
        external.inc connect.inc ewaldsum.inc fepsi.inc qqlist.inc nsix.inc \
        peboco.inc cell.inc conver.inc Atom_energy.f
	$(LDS) Atom_energy.f
analysis.o: control.inc coord.inc  inputdata.inc connect.inc gor.inc \
	system.inc cell.inc analysis.f
	$(LDS) analysis.f
atomtype.o: atomtype.f
	$(LDS) atomtype.f
anes.o: control.inc coord.inc coord2.inc ensemble.inc rosen.inc \
	inputdata.inc system.inc anes.f
	$(LDS) anes.f
boltz.o: boltz.f control.inc coord.inc system.inc poten.inc neigh.inc \
	cbmc.inc external.inc connect.inc ewaldsum.inc cell.inc \
	peboco.inc fepsi.inc qqlist.inc nsix.inc ipswpar.inc eepar.inc conver.inc boltz.f
	$(LDS) boltz.f
bondlength.o: bondlength.f
	$(LDS) bondlength.f
calcsolpar.o: calcsolpar.f control.inc inputdata.inc blkavg.inc coord.inc \
	ensemble.inc system.inc cell.inc calcsolpar.f 
	$(LDS) calcsolpar.f
calctor.o: control.inc fix.inc calctor.f
	$(LDS) calctor.f
charge.o: control.inc coord.inc poten.inc conver.inc ewaldsum.inc charge.f
	$(LDS) charge.f
close.o : close.f
	$(LDS) close.f
cone.o: conver.inc cone.f
	$(LDS) cone.f
coneangle.o: coneangle.f
	$(LDS) coneangle.f
config.o: control.inc coord.inc coord2.inc system.inc ensemble.inc cbmc.inc \
	rosen.inc inputdata.inc ewaldsum.inc poten.inc ipswpar.inc \
	eepar.inc config.f
	$(LDS) config.f
corp.o: control.inc coord.inc conver.inc poten.inc system.inc expsix.inc \
	merck.inc corp.f 
	$(LDS) corp.f
coru.o: control.inc coord.inc poten.inc conver.inc system.inc expsix.inc \
	merck.inc coru.f
	$(LDS) coru.f
coruz.o: control.inc conver.inc system.inc zeopoten.inc coruz.f
	$(LDS) coruz.f
ctrmas.o: control.inc coord.inc system.inc cell.inc ctrmas.f
	$(LDS) ctrmas.f
dipole.o: control.inc coord.inc coord2.inc ewaldsum.inc dipole.f
	$(LDS) dipole.f
dump.o: control.inc coord.inc poten.inc system.inc ensemble.inc cbmc.inc \
	conver.inc inpar.inc external.inc zeolite.inc \
	inputdata.inc blkavg.inc bnbsma.inc swtcmove.inc ewaldsum.inc \
	cell.inc gor.inc dump.f
	$(LDS) dump.f
energy.o: control.inc coord.inc system.inc neigh.inc poten.inc coord2.inc \
	external.inc connect.inc ewaldsum.inc fepsi.inc qqlist.inc nsix.inc 	peboco.inc cell.inc ipswpar.inc eepar.inc conver.inc garofalini.inc energy.f
	$(LDS) energy.f
erfunc.o: erfunc.f
	$(LDS) erfunc.f
exgrph.o: control.inc external.inc poten.inc exgrph.f
	$(LDS) exgrph.f
expand.o: control.inc coord.inc coord2.inc system.inc ensemble.inc \
	bnbsma.inc inputdata.inc expand.inc poten.inc expand.f
	$(LDS) expand.f
exsix.o: control.inc expsix.inc exsix.f
	$(LDS) exsix.f
exmuir.o: external.inc externalmuir.inc exmuir.f
	$(LDS) exmuir.f
exsami.o: external.inc exsami.f
	$(LDS) exsami.f
explct.o: control.inc coord.inc explct.f
	$(LDS) explct.f
exzeo.o: zeopoten.inc zeolite.inc control.inc grid.inc external.inc exzeo.f
	$(LDS) exzeo.f
exzeof.o: zeopoten.inc zeolite.inc control.inc exzeof.f
	$(LDS) exzeof.f
flucq.o: control.inc coord.inc coord2.inc ensemble.inc system.inc \
	inputdata.inc bnbsma.inc ewaldsum.inc fepsi.inc poten.inc flucq.f
	$(LDS) flucq.f
geometry.o: control.inc coord.inc connect.inc conver.inc \
	rosen.inc cbmc.inc geometry.f
	$(LDS) geometry.f
forcefield.o: zeolite.inc zeopoten.inc control.inc forcefield.f
	$(LDS) forcefield.f
inclus.o: control.inc coord.inc connect.inc inclus.f
	$(LDS) inclus.f
initia.o: control.inc coord.inc system.inc poten.inc conver.inc inpar.inc \
	ensemble.inc connect.inc cbmc.inc initia.f
	$(LDS) initia.f
linkcell.o: control.inc coord.inc system.inc linkcell.f
	$(LDS) linkcell.f
ljmuir.o: control.inc coord.inc poten.inc ljmuir.f
	$(LDS) ljmuir.f
ljpsur.o: external.inc ljpsur.f
	$(LDS) ljpsur.f
ljsami.o: external.inc ljsamipara.inc ljsami.f
	$(LDS) ljsami.f
matops.o: cell.inc coord.inc conver.inc system.inc matops.f
	$(LDS) matops.f
mbessel.o: mbessel.f
	$(LDS) mbessel.f
mmff.o: control.inc merck.inc mmff.f
	$(LDS) mmff.f
mimage.o: peboco.inc control.inc system.inc cell.inc mimage.f
	$(LDS) mimage.f
molsetup.o: control.inc coord.inc connect.inc molsetup.f
	$(LDS) molsetup.f
monola.o: control.inc coord.inc poten.inc system.inc ensemble.inc cbmc.inc \
	conver.inc inpar.inc external.inc zeolite.inc \
	inputdata.inc blkavg.inc bnbsma.inc swtcmove.inc ewaldsum.inc \
	cell.inc gor.inc ipswpar.inc eepar.inc monola.f
	$(LDS) monola.f
monper.o: control.inc coord.inc system.inc ensemble.inc cbmc.inc \
	inputdata.inc blkavg.inc bnbsma.inc connect.inc coord2.inc cell.inc \
	ipswpar.inc poten.inc monper.f
	$(LDS) monper.f
ninesix.o: control.inc nsix.inc ninesix.f
	$(LDS) ninesix.f
nrtab.o: control.inc coord.inc poten.inc connect.inc nrtab.inc nrtab.f
	$(LDS) nrtab.f
place.o: control.inc coord.inc cbmc.inc connect.inc fix.inc system.inc \
	rosen.inc place.f
	$(LDS) place.f
pgrid.o: pgrid.f
	$(LDS) pgrid.f
pressure.o: control.inc coord.inc system.inc poten.inc expsix.inc \
	merck.inc connect.inc ewaldsum.inc conver.inc \
	inputdata.inc fepsi.inc qqlist.inc ipswpar.inc pressure.f
	$(LDS) pressure.f
prvolume.o: control.inc coord.inc ensemble.inc system.inc \
	    inputdata.inc  bnbsma.inc ewaldsum.inc cell.inc garofalini.inc prvolume.f
	$(LDS) prvolume.f
random.o: random.f
	$(LDS) random.f
readdat.o: control.inc coord.inc cbmc.inc system.inc conver.inc poten.inc \
	inpar.inc external.inc externalmuir.inc zeolite.inc\
	nrtab.inc connect.inc inputdata.inc ewaldsum.inc \
	swtcmove.inc fepsi.inc expand.inc qqlist.inc clusterbias.inc \
	cell.inc gor.inc garofalini.inc conver.inc readdat.f
	$(LDS) readdat.f
rigfix.o: control.inc coord.inc cbmc.inc rosen.inc fix.inc connect.inc \
	system.inc rigfix.f 
	$(LDS) rigfix.f
rigrot.o: control.inc cbmc.inc coord.inc coord2.inc rosen.inc \
	system.inc ipswpar.inc rigrot.f
	$(LDS) rigrot.f
recip.o: control.inc coord.inc coord2.inc ewaldsum.inc poten.inc recip.f
	$(LDS) recip.f
recip_atom.o: control.inc coord.inc coord2.inc ewaldsum.inc poten.inc recip_atom.f
	 $(LDS) recip_atom.f
recippress.o: control.inc ewaldsum.inc coord.inc recippress.f
	$(LDS) recippress.f
recipsum.o: control.inc system.inc coord.inc ewaldsum.inc poten.inc \
	conver.inc cell.inc recipsum.f
	$(LDS) recipsum.f
rosenbluth.o: coord.inc system.inc ensemble.inc poten.inc conver.inc \
	cbmc.inc rosen.inc control.inc connect.inc fix.inc \
	ipswpar.inc rosenbluth.f
	$(LDS) rosenbluth.f
rotxyz.o: coord.inc coord2.inc ensemble.inc neigh2.inc system.inc \
	control.inc inputdata.inc bnbsma.inc fepsi.inc ipswpar.inc \
	eepar.inc rotxyz.f 
	$(LDS) rotxyz.f
rwztab.o: grid.inc zeolite.inc rwztab.f
	$(LDS) rwztab.f
safecbmc.o: coord.inc control.inc fix.inc connect.inc rosen.inc cbmc.inc \
	safecbmc.f
	$(LDS) safecbmc.f
safeschedule.o: coord.inc coord2.inc control.inc inputdata.inc rosen.inc \
	fix.inc connect.inc safeschedule.f
	$(LDS) safeschedule.f
setnn.o: coord.inc neigh2.inc neigh.inc system.inc control.inc setnn.f
	$(LDS) setnn.f
setpbc.o: system.inc peboco.inc control.inc setpbc.f
	$(LDS) setpbc.f
schedule.o: control.inc coord.inc cbmc.inc rosen.inc connect.inc \
	inputdata.inc schedule.f
	$(LDS) schedule.f
slitpore.o: poten.inc external.inc slitpore.f
	$(LDS) slitpore.f  	
sphere.o: sphere.f
	$(LDS) sphere.f
suijtab.o: control.inc coord.inc external.inc poten.inc system.inc \
	cell.inc merck.inc nsix.inc suijtab.f
	$(LDS) suijtab.f
sumup.o: coord.inc system.inc poten.inc conver.inc control.inc \
	external.inc zeolite.inc connect.inc inputdata.inc fepsi.inc\
        inputdata.inc qqlist.inc cell.inc nsix.inc peboco.inc ipswpar.inc \
	eepar.inc garofalini.inc sumup.f
	$(LDS) sumup.f
susami.o: external.inc ljsamipara.inc susami.f
	$(LDS) susami.f
suvibe.o: control.inc coord.inc connect.inc conver.inc contorsion.inc \
	suvibe.f
	$(LDS) suvibe.f
suzeo.o: zeopoten.inc zeolite.inc external.inc suzeo.f
	$(LDS) suzeo.f
swatch.o: coord.inc coord2.inc system.inc ensemble.inc cbmc.inc control.inc \
	inputdata.inc rosen.inc swtcmove.inc ewaldsum.inc cell.inc swatch.f
	$(LDS) swatch.f
swap.o: coord.inc coord2.inc system.inc ensemble.inc cbmc.inc rosen.inc \
	boltzmann.inc control.inc external.inc connect.inc fepsi.inc \
	inputdata.inc clusterbias.inc cell.inc eepar.inc swap.f
	$(LDS) swap.f
trim_swap.o: coord.inc coord2.inc system.inc ensemble.inc cbmc.inc rosen.inc \
        boltzmann.inc control.inc external.inc connect.inc fepsi.inc \
        inputdata.inc clusterbias.inc cell.inc trim_swap.f
	$(LDS) trim_swap.f
topmon.o: topmon.f
	$(LDS) topmon.f
traxyz.o: coord.inc coord2.inc ensemble.inc neigh2.inc system.inc \
	control.inc inputdata.inc bnbsma.inc fepsi.inc ipswpar.inc \
	eepar.inc garofalini.inc traxyz.f
	$(LDS) traxyz.f
update.o: coord.inc system.inc neigh.inc neigh2.inc control.inc blkavg.inc \
	update.f
	$(LDS) update.f
update1.o: coord.inc system.inc neigh.inc neigh2.inc control.inc blkavg.inc \
        update1.f
	$(LDS) update1.f
updnn.o: coord.inc system.inc neigh.inc neigh2.inc control.inc updnn.f
	$(LDS) updnn.f
virial.o: coord.inc system.inc poten.inc conver.inc control.inc \
	external.inc connect.inc eepar.inc virial.f
	$(LDS) virial.f
volume.o: coord.inc ensemble.inc system.inc control.inc \
	  inputdata.inc bnbsma.inc ewaldsum.inc cell.inc garofalini.inc neigh.inc volume.f
	$(LDS) volume.f
vtorso.o: conver.inc contorsion.inc vtorso.f
	$(LDS) vtorso.f
writepdb.o: control.inc zeolite.inc coord.inc writepdb.f
	$(LDS) writepdb.f
zeocoord.o: zeolite.inc zeopoten.inc zeocoord.f
	$(LDS) zeocoord.f
ztab.o: grid.inc zeolite.inc ztab.f
	$(LDS) ztab.f
ztest.o: grid.inc zeolite.inc ztest.f
	$(LDS) ztest.f
Intra_energy.o: control.inc coord.inc system.inc neigh.inc poten.inc coord2.inc \
	external.inc connect.inc ewaldsum.inc fepsi.inc qqlist.inc nsix.inc Intra_energy.f
	$(LDS) Intra_energy.f
deriv.o: control.inc system.inc coord.inc ipswpar.inc deriv.f
	$(LDS) deriv.f
ipswsetup.o: control.inc coord.inc ipswpar.inc inputdata.inc system.inc ipswsetup.f
	$(LDS) ipswsetup.f
eemove.o: coord.inc coord2.inc system.inc ensemble.inc cbmc.inc rosen.inc \
	boltzmann.inc control.inc external.inc connect.inc fepsi.inc \
	inputdata.inc clusterbias.inc cell.inc eepar.inc eemove.f
	$(LDS) eemove.f
eesetup.o: control.inc coord.inc eepar.inc inputdata.inc system.inc \
	nsix.inc poten.inc eesetup.f
	$(LDS) eesetup.f
ee_index_swap.o: control.inc coord.inc eepar.inc ee_index_swap.f
	$(LDS) ee_index_swap.f
ee_energy.o: control.inc coord.inc system.inc neigh.inc poten.inc coord2.inc \
	cell.inc external.inc connect.inc ewaldsum.inc fepsi.inc qqlist.inc \
	nsix.inc eepar.inc ee_energy.f
	$(LDS) ee_energy.f
ee_recip.o: control.inc coord.inc coord2.inc ewaldsum.inc poten.inc ee_recip.f
	$(LDS) ee_recip.f
v_elect_field.o: control.inc external.inc coord.inc v_elect_field.f
	$(LDS) v_elect_field.f
torsion.o:  torsion.inc conver.inc torsion.f
	$(LDS) torsion.f
garofalini.o: control.inc garofalini.inc poten.inc garofalini.f
	$(LDS) garofalini.f
garo3body.o: control.inc coord.inc garofalini.inc neigh.inc garo3body.f
	$(LDS) garo3body.f

##### remove object files
clean:
	rm -f *.o


