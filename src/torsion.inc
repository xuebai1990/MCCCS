! TORSION.INC

	integer::splpnts(10),nttor
	real(8)::deg(1500,2),tabtorso(1500,2),torderiv2(1500,2)
	real(8)::tordif(1500,2)

	common /torsion/ deg,tabtorso,torderiv2,tordif,splpnts,nttor
