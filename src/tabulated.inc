
      integer vibsplits(150), bendsplits(150), ntabvib, ntabbend
      integer vdWsplits(150,150), ntabvdW
      integer electsplits(150,150), ntabelect
      integer num_int_vib(150), num_int_bend(150), 
     &	num_int_vdW(150,150), num_int_elect(150,150)
      double precision vib(1500,150), tabvib(1500,150)
      double precision bend(1500,150), tabbend(1500,150)
      double precision rvdW(1500,150,150), tabvdW(1500,150,150)
      double precision relect(1500,150,150), tabelect(1500,150,150)
      double precision vibdiff(1500,150), benddiff(1500,150), 
     & 	vdWdiff(1500,150,150), electdiff(1500,150,150)
      logical L_vib_table, L_bend_table, L_vdW_table, 
     &  L_elect_table


      common /tabulated/ vib, tabvib, bend, tabbend, rvdW,
     & 	 tabvdW, relect, tabelect, vibdiff, benddiff, vdWdiff, 
     &   electdiff, electsplits, ntabelect, num_int_elect,
     & 	 vibsplits, bendsplits, vdWsplits, ntabvib, ntabbend,
     &   ntabvdW, num_int_vib, num_int_bend, num_int_vdW 
      common /ltab/ L_vib_table, L_bend_table, L_vdW_table,
     &   L_elect_table