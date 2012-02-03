      
      integer beadtyp,decode,tor_num,angle_num
      integer tor_1,tor_2,tor_3,tor_4,nhere
      double precision  aframe,nframe,ang_bin_size 		
      double precision angle_bin, angle_tot,tor_bin,tor_tot
      double precision gdefect,btrans,bg_plus,bg_minus
      double precision g_plus,g_minus,nskip,trans	
      double precision biganalhist,shell,combiganalhist,comshell
      double precision angle_avg,nnone,comnone,cmolec
      double precision angle_1,angle_2,angle_3    
      double precision  tor_bin_size   			
      double precision end_to_end,max_length,boxcom_rhoz,rhoz
      double precision bigboxcom_rhoz,bigrhoz,max_boxlz,boxmin
	
C This is being removed because of high memory requirement

      double precision pattern  
	
      dimension angle_bin(nbxmax,ntmax,angle_max,ang_bin_max)
      dimension angle_tot(nbxmax,ntmax,angle_max) 
      dimension angle_avg(nbxmax,ntmax,angle_max)
      dimension tor_bin(nbxmax,ntmax,tor_max,tor_bin_max)	
      dimension tor_tot(nbxmax,ntmax,tor_max)
      dimension g_minus(nbxmax,ntmax),g_plus(nbxmax,ntmax)
      dimension btrans(nbxmax,ntmax,tor_max) ,bg_plus(nbxmax,ntmax,
     &       tor_max),bg_minus(nbxmax,ntmax,tor_max)
      dimension gdefect(nbxmax,ntmax,tor_max+1)
      dimension biganalhist(nbxmax*ntmax*ntdifmx*ntmax*ntdifmx,nbinmx)
      dimension shell(nbxmax*ntmax*ntmax*ntdifmx*ntdifmx,nbinmx,2)
      dimension combiganalhist(nbxmax*ntmax*ntmax,nbinmx)
      dimension comshell(nbxmax*ntmax*ntmax,nbinmx,2)
      dimension nnone(nbxmax,ntmax*ntmax*ntdifmx*ntdifmx)
      dimension comnone(nbxmax,ntmax*ntmax)
      dimension cmolec(ntmax,nbxmax)
      dimension angle_num(ntmax)
      dimension angle_1(ntmax,angle_max),angle_2(ntmax,angle_max)
     &     ,angle_3(ntmax,angle_max)
      dimension tor_1(ntmax,tor_max),tor_2(ntmax,tor_max)
     &     ,tor_3(ntmax,tor_max),tor_4(ntmax,tor_max)     
      dimension beadtyp(ntdifmx)

      dimension pattern(nbxmax,ntmax,0:3**tor_max)      

      dimension decode(nntype),tor_num(ntmax),trans(nbxmax,ntmax)	
      dimension end_to_end(nbxmax,ntmax,nbinmax_ete),max_length(ntmax)
      dimension boxcom_rhoz(nbxmax,ntmax,-nbinmax_ete:nbinmax_ete)
      dimension rhoz(nbxmax,ntmax,nbinmax_ete),boxmin(nbxmax)
      dimension bigboxcom_rhoz(nbxmax,ntmax,-nbinmax_ete:nbinmax_ete)
      dimension bigrhoz(nbxmax,ntmax,nbinmax_ete),max_boxlz(nbxmax)		

      common /gora/ aframe, nframe,trans,tor_1,tor_2,tor_3,tor_4,
     & ang_bin_size
      common /gorb/ angle_bin, angle_tot,tor_bin,tor_tot,gdefect,
     &      btrans,bg_plus,bg_minus,g_plus,g_minus,nskip,tor_num,nhere
      common /gorc/ biganalhist,shell,combiganalhist,comshell ,
     &     angle_avg,nnone,comnone,cmolec,angle_num
      common /gord/ end_to_end,max_length, boxcom_rhoz,rhoz,beadtyp,
     & decode
      common /gore/ tor_bin_size,max_boxlz,boxmin	
      common /gorf/ bigboxcom_rhoz, bigrhoz
      common /gorg/ angle_1,angle_2,angle_3       

  
      common /gorf/ pattern
