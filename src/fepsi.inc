! FEPSI.INC

      logical::lfepsi
      real(8)::aslope,bslope,ashift,bshift,a0,b0

! ----- If lfepsi is true, the fluctuation of epsilon is used instead of
! ----- the fluctuation of the sigma.
 
        parameter (lfepsi = .false.)
! *** slope=0.3 a=3.05
!	parameter (a0 = 0.2003d0)
!	parameter (b0 = 1.3946d0)
!	parameter (aslope = 8.85d5)
!	parameter (bslope = 158.25d0)
!	parameter (ashift = 7.227d5)
!	parameter (bshift = 505.97d0)
! *** slope=0.3 a=2.90
!	parameter (a0 = 0.15561d0)
!	parameter (b0 = 1.2960d0)
!	parameter (aslope = 5.5475d5)
!	parameter (bslope = 131.23d0)
!	parameter (ashift = 4.121d5)
!	parameter (bshift = 381.96d0)
! *** slope=0.3 a=2.95
!	parameter (a0 = 0.16833d0)
!	parameter (b0 = 1.3299d0)
!	parameter (aslope = 6.5125d5)
!	parameter (bslope = 139.75d0)
!	parameter (ashift = 4.9975d5)
!	parameter (bshift = 419.81d0)
! *** slope=0.3 a=2.82
!	parameter (a0 = 0.13150d0)
!	parameter (b0 = 1.2574d0)
!	parameter (aslope = 4.2863d5)
!	parameter (bslope = 117.5d0)
!	parameter (ashift = 3.0208d5)
!	parameter (bshift = 323.6d0)
! *** slope=0.3 a=2.75
!	parameter (a0 = 0.11037d0)
!	parameter (b0 = 1.1942d0)
!	parameter (aslope = 3.4013d5)
!	parameter (bslope = 108d0)
!	parameter (ashift = 2.2868d5)
!	parameter (bshift = 285.04d0)
! *** slope=0.3 a=2.85
!	parameter (a0 = 0.14035d0)
!	parameter (b0 = 1.2596d0)
!	parameter (aslope = 4.7263d5)
!	parameter (bslope = 123.25d0)
!	parameter (ashift = 3.3978d5)
!	parameter (bshift = 347.63d0)
! *** slope=0.3 a=2.88
!	parameter (a0 = 0.14820d0)
!	parameter (b0 = 1.2689d0)
!	parameter (aslope = 5.2125d5)
!	parameter (bslope = 128.75d0)
!	parameter (ashift = 3.8220d5)
!	parameter (bshift = 371.29d0)
! *** slope=0.5 a=2.71
      parameter (a0 = -0.22818d0)
      parameter (b0 = 0.44662d0)
      parameter (aslope = 14.0738d5)
      parameter (bslope = 351.25d0)
      parameter (ashift = 3.5852d5)
      parameter (bshift = 358.94d0)
! *** slope=0.5 a=2.68
!	parameter (a0 = -0.23379d0)
!	parameter (b0 = 0.43630d0)
!	parameter (aslope = 12.78125d5)
!	parameter (bslope = 337.5d0)
!	parameter (ashift = 3.1904d5)
!	parameter (bshift = 337.86d0)
! *** slope=0.5 a=2.655
!	parameter (a0 = -0.23808d0)
!	parameter (b0 = 0.425507d0)
!	parameter (aslope = 11.7775d5)
!	parameter (bslope = 322.04d0)
!	parameter (ashift = 2.8902d5)
!	parameter (bshift = 326.875d0)
! *** slope=0.5 a=2.665
!	parameter (a0 = -0.23641d0)
!	parameter (b0 = 0.43409d0)
!	parameter (aslope = 12.171d5)
!	parameter (bslope = 330d0)
!	parameter (ashift = 3.0075d5)
!	parameter (bshift = 326.52d0)
! *** slope=0.5 a=2.65
!	parameter (a0 = -0.23903d0)
!	parameter (b0 = 0.42231d0)
!	parameter (aslope = 11.5875d5)
!	parameter (bslope = 325d0)
!	parameter (ashift = 2.8339d5)
!	parameter (bshift = 319.44d0)
! *** slope=0.5 a=2.60
!	parameter (a0 = -0.24796d0)
!	parameter (b0 = 0.40345d0)
!	parameter (aslope = 9.8238d5)
!	parameter (bslope = 304d0)
!	parameter (ashift = 2.3206d5)
!	parameter (bshift = 288.72d0)
! *** slope=0.5 a=2.55
!	parameter (a0 = -0.25703d0)
!	parameter (b0 = 0.38299d0)
!	parameter (aslope = 8.3075d5)
!	parameter (bslope = 284.38d0)
!	parameter (ashift = 1.8945d5)
!	parameter (bshift = 261.09d0)
! *** slope = 0.7 a=2.52
!	parameter (a0 = -0.39106d0)
!	parameter (b0 = 0.08431d0)
!	parameter (aslope = 25.22875d5)
!	parameter (bslope = 662.75d0)
!	parameter (ashift = 3.0689d5)
!	parameter (bshift = 335.43d0)
! *** slope = 0.7 a=2.51
!	parameter (a0 = -0.39233d0)
!	parameter (b0 = 0.08203d0)
!	parameter (aslope = 24.425d5)
!	parameter (bslope = 653.75d0)
!	parameter (ashift = 2.9499d5)
!	parameter (bshift = 328.60d0)
! *** slope = 0.7 a=2.50
!	parameter (a0 = -0.39357d0)
!	parameter (b0 = 0.07946d0)
!	parameter (aslope = 23.641d5)
!	parameter (bslope = 645d0)
!	parameter (ashift = 2.835d5)
!	parameter (bshift = 322.13d0)
! *** slope = 0.7 a=2.49
!	parameter (a0 = -0.394906d0)
!	parameter (b0 = 0.075678d0)
!	parameter (aslope = 22.88625d5)
!	parameter (bslope = 637.735d0)
!	parameter (ashift = 2.72469d5)
!	parameter (bshift = 316.22d0)
! *** slope = 0.7 a=2.48
!	parameter (a0 = -0.396175d0)
!	parameter (b0 = 0.072973d0)
!	parameter (aslope = 22.149125d5)
!	parameter (bslope = 629d0)
!	parameter (ashift = 2.61798d5)
!	parameter (bshift = 309.95d0)
! *** slope = 0.7 a=2.47
!	parameter (a0 = -0.39745d0)
!	parameter (b0 = 0.070157d0)
!	parameter (aslope = 21.4335d5)
!	parameter (bslope = 620.75d0)
!	parameter (ashift = 2.5151d5)
!	parameter (bshift = 306.89d0)
! *** slope = 0.7 a=2.45
!	parameter (a0 = -0.399983d0)
!	parameter (b0 = 0.064159d0)
!	parameter (aslope = 20.06425d5)
!	parameter (bslope = 604.75d0)
!	parameter (ashift = 2.32037d5)
!	parameter (bshift = 292.09d0)
! *** slope = 0.7 a=2.40
!	parameter (a0 = -0.40629d0)
!	parameter (b0 = 0.050088d0)
!	parameter (aslope = 16.9775d5)
!	parameter (bslope = 565.50d0)
!	parameter (ashift = 1.89213d5)
!	parameter (bshift = 263.93d0)
! *** slope = 0.9 a=2.30
!	parameter (a0 = -0.48332d0)
!	parameter (b0 = -0.12334d0)
!	parameter (aslope = 34.64125d5)
!	parameter (bslope = 1014.25d0)
!	parameter (ashift = 2.2815d5)
!	parameter (bshift = 294.26d0)
! *** slope = 0.9 a=2.32
!	parameter (a0 = -0.48137d0)
!	parameter (b0 = -0.11887d0)
!	parameter (aslope = 36.99125d5)
!	parameter (bslope = 1041.25d0)
!	parameter (ashift = 2.4744d5)
!	parameter (bshift = 306.23d0)
! *** slope = 0.3 a=2.85 
!	parameter (a0 = 0.0d0)
!	parameter (b0 = 0.0d0)
!	parameter (aslope = 3.0d5)
!	parameter (bslope = 0.0d0)
!	parameter (ashift = 8.0d5)
!	parameter (bshift = 1200.0d0)
      real(8)::favor(nmax),favor2(nmax)
      common /fdata/ favor,favor2




 
