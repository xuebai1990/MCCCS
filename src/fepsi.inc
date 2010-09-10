c FEPSI.INC
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff,
c John Stubbs, and Collin Wick and Ilja Siepmann
c
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License
c as published by the Free Software Foundation; either version 2
c of the License, or (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA  02111-1307, USA.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      logical::lfepsi
      real(8)::aslope,bslope,ashift,bshift,a0,b0

c ----- If lfepsi is true, the fluctuation of epsilon is used instead of
c ----- the fluctuation of the sigma.
 
        parameter (lfepsi = .false.)
c *** slope=0.3 a=3.05
c	parameter (a0 = 0.2003d0)
c	parameter (b0 = 1.3946d0)
c	parameter (aslope = 8.85d5)
c	parameter (bslope = 158.25d0)
c	parameter (ashift = 7.227d5)
c	parameter (bshift = 505.97d0)
c *** slope=0.3 a=2.90
c	parameter (a0 = 0.15561d0)
c	parameter (b0 = 1.2960d0)
c	parameter (aslope = 5.5475d5)
c	parameter (bslope = 131.23d0)
c	parameter (ashift = 4.121d5)
c	parameter (bshift = 381.96d0)
c *** slope=0.3 a=2.95
c	parameter (a0 = 0.16833d0)
c	parameter (b0 = 1.3299d0)
c	parameter (aslope = 6.5125d5)
c	parameter (bslope = 139.75d0)
c	parameter (ashift = 4.9975d5)
c	parameter (bshift = 419.81d0)
c *** slope=0.3 a=2.82
c	parameter (a0 = 0.13150d0)
c	parameter (b0 = 1.2574d0)
c	parameter (aslope = 4.2863d5)
c	parameter (bslope = 117.5d0)
c	parameter (ashift = 3.0208d5)
c	parameter (bshift = 323.6d0)
c *** slope=0.3 a=2.75
c	parameter (a0 = 0.11037d0)
c	parameter (b0 = 1.1942d0)
c	parameter (aslope = 3.4013d5)
c	parameter (bslope = 108d0)
c	parameter (ashift = 2.2868d5)
c	parameter (bshift = 285.04d0)
c *** slope=0.3 a=2.85
c	parameter (a0 = 0.14035d0)
c	parameter (b0 = 1.2596d0)
c	parameter (aslope = 4.7263d5)
c	parameter (bslope = 123.25d0)
c	parameter (ashift = 3.3978d5)
c	parameter (bshift = 347.63d0)
c *** slope=0.3 a=2.88
c	parameter (a0 = 0.14820d0)
c	parameter (b0 = 1.2689d0)
c	parameter (aslope = 5.2125d5)
c	parameter (bslope = 128.75d0)
c	parameter (ashift = 3.8220d5)
c	parameter (bshift = 371.29d0)
c *** slope=0.5 a=2.71
      parameter (a0 = -0.22818d0)
      parameter (b0 = 0.44662d0)
      parameter (aslope = 14.0738d5)
      parameter (bslope = 351.25d0)
      parameter (ashift = 3.5852d5)
      parameter (bshift = 358.94d0)
c *** slope=0.5 a=2.68
c	parameter (a0 = -0.23379d0)
c	parameter (b0 = 0.43630d0)
c	parameter (aslope = 12.78125d5)
c	parameter (bslope = 337.5d0)
c	parameter (ashift = 3.1904d5)
c	parameter (bshift = 337.86d0)
c *** slope=0.5 a=2.655
c	parameter (a0 = -0.23808d0)
c	parameter (b0 = 0.425507d0)
c	parameter (aslope = 11.7775d5)
c	parameter (bslope = 322.04d0)
c	parameter (ashift = 2.8902d5)
c	parameter (bshift = 326.875d0)
c *** slope=0.5 a=2.665
c	parameter (a0 = -0.23641d0)
c	parameter (b0 = 0.43409d0)
c	parameter (aslope = 12.171d5)
c	parameter (bslope = 330d0)
c	parameter (ashift = 3.0075d5)
c	parameter (bshift = 326.52d0)
c *** slope=0.5 a=2.65
c	parameter (a0 = -0.23903d0)
c	parameter (b0 = 0.42231d0)
c	parameter (aslope = 11.5875d5)
c	parameter (bslope = 325d0)
c	parameter (ashift = 2.8339d5)
c	parameter (bshift = 319.44d0)
c *** slope=0.5 a=2.60
c	parameter (a0 = -0.24796d0)
c	parameter (b0 = 0.40345d0)
c	parameter (aslope = 9.8238d5)
c	parameter (bslope = 304d0)
c	parameter (ashift = 2.3206d5)
c	parameter (bshift = 288.72d0)
c *** slope=0.5 a=2.55
c	parameter (a0 = -0.25703d0)
c	parameter (b0 = 0.38299d0)
c	parameter (aslope = 8.3075d5)
c	parameter (bslope = 284.38d0)
c	parameter (ashift = 1.8945d5)
c	parameter (bshift = 261.09d0)
c *** slope = 0.7 a=2.52
c	parameter (a0 = -0.39106d0)
c	parameter (b0 = 0.08431d0)
c	parameter (aslope = 25.22875d5)
c	parameter (bslope = 662.75d0)
c	parameter (ashift = 3.0689d5)
c	parameter (bshift = 335.43d0)
c *** slope = 0.7 a=2.51
c	parameter (a0 = -0.39233d0)
c	parameter (b0 = 0.08203d0)
c	parameter (aslope = 24.425d5)
c	parameter (bslope = 653.75d0)
c	parameter (ashift = 2.9499d5)
c	parameter (bshift = 328.60d0)
c *** slope = 0.7 a=2.50
c	parameter (a0 = -0.39357d0)
c	parameter (b0 = 0.07946d0)
c	parameter (aslope = 23.641d5)
c	parameter (bslope = 645d0)
c	parameter (ashift = 2.835d5)
c	parameter (bshift = 322.13d0)
c *** slope = 0.7 a=2.49
c	parameter (a0 = -0.394906d0)
c	parameter (b0 = 0.075678d0)
c	parameter (aslope = 22.88625d5)
c	parameter (bslope = 637.735d0)
c	parameter (ashift = 2.72469d5)
c	parameter (bshift = 316.22d0)
c *** slope = 0.7 a=2.48
c	parameter (a0 = -0.396175d0)
c	parameter (b0 = 0.072973d0)
c	parameter (aslope = 22.149125d5)
c	parameter (bslope = 629d0)
c	parameter (ashift = 2.61798d5)
c	parameter (bshift = 309.95d0)
c *** slope = 0.7 a=2.47
c	parameter (a0 = -0.39745d0)
c	parameter (b0 = 0.070157d0)
c	parameter (aslope = 21.4335d5)
c	parameter (bslope = 620.75d0)
c	parameter (ashift = 2.5151d5)
c	parameter (bshift = 306.89d0)
c *** slope = 0.7 a=2.45
c	parameter (a0 = -0.399983d0)
c	parameter (b0 = 0.064159d0)
c	parameter (aslope = 20.06425d5)
c	parameter (bslope = 604.75d0)
c	parameter (ashift = 2.32037d5)
c	parameter (bshift = 292.09d0)
c *** slope = 0.7 a=2.40
c	parameter (a0 = -0.40629d0)
c	parameter (b0 = 0.050088d0)
c	parameter (aslope = 16.9775d5)
c	parameter (bslope = 565.50d0)
c	parameter (ashift = 1.89213d5)
c	parameter (bshift = 263.93d0)
c *** slope = 0.9 a=2.30
c	parameter (a0 = -0.48332d0)
c	parameter (b0 = -0.12334d0)
c	parameter (aslope = 34.64125d5)
c	parameter (bslope = 1014.25d0)
c	parameter (ashift = 2.2815d5)
c	parameter (bshift = 294.26d0)
c *** slope = 0.9 a=2.32
c	parameter (a0 = -0.48137d0)
c	parameter (b0 = -0.11887d0)
c	parameter (aslope = 36.99125d5)
c	parameter (bslope = 1041.25d0)
c	parameter (ashift = 2.4744d5)
c	parameter (bshift = 306.23d0)
c *** slope = 0.3 a=2.85 
c	parameter (a0 = 0.0d0)
c	parameter (b0 = 0.0d0)
c	parameter (aslope = 3.0d5)
c	parameter (bslope = 0.0d0)
c	parameter (ashift = 8.0d5)
c	parameter (bshift = 1200.0d0)
      real(8)::favor(nmax),favor2(nmax)
      common /fdata/ favor,favor2




 
