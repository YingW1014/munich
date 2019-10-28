C------------------------------------------------------------------------
C     Copyright (C) 2001-2008, ENPC - INRIA - EDF R&D
C     
C     This file is part of the air quality modeling system Polyphemus.
C    
C     Polyphemus is developed in the INRIA - ENPC joint project-team
C     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
C    
C     Polyphemus is free software; you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published
C     by the Free Software Foundation; either version 2 of the License,
C     or (at your option) any later version.
C     
C     Polyphemus is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
C     General Public License for more details.
C     
C     For more information, visit the Polyphemus web site:
C     http://cerea.enpc.fr/polyphemus/
C
C------------------------------------------------------------------------

      SUBROUTINE LU_decompose_leighton (ns,M)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes an LU factorization of M.
C	  This routine is automatically generated by LU.cpp
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C
C	  nesp: matrix row number from the chemical species number
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     M: on entry, an invertible matrix. On exit, an LU factorization
C     of M.
C     
C     -- OUTPUT VARIABLES
C     
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C------------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C     
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C------------------------------------------------------------------------

      IMPLICIT NONE


      INTEGER ns
      DOUBLE PRECISION M(ns,ns)
      DOUBLE PRECISION temp


















































C     Upper part.
      M(49, 52) = M(49, 52) / M(49, 49)

C     Upper part.
      M(50, 51) = M(50, 51) / M(50, 50)
C     Upper part.
      temp = M(50, 49) * M(49, 52)
      M(50, 52) = ( M(50, 52) - temp ) / M(50, 50)

C     Lower part.
      temp = M(51, 50) * M(50, 51)
      M(51, 51) = M(51, 51) - temp
C     Lower part.
      temp = M(52, 50) * M(50, 51)
      M(52, 51) = M(52, 51) - temp
C     Upper part.
      temp = M(51, 50) * M(50, 52)
      M(51, 52) = ( M(51, 52) - temp ) / M(51, 51)

C     Lower part.
      temp = M(52, 50) * M(50, 52)
      temp = temp + M(52, 51) * M(51, 52)
      M(52, 52) = M(52, 52) - temp


      END
