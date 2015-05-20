#DEFINE _CHECK_DIPOLE_MOMMAP 0
#DEFINE _APPLY_CUTS 0

  MODULE ModDipoles_eeTTB
    use ModTopdecay
    implicit none


    double precision, private, parameter :: NCol=3d0
    double precision, private, parameter :: TR=1d0/2d0
    double precision, private, parameter :: CA=2d0*TR*NCol
    double precision, private, parameter :: CF=TR*(NCol**2-1d0)/NCol

    type :: Dipole
       double precision DipoleValue
       double precision MomTd(0:3,1:4)
       double precision Mass2Td(1:4)
       double precision alphaCut
    end type Dipole

    type(Dipole),private :: TheDipoles(1:4)

    double precision, private :: yRndDK(1:8),xFrac,r_sc
    integer, private :: NBin(1:20)
    
    logical, parameter :: invert_alphaCut = .false.

  contains



    
    SUBROUTINE Dipoles_eettb(nDipole,MomExt,MomExtTd,Dipole)! global norm:   4d0*Pi*alpha_s
      use ModParameters
      use ModKinematics
      use ModMisc
      implicit none
      integer :: nDipole,a,i,b,j,k
      real(8) :: MomExt(1:4,1:5),MomExtTd(1:4,1:4),Q(1:4),QTd(1:4),KSum(1:4)
      real(8) :: sab,sai,sbi,sij,sik,skj,x,v,y,yp,z,Q2,mu2,mu,MomFac1,MomFac2,MomFac3
      real(8) :: Dipole
      
      
      if(nDipole.eq.1) then
         i=4; j=3; k=5! final-final
         
         sij = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,j))
         sik = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,k))
         skj = 2d0*(MomExt(1:4,k)).dot.(MomExt(1:4,j))
         Q2 = 2d0*m_Top**2+ sij+sik+skj
         mu2 = m_Top**2/dabs(Q2)
         mu = dsqrt(mu2)
         y = sij/(sij+skj+sik)
         z = sik/(skj+sik)
         v = (1d0-y)*dsqrt((1d0-4d0*mu2)/( (2d0*mu2+(1d0-2d0*mu2)*(1d0-y))**2-4d0*mu2))
         yp = 1d0-2d0*mu*(1d0-mu)/(1d0-2d0*mu2)
         if( .not. invert_alphacut ) then
              if( y.gt.alpha_ff*yp) then
                  Dipole = (0d0,0d0)
                  return
              endif
         else
              if( y.lt.alpha_ff*yp) then
                  Dipole = (0d0,0d0)
                  return
              endif
         endif

         MomFac1 = dsqrt(((sij+sik+skj)**2-4d0*m_Top**4)/((sik+skj)**2-4d0*m_Top**2*(m_Top**2+sij)))
         MomFac2 = (0.5d0*(sik+skj)+m_Top**2)/Q2
         MomFac3 = MomFac2 + 0.5d0*sij/Q2
         
         Q(1:4) = MomExt(1:4,i)+MomExt(1:4,j)+MomExt(1:4,k)
         MomExtTd(1:4,k-1) = MomFac1*(MomExt(1:4,k)-MomFac2*Q(1:4)) + MomFac3*Q(1:4)
         MomExtTd(1:4,i-1) = Q(1:4) - MomExtTd(1:4,k-1)
         MomExtTd(1:4,1:2) = MomExt(1:4,1:2)
         
         Dipole = -1d0/sij * (2d0/(1d0-z+y*z)-v*(1d0+z+2d0*m_Top**2/sij))


      elseif(nDipole.eq.2) then
         i=5; j=3; k=4! final-final

         sij = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,j))
         sik = 2d0*(MomExt(1:4,i)).dot.(MomExt(1:4,k))
         skj = 2d0*(MomExt(1:4,k)).dot.(MomExt(1:4,j))
         Q2 = 2d0*m_Top**2 + sij+sik+skj
         mu2 = m_Top**2/dabs(Q2)
         mu = dsqrt(mu2)
         y = sij/(sij+skj+sik)
         z = sik/(skj+sik)
         v = (1d0-y)*dsqrt((1d0-4d0*mu2)/( (2d0*mu2+(1d0-2d0*mu2)*(1d0-y))**2-4d0*mu2))
         yp = 1d0-2d0*mu*(1d0-mu)/(1d0-2d0*mu2)
         if( .not. invert_alphacut ) then
              if( y.gt.alpha_ff*yp) then
                  Dipole = (0d0,0d0)
                  return
              endif
         else
              if( y.lt.alpha_ff*yp) then
                  Dipole = (0d0,0d0)
                  return
              endif
         endif
         MomFac1 = dsqrt(((sij+sik+skj)**2-4d0*m_Top**4)/((sik+skj)**2-4d0*m_Top**2*(m_Top**2+sij)))
         MomFac2 = (0.5d0*(sik+skj)+m_Top**2)/Q2
         MomFac3 = MomFac2 + 0.5d0*sij/Q2
         
         Q(1:4) = MomExt(1:4,i)+MomExt(1:4,j)+MomExt(1:4,k)
         MomExtTd(1:4,k-1) = MomFac1*(MomExt(1:4,k)-MomFac2*Q(1:4)) + MomFac3*Q(1:4)
         MomExtTd(1:4,i-1) = Q(1:4) - MomExtTd(1:4,k-1)
         MomExtTd(1:4,1:2) = MomExt(1:4,1:2)
         

         Dipole = -1d0/sij * (2d0/(1d0-z+y*z)-v*(1d0+z+2d0*m_Top**2/sij))

      endif

    END SUBROUTINE
    





  END MODULE ModDipoles_eeTTB
