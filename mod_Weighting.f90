  module ModWeighting
    implicit none
    
    integer,private,parameter :: NumMaxHisto=45
    
    contains

      subroutine ReadWeights(useHisto,filename)
! reads in the weights (=kfactors) from file filename, for histograms given in useHisto
! puts weights into global Weights         
        use ModParameters
        implicit none
      integer :: useHisto(1:NumMaxHisto)
      character :: filename*(80)
      character(len=*),parameter :: fmt1 = "(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)"
      real(8) :: BinVal=-1d-99,Val=-1d-99,Error=-1d-99
      integer :: NHisto=-999999,Hits=-999999
      integer :: BinCount,i,j
      character :: dummy*(1)



! if we use an average, use this
      Weights=-999d0
! if we use a geometric average, use this
!      kfactor=1d0      

      i=1
      j=0
      
! open the k-factor file
      do while (useHisto(i).ne.0)
         open(unit=1010,file=trim(filename),form='formatted',access='sequential')
         BinCount=0
         do while(.not.eof(1010))
            read(unit=1010,fmt="(A)") dummy
            if(dummy(1:1).eq."#") cycle
            backspace(unit=1010) ! go to the beginning of the line
            read(unit=1010,fmt=fmt1) NHisto,dummy,BinVal,dummy,Val,dummy,Error,dummy,Hits,dummy
                      
            if(NHisto.ne.useHisto(i)) cycle          
            BinCount=BinCount+1

            Weights(i,BinCount)=Val

         enddo
         close(1010)
       i=i+1
    enddo

  end subroutine ReadWeights

    subroutine ReWeight(ds_in,NBin,useHisto,ds_out)
! takes input ds_in and uses k-factors stored in Weights to reweight
! NBin contains info about bins that this ph space point lives in --
! this is used to pick the reweighting factors
! output ds_out
      use ModParameters
      implicit none 
      integer :: useHisto(1:NumMaxHisto),NBin(1:NumMaxHisto)
      complex(8) :: ds_in,ds_out
      integer    :: i,div
      real(8)    :: TotalWeight,CurrentWeight

      div=0
      TotalWeight=0d0
      do i=1,size(Weights,1)
         CurrentWeight=Weights(i,NBin(useHisto(i)))              
         

         if ( abs(CurrentWeight+999d0) .lt. 1d-10) then
           div=div-1 
         else
! this is NOT the ideal way to combine: better use to (NLO1+NLO2+...)/(LO1+LO2+...)
            TotalWeight = TotalWeight + CurrentWeight
         endif
      enddo
      div=div+size(Weights,1)
      if (div .ne. 0) then
         TotalWeight=TotalWeight/div
      else
         TotalWeight=1d0
      endif

      ds_out=ds_in * TotalWeight


    end subroutine ReWeight


  end module ModWeighting
