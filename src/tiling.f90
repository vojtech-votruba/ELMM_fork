module Tiling
!computation of tile sizes
  use Kinds, only: int8,knd
 !$ use OMP_LIB

  implicit none

  integer :: tilesize = 2**16 !size of loop tile in bytes
                              !should be close to the size of the L1 cache

  integer,dimension(8),protected :: tilenx,tileny,tilenz
  !index is the number of real arrays needed in the tiled loop

  contains

    subroutine InitTiles(Prnx,Prny,Prnz)
#ifdef PAR
      use custom_par
#endif
      integer,intent(in) :: Prnx,Prny,Prnz
      integer :: omp_threads = 1
      integer :: narrays
      integer :: bnx,bny,bnz,bn
      integer :: bnzl,bnzu,bnzo
      real    :: rbn
      integer :: i

      !$omp parallel default(private) shared(tilenx,tileny,tilenz,tilesize,Prnx,Prny,Prnz,omp_threads)
      !$omp single
      !$ omp_threads = omp_get_num_threads()
      !$omp end single

      !$omp do
      do narrays=1,size(tilenx)

        bn = int( tilesize/(storage_size(1._knd)/storage_size('a')) )

        bn = bn/narrays !four arrays U,V,W,Visc

        rbn = bn**(1./3.)/(Prnx*Prny*Prnz)**(1./3.) !geometric mean to scale the tiling box

        bnx = max(nint(Prnx*rbn),1)  !first try for tiling block size
        bny = max(nint(Prny*rbn),1)
        bnz = max(nint(Prnz*rbn),1)

        if (mod(Prnz/bnz,omp_threads)/=0) then
          bnzo = bnz !original

          bnzu = bnz !upper try
          i = 0
          do
            i = i + 1

            if (mod(Prnz/(bnz+i),omp_threads)==0) then
              bnzu = bnz + i
              exit
            end if

            if (i>=bnz) exit
          end do

          bnzl = bnz !lower try
          if (bnzl>1) then
            i = 0
            do
              i = i - 1
              if (mod(Prnz/(bnz+i),omp_threads)==0) then
                bnzl = bnz + i
                exit
              end if

              if (i<=-bnz/2) exit
            end do
          end if

          if (bnzu/=bnz .and. bnzl/=bnz) then  !choose the closer one
            if (abs(bnzu-bnz)<=abs(bnzl-bnz)) then
              bnz = bnzu
            else
              bnz = bnzl
            end if
          else if (bnzu/=bnz) then
            bnz = bnzu
          else if (bnzl/=bnz) then
            bnz = bnzl
          end if

          rbn = sqrt(real(bn)/bnz)/sqrt(real(Prnx*Prny)) !geometric mean to scale the tiling box

          bnx = nint(Prnx*rbn)  !first try for tiling block size
          bny = nint(Prny*rbn)

        end if


        tilenx(narrays) = bnx
        tileny(narrays) = bny
        tilenz(narrays) = bnz

      end do
      !$omp end do

      !$omp workshare
      where (tilenx<1) tilenx=1
      where (tileny<1) tileny=1
      where (tilenz<1) tilenz=1
      !$omp end workshare

      !$omp end parallel

    end subroutine InitTiles

end module Tiling
