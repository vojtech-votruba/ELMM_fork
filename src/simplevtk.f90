module vtkarray
 !Simple module to output arrays for visualization. No physical coordinates are used, only the position in the array.
 !Mostly only for debugging.

  use iso_fortran_env, only: real32, real64
  use Endianness, only: BigEnd

  implicit none

  interface VtkArrayAscii
    module procedure SVtkArrayAscii
    module procedure DVtkArrayAscii
  end interface

  interface VtkArrayBin
    module procedure SVtkArrayBin
    module procedure DVtkArrayBin
  end interface
  
  character, parameter :: lf = achar(10)

contains

  subroutine SVtkArrayAscii(fname, A)
    character(len=*), intent(in)  :: fname
    real(real32), intent(in)      :: A(:,:,:)
    integer                       :: nx,ny,nz
    integer                       :: i
    integer                       :: unit
    character(len=40)             :: str

    nx = ubound(A,1)
    ny = ubound(A,2)
    nz = ubound(A,3)

    open(newunit=unit,file=fname,status="replace",action="write")
    write(unit,"(A)") "# vtk DataFile Version 2.0"
    write(unit,"(A)") "CLMM output file"
    write(unit,"(A)") "ASCII"
    write(unit,"(A)") "DATASET RECTILINEAR_GRID"
    str="DIMENSIONS"
    write(str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
    write(unit,"(A)") trim(str)
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"float"
    write(unit,"(A)") str
    write(unit,*) (i, i=1,nx)
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"float"
    write(unit,"(A)") trim(str)
    write(unit,*) (i, i=1,ny)
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"float"
    write(unit,"(A)") trim(str)
    write(unit,*) (i, i=1,nz)
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(unit,"(A)") trim(str)


    write(unit,"(A)") "SCALARS array float"
    write(unit,"(A)") "LOOKUP_TABLE default"

    write(unit,'(*(g0,/))') A(1:nx,1:ny,1:nz)

    write(unit,*)
    close(unit)

  end subroutine SVtkArrayAscii

  subroutine DVtkArrayAscii(fname, A)
    character(len=*), intent(in)  :: fname
    real(real64), intent(in)      :: A(:,:,:)
    integer                       :: nx,ny,nz
    integer                       :: i
    integer                       :: unit
    character(len=40)             :: str

    nx = ubound(A,1)
    ny = ubound(A,2)
    nz = ubound(A,3)

    open(newunit=unit,file=fname,status="replace",action="write")
    write(unit,"(A)") "# vtk DataFile Version 2.0"
    write(unit,"(A)") "CLMM output file"
    write(unit,"(A)") "ASCII"
    write(unit,"(A)") "DATASET RECTILINEAR_GRID"
    str="DIMENSIONS"
    write(str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
    write(unit,"(A)") trim(str)
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"double"
    write(unit,"(A)") str
    write(unit,*) (i, i=1,nx)
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"double"
    write(unit,"(A)") trim(str)
    write(unit,*) (i, i=1,ny)
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"double"
    write(unit,"(A)") trim(str)
    write(unit,*) (i, i=1,nz)
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(unit,"(A)") trim(str)


    write(unit,"(A)") "SCALARS array double"
    write(unit,"(A)") "LOOKUP_TABLE default"

    write(unit,'(*(g0,/))') A(1:nx,1:ny,1:nz)

    write(unit,*)
    close(unit)

  end subroutine DVtkArrayAscii
  
  subroutine SVtkArrayBin(fname, A, offsets)
    character(len=*), intent(in)  :: fname
    real(real32), intent(in)      :: A(:,:,:)
    integer, optional, intent(in) :: offsets(3)
    integer                       :: nx,ny,nz
    integer                       :: i
    integer                       :: unit
    integer                       :: offs(3)
    character(len=70)             :: str

    nx = ubound(A,1)
    ny = ubound(A,2)
    nz = ubound(A,3)
    
    if (present(offsets)) then
      offs = offsets
    else
      offs = 0
    end if

    open(newunit=unit,file=fname,access='stream',form='unformatted',status="replace",action="write")
    write(unit) "# vtk DataFile Version 2.0", lf
    write(unit) "CLMM output file", lf
    write(unit) "BINARY", lf
    write(unit) "DATASET RECTILINEAR_GRID", lf
    str="DIMENSIONS"
    write(str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
    write(unit) str, lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(1), i=1,nx)], real32)), lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(2), i=1,ny)], real32)), lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(3), i=1,nz)], real32)), lf
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(unit) str, lf


    write(unit) "SCALARS array float", lf
    write(unit) "LOOKUP_TABLE default", lf

    write(unit) Bigend(A(1:nx,1:ny,1:nz)), lf

    close(unit)

  end subroutine SVtkArrayBin

  subroutine DVtkArrayBin(fname, A, offsets)
    character(len=*), intent(in)  :: fname
    real(real64), intent(in)      :: A(:,:,:)
    integer, optional, intent(in) :: offsets(3)
    integer                       :: nx,ny,nz
    integer                       :: i
    integer                       :: unit
    integer                       :: offs(3)
    character(len=70)             :: str

    nx = ubound(A,1)
    ny = ubound(A,2)
    nz = ubound(A,3)
    
    if (present(offsets)) then
      offs = offsets
    else
      offs = 0
    end if

    open(newunit=unit,file=fname,access='stream',form='unformatted',status="replace",action="write")
    write(unit) "# vtk DataFile Version 2.0", lf
    write(unit) "CLMM output file", lf
    write(unit) "BINARY", lf
    write(unit) "DATASET RECTILINEAR_GRID", lf
    str="DIMENSIONS"
    write(str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
    write(unit) str, lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(1), i=1,nx)], real32)), lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(2), i=1,ny)], real32)), lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(3), i=1,nz)], real32)), lf
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(unit) str, lf


    write(unit) "SCALARS array float", lf
    write(unit) "LOOKUP_TABLE default", lf

    write(unit) BigEnd(real(A(1:nx,1:ny,1:nz), real32)), lf

    close(unit)

  end subroutine DVtkArrayBin
  
end module vtkarray


