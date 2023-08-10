module Kinds
  use iso_fortran_env
#ifdef DPREC  
  integer,parameter :: rp = real64
#else
  integer,parameter :: rp = real32
#endif
end module

module Endianness

  use iso_fortran_env
  
  implicit none

  private
  public :: GetEndianness, BigEnd, SwapB, littleendian

  logical, save :: littleendian = .true.

  interface BigEnd
    module procedure BigEnd16
    module procedure BigEnd32
    module procedure BigEnd32_a1
    module procedure BigEnd32_a2
    module procedure BigEnd32_a3
    module procedure BigEnd32_a4
    module procedure BigEnd64
    module procedure BigEnd64_a1
    module procedure BigEnd64_a2
    module procedure BigEnd64_a3
    module procedure BigEnd64_a4
  end interface

  interface SwapB
    module procedure SwapB32
    module procedure SwapB64
  end interface

  contains

    subroutine GetEndianness
      character(4) :: bytes !may not work on some processors

      bytes = transfer(1_int32,bytes)
      if (ichar(bytes(4:4))==1) then
        littleendian=.false.
      else
        littleendian=.true.
      endif
    end subroutine GetEndianness

   
    elemental function BigEnd16(x) result(res)
      integer(int16) :: res
      integer(int16),intent(in)::x
      character(2) :: bytes
      
      if (.not.littleendian) then
        res = x
      else
        bytes = transfer(x,bytes)
        res = ichar(bytes(2:2),int16)
        res = ior( ishft(ichar(bytes(1:1),int16),8), res )
      endif
    end function
    
    function BigEnd32(x) result(res)
      real(real32),intent(in) :: x
      real(real32) :: res
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a1(x) result(res)
      real(real32),intent(in) :: x(:)
      real(real32) :: res(size(x,1))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a2(x) result(res)
      real(real32),intent(in) :: x(:,:)
      real(real32) :: res(size(x,1),size(x,2))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a3(x) result(res)
      real(real32),intent(in) :: x(:,:,:)
      real(real32) :: res(size(x,1),size(x,2),size(x,3))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a4(x) result(res)
      real(real32),intent(in) :: x(:,:,:,:)
      real(real32) :: res(size(x,1),size(x,2),size(x,3),size(x,4))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64(x) result(res)
      real(real64),intent(in) :: x
      real(real64) :: res
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a1(x) result(res)
      real(real64),intent(in) :: x(:)
      real(real64) :: res(size(x,1))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a2(x) result(res)
      real(real64),intent(in) :: x(:,:)
      real(real64) :: res(size(x,1),size(x,2))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a3(x) result(res)
      real(real64),intent(in) :: x(:,:,:)
      real(real64) :: res(size(x,1),size(x,2),size(x,3))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a4(x) result(res)
      real(real64),intent(in) :: x(:,:,:,:)
      real(real64) :: res(size(x,1),size(x,2),size(x,3),size(x,4))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    elemental function SwapB32(x) result(res)
      real(real32) :: res
      real(real32),intent(in) :: x
      character(4) :: bytes
      integer(int32) :: t
      real(real32) :: rbytes, rt
      !nicer looking TRANSFER is problematic to optimize by ifort
      ! gfortran would be OK with that
      equivalence (rbytes, bytes)
      equivalence (t, rt)
      
      rbytes = x

      t = ichar(bytes(4:4),int32)

      t = ior( ishftc(ichar(bytes(3:3),int32),8),  t )

      t = ior( ishftc(ichar(bytes(2:2),int32),16), t )

      t = ior( ishftc(ichar(bytes(1:1),int32),24), t )

      res = rt
        
    end function
    
    elemental function SwapB64(x) result(res)
      real(real64) :: res
      real(real64),intent(in) :: x
      character(8) :: bytes
      integer(int64) :: t
      real(real64) :: rbytes, rt
      equivalence (rbytes, bytes)
      equivalence (t, rt)
      
      rbytes = x

      t = ichar(bytes(8:8),int64)

      t = ior( ishftc(ichar(bytes(7:7),int64),8),  t )

      t = ior( ishftc(ichar(bytes(6:6),int64),16), t )

      t = ior( ishftc(ichar(bytes(5:5),int64),24), t )

      t = ior( ishftc(ichar(bytes(4:4),int64),32), t )

      t = ior( ishftc(ichar(bytes(3:3),int64),40), t )

      t = ior( ishftc(ichar(bytes(2:2),int64),48), t )

      t = ior( ishftc(ichar(bytes(1:1),int64),56), t )

      res = rt

    end function

end module Endianness




module Types
  use Kinds
  
  implicit none

  type grid
    integer :: nx,ny,nz
    integer :: offx,offy,offz
    integer :: unit
    real(rp),allocatable :: x(:),y(:),z(:)
    character(100) :: fname
  contains
    procedure :: read_header
    procedure :: read_title
    procedure :: read_scalar
    procedure :: read_vector
    procedure :: skip_var
  end type
  
  integer, parameter :: SCALAR = 7, VECTOR = 8
  character,parameter :: lf = achar(10)

contains
  

  subroutine read_header(g)
    use Endianness
    class(grid),intent(inout) :: g
    integer :: io
    real(real32), allocatable :: tmp(:)

    open(newunit=g%unit,file=g%fname,access="stream",form="unformatted",status="old",action="read",iostat=io)

    if (io/=0) then
      g%nx = 0
      g%ny = 0
      g%nz = 0
      return
    end if

    call skip_to("X_COORDINATES", io)
    if (io/=0) then
      write(*,*) "Error, cannot find X_COORDINATES in '",g%fname,"'"
      stop
    end if
    call get_number(g%nx)
    call skip_line
    
    allocate(tmp(g%nx))
    read(g%unit,iostat=io) tmp
    g%x=BigEnd(tmp)
    deallocate(tmp)

    call skip_to("Y_COORDINATES", io)
    if (io/=0) then
      write(*,*) "Error, cannot find Y_COORDINATES in '",g%fname,"'"
      stop
    end if
    call get_number(g%ny)
    call skip_line
    
    allocate(tmp(g%ny))
    read(g%unit,iostat=io) tmp
    g%y=BigEnd(tmp)
    deallocate(tmp)

    call skip_to("Z_COORDINATES", io)
    if (io/=0) then
      write(*,*) "Error, cannot find Z_COORDINATES in '",g%fname,"'"
      stop
    end if
    call get_number(g%nz)
    call skip_line
    
    allocate(tmp(g%nz))
    read(g%unit,iostat=io) tmp
    g%z=BigEnd(tmp)
    deallocate(tmp)
    
    call skip_line
    call skip_line
    
  contains
    subroutine skip_line
      character :: ch
      do
        read(g%unit) ch
        if (ch==new_line("a")) return
      end do
    end subroutine

    subroutine get_number(n)
      integer, intent(out) :: n
      character :: ch
      integer :: io
      character(:), allocatable :: num_str
      
      do
        read(g%unit, iostat=io) ch
        if (io/=0) then
          write(*,*) "Error reading from '",g%fname,"'."
          stop 3
        end if
        if (ch/=' ') exit
      end do
      
      num_str=ch
      do
        read(g%unit) ch
        if (ch<'0' .or. ch>'9') exit
        num_str = num_str // ch
      end do

      read(num_str,*) n
    end subroutine

    subroutine skip_to(str, stat)
      character(*), intent(in) :: str
      integer, intent(out) :: stat
      character :: ch
      integer :: io

      do
        read(g%unit, iostat=io) ch

        if (io/=0) then
          stat = 1
          return
        end if

        if (ch==str(1:1)) then
          call check(str(2:), stat)
          if (stat == 0) return
        end if

      end do
    end subroutine
    
    subroutine check(str, stat)
      character(*), intent(in) :: str
      integer, intent(out) :: stat
      character :: ch
      integer :: i, io

      stat = 1
      i = 0

      do
        i = i + 1

        read(g%unit, iostat=io) ch

        if (io/=0) return

        if (ch/=str(i:i)) return

        if (i==len(str)) then
          stat = 0
          return
        end if
      end do
    end subroutine
   
  end subroutine
  
  subroutine read_title(g,status,title)
    use iso_fortran_env
    class(grid),intent(in) :: g
    integer :: status
    character(:),allocatable,intent(out) :: title
    character(7) :: vtype
    character :: ch
    integer :: io, n, lf_pos
    character(256) :: msg
  
    do
      vtype=""
      read(g%unit,iostat=io,iomsg=msg) vtype
      
      if (io==iostat_end) then
        status = 0
        return
      else if (io/=0) then
        write(*,'(*(g0))') io, trim(msg), "  ", g%fname, " ", "'",vtype,"'"
        stop "Error reading title."
      end if
      
      if (vtype=='VECTORS') then
      
        status = VECTOR
        title = vtype
        do
          read(g%unit) ch
          title = title // ch
          if (ch==lf) exit
        end do
        
        exit
        
      else if (vtype=='SCALARS') then
      
        status = SCALAR
        title = vtype
        n = 0
        do
          read(g%unit) ch
          title = title // ch
          if (ch==lf) n = n + 1
          if (n==2) exit
        end do
        
        exit
        
      else
        !if our vtype buffer contains the end-of-line LF character, move back just behind it
        !otherwise skip to the nearest one
        lf_pos = scan(vtype, achar(10),.true.)
        if (lf_pos>0) then
          call stream_back(len(vtype)-lf_pos+1)
        else
          call skip_line
        end if
        
      end if
    end do
    
  contains
  
    subroutine skip_line
      character :: ch
      do
        read(g%unit) ch
        if (ch==new_line("a")) return
      end do
    end subroutine
    
    subroutine stream_back(n)
      integer, intent(in) :: n
      integer :: pos
      character :: dummy
      inquire(unit=g%unit,pos=pos)
      read(g%unit,pos=pos-n) dummy
    end subroutine

  end subroutine
  
  subroutine read_scalar(g,buf)
    use Endianness
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,:,:)
    real(real32), allocatable :: tmp(:,:,:) 
    character :: ch
    
    if (rp==real32) then
      read(g%unit) buf(1:g%nx,1:g%ny,1:g%nz)
      buf(1:g%nx,1:g%ny,1:g%nz) = BigEnd(buf(1:g%nx,1:g%ny,1:g%nz))
    else
      allocate(tmp(1:g%nx,1:g%ny,1:g%nz))
      read(g%unit) tmp
      buf(1:g%nx,1:g%ny,1:g%nz) = BigEnd(tmp)
    end if
    read(g%unit) ch
  end subroutine
  
  subroutine read_vector(g,buf)
    use Endianness
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,:,:,:)
    real(real32), allocatable :: tmp(:,:,:,:) 
    character :: ch

   if (rp==real32) then
      read(g%unit) buf(:,1:g%nx,1:g%ny,1:g%nz)
      buf(:,1:g%nx,1:g%ny,1:g%nz) = BigEnd(buf(:,1:g%nx,1:g%ny,1:g%nz))
    else
      allocate(tmp(3,1:g%nx,1:g%ny,1:g%nz))
      read(g%unit) tmp
      buf(:,1:g%nx,1:g%ny,1:g%nz) = BigEnd(tmp)
    end if
    read(g%unit) ch
  end subroutine
  
  subroutine skip_var(g,stat)
    use iso_c_binding
    class(grid),intent(in) :: g
    integer, intent(in) :: stat
    integer(c_size_t) :: nbytes, fpos
    nbytes = c_sizeof(1_real32) * &
            int(g%nx, c_size_t) * &
            int(g%ny, c_size_t) * &
            int(g%nz, c_size_t)
    if (stat==VECTOR) nbytes = nbytes *3
    nbytes = nbytes + 1 !EOL
    
    inquire(g%unit, pos = fpos)
    read(g%unit, pos=fpos + nbytes)
  end subroutine
  
end module Types






program spectrum_1D
  use fftw3
  use Kinds
  use Types
  use Endianness
  
  implicit none
  
  integer :: nx, ny, nz
  real(rp), allocatable :: sc(:,:,:),vec(:,:,:,:)

  real(rp), allocatable, dimension(:) :: kappas, spectrum_1comp
  real(rp), allocatable, dimension(:,:) :: spectrum
  real(rp) :: variance, variances(3)
 
  character(:), allocatable :: title
  integer :: status
  
  character(:), allocatable :: in_name, out_name, field, arg, field_title
  integer :: iarg, arg_len
  
  integer :: n, iu, j
  
  type(grid) :: in_g
  
  integer :: dir = 1
  
  integer :: xrange(2)=[-1,huge(1)], &
             yrange(2)=[-1,huge(1)], &
             zrange(2)=[-1,huge(1)]
  
  call GetEndianness

  if (command_argument_count()<2) then
    write(*,*) 'Example usage: 1D_spectrum U.vtk spectrum.txt [dir=1 field=U xrange=1,10 yrange=1,10 zrange=1,10]'
    stop 1
  end if
  
  
  call get_command_argument(1, length=arg_len)
  allocate(character(arg_len) :: in_name)
  call get_command_argument(1, value=in_name)
 
  call get_command_argument(2, length=arg_len)
  allocate(character(arg_len) :: out_name)
  call get_command_argument(2, value=out_name)
 
  do iarg = 3, command_argument_count()
    call get_command_argument(iarg, length=arg_len)
    allocate(character(arg_len) :: arg)
    call get_command_argument(iarg, value=arg)
    
    !NOTE:quick and dirty, but working
    if (arg(1:4)=="dir=") then
      read(arg(5:),*) dir
    else if (arg(1:7)=="xrange=") then
      read(arg(8:),*) xrange
    else if (arg(1:7)=="yrange=") then
      read(arg(8:),*) yrange
    else if (arg(1:7)=="zrange=") then
      read(arg(8:),*) zrange
    else if (arg(1:6)=="field=") then
      field = arg(7:)
    else
      write(*,*) "Error, unrecognized command line argument:","'",arg,"'"
      stop 1
    end if
    deallocate(arg)
  end do
  
  if (.not.allocated(field)) allocate(character(0) :: field)
    
  
  in_g%fname = in_name
  call in_g%read_header
  
  
  nx = in_g%nx
  ny = in_g%ny
  nz = in_g%nz
  write(*,*) in_g%nx, in_g%ny, in_g%nz
  
  allocate(sc(1:nx,1:ny,1:nz))
  allocate(vec(1:3,1:nx,1:ny,1:nz))
  
  do
    call in_g%read_title(status,title)
    
    if (status/=SCALAR .and. status/=VECTOR) exit

    field_title = get_field_name(title)

    if (len(field)<=0.or.field_title==field) then
    
      if (status==SCALAR) then
        call get_buffer(status,sc,vec)

        call do_transform(sc, dir, kappas, spectrum_1comp, variance)
        
        allocate(spectrum(size(spectrum_1comp),1))
        spectrum(:,1) = spectrum_1comp
        
        n = size(kappas)
        
        open(newunit=iu,file=out_name)
        do j = 1, n
          write(iu,*) kappas(j), spectrum(j,1)
        end do 
        close(iu)
        
        write(*,*) "average 1D variance:", variance
        write(*,*) "spectrum integral:", sum(spectrum)/n *kappas(n)
        write(*,*) "full 3D array variance:", array_variance(sc)

      else if (status==VECTOR) then
        call get_buffer(status,sc,vec)
        
        sc = vec(1,:,:,:)
        call do_transform(sc, dir, kappas, spectrum_1comp, variances(1))        
        
        allocate(spectrum(size(spectrum_1comp),3))
        
        spectrum(:,1) = spectrum_1comp
        sc = vec(2,:,:,:)
        call do_transform(sc, dir, kappas, spectrum_1comp, variances(2))        
        spectrum(:,2) = spectrum_1comp
        sc = vec(3,:,:,:)
        call do_transform(sc, dir, kappas, spectrum_1comp, variances(3))        
        spectrum(:,3) = spectrum_1comp
        
        
        n = size(kappas)
        
        open(newunit=iu,file=out_name)
        do j = 1, n
          write(iu,*) kappas(j), spectrum(j,:), sum(spectrum(j,:))/2
        end do 
        close(iu)
        
        write(*,*) "average 1D variances:", variances
        write(*,*) "spectrum integrals:", sum(spectrum,dim=1)/n *kappas(n)
        write(*,*) "<E11 + E22 + E33>/2:", sum(variances)/2
        forall(j=1:3) variances(j) = array_variance(vec(j,:,:,:))
        write(*,*) "full 3D variances and TKE:", variances, sum(variances)/2
      end if
      exit
    else
      call in_g%skip_var(status)
    end if
  end do
  
contains

  subroutine do_transform(u3d, dir, kappas, spectrum, u_var)
    real(rp), intent(in) :: u3d(:,:,:)
    integer, intent(in) :: dir    
    real(rp), allocatable, intent(out) :: kappas(:), spectrum(:)
    real(rp), intent(out) :: u_var
    real(rp), allocatable :: u(:), sp(:), sp_avg(:)
    complex(rp), allocatable :: u_hat(:)
    integer :: n, nj, nk
    integer :: mini, minj, mink, maxi, maxj, maxk
    type(c_ptr), save :: forw = c_null_ptr
    integer :: i, j, k
    real(rp) :: dx, lx, u_m
    real(rp), parameter :: pi = acos(-1.0_rp)
    
    n = size(u3d,dir)
    
    select case (dir)
      case (1)
        nj = size(u3d, 2)
        nk = size(u3d, 3)

        mini = max(1,xrange(1))
        maxi = min(n,xrange(2))

        minj = max(1,yrange(1))
        maxj = min(nj,yrange(2))

        mink = max(1,zrange(1))
        maxk = min(nk,zrange(2))

      case (2)
        nj = size(u3d, 1)
        nk = size(u3d, 3)

        mini = max(1,yrange(1))
        maxi = min(n,yrange(2))

        minj = max(1,xrange(1))
        maxj = min(nj,xrange(2))

        mink = max(1,zrange(1))
        maxk = min(nk,zrange(2))

      case (3)
        nj = size(u3d, 1)
        nk = size(u3d, 2)

        mini = max(1,zrange(1))
        maxi = min(n,zrange(2))

        minj = max(1,xrange(1))
        maxj = min(nj,xrange(2))

        mink = max(1,yrange(1))
        maxk = min(nk,yrange(2))

    end select
    
    n = maxi - mini + 1
    nj = maxj - minj + 1
    nk = maxk - mink + 1
    
    allocate(u(n))
    allocate(u_hat(0:n/2))
    allocate(sp(n/2))
    allocate(sp_avg(n/2))
    
    sp_avg = 0
    u_var = 0
    
    if (.not. c_associated(forw)) &
      forw = fftw_plan_gen(size(u), &
                  u, u_hat, FFTW_UNALIGNED+FFTW_ESTIMATE)
                
    do k = mink, maxk
      do j = minj, maxj
        select case (dir)
          case (1)
            u = u3d(mini:maxi,j,k)
          case (2)
            u = u3d(j,mini:maxi,k)
          case (3)
            u = u3d(j,k,mini:maxi)
        end select
        
#ifdef DPREC
        call fftw_execute_dft_r2c(forw, u, u_hat)
#else
        call fftwf_execute_dft_r2c(forw, u, u_hat)
#endif
        sp = real(u_hat(1:) * conjg(u_hat(1:)))
        sp = sp / n
        
        sp_avg = sp_avg + sp
        
        u_m = sum(u)/n
        do i = 1, n
          u_var = u_var + (u(i)-u_m)**2
        end do
      end do
    end do
    
    sp_avg = sp_avg / (nj*nk)
    
    u_var = u_var / (n*nj*nk)

    select case (dir)
      case(1)
        dx = (in_g%x(n)-in_g%x(1))/(n-1)
        lx = dx * n
      case(2)
        dx = (in_g%y(n)-in_g%y(1))/(n-1)
        lx = dx * n
      case(3)
        dx = (in_g%z(n)-in_g%z(1))/(n-1)
        lx = dx * n
    end select
    
    spectrum = dx * sp_avg / pi
    
    kappas = [(2*pi*j/lx, j = 1, size(sp_avg))]
        
  end subroutine
  
  subroutine get_buffer(status,sc,vec)
    integer,intent(in) :: status
    real(rp),intent(out) :: sc(:,:,:),vec(:,:,:,:)

    if (status==SCALAR) then
      call in_g%read_scalar(sc)
    else
      call in_g%read_vector(vec)
    end if
  end subroutine
  
  function get_field_name(line) result(res)
    character(:), allocatable :: res
    character(*), intent(in) :: line
    integer :: s1, s2
    s1 = scan(line,' ')
    s2 = scan(line(s1+1:),' ') + s1

    res = line(s1+1:s2-1)
  end function
  
  pure function array_variance(a) result(var)
    real(rp), intent(in) :: a(:,:,:)
    real(rp) :: am, var
    integer :: n
    n = size(a)
    am = sum(a) / n
    var = sum((a-am)**2)/n
  end function
end program
