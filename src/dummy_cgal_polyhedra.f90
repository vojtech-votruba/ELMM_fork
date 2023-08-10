module CGAL_Polyhedra
  use iso_c_binding
  use stop_procedures
  
  implicit none
  
  private
  
  public cgal_polyhedron_read, &
         cgal_polyhedron_closest, &
         cgal_polyhedron_inside, &
         cgal_polyhedron_intersects_ray, &
         cgal_polyhedron_bbox, &
         cgal_polyhedron_finalize
  
  type, bind(C) :: d3
    real(c_double) :: x, y, z
  end type
  
  interface cgal_polyhedron_closest
    module procedure cgal_polyhedron_closest_s
    module procedure cgal_polyhedron_closest_d
  end interface
  
  interface cgal_polyhedron_inside
    module procedure cgal_polyhedron_inside_s
    module procedure cgal_polyhedron_inside_d
  end interface
  
  interface cgal_polyhedron_intersects_ray
    module procedure cgal_polyhedron_intersects_ray_s
    module procedure cgal_polyhedron_intersects_ray_d
  end interface
  
  interface cgal_polyhedron_bbox
    module procedure cgal_polyhedron_bbox_s
    module procedure cgal_polyhedron_bbox_d
  end interface
  
contains


  subroutine cgal_polyhedron_read(ptree, fname)
    type(c_ptr),intent(out) :: ptree
    character(*),intent(in) :: fname

    call error_stop("Dummy CGAL!")
  end subroutine

  subroutine cgal_polyhedron_closest_s(ptree, xq,yq,zq, xn,yn,zn)
    type(c_ptr),intent(in) :: ptree
    real(c_float),intent(in)  :: xq,yq,zq
    real(c_float),intent(out) :: xn,yn,zn
    type(d3) :: query, near

    call error_stop("Dummy CGAL!")
  end subroutine
  
  subroutine cgal_polyhedron_closest_d(ptree, xq,yq,zq, xn,yn,zn)
    type(c_ptr),intent(in) :: ptree
    real(c_double),intent(in)  :: xq,yq,zq
    real(c_double),intent(out) :: xn,yn,zn
    type(d3) :: query, near

    call error_stop("Dummy CGAL!")
  end subroutine
  
  subroutine cgal_polyhedron_bbox_s(ptree, xmin,ymin,zmin, xmax,ymax,zmax)
    type(c_ptr),intent(in) :: ptree
    real(c_float),intent(out)  :: xmin,ymin,zmin
    real(c_float),intent(out) :: xmax,ymax,zmax
    type(d3) :: min, max

    call error_stop("Dummy CGAL!")
  end subroutine
  
  subroutine cgal_polyhedron_bbox_d(ptree, xmin,ymin,zmin, xmax,ymax,zmax)
    type(c_ptr),intent(in) :: ptree
    real(c_double),intent(out)  :: xmin,ymin,zmin
    real(c_double),intent(out) :: xmax,ymax,zmax
    type(d3) :: min, max

    call error_stop("Dummy CGAL!")
  end subroutine
  
  function cgal_polyhedron_inside_s(ptree, xq,yq,zq, xr,yr,zr) result(res)
    logical :: res
    type(c_ptr),intent(in) :: ptree
    real(c_float),intent(in) :: xq,yq,zq, xr,yr,zr
    type(d3) :: query, ref
    
    res = .false.
  end function

  function cgal_polyhedron_inside_d(ptree, xq,yq,zq, xr,yr,zr) result(res)
    logical :: res
    type(c_ptr),intent(in) :: ptree
    real(c_double),intent(in) :: xq,yq,zq, xr,yr,zr
    type(d3) :: query, ref
 
    res = .false.
  end function
  
  function cgal_polyhedron_intersects_ray_s(ptree, xc,yc,zc, a,b,c) result(res)
    logical :: res
    type(c_ptr),intent(in) :: ptree
    real(c_float),intent(in) :: xc,yc,zc, a,b,c
    
    res = .false.
  end function

  function cgal_polyhedron_intersects_ray_d(ptree, xc,yc,zc, a,b,c) result(res)
    logical :: res
    type(c_ptr),intent(in) :: ptree
    real(c_double),intent(in) :: xc,yc,zc, a,b,c
    
    res = .false.
  end function

  subroutine cgal_polyhedron_finalize(ptree)
    type(c_ptr),intent(inout) :: ptree
    
    call error_stop("Dummy CGAL!")
  end subroutine

  
end module CGAL_Polyhedra
