module MODEL

  !-----------------------------------------------------------------------------
  ! MODEL 
  ! A Poisson-Bolzmann equation solver
  ! Navier-Stokes Equation and Continuity equation
  ! \nabla pressure = -f = - \rho (u'_x^2 + 2 u'_y v'_x + v'_y^2)
  !-----------------------------------------------------------------------------

  use ESMF
  
  implicit none
  
  private

  integer, parameter     :: nx=400, ny=400
  real*8, pointer        :: pressure(:,:), new_p(:,:), f(:,:)
  type(ESMF_Grid)        :: model_grid
  
  public pressure
  public model_initialize
  public model_run
  public model_finalize
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine model_initialize(rc)
    integer, intent(out)          :: rc

    type(ESMF_Field)              :: field
    
    rc = ESMF_SUCCESS

    ! Create Model Grid
    model_grid =  ESMF_GridCreate1PeriDimUfrm(maxIndex=(/nx,ny/), &
      minCornerCoord=(/0._ESMF_KIND_R8, 0._ESMF_KIND_R8/), &
      maxCornerCoord=(/200._ESMF_KIND_R8, 200._ESMF_KIND_R8/), &
      coordSys=ESMF_COORDSYS_CART, &
      staggerLocList=(/ESMF_STAGGERLOC_CENTER/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out

    ! Create Model pressure Field, allocate memory with halo
    field = ESMF_FieldCreate(model_grid, typekind=ESMF_TYPEKIND_R8, &
      totalLWidth=(/1,1/), totalUWidth=(/1,1/), &
      name="pressure", &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out

    ! Retrieve pressure fortran array pointer from pressure field.
    call ESMF_FieldGet(field, farrayPtr=pressure, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out

    allocate(new_p(lbound(pressure,1):ubound(pressure,1), &
                   lbound(pressure,2):ubound(pressure,2)))
    allocate(f    (lbound(pressure,1):ubound(pressure,1), &
                   lbound(pressure,2):ubound(pressure,2)))

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine model_run

    integer                         :: i, j
    

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine model_finalize

    ! ESMF_Fields will be garbage collected

    ! Manually delete user allocated memory
    deallocate(new_p, f)

  end subroutine
  
  !-----------------------------------------------------------------------------

end module

