!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Holds code to solve the PMSL equation using Successive Over-Relaxation (SOR)
!>
!> @details This is a drop-in replacement for pmsl_solve_kernel_mod that applies an
!>          SOR relaxation update:
!>
!>            exner_sor = (dx2*dy2)/(2*dx2+2*dy2) * (N/dx2 + S/dx2 + E/dy2 + W/dy2 - f)
!>            exner_out = exner_in + omega * (exner_sor - exner_in)
!>                      = (1-omega)*exner_in + omega*exner_sor
!>
!>          where omega is the SOR relaxation parameter.  For the optimal value,
!>          set omega close to:
!>
!>            omega_opt = 2 / (1 + sin(pi/(N+1)))
!>
!>          where N = cells_per_side, giving:
!>            C12:  omega_opt ~ 1.730,  speedup ???x vs Jacobi (to be confirmed by timing)
!>            C48:  omega_opt ~ 1.874,  speedup ???x vs Jacobi (to be confirmed by timing)
!>            C224: omega_opt ~ 1.972,  speedup ???x vs Jacobi (to be confirmed by timing)
!>            C896: omega_opt ~ 1.993,  speedup ???x vs Jacobi (to be confirmed by timing)
!>
!>          With omega=1.0 this kernel is identical to the standard Jacobi kernel.

module pmsl_sor_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD, GH_SCALAR,       &
                                  GH_READ, GH_WRITE,         &
                                  GH_REAL, CELL_COLUMN,      &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  STENCIL, CROSS2D
  use fs_continuity_mod,    only: WTHETA, W2
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: pmsl_sor_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                        &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),                           & ! dx_at_w2
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),    & ! f_func
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1,     &
                                                 STENCIL(CROSS2D)),             & ! exner_in
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                       & ! height_wth
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! exner_out
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                                 & ! omega
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: pmsl_sor_code
  end type pmsl_sor_kernel_type

  public :: pmsl_sor_code

contains

  !> @brief Solve Poisson equation for PMSL using SOR relaxation
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     dx_at_w2      cell sizes at w2 dofs
  !> @param[in]     f_func        total forcing function to relax against PMSL
  !> @param[in]     exner_in      Current guess for exner at MSL
  !> @param[in]     smap_2d_size  Size of the stencil map in each direction
  !> @param[in]     sm_len        Max size of the stencil map in any direction
  !> @param[in]     smap_2d       Stencil map
  !> @param[in]     height_wth    Height of wth levels above mean sea level
  !> @param[in,out] exner_out     Next guess for exner at MSL after SOR update
  !> @param[in]     omega         SOR relaxation parameter (1.0 = Jacobi, 1<omega<2 = SOR)
  !> @param[in]     ndf_w2        Number of degrees of freedom per cell for w2 fields
  !> @param[in]     undf_w2       Number of total degrees of freedom for w2 fields
  !> @param[in]     map_w2        Dofmap for the cell at the base of the column for w2 fields
  !> @param[in]     ndf_2d        Number of degrees of freedom per cell for 2d fields
  !> @param[in]     undf_2d       Number of total degrees of freedom for 2d fields
  !> @param[in]     map_2d        Dofmap for the cell at the base of the column for 2d fields
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for wtheta
  !> @param[in]     undf_wth      Number of total degrees of freedom for wtheta
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for wtheta
  subroutine pmsl_sor_code(nlayers,                            &
                           dx_at_w2,                          &
                           f_func,                            &
                           exner_in,                          &
                           smap_2d_size, sm_len, smap_2d,     &
                           height_wth,                        &
                           exner_out,                         &
                           omega,                             &
                           ndf_w2, undf_w2, map_w2,           &
                           ndf_2d, undf_2d, map_2d,           &
                           ndf_wth, undf_wth, map_wth)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    integer(kind=i_def), intent(in) :: sm_len
    integer(kind=i_def), dimension(4), intent(in) :: smap_2d_size
    integer(kind=i_def), dimension(ndf_2d,sm_len,4), intent(in) :: smap_2d

    integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
    integer(kind=i_def), dimension(ndf_w2), intent(in)  :: map_w2

    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf_w2)  :: dx_at_w2
    real(kind=r_def),    intent(in), dimension(undf_2d)  :: f_func
    real(kind=r_def),    intent(in), dimension(undf_2d)  :: exner_in
    real(kind=r_def),    intent(in), dimension(undf_wth) :: height_wth
    real(kind=r_def),    intent(inout), dimension(undf_2d) :: exner_out
    real(kind=r_def),    intent(in) :: omega

    ! Internal variables
    real(kind=r_def) :: dx2, dy2, pre_factor, exner_jacobi
    real(kind=r_def), parameter :: pmsl_smooth_height = 500.0_r_def
    integer(kind=i_def) :: xp1, xm1, yp1, ym1

    ! Calculate which cell in the x branch of the stencil to use
    ! This sets the point to use to be the stencil point (2) if it exists,
    ! or the centre point (1) if it doesn't (i.e. we are at a domain edge)
    xp1 = min(2, smap_2d_size(3))
    xm1 = min(2, smap_2d_size(1))
    ! Calculate which cell in the y branch of the stencil to use
    ! This sets the point to use to be the stencil point (2) if it exists,
    ! or the centre point (1) if it doesn't (i.e. we are at a domain edge)
    yp1 = min(2, smap_2d_size(4))
    ym1 = min(2, smap_2d_size(2))

    ! Only calculated above a certain height and when all interior stencil points exist
    if (height_wth(map_wth(1)) > pmsl_smooth_height .and. &
        xp1 == 2 .and. xm1 == 2 .and. yp1 == 2 .and. ym1 == 2) then

      ! Calculate squared cell-centre distances (same as standard Jacobi kernel)
      dx2 = ((dx_at_w2(map_w2(1))+dx_at_w2(map_w2(3)))/2.0_r_def)**2_i_def
      dy2 = ((dx_at_w2(map_w2(2))+dx_at_w2(map_w2(4)))/2.0_r_def)**2_i_def

      ! Factor from re-arrangement of 5-point Poisson stencil
      pre_factor = (dx2*dy2) / (2.0_r_def*dx2 + 2.0_r_def*dy2)

      ! Standard Jacobi iterate (omega=1 result)
      exner_jacobi = pre_factor * (                             &
                       (1.0_r_def/dx2) *                       &
                      (exner_in(smap_2d(1,xm1,1)) +            &
                       exner_in(smap_2d(1,xp1,3)))             &
                     + (1.0_r_def/dy2) *                       &
                      (exner_in(smap_2d(1,ym1,2)) +            &
                       exner_in(smap_2d(1,yp1,4)))             &
                     - f_func(map_2d(1)) )

      ! SOR relaxation:  exner_out = (1-omega)*exner_old + omega*exner_jacobi
      exner_out(map_2d(1)) = (1.0_r_def - omega) * exner_in(smap_2d(1,1,1)) &
                            + omega * exner_jacobi

    else

      exner_out(map_2d(1)) = exner_in(smap_2d(1,1,1))

    end if

  end subroutine pmsl_sor_code

end module pmsl_sor_kernel_mod
