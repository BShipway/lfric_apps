!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Two-step Chebyshev semi-iterative combination kernel for PMSL
!>
!> @details Computes the combination step of the Golub-Varga (1961) two-step
!>          Chebyshev semi-iterative method:
!>
!>            x^{k+1} = gamma * J(x^k) + (1-gamma) * x^{k-1}
!>
!>          where J(x^k) is the standard Jacobi iterate and x^{k-1} = exner_prev.
!>          The Chebyshev weight recurrence is:
!>            gamma_1 = 1.0
!>            gamma_{k+1} = 1 / (1 - rho_J^2 * gamma_k / 4)
!>
!>          This formulation is stable even for gamma > 1 because the
!>          extrapolation uses x^{k-1}, providing damping that is absent
!>          in the equivalent single-step omega-Jacobi method.  The single-step
!>          stability limit omega < 2/(1+rho_J) ~ 1.001 would be violated by
!>          gamma_1 ~ 1.99 at large N; single-step is therefore unusable.
!>
!>          The algorithm module handles step 1 separately (pure Jacobi via
!>          pmsl_sor_kernel_mod with omega=1.0); this kernel is called for
!>          steps k = 2, 3, ..., n_iter_pmsl.

module pmsl_chebyshev_kernel_mod

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
  type, public, extends(kernel_type) :: pmsl_chebyshev_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                        &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),                           & ! dx_at_w2
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),    & ! f_func
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1,     &
                                                 STENCIL(CROSS2D)),             & ! exner_in  (x^k)
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                       & ! height_wth
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! exner_out (x^{k+1})
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),    & ! exner_prev (x^{k-1})
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                                 & ! gamma_k
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: pmsl_chebyshev_code
  end type pmsl_chebyshev_kernel_type

  public :: pmsl_chebyshev_code

contains

  !> @brief Two-step Chebyshev combination: x^{k+1} = gamma*J(x^k) + (1-gamma)*x^{k-1}
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     dx_at_w2      cell sizes at w2 dofs
  !> @param[in]     f_func        total forcing function for the PMSL Poisson equation
  !> @param[in]     exner_in      Current iterate x^k for exner at MSL
  !> @param[in]     smap_2d_size  Size of the stencil map in each direction
  !> @param[in]     sm_len        Max size of the stencil map in any direction
  !> @param[in]     smap_2d       Stencil map
  !> @param[in]     height_wth    Height of wth levels above mean sea level
  !> @param[in,out] exner_out     Next iterate x^{k+1}: gamma*J(x^k) + (1-gamma)*x^{k-1}
  !> @param[in]     exner_prev    Previous iterate x^{k-1}
  !> @param[in]     gamma_k       Chebyshev weight for this step
  !> @param[in]     ndf_w2        Number of degrees of freedom per cell for w2 fields
  !> @param[in]     undf_w2       Number of total degrees of freedom for w2 fields
  !> @param[in]     map_w2        Dofmap for the cell at the base of the column for w2 fields
  !> @param[in]     ndf_2d        Number of degrees of freedom per cell for 2d fields
  !> @param[in]     undf_2d       Number of total degrees of freedom for 2d fields
  !> @param[in]     map_2d        Dofmap for the cell at the base of the column for 2d fields
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for wtheta
  !> @param[in]     undf_wth      Number of total degrees of freedom for wtheta
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for wtheta
  subroutine pmsl_chebyshev_code(nlayers,                            &
                                 dx_at_w2,                          &
                                 f_func,                            &
                                 exner_in,                          &
                                 smap_2d_size, sm_len, smap_2d,     &
                                 height_wth,                        &
                                 exner_out,                         &
                                 exner_prev,                        &
                                 gamma_k,                           &
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
    real(kind=r_def),    intent(in), dimension(undf_2d)  :: exner_prev
    real(kind=r_def),    intent(in) :: gamma_k

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

      ! Standard Jacobi iterate J(x^k)
      exner_jacobi = pre_factor * (                             &
                       (1.0_r_def/dx2) *                       &
                      (exner_in(smap_2d(1,xm1,1)) +            &
                       exner_in(smap_2d(1,xp1,3)))             &
                     + (1.0_r_def/dy2) *                       &
                      (exner_in(smap_2d(1,ym1,2)) +            &
                       exner_in(smap_2d(1,yp1,4)))             &
                     - f_func(map_2d(1)) )

      ! Two-step Chebyshev combination:
      !   x^{k+1} = gamma * J(x^k) + (1-gamma) * x^{k-1}
      exner_out(map_2d(1)) = gamma_k * exner_jacobi &
                           + (1.0_r_def - gamma_k) * exner_prev(map_2d(1))

    else

      exner_out(map_2d(1)) = exner_in(smap_2d(1,1,1))

    end if

  end subroutine pmsl_chebyshev_code

end module pmsl_chebyshev_kernel_mod
