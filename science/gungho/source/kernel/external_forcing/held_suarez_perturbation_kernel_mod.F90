!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Perturbs the Held-Suarez theta increment by a small random factor.
!>
!> The perturbation was part of held_suarez_fv_kernel_mod; it lives in a
!> kernel of its own so that the forcing kernel draws no random numbers and is
!> invoked only when hs_random is set.
!>
module held_suarez_perturbation_kernel_mod

  use argument_mod,  only: arg_type,          &
                           GH_FIELD, GH_REAL, &
                           GH_READWRITE, CELL_COLUMN
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type
  use fs_continuity_mod, only: Wtheta

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: held_suarez_perturbation_kernel_type
    private
    type(arg_type) :: meta_args(1) = (/                  &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta) &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: held_suarez_perturbation_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: held_suarez_perturbation_code

contains

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in,out] dtheta Real array, theta increment data
!! @param[in] ndf_wth The number of degrees of freedom per cell for wth
!! @param[in] undf_wth The number of unique degrees of freedom for wth
!! @param[in] map_wth Integer array holding the dofmap for the cell at the
!>            base of the column for wth
subroutine held_suarez_perturbation_code(nlayers, dtheta,              &
                                         ndf_wth, undf_wth, map_wth)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth

  real(kind=r_def), dimension(undf_wth), intent(inout) :: dtheta

  integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth

  ! Internal variables
  real(kind=r_def) :: pert(nlayers+1)

  call random_number(pert)

  dtheta(map_wth(1):map_wth(1)+nlayers) = dtheta(map_wth(1):map_wth(1)+nlayers) &
     *(1.0_r_def + (pert-0.5_r_def)*.0001_r_def)

end subroutine held_suarez_perturbation_code

end module held_suarez_perturbation_kernel_mod
