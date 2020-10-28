!-------------------------------------------------------------------------------
! adsj C wrapper
!-------------------------------------------------------------------------------
!
! MODULE: adsj_c_wrapper [external members: `_adsj`]
!
! DESCRIPTION:
!> Hooks used in the nnlo-bridge code to interface to APPLgrid & fastNLO
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
module adsj_c_wrapper
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none
  public

contains

  !-----------------------------------------------------------------------------
  !> @brief
  !> hypergeometric Gauss 2F1 for real values
  !
  !> @param[in] a
  !> @param[in] b
  !> @param[in] c
  !> @param[in] z
  !> @return 2F1
  !-----------------------------------------------------------------------------
  function r2F1_adsj(a,b,c,z) bind(C)
    real(c_double) :: r2F1_adsj
    real(c_double), intent(in) :: a,b,c,z
    complex*16 :: f21
    external f21
    r2F1_adsj = real( f21(complex(a,0d0),complex(b,0d0),complex(c,0d0),complex(z,0d0)) )
  end function r2F1_adsj

  !-----------------------------------------------------------------------------
  !> @brief
  !> Appell F1 for real values
  !
  !> @param[in] a
  !> @param[in] b1
  !> @param[in] b2
  !> @param[in] c
  !> @param[in] x
  !> @param[in] y
  !> @return F1
  !-----------------------------------------------------------------------------
  function rF1_adsj(a,b1,b2,c,x,y) bind(C)
    real(c_double) :: rF1_adsj
    real(c_double), intent(in) :: a,b1,b2,c,x,y
    complex*16 :: f1
    external f1
    rF1_adsj = real( f1(complex(a,0d0),complex(b1,0d0),complex(b2,0d0),complex(c,0d0),complex(x,0d0),complex(y,0d0)) )
  end function rF1_adsj

end module adsj_c_wrapper
