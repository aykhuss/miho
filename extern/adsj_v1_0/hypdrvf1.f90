    !-------------------------------------------------------------------
    ! subroutine: hypdrvf1
    !
    ! package   : F1
    !
    ! Language  : Fortran 90
    !
    ! author    : F. Colavecchia (flavioc@lanl.gov)
    !                            (flavioc@cab.cnea.gov.ar)
    !
    ! date      : 3/26/97      version: 0.1
    ! revision  : 6/25/02      version: 1.0
    !
    ! purpose   :  Giving the function y(1) and its derivatives y(2),y(3)
    !			         compute the derivatives dyds at z=z0+s*dz, see eq. (33)
    !              of paper CPC 138 (2001) 29.
    !
    ! input     :    s  -> real variable
    !               y(3)-> function, 1st and 2nd derivatives at z
    !
    ! output    :  dyds(3)  ->  function, 1st and 2nd derivatives at z=z0+s*dz
    !
    ! modules   :
    !
    !
    ! common    :  /hypgf1/ aa,bb,bbp,cc,xx,yy,z0,dz
    !              contains the parameters (aa,bb,bbp,cc) 
    !              and variables (xx,yy) of the F1 function.
    !
    ! notes     :  most of the code is f77, few things come from f90.
    !
    !-------------------------------------------------------------------
    subroutine hypdrvf1(s,y,dyds)
    double precision s
    complex*16 y(3),dyds(3),aa,bb,bbp,cc,z0,dz,z,uu,vv,xx,yy
    complex*16 f0,f1,f2,f3
    common /hypgf1/ aa,bb,bbp,cc,xx,yy,z0,dz

    uu = xx
    vv = yy
    z=z0+s*dz
    f0 = aa*z*(bb*uu*(1 - cc + (1 + aa)*z*vv) + (1 - cc + (1 + aa)*z*uu)*vv*bbp)
    f1 = z*(cc**2 - cc*(1 + z*(uu + aa*uu + bb*uu + vv + aa*vv)) +             &
 	       z*uu*(aa**2*z*vv + 2*(1 + bb)*z*vv + aa*(-bb + 3*z*vv + 2*bb*z*vv)) + &
         z*(-aa - cc + 2*z*uu + 2*aa*z*uu)*vv*bbp)
    f2 = z**2*(-(cc*(-2 + z*(uu + vv))) -                                    & 
	       z*((2 + aa)*vv + uu*(2 + aa + bb - 4*z*vv - 2*aa*z*vv - bb*z*vv)) + &
         z*(-1 + z*uu)*vv*bbp)
    f3 = z**3*(-1 + z*uu)*(-1 + z*vv)
    dyds(1)=y(2)*dz
    dyds(2)=y(3)*dz
    dyds(3)=-dz/f3*(f2*y(3)+f1*y(2)+f0*y(1))
    return
    end
