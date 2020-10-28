    !-------------------------------------------------------------------
    ! program   : f1
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
    ! purpose   :  Driver for testing the hypergeometric F1 routines
    !              Computes the F1 function in a square lattice of
    !              nmax points from -range to range in each variable
    !              Tests give data results as shown in the paper CPC 138 (2001) 29,
    !              in the following refered as paper I. 
    !
    !              The program asks for a test, 1-4. 
    !              Test 1: read the input parameters and variables from file
    !                      f1.data and computes the f1 function for those values.
    !              Test 2-4: Generate the Tables 2-4 from paper 1.
    !
    !
    ! modules   :
    !
    !
    ! common    :
    !
    !
    ! notes     :  most of the code is f77, few things come from f90.
    !
    !
    !-------------------------------------------------------------------
    program driver
    implicit none

    complex*16 ca,cb,cbp,cc,ccp,c1,ci,cw,cx,cy,cz
    complex*16 cf1,cf1p,f1,f2,g2,f21,caux,cf1pp
    complex*16 f1bnl,f1inte,hypgeof1,hypgeo,f1conv
    real*8 ar,ai,br,bi,bpr,bpi,ccr,cci,ccpr,ccpi
    real*8 test,x,y,step,range,delx
    real*8 eps,zero
    integer*4 n,nmax,ni,nj,randomtest,iresult
    integer ic,flag_view,itest
    external f1,f21_new
    !
    !  Some Constants
    !
    c1=(1.d0,0.d0)
    ci=(0.d0,1.d0)
    zero = 1d-12
!    goto 111
    write(*,*) '------------------------------------------'
    write(*,*) '             Select the test'
    write(*,*) '------------------------------------------'
    write(*,*) 'Test 1  : read parameters and variables'
    write(*,*) '          from f1.data   '
    write(*,*) 'Test 2-4: Generate Tables 2-4 of '
    write(*,*) '          Comp. Phys. Comm. 138 (2001) 29'
    write(*,*) '------------------------------------------'
    do while(itest.gt.4.or.itest.le.0)
        read(*,*) itest
    end do
    
!111    itest=4
    flag_view = 1
    !
    !  by default the results do not show in the screen
    !  uncomment next line to display the computation in the screen
    !flag_view = 1
    !  
    !    Tables
    !
    select case(itest)
        case(1) 
            !
            !  Write your data in file f1.data for
            !  computation of f1
            !
            open(unit=100,file='f1.data')
            read(100,*) ca
            read(100,*) cb
            read(100,*) cbp
            read(100,*) cc
            read(100,*) x
            read(100,*) y
            close(100)           
!            ca = -0.5d0
!            cb = 2d0
!            cbp = 1d0
!            cc = cb+cbp            
!            x = 2.5d0
!            y = 1.5d0            
                cf1   = f1bnl(ca,cb,cbp,cc,x,y)                           ! f1 by series
            !cf1 = f1(ca,cb,cbp,cc,x,y)
            call writef1(6,ca,cb,cbp,cc,x*(1d0,0d0),y*(1d0,0d0),cf1)
            !write(*,*) '|f1| = ',cdabs(cf1)
        case(2)
            !
            !  Table 2 of the paper I: Comparison between the 
            !  function in the convergence region and the exact
            !  result. See paper for a discussion.
            !
            !  grid data
            range = 0.95d0
            nmax = 6
            step = 2*range/(nmax-1)
            !  Example data
            ca = 1d0
            cb = 2+ci
            cbp = 1.5d0-0.5*ci
            cc = ca
            open(unit=10,file='table2.dat',status='unknown')
            write(10,*) 'Table 2 of paper 1, corrected'
            write(*,*) 'Table 2 of paper 1, corrected'
            !   Loop
            delx = 0.42d0
            do ni=0,nmax-2
              x = -range  + delx*ni
              do nj=0,nmax-1
                y = -range + 0.38d0*nj
                cf1p  = hypgeof1(ca,cb,cbp,cc,x*(1d0,0d0),y*(1d0,0d0))  !  ODE Integration
                cf1pp = (1d0-x)**(-cb)*(1d0-y)**(-cbp)                  !   Exact result
                write(10,'(2f7.2,1x,2f26.16,1x,1es8.1)') x,y,cdabs(cf1p),cdabs(cf1pp),         &
                dabs((cdabs(cf1p)-cdabs(cf1pp))/cdabs(cf1pp))
                if(flag_view.eq.1) then
                  write(*,'(2f7.2,1x,2f26.16,1x,1es8.1)') x,y,cdabs(cf1p),cdabs(cf1pp),         &
                  dabs((cdabs(cf1p)-cdabs(cf1pp))/cdabs(cf1pp))
                endif
                call flush(10)
              end do  ! nj loop
            end do    ! ni loop
            close(10)
        case(3)
            !
            !  Table 3 of the paper I: Comparison between the 
            !  functions in the convergence region, computed with ODE integration
            !  and the series function, and the exact
            !  result. See paper for a discussion.
            !
            !  grid data
            range = 0.95d0
            nmax = 6
            step = 2*range/(nmax-1)
            !  Example data
            ca = 1d0
            cb = 3+ci
            cbp = 2d0-0.5*ci
            cc = cb+cbp
            !   Loop
            open(unit=10,file='table3.dat',status='unknown')
            write(10,*) 'Table 3 of paper 1'
            write(*,*) 'Table 3 of paper 1'
            do ni=0,nmax-1
              x = -range + step*ni
              do nj=0,nmax-1
                y = -range + step*nj
                cf1   = f1bnl(ca,cb,cbp,cc,x,y)                           ! f1 by series
                cf1p  = hypgeof1(ca,cb,cbp,cc,x*(1d0,0d0),y*(1d0,0d0))    ! f1 by ODE integration
                cf1pp = (1d0-y)**(-ca)*f21(ca,cb,cbp+cb,(y-x)/(y-1)*c1)   ! 2F1 reduction
                write(10,'(2f7.2,1x,3f18.8,1x,2es8.1)') x,y,cdabs(cf1),cdabs(cf1p),cdabs(cf1pp),         &
                dabs((cdabs(cf1)-cdabs(cf1pp))/cdabs(cf1pp)),dabs((cdabs(cf1p)-cdabs(cf1pp))/cdabs(cf1pp))
                if(flag_view.eq.1) then
                  write(*,'(2f7.2,1x,3f18.8,1x,2es8.1)') x,y,cdabs(cf1),cdabs(cf1p),cdabs(cf1pp),         &
                  dabs((cdabs(cf1)-cdabs(cf1pp))/cdabs(cf1pp)),dabs((cdabs(cf1p)-cdabs(cf1pp))/cdabs(cf1pp))
                endif
                call flush(10)
              end do  ! nj loop
            end do    ! ni loop
            close(10)
      case(4)
        !  grid data
        range = 3.5d0
        nmax = 8
        step = 2*range/(nmax-1)
        !  Example data
        ca = -0.5d0
        cb = 2d0
        cbp = 1d0
        cc = cb+cbp
        !   Loop
        open(unit=10,file='table4.dat',status='unknown')
        write(10,*) 'Table 4 of paper 1'
        write(*,*)  'Table 4 of paper 1'
        do ni=0,nmax-1
          x = -range + step*ni
          do nj=0,nmax-1
            y = -range + step*nj
            cf1 = f1(ca,cb,cbp,cc,x,y)        ! F1 function
            cf1p= (1d0-y)**(-ca)*f21(ca,cb,cbp+cb,(y-x)/(y-1)*c1)  ! 2F1 reduction
            cz = (y-x)/(y-1)*c1          ! Exact case
            cf1pp = 4/15.0*(2-(2+3*cz)*(cdsqrt(1-cz))**3)*(1-y)**(-ca)/cz**2
            if(cdabs(cz).lt.zero) then      ! Limiting case of the exact expression
              cf1pp = (1-y)**(-ca)      ! for z->0
            end if
            write(10,'(2f7.2,1x,3f18.8,1x,2es8.1)') x,y,cdabs(cf1),cdabs(cf1p),cdabs(cf1pp),         &
            dabs((cdabs(cf1)-cdabs(cf1pp))/cdabs(cf1pp)),dabs((cdabs(cf1p)-cdabs(cf1pp))/cdabs(cf1pp))
            call flush(10)
            if(flag_view.eq.1) then
              write(*,'(2f7.2,1x,3f18.8,1x,2es8.1)') x,y,cdabs(cf1),cdabs(cf1p),cdabs(cf1pp),         &
              dabs((cdabs(cf1)-cdabs(cf1pp))/cdabs(cf1pp)),dabs((cdabs(cf1p)-cdabs(cf1pp))/cdabs(cf1pp))
            endif
          end do  ! nj loop
        end do    ! ni loop
        close(10)
      case default
        write(*,*)'Unimplemented test '
        stop
    end select

5353   format(i1)
5454  format(a7,2f10.5)
5555  format(f10.5)
5656  format(a7,2e10.5)
5757  format(a7,e10.5)


end

