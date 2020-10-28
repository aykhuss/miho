!  rkf45.f90  02 June 2000
!
      subroutine rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
!
!     fehlberg fourth-fifth order runge-kutta method
!
!     written by h.a.watts and l.f.shampine
!                   sandia laboratories
!                  albuquerque,new mexico
!
!    rkf45 is primarily designed to solve non-stiff and mildly stiff
!    differential equations when derivative evaluations are inexpensive.
!    rkf45 should generally not be used when the user is demanding
!    high accuracy.
!
! abstract
!
!    subroutine  rkf45  integrates a system of neqn first order
!    ordinary differential equations of the form
!             dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
!              where the y(i) are given at t .
!    typically the subroutine is used to integrate from t to tout but it
!    can be used as a one-step integrator to advance the solution a
!    single step in the direction of tout.  on return the parameters in
!    the call list are set for continuing the integration. the user has
!    only to call rkf45 again (and perhaps define a new value for tout).
!    actually, rkf45 is an interfacing routine which calls subroutine
!    rkfs for the solution.  rkfs in turn calls subroutine  fehl which
!    computes an approximate solution over one step.
!
!    rkf45  uses the runge-kutta-fehlberg (4,5)  method described
!    in the reference
!    e.fehlberg , low-order classical runge-kutta formulas with stepsize
!                 control , nasa tr r-315
!
!    the performance of rkf45 is illustrated in the reference
!    l.f.shampine,h.a.watts,s.davenport, solving non-stiff ordinary
!                 differential equations-the state of the art ,
!                 sandia laboratories report sand75-0182 ,
!                 to appear in siam review.
!
!
!    the parameters represent-
!      f -- subroutine f(t,y,yp) to evaluate derivatives yp(i)=dy(i)/dt
!      neqn -- number of equations to be integrated
!      y(*) -- solution vector at t
!      t -- independent variable
!      tout -- output point at which solution is desired
!      relerr,abserr -- relative and absolute error tolerances for local
!            error test. at each step the code requires that
!                 abs(local error) .le. relerr*abs(y) + abserr
!            for each component of the local error and solution vectors
!      iflag -- indicator for status of integration
!      work(*) -- array to hold information internal to rkf45 which is
!            necessary for subsequent calls. must be dimensioned
!            at least  3+6*neqn
!      iwork(*) -- integer array used to hold information internal to
!            rkf45 which is necessary for subsequent calls. must be
!            dimensioned at least  5
!
!
!  first call to rkf45
!
!    the user must provide storage in his calling program for the arrays
!    in the call list  -      y(neqn) , work(3+6*neqn) , iwork(5)  ,
!    declare f in an external statement, supply subroutine f(t,y,yp) and
!    initialize the following parameters-
!
!      neqn -- number of equations to be integrated.  (neqn .ge. 1)
!      y(*) -- vector of initial conditions
!      t -- starting point of integration , must be a variable
!      tout -- output point at which solution is desired.
!            t=tout is allowed on the first call only, in which case
!            rkf45 returns with iflag=2 if continuation is possible.
!      relerr,abserr -- relative and absolute local error tolerances
!            which must be non-negative. relerr must be a variable while
!            abserr may be a constant. the code should normally not be
!            used with relative error control smaller than about 1.e-8 .
!            to avoid limiting precision difficulties the code requires
!            relerr to be larger than an internally computed relative
!            error parameter which is machine dependent. in particular,
!            pure absolute error is not permitted. if a smaller than
!            allowable value of relerr is attempted, rkf45 increases
!            relerr appropriately and returns control to the user before
!            continuing the integration.
!      iflag -- +1,-1  indicator to initialize the code for each new
!            problem. normal input is +1. the user should set iflag=-1
!            only when one-step integrator control is essential. in this
!            case, rkf45 attempts to advance the solution a single step
!            in the direction of tout each time it is called. since this
!            mode of operation results in extra computing overhead, it
!            should be avoided unless needed.
!
!
!  output from rkf45
!
!      y(*) -- solution at t
!      t -- last point reached in integration.
!      iflag = 2 -- integration reached tout. indicates successful retur
!                   and is the normal mode for continuing integration.
!            =-2 -- a single successful step in the direction of tout
!                   has been taken. normal mode for continuing
!                   integration one step at a time.
!            = 3 -- integration was not completed because relative error
!                   tolerance was too small. relerr has been increased
!                   appropriately for continuing.
!            = 4 -- integration was not completed because more than
!                   3000 derivative evaluations were needed. this
!                   is approximately 500 steps.
!            = 5 -- integration was not completed because solution
!                   vanished making a pure relative error test
!                   impossible. must use non-zero abserr to continue.
!                   using the one-step integration mode for one step
!                   is a good way to proceed.
!            = 6 -- integration was not completed because requested
!                   accuracy could not be achieved using smallest
!                   allowable stepsize. user must increase the error
!                   tolerance before continued integration can be
!                   attempted.
!            = 7 -- it is likely that rkf45 is inefficient for solving
!                   this problem. too much output is restricting the
!                   natural stepsize choice. use the one-step integrator
!                   mode.
!            = 8 -- invalid input parameters
!                   this indicator occurs if any of the following is
!                   satisfied -   neqn .le. 0
!                                 t=tout  and  iflag .ne. +1 or -1
!                                 relerr or abserr .lt. 0.
!                                 iflag .eq. 0  or  .lt. -2  or  .gt. 8
!      work(*),iwork(*) -- information which is usually of no interest
!                   to the user but necessary for subsequent calls.
!                   work(1),...,work(neqn) contain the first derivatives
!                   of the solution vector y at t. work(neqn+1) contains
!                   the stepsize h to be attempted on the next step.
!                   iwork(1) contains the derivative evaluation counter.
!
!
!  subsequent calls to rkf45
!
!    subroutine rkf45 returns with all information needed to continue
!    the integration. if the integration reached tout, the user need onl
!    define a new tout and call rkf45 again. in the one-step integrator
!    mode (iflag=-2) the user must keep in mind that each step taken is
!    in the direction of the current tout. upon reaching tout (indicated
!    by changing iflag to 2),the user must then define a new tout and
!    reset iflag to -2 to continue in the one-step integrator mode.
!
!    if the integration was not completed but the user still wants to
!    continue (iflag=3,4 cases), he just calls rkf45 again. with iflag=3
!    the relerr parameter has been adjusted appropriately for continuing
!    the integration. in the case of iflag=4 the function counter will
!    be reset to 0 and another 3000 function evaluations are allowed.
!
!    however,in the case iflag=5, the user must first alter the error
!    criterion to use a positive value of abserr before integration can
!    proceed. if he does not,execution is terminated.
!
!    also,in the case iflag=6, it is necessary for the user to reset
!    iflag to 2 (or -2 when the one-step integration mode is being used)
!    as well as increasing either abserr,relerr or both before the
!    integration can be continued. if this is not done, execution will
!    be terminated. the occurrence of iflag=6 indicates a trouble spot
!    (solution is changing rapidly,singularity may be present) and it
!    often is inadvisable to continue.
!
!    if iflag=7 is encountered, the user should use the one-step
!    integration mode with the stepsize determined by the code or
!    consider switching to the adams codes de/step,intrp. if the user
!    insists upon continuing the integration with rkf45, he must reset
!    iflag to 2 before calling rkf45 again. otherwise,execution will be
!    terminated.
!
!    if iflag=8 is obtained, integration can not be continued unless
!    the invalid input parameters are corrected.
!
!    it should be noted that the arrays work,iwork contain information
!    required for subsequent integration. accordingly, work and iwork
!    should not be altered.
!
!
      integer neqn,iflag,iwork(5)
      double precision y(neqn),t,tout,relerr,abserr,work(1)
!
      external f
!
      integer k1,k2,k3,k4,k5,k6,k1m
!
!
!     compute indices for the splitting of the work array
!
      k1m=neqn+1
      k1=k1m+1
      k2=k1+neqn
      k3=k2+neqn
      k4=k3+neqn
      k5=k4+neqn
      k6=k5+neqn
!
!     this interfacing routine merely relieves the user of a long
!     calling list via the splitting apart of two working storage
!     arrays. if this is not compatible with the users compiler,
!     he must use rkfs directly.
!
      call rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,work(1),work(k1m),
     1          work(k1),work(k2),work(k3),work(k4),work(k5),work(k6),
     2          work(k6+1),iwork(1),iwork(2),iwork(3),iwork(4),iwork(5))
!
      return
      end
      subroutine rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,yp,h,f1,f2,f3,
     1                f4,f5,savre,savae,nfe,kop,init,jflag,kflag)
!
!     fehlberg fourth-fifth order runge-kutta method
!
!
!     rkfs integrates a system of first order ordinary differential
!     equations as described in the comments for rkf45 .
!     the arrays yp,f1,f2,f3,f4,and f5 (of dimension at least neqn) and
!     the variables h,savre,savae,nfe,kop,init,jflag,and kflag are used
!     internally by the code and appear in the call list to eliminate
!     local retention of variables between calls. accordingly, they
!     should not be altered. items of possible interest are
!         yp - derivative of solution vector at t
!         h  - an appropriate stepsize to be used for the next step
!         nfe- counter on the number of derivative function evaluations
!
!
      logical hfaild,output
!
      integer  neqn,iflag,nfe,kop,init,jflag,kflag
      double precision  y(neqn),t,tout,relerr,abserr,h,yp(neqn),
     1  f1(neqn),f2(neqn),f3(neqn),f4(neqn),f5(neqn),savre,
     2  savae
!
      external f
!
      double precision  a,ae,dt,ee,eeoet,esttol,et,hmin,remin,rer,s,
     1  scale,tol,toln,twoeps,u26,ypk
!
      integer  k,maxnfe,mflag
!
      double precision  dabs,dmax1,dmin1,dsign,d1mach
!
!  remin is the minimum acceptable value of relerr.  attempts
!  to obtain higher accuracy with this subroutine are usually
!  very expensive and often unsuccessful.
!
      data remin/1.d-12/
!
!
!     the expense is controlled by restricting the number
!     of function evaluations to be approximately maxnfe.
!     as set, this corresponds to about 500 steps.
!
      data maxnfe/3000/
!
!   here two constants emboding the machine epsilon is present
!   twoesp is set to twice the machine epsilon while u26 is set
!   to 26 times the machine epsilon
!
!     data twoeps, u26/4.4d-16, 5.72d-15/                               ***
      twoeps = 2.*d1mach(4)                                             ***
      u26 = 13.*twoeps                                                  ***
!
!
!     check input parameters
!
!
      if (neqn .lt. 1) go to 10
      if ((relerr .lt. 0.0d0)  .or.  (abserr .lt. 0.0d0)) go to 10
      mflag=iabs(iflag)
      if ((mflag .ge. 1)  .and.  (mflag .le. 8)) go to 20
!
!     invalid input
   10 iflag=8
      return
!
!     is this the first call
   20 if (mflag .eq. 1) go to 50
!
!     check continuation possibilities
!
      if ((t .eq. tout) .and. (kflag .ne. 3)) go to 10
      if (mflag .ne. 2) go to 25
!
!     iflag = +2 or -2
      if (kflag .eq. 3) go to 45
      if (init .eq. 0) go to 45
      if (kflag .eq. 4) go to 40
      if ((kflag .eq. 5)  .and.  (abserr .eq. 0.0d0)) go to 30
      if ((kflag .eq. 6)  .and.  (relerr .le. savre)  .and.
     1    (abserr .le. savae)) go to 30
      go to 50
!
!     iflag = 3,4,5,6,7 or 8
   25 if (iflag .eq. 3) go to 45
      if (iflag .eq. 4) go to 40
      if ((iflag .eq. 5) .and. (abserr .gt. 0.0d0)) go to 45
!
!     integration cannot be continued since user did not respond to
!     the instructions pertaining to iflag=5,6,7 or 8
   30 stop
!
!     reset function evaluation counter
   40 nfe=0
      if (mflag .eq. 2) go to 50
!
!     reset flag value from previous call
   45 iflag=jflag
      if (kflag .eq. 3) mflag=iabs(iflag)
!
!     save input iflag and set continuation flag value for subsequent
!     input checking
   50 jflag=iflag
      kflag=0
!
!     save relerr and abserr for checking input on subsequent calls
      savre=relerr
      savae=abserr
!
!     restrict relative error tolerance to be at least as large as
!     2*eps+remin to avoid limiting precision difficulties arising
!     from impossible accuracy requests
!
      rer=twoeps+remin
      if (relerr .ge. rer) go to 55
!
!     relative error tolerance too small
      relerr=rer
      iflag=3
      kflag=3
      return
!
   55 dt=tout-t
!
      if (mflag .eq. 1) go to 60
      if (init .eq. 0) go to 65
      go to 80
!
!     initialization --
!                       set initialization completion indicator,init
!                       set indicator for too many output points,kop
!                       evaluate initial derivatives
!                       set counter for function evaluations,nfe
!                       evaluate initial derivatives
!                       set counter for function evaluations,nfe
!                       estimate starting stepsize
!
   60 init=0
      kop=0
!
      a=t
      call f(a,y,yp)
      nfe=1
      if (t .ne. tout) go to 65
      iflag=2
      return
!
!
   65 init=1
      h=dabs(dt)
      toln=0.
      do 70 k=1,neqn
        tol=relerr*dabs(y(k))+abserr
        if (tol .le. 0.) go to 70
        toln=tol
        ypk=dabs(yp(k))
        if (ypk*h**5 .gt. tol) h=(tol/ypk)**0.2d0
   70 continue
      if (toln .le. 0.0d0) h=0.0d0
      h=dmax1(h,u26*dmax1(dabs(t),dabs(dt)))
      jflag=isign(2,iflag)
!
!
!     set stepsize for integration in the direction from t to tout
!
   80 h=dsign(h,dt)
!
!     test to see if rkf45 is being severely impacted by too many
!     output points
!
      if (dabs(h) .ge. 2.0d0*dabs(dt)) kop=kop+1
      if (kop .ne. 100) go to 85
!
!     unnecessary frequency of output
      kop=0
      iflag=7
      return
!
   85 if (dabs(dt) .gt. u26*dabs(t)) go to 95
!
!     if too close to output point,extrapolate and return
!
      do 90 k=1,neqn
   90   y(k)=y(k)+dt*yp(k)
      a=tout
      call f(a,y,yp)
      nfe=nfe+1
      go to 300
!
!
!     initialize output point indicator
!
   95 output= .false.
!
!     to avoid premature underflow in the error tolerance function,
!     scale the error tolerances
!
      scale=2.0d0/relerr
      ae=scale*abserr
!
!
!     step by step integration
!
  100 hfaild= .false.
!
!     set smallest allowable stepsize
!
      hmin=u26*dabs(t)
!
!     adjust stepsize if necessary to hit the output point.
!     look ahead two steps to avoid drastic changes in the stepsize and
!     thus lessen the impact of output points on the code.
!
      dt=tout-t
      if (dabs(dt) .ge. 2.0d0*dabs(h)) go to 200
      if (dabs(dt) .gt. dabs(h)) go to 150
!
!     the next successful step will complete the integration to the
!     output point
!
      output= .true.
      h=dt
      go to 200
!
  150 h=0.5d0*dt
!
!
!
!     core integrator for taking a single step
!
!     the tolerances have been scaled to avoid premature underflow in
!     computing the error tolerance function et.
!     to avoid problems with zero crossings,relative error is measured
!     using the average of the magnitudes of the solution at the
!     beginning and end of a step.
!     the error estimate formula has been grouped to control loss of
!     significance.
!     to distinguish the various arguments, h is not permitted
!     to become smaller than 26 units of roundoff in t.
!     practical limits on the change in the stepsize are enforced to
!     smooth the stepsize selection process and to avoid excessive
!     chattering on problems having discontinuities.
!     to prevent unnecessary failures, the code uses 9/10 the stepsize
!     it estimates will succeed.
!     after a step failure, the stepsize is not allowed to increase for
!     the next attempted step. this makes the code more efficient on
!     problems having discontinuities and more effective in general
!     since local extrapolation is being used and extra caution seems
!     warranted.
!
!
!     test number of derivative function evaluations.
!     if okay,try to advance the integration from t to t+h
!
  200 if (nfe .le. maxnfe) go to 220
!
!     too much work
      iflag=4
      kflag=4
      return
!
!     advance an approximate solution over one step of length h
!
  220 call fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,f1)
      nfe=nfe+5
!
!     compute and test allowable tolerances versus local error estimates
!     and remove scaling of tolerances. note that relative error is
!     measured with respect to the average of the magnitudes of the
!     solution at the beginning and end of the step.
!
      eeoet=0.0d0
      do 250 k=1,neqn
        et=dabs(y(k))+dabs(f1(k))+ae
        if (et .gt. 0.0d0) go to 240
!
!       inappropriate error tolerance
        iflag=5
        return
!
  240   ee=dabs((-2090.0d0*yp(k)+(21970.0d0*f3(k)-15048.0d0*f4(k)))+
     1                        (22528.0d0*f2(k)-27360.0d0*f5(k)))
  250   eeoet=dmax1(eeoet,ee/et)
!
      esttol=dabs(h)*eeoet*scale/752400.0d0
!
      if (esttol .le. 1.0d0) go to 260
!
!
!     unsuccessful step
!                       reduce the stepsize , try again
!                       the decrease is limited to a factor of 1/10
!
      hfaild= .true.
      output= .false.
      s=0.1d0
      if (esttol .lt. 59049.0d0) s=0.9d0/esttol**0.2d0
      h=s*h
      if (dabs(h) .gt. hmin) go to 200
!
!     requested error unattainable at smallest allowable stepsize
      iflag=6
      kflag=6
      return
!
!
!     successful step
!                        store solution at t+h
!                        and evaluate derivatives there
!
  260 t=t+h
      do 270 k=1,neqn
  270   y(k)=f1(k)
      a=t
      call f(a,y,yp)
      nfe=nfe+1
!
!
!                       choose next stepsize
!                       the increase is limited to a factor of 5
!                       if step failure has just occurred, next
!                          stepsize is not allowed to increase
!
      s=5.0d0
      if (esttol .gt. 1.889568d-4) s=0.9d0/esttol**0.2d0
      if (hfaild) s=dmin1(s,1.0d0)
      h=dsign(dmax1(s*dabs(h),hmin),h)
!
!     end of core integrator
!
!
!     should we take another step
!
      if (output) go to 300
      if (iflag .gt. 0) go to 100
!
!
!     integration successfully completed
!
!     one-step mode
      iflag=-2
      return
!
!     interval mode
  300 t=tout
      iflag=2
      return
!
      end
      subroutine fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,s)
!
!     fehlberg fourth-fifth order runge-kutta method
!
!    fehl integrates a system of neqn first order
!    ordinary differential equations of the form
!             dy(i)/dt=f(t,y(1),---,y(neqn))
!    where the initial values y(i) and the initial derivatives
!    yp(i) are specified at the starting point t. fehl advances
!    the solution over the fixed step h and returns
!    the fifth order (sixth order accurate locally) solution
!    approximation at t+h in array s(i).
!    f1,---,f5 are arrays of dimension neqn which are needed
!    for internal storage.
!    the formulas have been grouped to control loss of significance.
!    fehl should be called with an h not smaller than 13 units of
!    roundoff in t so that the various independent arguments can be
!    distinguished.
!
!
      integer  neqn
      double precision  y(neqn),t,h,yp(neqn),f1(neqn),f2(neqn),
     1  f3(neqn),f4(neqn),f5(neqn),s(neqn)
!
      double precision  ch
      integer  k
!
      ch=h/4.0d0
      do 221 k=1,neqn
  221   f5(k)=y(k)+ch*yp(k)
      call f(t+ch,f5,f1)
!
      ch=3.0d0*h/32.0d0
      do 222 k=1,neqn
  222   f5(k)=y(k)+ch*(yp(k)+3.0d0*f1(k))
      call f(t+3.0d0*h/8.0d0,f5,f2)
!
      ch=h/2197.0d0
      do 223 k=1,neqn
  223   f5(k)=y(k)+ch*(1932.0d0*yp(k)+(7296.0d0*f2(k)-7200.0d0*f1(k)))
      call f(t+12.0d0*h/13.0d0,f5,f3)
!
      ch=h/4104.0d0
      do 224 k=1,neqn
  224   f5(k)=y(k)+ch*((8341.0d0*yp(k)-845.0d0*f3(k))+
     1                            (29440.0d0*f2(k)-32832.0d0*f1(k)))
      call f(t+h,f5,f4)
!
      ch=h/20520.0d0
      do 225 k=1,neqn
  225   f1(k)=y(k)+ch*((-6080.0d0*yp(k)+(9295.0d0*f3(k)-
     1         5643.0d0*f4(k)))+(41040.0d0*f1(k)-28352.0d0*f2(k)))
      call f(t+h/2.0d0,f1,f5)
!
!     compute approximate solution at t+h
!
      ch=h/7618050.0d0
      do 230 k=1,neqn
  230   s(k)=y(k)+ch*((902880.0d0*yp(k)+(3855735.0d0*f3(k)-
     1        1371249.0d0*f4(k)))+(3953664.0d0*f2(k)+
     2        277020.0d0*f5(k)))
!
      return
      end

