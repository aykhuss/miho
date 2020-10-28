c***********************************************************************
c
c   program :  complex hypergeometric function
c
c   notation :  F(a,b;c;z)
c
c   usage:   The program has been written to be used as a function
c            where the user supplies the hypergeometric parameters
c            (a,b, and c), the independent variable (z).
c
c   warning:  This code is under construction. Warning messages
c             will appear in regions where the code has not yet
c             been completed.
c
c   written by:
c
c        Robert C. Forrey
c        Institute for Theoretical Atomic and Molecular Physics
c        Harvard-Smithsonian Center for Astrophysics
c        60 Garden Street Cambridge, MA 02138
c        rforrey@cfa.harvard.edu
c
c***********************************************************************
c
c  function name      - chyp(a,b,c,z)
c
c  computation
c  performed          - calculates the hypergeometric function.
c
c  arguments       z  - the independent variable of the hypergeometric
c                       function.
c
c               a,b,c - parameters of the hypergeometric function.
c
c  precision          - complex*16
c
c  language           - fortran
c
c***********************************************************************

      complex*16 function chyp(a,b,c,z)

      implicit none
      real*8  zero,one,pi,tmp1,tmp2,x,y,r,phi
      parameter (zero=0.d0,one=1.d0)
      integer flag,n
      complex*16  i,a,b,c,z,w,f1,f2,cgamm,csum,
     #            ctmp,ctmp1,ctmp2,ctmp3,ctmp4,ogamm

      pi=dacos(-1.d0)
      i=dcmplx(0.d0,1.d0)

c  Set error messages

      tmp1=dreal(a-b)
      tmp2=dimag(a-b)
      n=nint(tmp1)
      tmp1=tmp1+tmp2-n
!      if(tmp1.eq.zero)then
!        write(*,*)'transformation error in chyp'
!        return
!      endif

      tmp1=dreal(c-a-b)
      tmp2=dimag(c-a-b)
      n=nint(tmp1)
      tmp1=tmp1+tmp2-n
!      if(tmp1.eq.zero)then
!        write(*,*)'transformation error in chyp'
!        return
!      endif

c  Handle the special cases when z=0 or z=1

      tmp1=cdabs(z)
      if (tmp1.eq.zero) then
        chyp=dcmplx(one)
        return
      endif

      ctmp=dcmplx(one)
      if (z.eq.ctmp) then
c       chyp=cgamm(c)*cgamm(c-a-b)/cgamm(c-a)/cgamm(c-b)
        chyp=cgamm(c)*cgamm(c-a-b)*ogamm(c-a)*ogamm(c-b)
        return
      endif

c  Transform to a new variable w whose magnitude is the 
c  lowest possible value that lies between 0 and 1.

      tmp1=cdabs(z)
      flag=3

      w=1/(1-z)
      tmp2=cdabs(w)
      if(tmp2 .lt. tmp1) then
        tmp1=tmp2
        flag=1
      endif

      w=z/(z-1)
      tmp2=cdabs(w)
      if(tmp2 .lt. tmp1) then
        tmp1=tmp2
        flag=2
      endif

      w=1-z
      tmp2=cdabs(w)
      if(tmp2 .lt. tmp1) then
        tmp1=tmp2
        flag=4
      endif

      w=1-1/z
      tmp2=cdabs(w)
      if(tmp2 .lt. tmp1) then
        tmp1=tmp2
        flag=5
      endif

      w=1/z
      tmp2=cdabs(w)
      if(tmp2 .lt. tmp1) then
        tmp1=tmp2
        flag=6
      endif

      if(tmp1.ge.one)then
        write(*,*)'error in transformation'
        return
      endif

      write(20,*)'flag=',flag

c  Compute the hypergeometric function in z via the transformation
c  theory and the series representation in w.

      if (flag .eq. 1)then

        w=1/(1-z)
        f1=csum(a,c-b,a-b+1,w)
        f2=csum(b,c-a,b-a+1,w)
        x=dreal(1-z)
        y=dimag(1-z)
        if(x.gt.zero .and. y.ge.zero)then
          phi=datan(dabs(y/x))
        elseif(x.lt.zero .and. y.gt.zero)then
          phi=pi-datan(dabs(y/x))
        elseif(x.lt.zero .and. y.lt.zero)then
          phi=-pi+datan(dabs(y/x))
        elseif(x.gt.zero .and. y.lt.zero)then
          phi=-datan(dabs(y/x))
        elseif(x.eq.zero .and. y.gt.zero)then
          phi=pi/2
        elseif(x.eq.zero .and. y.lt.zero)then
          phi=-pi/2
        else
          write(*,*)'error in transformation 1'
          return
        endif
        r=dsqrt(x**2+y**2)
        ctmp=r*cdexp(i*phi)
        ctmp1=ctmp**(-a)
        ctmp2=ctmp**(-b)
c       chyp=f1*ctmp1*cgamm(c)*cgamm(b-a)/cgamm(b)/cgamm(c-a)
c    #      +f2*ctmp2*cgamm(c)*cgamm(a-b)/cgamm(a)/cgamm(c-b)
        chyp=f1*ctmp1*cgamm(c)*cgamm(b-a)*ogamm(b)*ogamm(c-a)
     #      +f2*ctmp2*cgamm(c)*cgamm(a-b)*ogamm(a)*ogamm(c-b)
c       write(*,'(2(2 e13.6,2x))')ctmp1,ctmp1-w**a
c       write(*,'(2(2 e13.6,2x))')ctmp2,ctmp2-w**b

      elseif (flag .eq. 2)then

        w=z/(z-1)
        chyp=(1-w)**a*csum(a,c-b,c,w)

      elseif (flag .eq. 3)then

        w=z
        chyp=csum(a,b,c,w)

      elseif (flag .eq. 4)then

        w=1-z
        f1=csum(a,b,a+b-c+1,w)
        f2=csum(c-a,c-b,c-a-b+1,w)
        x=dreal(w)
        y=dimag(w)
        if(x.gt.zero .and. y.ge.zero)then
          phi=datan(dabs(y/x))
        elseif(x.lt.zero .and. y.gt.zero)then
          phi=pi-datan(dabs(y/x))
        elseif(x.lt.zero .and. y.lt.zero)then
          phi=-pi+datan(dabs(y/x))
        elseif(x.gt.zero .and. y.lt.zero)then
          phi=-datan(dabs(y/x))
        elseif(x.eq.zero .and. y.gt.zero)then
          phi=pi/2
        elseif(x.eq.zero .and. y.lt.zero)then
          phi=-pi/2
        else
          write(*,*)'error in transformation 4'
          return
        endif
        r=dsqrt(x**2+y**2)
        ctmp=r*cdexp(i*phi)
        ctmp2=ctmp**(c-a-b)
c       chyp=f1*cgamm(c)*cgamm(c-a-b)/cgamm(c-a)/cgamm(c-b)
c    #      +f2*ctmp2*cgamm(c)*cgamm(a+b-c)/cgamm(a)/cgamm(b)
        chyp=f1*cgamm(c)*cgamm(c-a-b)*ogamm(c-a)*ogamm(c-b)
     #      +f2*ctmp2*cgamm(c)*cgamm(a+b-c)*ogamm(a)*ogamm(b)
c       write(*,'(2(2 e13.6,2x))')ctmp2,ctmp2-w**(c-a-b)

      elseif (flag .eq. 5)then

        w=1-1/z
        f1=csum(a,a-c+1,a+b-c+1,w)
        f2=csum(c-a,1-a,c-a-b+1,w)
        x=dreal(z)
        y=dimag(z)
        if(x.gt.zero .and. y.ge.zero)then
          phi=datan(dabs(y/x))
        elseif(x.lt.zero .and. y.gt.zero)then
          phi=pi-datan(dabs(y/x))
        elseif(x.lt.zero .and. y.lt.zero)then
          phi=-pi+datan(dabs(y/x))
        elseif(x.gt.zero .and. y.lt.zero)then
          phi=-datan(dabs(y/x))
        elseif(x.eq.zero .and. y.gt.zero)then
          phi=pi/2
        elseif(x.eq.zero .and. y.lt.zero)then
          phi=-pi/2
        else
          write(*,*)'error in transformation 5'
          return
        endif
        r=dsqrt(x**2+y**2)
        ctmp=r*cdexp(i*phi)
        ctmp1=ctmp**(-a)
        ctmp2=ctmp**(a-c)
        x=dreal(1-z)
        y=dimag(1-z)
        if(x.gt.zero .and. y.ge.zero)then
          phi=datan(dabs(y/x))
        elseif(x.lt.zero .and. y.gt.zero)then
          phi=pi-datan(dabs(y/x))
        elseif(x.lt.zero .and. y.le.zero)then
          phi=-pi+datan(dabs(y/x))
        elseif(x.gt.zero .and. y.lt.zero)then
          phi=-datan(dabs(y/x))
        elseif(x.eq.zero .and. y.gt.zero)then
          phi=pi/2
        elseif(x.eq.zero .and. y.lt.zero)then
          phi=-pi/2
        else
          write(*,*)'error in transformation 5'
          return
        endif
        r=dsqrt(x**2+y**2)
        ctmp=r*cdexp(i*phi)
        ctmp2=ctmp2*ctmp**(c-a-b)
c       chyp=f1*ctmp1*cgamm(c)*cgamm(c-a-b)/cgamm(c-a)/cgamm(c-b)
c    #      +f2*ctmp2*cgamm(c)*cgamm(a+b-c)/cgamm(a)/cgamm(b)
        chyp=f1*ctmp1*cgamm(c)*cgamm(c-a-b)*ogamm(c-a)*ogamm(c-b)
     #      +f2*ctmp2*cgamm(c)*cgamm(a+b-c)*ogamm(a)*ogamm(b)
c       write(*,'(2(2 e13.6,2x))')ctmp1,ctmp1-z**(-a)
c       write(*,'(2(2 e13.6,2x))')ctmp2,ctmp2-z**(a-c)*(1-z)**(c-a-b)

      elseif (flag .eq. 6)then

        w=1/z
        f1=csum(a,a-c+1,a-b+1,w)
        f2=csum(b-c+1,b,b-a+1,w)
        x=dreal(-z)
        y=dimag(-z)
        if(x.gt.zero .and. y.ge.zero)then
          phi=datan(dabs(y/x))
        elseif(x.lt.zero .and. y.gt.zero)then
          phi=pi-datan(dabs(y/x))
        elseif(x.lt.zero .and. y.le.zero)then
          phi=-pi+datan(dabs(y/x))
        elseif(x.gt.zero .and. y.lt.zero)then
          phi=-datan(dabs(y/x))
        elseif(x.eq.zero .and. y.gt.zero)then
          phi=pi/2
        elseif(x.eq.zero .and. y.lt.zero)then
          phi=-pi/2
        else
          write(*,*)'error in transformation 6'
          return
        endif
        r=dsqrt(x**2+y**2)
        ctmp=r*cdexp(i*phi)
        ctmp1=ctmp**(-a)
        ctmp2=ctmp**(-b)
c       chyp=f1*ctmp1*cgamm(c)*cgamm(b-a)/cgamm(b)/cgamm(c-a)
c    #      +f2*ctmp2*cgamm(c)*cgamm(a-b)/cgamm(a)/cgamm(c-b)
        chyp=f1*ctmp1*cgamm(c)*cgamm(b-a)*ogamm(b)*ogamm(c-a)
     #      +f2*ctmp2*cgamm(c)*cgamm(a-b)*ogamm(a)*ogamm(c-b)
c       write(*,'(2(2 e13.6,2x))')ctmp1,ctmp1-(-z)**(-a)
c       write(*,'(2(2 e13.6,2x))')ctmp2,ctmp2-(-z)**(-b)

      endif

c     write(20,'(i2,2x,2(e13.6,2x))')flag,cdabs(z),cdabs(w)

      return
      end

c*******************************************************************
c
c234567
      complex*16 function csum(a,b,c,z)
      implicit none
      integer n
      real*8 test,eps
      complex*16 a,b,c,z,num,den,tmp,sum

      eps=1.d-14

      tmp=dcmplx(1.d0)
      sum=dcmplx(0.d0)
      do 10 n=1,1000
        sum=sum+tmp
        test=cdabs(tmp/sum)
        if(test.gt.eps)then
          num=(a+n-1)*(b+n-1)
          den=(c+n-1)*n
          tmp=tmp*num*z/den
        else
          go to 20
        endif
  10  continue
      write(*,*)'warning in csum: sum not completely converged,
     #|z|=',cdabs(z)
  20  continue

      csum=sum

      return
      end

c*******************************************************************
c
c  function name      - ogamm
c
c  computation
c  performed          - calculates one over the gamma function
c
c  usage              - ogamm(z)
c
c  argument         z - any complex number
c
c  precision          - double
c
c  language           - fortran
c
c*******************************************************************

      complex*16 function ogamm(z)                                       
      implicit none
      real*8 one
      parameter (one=1.d0)
      complex*16 z,ztmp,coeff,cgamm

            ztmp=z
            coeff=dcmplx(one)
  1         if (dreal(ztmp).lt.one) then
              coeff=coeff*ztmp
              ztmp=ztmp+one
              go to 1
            endif

            ogamm=coeff/cgamm(ztmp)

            return
            end

c*******************************************************************
c                                                                   
c  function name      - cgamm                                         
c                                                                      
c  computation                                                        
c  performed          - calculates the complex gamma function   
c                                                                     
c  usage              - cgamm(z)                                         
c                                                                     
c  argument         z - any complex number 
c                                                                     
c  precision          - double                                        
c                                                                    
c  language           - fortran                                      
c                                                                    
c*******************************************************************

      complex*16 function cgamm(z)
      implicit none
      integer k
      real*8  c(0:15),pi,one
      parameter (one=1.d0)
      complex*16  z,ztmp,sum,tmp,const

      pi=acos(-1.d0)
      call coeff(c,15)

      ztmp=z
      const=dcmplx(one)
  1   if (dreal(ztmp).lt.one) then
        const=const*ztmp
        ztmp=ztmp+one
        go to 1
      endif

      sum=cmplx(0.d0)
      tmp=cmplx(1.d0)
      do 10 k=0,15
        sum=sum+c(k)*tmp
        tmp=tmp*(ztmp-k-1)/(ztmp+k)
   10 continue

      cgamm=sum*dsqrt(2*pi)*(ztmp+9.d0/2)**(ztmp-.5d0)
     #         *exp(-ztmp-9.d0/2)/const

      return
      end
                                                                   
c*******************************************************************
c                                                                  
c  subroutine name    - coeff                                      
c                                                                  
c  computation                                                     
c  performed          - tabulates the gamma function coefficients  
c                       (see luke vol.2 p.304)
c                                                                  
c  usage              - call coeff(c,n)                            
c                                                                 
c  arguments       c  - the array (output) which contains the     
c                       coefficients.                       
c                                                                 
c                  n  - the dimension (input) of the array 'c'.   
c                                                                
c  precision          - double                                
c
c  language           - fortran 77                                
c                                                                 
c*******************************************************************
                                                                  
      subroutine coeff(c,n)                                         
                                                                  
      real*8  c(0:n)                                              
                                                                   
      c(0) = 41.624436916439068d0                
      c(1) =-51.224241022374774d0                 
      c(2) = 11.338755813488977d0                
      c(3) = -0.747732687772388d0               
      c(4) =  0.008782877493061d0               
      c(5) = -0.000001899030264d0                
      c(6) =  0.000000001946335d0                
      c(7) = -0.000000000199345d0                 
      c(8) =  0.000000000008433d0                
      c(9) =  0.000000000001486d0               
      c(10)= -0.000000000000806d0               
      c(11)=  0.000000000000293d0                
      c(12)= -0.000000000000102d0              
      c(13)=  0.000000000000037d0              
      c(14)= -0.000000000000014d0               
      c(15)=  0.000000000000006d0                 

      return                                                      
      end                                                         
