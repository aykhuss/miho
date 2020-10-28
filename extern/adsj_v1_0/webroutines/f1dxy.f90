!----------------------------------------------------------------------
!	function df1dxy										20/2/97
!	______________
!
!	by:		f. colavecchia
!		
!			calculates the derivative of order (m,n) of
!			the hypergeometric function 
!			f1(a,b,b',c,x,y) 
!
!	input:	
!			complex*16 a,b,bp,c,x,y
!			integer*4 m,n
!			
!	output: complex*16 f1dxy
!
!	limitations: 
!
!
!
!--------------------------------------------------------------------
!
complex*16 function f1dxy(ca,cb,cbp,cc,cx,cy,m,n)

implicit none

integer*4 m,n
complex*16 ca,cb,cbp,cc,cx,cy
complex*16 capoch,cbpoch,cbppoch,ccpoch,caux
complex*16 f1,pochhammer

!
!	Calculating pochhammer symbols
!
capoch = Pochhammer(ca,n+m)
cbpoch = Pochhammer(cb,m)
cbppoch= Pochhammer(cbp,n)
ccpoch = Pochhammer(cc,n+m)

caux = f1(ca+m+n,cb+m,cbp+n,cc+m+n,dreal(cx),dreal(cy))
f1dxy = capoch*cbpoch*cbppoch/ccpoch*caux

return

end

