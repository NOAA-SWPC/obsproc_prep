!
!                                 ******************************************
!                                 *             MODULE pspl                *
!                                 *  R. J. Purser, NOAA/NCEP/EMC, May 2014 * 
!                                 *  jim.purser@noaa.gov                   *
!                                 *                                        *
!                                 ******************************************
!
! A collection of handy spline routines
!
! Last modified:
!  Keyser (2014-12-12) - print written to unit 41 rather than stdout (for use in
!                        prepobs_prepacqc program - limits amount of stdout)
! 
! DIRECT DEPENDENCIES:
! Modules:   pietc, pkind, pmat2
! Libraries: pmat
!
!
!=============================================================================
module pspl
!=============================================================================
use pkind, only: dp 
use pietc, only: T,F,u0,o2,u1 ! True, False, 0., .5,  1.
implicit none
private
public:: expm,expmm,coshm,sinhm,coshmm,xcms,                bnewton, &
     fit_tspline,eval_tspline,int_tspline,eval_itspline,             &
     fit_uspline,eval_uspline,int_uspline,eval_iuspline,             &
     best_slalom,count_gates,set_gates,set_posts,                    &
     count_routes,list_routes,next_route,slalom_spline,convertd,convertd_back

integer, parameter:: ihu=1025      ! "Huge" integer parameter
real(dp),parameter:: hu=huge(hu)/2 ! "Huge" real parameter
real(dp),parameter:: eps=epsilon(eps),heps=.01,geps=1.e-10 ! Small parameters

interface expm;          module procedure expm;                  end interface
interface expmm;         module procedure expmm;                 end interface
interface coshm;         module procedure coshm;                 end interface
interface sinhm;         module procedure sinhm;                 end interface
interface coshmm;        module procedure coshmm;                end interface
interface xcms;          module procedure xcms;                  end interface
interface bnewton;       module procedure tbnewton,ubnewton;     end interface
interface fit_tspline;   module procedure fit_gtspline,fit_tspline
                                                                 end interface
interface eval_tspline
   module procedure eval_tspline,eval_tsplined,eval_tsplinedd,eval_tsplineddd
                                                                 end interface
interface int_tspline;   module procedure int_tspline;           end interface
interface eval_itspline; module procedure eval_itspline;         end interface
interface fit_uspline;   module procedure fit_guspline,fit_uspline
                                                                 end interface
interface eval_uspline
   module procedure eval_uspline,eval_usplined,eval_usplinedd,eval_usplineddd
                                                                 end interface
interface int_uspline;   module procedure int_uspline;           end interface
interface eval_iuspline; module procedure eval_iuspline;         end interface
interface best_slalom;   module procedure best_tslalom,best_uslalom
                                                                 end interface
interface count_gates;   module procedure count_gates;           end interface
interface set_gates;     module procedure set_gates;             end interface
interface set_posts;     module procedure set_posts;             end interface
interface count_routes;  module procedure count_routes;          end interface
interface list_routes;   module procedure list_routes;           end interface
interface next_route;    module procedure next_route;            end interface
interface slalom_spline; module procedure slalom_tspline,slalom_uspline
                                                                 end interface
interface convertd;      module procedure convertd;              end interface
interface convertd_back; module procedure convertd_back;         end interface

contains

!=============================================================================
function expm(x) result(e)!                                             [expm]
!=============================================================================
! exp(x)-1 (approximately x for small x)
! = I^(1)exp(x), where I^(p) is the integral iterated p times
real(dp),intent(in ):: x
real(dp)            :: e
!-----------------------------------------------------------------------------
real(dp):: p
integer :: i
!=============================================================================
if(abs(x)>o2)then
   e=exp(x)-u1
else
   p=x; e=p
   do i=2,19; p=p*x/i; e=e+p; if(abs(p)<=abs(e*eps))return; enddo
endif
end function expm

!=============================================================================
function expmm(x) result(e)!                                           [expmm]
!=============================================================================
! exp(x)-1-x (approximately x^2/2 for small x)
! = I^(2)exp(x), where I^(p) is the integral iterated p times
real(dp),intent(in ):: x
real(dp)            :: e
!-----------------------------------------------------------------------------
real(dp):: p
integer :: i
!=============================================================================
if(abs(x)>o2)then
   e=exp(x)-u1-x
else
   p=x*x*o2; e=p
   do i=3,25; p=p*x/i; e=e+p; if(abs(p)<=abs(e*eps))return; enddo
endif
end function expmm

!=============================================================================
function coshm(x) result(c)!                                           [coshm]
!=============================================================================
! cosh(x)-1  (approximately x**2/2 for small x)
! =I^(2)cosh(x), where I^(p) is the integral iterated p times
real(dp),intent(in ):: x
real(dp)            :: c
!-----------------------------------------------------------------------------
c=2*sinh(x*o2)**2
end function coshm

!=============================================================================
function sinhm(x) result(s)!                                           [sinhm]
!=============================================================================
! sinh(x)-x (approximately x**3/6 for small x)
! =I^(3)cosh(x), where I^(p) is the integral iterated p times
real(dp),intent(in ):: x
real(dp)            :: s
!-----------------------------------------------------------------------------
real(dp):: p,xx
integer :: i
!=============================================================================
if(abs(x)>o2)then
   s=sinh(x)-x
else
   p=x**3/6;  s=p;  xx=x*x
   do i=5,19,2; p=p*xx/(i*(i-1)); s=s+p; if(abs(p)<=abs(s*eps))return; enddo
endif
end function sinhm

!=============================================================================
function coshmm(x) result(c)!                                         [coshmm]
!=============================================================================
! cosh(x)-1-x^2/2  (approximately x**4/24 for small x)
! =I^(4)cosh(x), where I^(p) is the integral iterated p times
real(dp),intent(in ):: x
real(dp)            :: c
!-----------------------------------------------------------------------------
real(dp)            :: xh
!=============================================================================
xh=x*o2
c=sinhm(xh)*(2*sinh(xh)+x)
end function coshmm

!=============================================================================
function xcms(x) result(e)!                                             [xcms]
!=============================================================================
real(dp),intent(in ):: x
real(dp)            :: e
!-----------------------------------------------------------------------------
real(dp):: p,xx
integer :: i,i2
!=============================================================================
! x*coshm(x)-sinhm(x) (approximately x**3/3 for small x)
if(abs(x)>o2)then
   e=x*coshm(x)-sinhm(x)
else
   p=x**3/3;  e=p;  xx=x*x
   do i=2,15
      i2=i*2; p=p*xx/(i2*(i2+1)); e=e+i*p; if(abs(p)<=abs(e*eps))return
   enddo
endif
end function xcms

!=============================================================================
subroutine tbnewton(nh,m,bigT,halfgate,ts,hs,tp,p,q, te,dhdt,FF)!    [bnewton]
!=============================================================================
! Perform a "bounded Newton" iteration to estimate the vertical velocity,
! dh/dt, as the trajectory passes through the nh "gates", each at height
! hs(i) and centered at time ts(i) with the halfwidth of the gate equal to
! halfgate. The characteristic timescale of the fitted spline is bigT.
! the time-nodes rescaled by this bigT are at the m points tp, and the 
! corresponding "p" and "q" coefficients of the tensioned spline are p and q
! respectively. 
! The output is an array of n estimates, dh/dt, at each gate's hs. This
! involves estimating the time t of passage through each gate, and is done
! generically by Newton's method, except that the newton increments are
! bounded to be within a range "gate" = 2*halfgate, so that we eliminate
! wild excursions when dh/dt is actually very small (or vanishes). When such
! an excusion is detected, the returned value of dh/dt is simply that
! evaluated at ts(i), and no further attempt at Newton refinement is made
! at this i. (The vertical motion in such cases is essentially negligible
! in any case, and very likely is multivalued as the trajectory wavers about
! this gate's value of hs(i).)
!=============================================================================
integer,               intent(in ):: nh,m
real(dp),              intent(in ):: bigT,halfgate
real(dp),dimension(nh),intent(in ):: ts,hs
real(dp),dimension(m) ,intent(in ):: tp,p,q
real(dp),dimension(nh),intent(out):: dhdt, te
logical,               intent(out):: FF
!-----------------------------------------------------------------------------
integer,parameter    :: nit=12
real(dp),dimension(m):: tr
real(dp)             :: gate,tee,he,hac,dhadt,dh,dt
integer              :: i,it
!=============================================================================
gate=2*halfgate/bigT
tr=tp/bigT
do i=1,nh
   tee=ts(i)/bigT
   he=hs(i)
!  Use Newton iterations to estimate the rescaled time, tee, at which the 
!  height is he
   it = 1
   do while (it <= nit)
      call eval_tspline(m,tr,p,q, tee,hac,dhadt)
      if(it==1)dhdt(i)=dhadt/bigT
      if(dhadt==u0)exit
      dh=hac-he
      dt=-dh/dhadt
      if(abs(dt)>gate)then
         write(41,*) 'WARNING! In tbnewton; i,it,dt/gate = ',i,it,dt/gate
         exit
      endif
      if(abs(dh)<heps)then
         dhdt(i)=dhadt/bigT
         exit
      endif
      tee=tee+dt
      it = it + 1
   enddo
   FF=(it>nit)
   if(FF)then
      write(41,'("In tbnewton; Newton iterations seem not to be")')
      write(41,'("converging at i=",i3)'),i
      write(41,'("tee,he,hac,heps,dhadt:",5(1x,e11.4))'),tee,he,hac,heps,dhadt
   endif
   te(i) = tee
enddo
end subroutine tbnewton

!=============================================================================
subroutine ubnewton(nh,m,halfgate,ts,hs,tp,p,q, te,dhdt,FF)!         [bnewton]
!=============================================================================
! Like tbnewton, but for the case of untensioned (i.e., cubic) splines
!=============================================================================
integer,               intent(in ):: nh,m
real(dp),              intent(in ):: halfgate
real(dp),dimension(nh),intent(in ):: ts,hs
real(dp),dimension(m) ,intent(in ):: tp,p,q
real(dp),dimension(nh),intent(out):: dhdt, te
logical,               intent(out):: FF
!-----------------------------------------------------------------------------
integer,parameter    :: nit=12
real(dp),dimension(m):: tr
real(dp)             :: gate,tee,he,hac,dhadt,dh,dt
integer              :: i,it
!=============================================================================
gate=2*halfgate
tr=tp
do i=1,nh
   tee=ts(i)
   he=hs(i)
!  Use Newton iterations to estimate the rescaled time, tee, at which the 
!  height is he
   it = 1
   do while (it <= nit)
      call eval_uspline(m,tr,p,q, tee,hac,dhadt)
      if(it==1)dhdt(i)=dhadt
      if(dhadt==u0)exit
      dh=hac-he
      dt=-dh/dhadt
      if(abs(dt)>gate)then
         write(41,*) 'WARNING! In ubnewton; i,it,dt/gate = ',i,it,dt/gate
         exit
      endif
      if(abs(dh)<heps)then
         dhdt(i)=dhadt
         exit
      endif
      tee=tee+dt
      it = it + 1
   enddo
   FF=(it>nit)
   if(FF)then
      write(41,'("In ubnewton; Newton iterations seem not to be")')
      write(41,'("converging at i=",i3)'),i
      write(41,'("tee,he,hac,heps,dhadt:",5(1x,e11.4))'),tee,he,hac,heps,dhadt
   endif
   te(i) = tee
enddo
end subroutine ubnewton

!============================================================================
subroutine fit_gtspline(n,xs,ys,on,q,j,yac,en,FF)!              [fit_tspline]
!============================================================================
! Fit the gappy tensioned spline, where only those nodes flagged "on"
! are effective in the fitting procedure. Owing to the fact that, where
! constraints are not "on" the spline will generally depart from ys, the
! actual y (yac) is returned for all nodes, regardless of the partial
! duplication with the given ys. In other respects, this is just
! like fit_tspline.
!============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,ys
logical, dimension(n),intent(in ):: on
real(dp),dimension(n),intent(out):: q,j,yac
real(dp),             intent(out):: en
logical,              intent(out):: FF
!----------------------------------------------------------------------------
real(dp),dimension(n):: xa,ya,qa,ja
integer              :: i,k,m
!============================================================================
m=0
do i=1,n
   if(on(i))then; m=m+1; xa(m)=xs(i); ya(m)=ys(i); endif
enddo
call fit_tspline(m,xa(1:m),ya(1:m),qa(1:m),ja(1:m),en,FF)
if(FF)then
   write(41,*) 'In fit_gtspline; failure flag raised at call to fit_tspline'
   return
endif
k=0
do i=1,n
   if(on(i))then
      k=k+1
      q(i)=qa(k)
      j(i)=ja(k)
      yac(i)=ys(i)
   else
      call eval_tsplined(m,xa(1:m),ya(1:m),qa(1:m),xs(i), yac(i),q(i))
      j(i)=0
   endif
enddo
end subroutine fit_gtspline

!============================================================================
subroutine fit_tspline(n,xs,p,q,j,en,FF)!                       [fit_tspline]
!============================================================================
! Solve for the coefficients, the 3rd-derivative jumps, and the energy,
! of the standardized tensioned spline passing through the n nodes at (xs,p).
!
! The value and successive derivatives on the immediate positive side of
! each node, xs(i), are to be found as p(i), q(i), r(i), s(i), with j(i)
! being the discontinuity of 3rd-derivative s between the negative and positive
! side of the node (value itself, and other derivatives, remaining continuous).
! In addition, p(0), q(0), r(0) and s(0) are the value and derivatives on the 
! immediate negative side of xs(1). The spline solution minimizes elastic
! and tensional energy, en, defined as the integral dx of half the sum of the
! squared first and second derivatives over the whole line. Euler-Lagrange
! implies the solution is expressible in each segment between or beyond nodes:
!       y(x') = p + q*x' + r*coshm(x') + s*sinhm(x')
! where x' = x-xs(i) is the local coordinate relative to the relevant node
! (the node at the left of the segent except that, implicitly, we take
! xs(0)===xs(1), and the two functions, coshm and sinhm, are defined:
! coshm(x) = cosh(x)-1
! sinhm(x) = sinh(x)-x.
! The solution in segment 0, i.e., x< xs(1), must exponentially decay towards
! a constant as x--> -infinity, while that for segment n must likewise decay
! as x--> +infinity, in order that energy remains finite. Thus, q(0)=r(0)=s(0)
! and q(n)=-r(n)=s(n) always. Solutions in these infinite end segments are
! therefore expressible in terms only of p(0),q(0) for segment 0 and in terms
! only of p(n), q(n) for segment n and is linear in these coefficients.
! Between consecutive nodes (segments 0<i<n) the solution y(x') is expressible
! in terms only of p(i),q(i),p(i+1),q(i+1) and is linear in these coefficients.
! Being a quadratic functional, the total energy is therefore expressible
! in vector-matrix quadratic form:
! En = (1/2)*p^T*PP*p + q^T*QP*p + (1/2)*q^T^QQ*q
! where p and q are here taken as the column vectors of all their discrete
! values, i=0,..n, and PP, QP and QQ are certain tri-diagonal matrices. 
! The spline solution is found as the energy-minimizing solution to:
! QP*p + QQ*q = 0
! where the p are known (the y(i)) and q is the vector of unknowns. Having 
! solved for q, we can immediately deduce the r and s, and hence the jump, j.
! Finally, we can also return the value of the energy as well.
! The use of a tri-diagonal solver, though seemingly more complicated than
! naive "shooting" alternatives, is very much better conditioned numerically,
! and will succeed in very long data series where the naive methods invariably
! fail.
! In practice, owing to the formal symmetries of the energy in each interval,
! we need only consider the change, p(i+1)-p(i), in each interval, and
! this "odd-symmetry" part of vector, p, only couples with the corresponding
! symmetry in the q (which is the part, (q(i)+q(i+1))), so two of the
! tridiagonals actually reduce to diagonals, simplifying the algebra.
!=============================================================================
use pmat2, only: ldltb, ltdlbv
integer,                intent(in ):: n
real(dp),dimension(  n),intent(in ):: xs,p
real(dp),dimension(  n),intent(out):: q,j
real(dp),               intent(out):: en
logical,                intent(out):: FF
!----------------------------------------------------------------------------
integer                   :: i,ip
real(dp)                  :: x,ch,sh,sa,sb,sap,ccc,xcmsx2,egg,ehh
real(dp),dimension(n-1)   :: difp,sumq,cpp,cqp
real(dp),dimension(n,-1:0):: qq ! <- Tridiagonal, stored as rows of nonupper
!=============================================================================
FF=F
if(n<1)stop 'In fit_tspline; size of data array must be positive'
if(n==1)then; q=0; j=0; en=0; return; endif
! apply a strict monotonicity check on the xs:
do i=2,n
   if(xs(i-1)>=xs(i)) then
      FF=T
      write(41,*) 'In fit_tspline; xs data must increase strictly monotonically'
      return
   end if
enddo
! Initialize tri-diagonal kernels for the energy definition:
qq=0 ! <- initialize symmetric tridiagonal, kernel for q^T*QQ*q
     ! where "q" are the dp/dx at each node.
! The coefficients in the quadratic form defining the spline energy also
! include terms involving factors (p(ip)-p(i))*(q(i)+q(ip)) and
! (p(ip)-p(i))*(p(ip)-p(i)), but these can be dealt with using, respectively,
! the  matrices cqp and cpp which are simply diagonal. It is the symmetries
! in the defiition of energy that allow this simplification.

! Loop over the intervals bounded by consecutive nodes:
do i=1,n-1
   ip=i+1
   difp(i)=p(ip)-p(i)
   x=(xs(ip)-xs(i))*o2 ! Halfwidth of interval
   ch=cosh(x);  sh=sinh(x)
   xcmsx2=xcms(x)*2
! egg relates to the odd-g-basis function's energy integral coefficient
! ehh relates to the even-g-basis function's energy integral coefficient
   egg=x*sh/xcmsx2; ehh=ch/(2*sh)
! ccc is the coefficient of energy integral coupling g(i)*g(i) and g(ip)*g(ip)
   ccc=egg+ehh
   cpp(i)=ch/xcmsx2          ! Energy coefficient for difp(i)*difp(i)...
   cqp(i)=-difp(i)*sh/xcmsx2 ! ..and for difp(i)*sumq(i)
   qq(i,0)=qq(i,0)+ccc; qq(ip,-1)=qq(ip,-1)+egg-ehh; qq(ip,0)=qq(ip,0)+ccc
enddo
! Add the exterior energy contributions to qq at both ends:
qq(1,0)=qq(1,0)+1
qq(n,0)=qq(n,0)+1

! Temporarily, q is made the vector of forcings in the tridiagonal linear
! system from which the final spline coefficients, q, are solved in place.
q(1:n-1)=-cqp; q(n)=0
q(2:n)=q(2:n)-cqp

! The following 2 lines solve the tridiagonal system for q:
call ldltb(n,1,qq)    ! <- Decompose qq into factors, L*(1/D)*L^T, L=lower
call ltdlbv(n,1,qq,q) ! <- Back-substitute, thus solving for q
sumq=q(1:n-1)+q(2:n)   ! <-pairwise sums of derivatives, q:

! The minimizing energy can now be evaluated as a sum of only 2 terms:
en=o2*(dot_product(difp**2,cpp)+dot_product(sumq,cqp))

! Finally, evaluate the 3rd-derivative "jumps", j, at each node:
! Here, sb is the 3rd-derivative at the right end, sa that at the left end,
! of whichever interval is under consideration, but for interior intervals,
! sa = sap+q(i) and sb=sap+q(i+1).
sb=q(1)
do i=1,n-1
   ip=i+1
   x=o2*(xs(ip)-xs(i))
   xcmsx2=xcms(x)*2
   ch=cosh(x);  sh=sinh(x)
   sap=(sh*sumq(i)-ch*difp(i))/xcmsx2
   sa=sap+q(i)
   j(i)=sa-sb
   sb  =sap+q(ip)
enddo
j(n)=q(n)-sb ! Final "sa" is just q(n) for the right exterior
end subroutine fit_tspline
  
!=============================================================================
subroutine int_tspline(n,xs,p,q, m)!                             [int_tspline]
!=============================================================================
! Take the sets of n parameters p and q of the tensioned spline
! and return the values of its integral at the n-1 interval midpoints, and
! the value at the last node, assuming that the integral at the first node
! is set to zero.
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),dimension(n),intent(out):: m
!-----------------------------------------------------------------------------
real(dp):: a,b,c,d,e,t2,x,pa,pd,qa,qd,shx,chmx,shmx,chmmx,xcmsx
integer :: i,ip
!=============================================================================
! e is the running integral as we loop over successive nodes, so it starts out
! zero at the first node:
e=u0
! Loop over intervals:
do i=1,n-1
   ip=i+1
   x=(xs(ip)-xs(i))*o2 !<- interval half-width
   t2=x*x*o2
   shx  =sinh  (x)
   chmx =coshm (x)
   shmx =sinhm (x)
   chmmx=coshmm(x)
   xcmsx=xcms  (x)
   pa=(p(ip)+p(i))*o2
   pd=(p(ip)-p(i))*o2/x
   qa=(q(ip)+q(i))*o2
   qd=(q(ip)-q(i))*o2/shx
! a,b,c,d are analogous to the Taylor coefficients of a cubic about the 
! interval midpoint, but more precisely, c and d relate to basis functions
! coshm and sinhm (instead of x**2/2 and x**3/6 for the perfect cubic).
   c=qd
   a=pa-c*chmx
   d=(qa-pd)*x/xcmsx
   b=qa-d*chmx
   m(i)=e+a*x -b*t2 +c*shmx -d*chmmx
   e=e+2*(a*x+c*shmx)
enddo
m(n)=e
end subroutine int_tspline

!============================================================================
subroutine fit_guspline(n,xs,ys,on,q,j,yac,en,FF)!              [fit_uspline]
!============================================================================
! Fit the gappy untensioned spline, where only those nodes flagged "on"
! are effective in the fitting procedure. Owing to the fact that, where
! constraints are not "on" the spline will generally depart from ys, the
! actual y (yac) is returned for all nodes, regardless of the partial
! duplication with the given ys. In other respects, this is just
! like fit_tspline.
!============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,ys
logical, dimension(n),intent(in ):: on
real(dp),dimension(n),intent(out):: q,j,yac
real(dp),             intent(out):: en
logical,              intent(out):: FF
!----------------------------------------------------------------------------
real(dp),dimension(n):: xa,ya,qa,ja
integer              :: i,k,m
!============================================================================
m=0
do i=1,n
   if(on(i))then; m=m+1; xa(m)=xs(i); ya(m)=ys(i); endif
enddo
call fit_uspline(m,xa(1:m),ya(1:m),qa(1:m),ja(1:m),en,FF)
if(FF)then
   write(41,*) 'In fit_guspline; failure flag raised at call to fit_uspline'
   return
endif
k=0
do i=1,n
   if(on(i))then
      k=k+1
      q(i)=qa(k)
      j(i)=ja(k)
      yac(i)=ys(i)
   else
      call eval_usplined(m,xa(1:m),ya(1:m),qa(1:m),xs(i), yac(i),q(i))
      j(i)=0
   endif
enddo
end subroutine fit_guspline

!============================================================================
subroutine fit_uspline(n,xs,p,q,j,en,FF)!                        [fit_uspline]
!============================================================================
! Solve for the coefficients, the 3rd-derivative jumps, and the energy,
! of the untensioned (cubic) spline passing through the n nodes at (xs,p).
!
! The algorithm follows the pattern given in fit_tspline, except that the
! hyperbolic functions are all replaced by their asymptotic (x--> 0) limiting
! forms. These limiting forms are as follows:
! cosh(x) --> 1
! sinh(x) --> x
! coshm(x) --> x**2/2
! sinhm(x) --> x**3/6
! xcms(x)  --> x**3/3
!=============================================================================
use pietc, only: o3
use pmat2, only: ldltb, ltdlbv
integer,                intent(in ):: n
real(dp),dimension(  n),intent(in ):: xs,p
real(dp),dimension(  n),intent(out):: q,j
real(dp),               intent(out):: en
logical,                intent(out):: FF
!----------------------------------------------------------------------------
integer                   :: i,ip
real(dp)                  :: x,x2,sa,sb,ccc,xcmsx2
real(dp),dimension(n-1)   :: difp,sumq,cpp,cqp
real(dp),dimension(n,-1:0):: qq ! <- Tridiagonal, stored as rows of nonupper
!=============================================================================
FF=F
if(n<1)stop 'In fit_uspline; size of data array must be positive'
if(n==1)then; q=0; j=0; en=0; return; endif
! apply a strict monotonicity check on the xs:
do i=2,n
   if(xs(i-1)>=xs(i)) then
      FF=T
      write(41,*) 'In fit_uspline; xs data must increase strictly monotonically'
      return
   end if
enddo
! Initialize tri-diagonal kernels for the energy definition:
qq=0 ! <- initialize symmetric tridiagonal, kernel for q^T*QQ*q
     ! where "q" are the dp/dx at each node.
! The coefficients in the quadratic form defining the spline energy also
! include terms involving factors (p(ip)-p(i))*(q(i)+q(ip)) and
! (p(ip)-p(i))*(p(ip)-p(i)), but these can be dealt with using, respectively,
! the  matrices cqp and cpp which are simply diagonal. It is the symmetries
! in the defiition of energy that allow this simplification.

! Loop over the intervals bounded by consecutive nodes:
do i=1,n-1
   ip=i+1
   difp(i)=p(ip)-p(i)
   x2=xs(ip)-xs(i); x=o2*x2! Width, and halfwidth of interval
   xcmsx2=o3*x**3*2

! ccc is the coefficient of energy integral coupling g(i)*g(i) and g(ip)*g(ip)
   ccc=2/x
   cpp(i)=u1/xcmsx2          ! Energy coefficient for difp(i)*difp(i)...
   cqp(i)=-difp(i)*x/xcmsx2 ! ..and for difp(i)*sumq(i)
   qq(i,0)=qq(i,0)+ccc; qq(ip,-1)=qq(ip,-1)+1/x; qq(ip,0)=qq(ip,0)+ccc
enddo
! There is NO exterior energy contributions to qq at both ends:

! Temporarily, q is made the vector of forcings in the tridiagonal linear
! system from which the final spline coefficients, q, are solved in place.
q(1:n-1)=-cqp; q(n)=0
q(2:n)=q(2:n)-cqp

! The following 2 lines solve the tridiagonal system for q:
call ldltb(n,1,qq)    ! <- Decompose qq into factors, L*(1/D)*L^T, L=lower
call ltdlbv(n,1,qq,q) ! <- Back-substitute, thus solving for q
sumq=q(1:n-1)+q(2:n)  ! <-pairwise sums of derivatives, q:

! The minimizing energy can now be evaluated as a sum of only 2 terms:
en=o2*(dot_product(difp**2,cpp)+dot_product(sumq,cqp))

! Finally, evaluate the 3rd-derivative "jumps", j, at each node:
! Here, sb and sa are the 3rd-derivatives in consecutive intervals
sb=0
do i=1,n-1
   ip=i+1
   x=o2*(xs(ip)-xs(i))
   xcmsx2=o3*x**3*2
   sa=(x*sumq(i)-difp(i))/xcmsx2
   j(i)=sa-sb
   sb  =sa
enddo
j(n)=-sb ! Final "sa" is just 0 for the right exterior
end subroutine fit_uspline
  
!=============================================================================
subroutine int_uspline(n,xs,p,q, m)!                             [int_uspline]
!=============================================================================
! Take the sets of n parameters p and q of the untensioned cubic spline
! and return the values of its integral at the n-1 interval midpoints, and
! the value at the last node, assuming that the integral at the first node
! is set to zero.
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),dimension(n),intent(out):: m
!-----------------------------------------------------------------------------
real(dp),parameter:: u3o2=3*o2
real(dp):: a,b,c,d,e,t2,t3,t4,x,pa,pd,qa,qd
integer :: i,ip
!=============================================================================
! e is the running integral as we loop over successive nodes, so it starts out
! zero at the first node:
e=u0
! Loop over intervals:
do i=1,n-1
   ip=i+1
   x=(xs(ip)-xs(i))*o2 !<- interval half-width
   t2=x*x/2
   t3=t2*x/3
   t4=t3*x/4
   pa=(p(ip)+p(i))*o2
   pd=(p(ip)-p(i))*o2/x
   qa=(q(ip)+q(i))*o2
   qd=(q(ip)-q(i))*o2/x
! a,b,c,d are the Taylor coefficients of the cubic about the interval midpoint:
   c=qd
   a=pa-c*t2
   d=(qa-pd)*u3o2/t2
   b=qa-d*t2
   m(i)=e+a*x-b*t2+c*t3-d*t4
   e=e+2*(a*x+c*t3)
enddo
m(n)=e
end subroutine int_uspline

!=============================================================================
subroutine eval_tspline(n,xs,p,q, x,y)!                         [eval_tspline]
!=============================================================================
! Assuming the 1st derivatives, q, are correctly given at the n nodes, xs,
! of the standardized tensioned spline, where p are the nodal values, 
! evaluate the spline function y at the location x.
! First find the nonvanishing interval in which x resides, then expand
! y using basis functions implied by the interval-end values of p and q
! using the interval midpoint as local origin when x is interior, or the
! single interval endpoint when it is exterior.
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y
!----------------------------------------------------------------------------
integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm
!============================================================================
if(x<=xs(1))then; xr=x-xs(1); y=p(1)+q(1)*expm( xr); return; endif
if(x>=xs(n))then; xr=x-xs(n); y=p(n)-q(n)*expm(-xr); return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle ! <- only consider intervals of positive width
   if(xs(ib)>=x)exit ! <- exit once finite interval straddling x is found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2   ! <- halfwidth of interval
xr=x-xs(ia)-xh          ! <- x relative to interval midpoint
pm=(p(ia)+p(ib))*o2     ! average of end values
qm=(p(ib)-p(ia))/(xh*2) ! average gradient
qah=q(ia)*o2; qbh=q(ib)*o2
qxh=qah+qbh-qm ! Half the total excess q at interval ends
qdh=qbh-qah    ! Half the difference of q at interval ends
shh=sinh(xh);   chh=cosh(xh)
sh =sinh(xr);    ch=cosh(xr)
shm=sinhm(xr);  chm=coshm(xr)
shhm=sinhm(xh); chhm=coshm(xh)
xcmsh=xcms(xh)
qdh=qdh/shh; qxh=qxh/xcmsh ! <- rescale qdh, qxh
y=pm+xr*qm +qdh*(chm-chhm) +  qxh*(xh*shm-xr*shhm)
end subroutine eval_tspline
   
!=============================================================================
subroutine eval_tsplined(n,xs,p,q, x,y,dydx)!                   [eval_tspline]
!=============================================================================
! Like eval_tspline, but also return the derivative dy/dx
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx
!----------------------------------------------------------------------------
integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm,&
           qemxr
!============================================================================
if(x<=xs(1))then
   xr=x-xs(1); qemxr=q(1)*expm( xr); y=p(1)+qemxr; dydx=qemxr+q(1); return
endif
if(x>=xs(n))then
   xr=x-xs(n); qemxr=q(n)*expm(-xr); y=p(n)-qemxr; dydx=qemxr+q(n); return
endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle ! <- only consider intervals of positive width
   if(xs(ib)>=x)exit    ! <- exit once finite interval straddling x is found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2 ! <- halfwidth of interval
xr=x-xs(ia)-xh          ! <- x relative to interval midpoint
pm=(p(ia)+p(ib))*o2   ! average of end values
qm=(p(ib)-p(ia))/(xh*2) ! average gradient
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm ! Half the total excess q at interval ends
qdh=qbh-qah    ! Half the difference of q at interval ends
shh=sinh(xh);   chh=cosh(xh)
sh =sinh(xr);   ch=cosh(xr)
shm=sinhm(xr);  chm=coshm(xr)
shhm=sinhm(xh); chhm=coshm(xh)
xcmsh=xcms(xh)
qdh=qdh/shh; qxh=qxh/xcmsh ! <- rescale qdh, qxh
y=pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx=qm+qdh*sh +qxh*(xh*chm-shhm)
end subroutine eval_tsplined

!=============================================================================
subroutine eval_tsplinedd(n,xs,p,q, x,y,dydx,ddydxx)!           [eval_tspline]
!=============================================================================
! Like eval_tspline, but also return the derivative dy/dx
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx,ddydxx
!----------------------------------------------------------------------------
integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm,&
           qemxr
!============================================================================
if(x<=xs(1))then
   xr=x-xs(1); qemxr=q(1)*expm( xr); y=p(1)+qemxr; dydx=qemxr+q(1)
   ddydxx=dydx; return
endif
if(x>=xs(n))then
   xr=x-xs(n); qemxr=q(n)*expm(-xr); y=p(n)-qemxr; dydx=qemxr+q(n)
   ddydxx=-dydx; return
endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle ! <- only consider intervals of positive width
   if(xs(ib)>=x)exit    ! <- exit once finite interval straddling x is found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2 ! <- halfwidth of interval
xr=x-xs(ia)-xh          ! <- x relative to interval midpoint
pm=(p(ia)+p(ib))*o2   ! average of end values
qm=(p(ib)-p(ia))/(xh*2) ! average gradient
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm ! Half the total excess q at interval ends
qdh=qbh-qah    ! Half the difference of q at interval ends
shh=sinh(xh);   chh=cosh(xh)
sh =sinh(xr);   ch=cosh(xr)
shm=sinhm(xr);  chm=coshm(xr)
shhm=sinhm(xh); chhm=coshm(xh)
xcmsh=xcms(xh)
qdh=qdh/shh; qxh=qxh/xcmsh ! <- rescale qdh, qxh
y=pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx=qm+qdh*sh +qxh*(xh*chm-shhm)
ddydxx=qdh*ch +qxh*xh*sh
end subroutine eval_tsplinedd

!=============================================================================
subroutine eval_tsplineddd(n,xs,p,q, x,y,dydx,ddydxx,dddydxxx)! [eval_tspline]
!=============================================================================
! Like eval_tspline, but also return the derivative dy/dx
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx,ddydxx,dddydxxx
!----------------------------------------------------------------------------
integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm,&
           qemxr
!============================================================================
if(x<=xs(1))then
   xr=x-xs(1); qemxr=q(1)*expm( xr); y=p(1)+qemxr; dydx=qemxr+q(1)
   ddydxx=dydx; dddydxxx=dydx; return
endif
if(x>=xs(n))then
   xr=x-xs(n); qemxr=q(n)*expm(-xr); y=p(n)-qemxr; dydx=qemxr+q(n)
   ddydxx=-dydx; dddydxxx=dydx; return
endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle ! <- only consider intervals of positive width
   if(xs(ib)>=x)exit    ! <- exit once finite interval straddling x is found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2 ! <- halfwidth of interval
xr=x-xs(ia)-xh          ! <- x relative to interval midpoint
pm=(p(ia)+p(ib))*o2   ! average of end values
qm=(p(ib)-p(ia))/(xh*2) ! average gradient
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm ! Half the total excess q at interval ends
qdh=qbh-qah    ! Half the difference of q at interval ends
shh=sinh(xh);   chh=cosh(xh)
sh =sinh(xr);   ch=cosh(xr)
shm=sinhm(xr);  chm=coshm(xr)
shhm=sinhm(xh); chhm=coshm(xh)
xcmsh=xcms(xh)
qdh=qdh/shh; qxh=qxh/xcmsh ! <- rescale qdh, qxh
y       =pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx    =      qm +qdh*sh          +  qxh*(xh*chm-   shhm)
ddydxx  =          qdh*ch          +  qxh* xh*sh
dddydxxx=          qdh*sh          +  qxh* xh*ch
end subroutine eval_tsplineddd

!=============================================================================
subroutine eval_itspline(n,xs, p,q,m,  x,y)!                   [eval_itspline]
!=============================================================================
! Evaluate the integrated tension spline at x, returning the value, y.
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q,m
real(dp),             intent(in ):: x
real(dp),             intent(out):: y
!-----------------------------------------------------------------------------
real(dp):: a,b,c,d,t2,xh,shx,chmx,shmx,chmmx,xcmsx,xr,pa,pd,qa,qd
integer :: ia,ib
!=============================================================================
if(x<=xs(1))then; xr=x-xs(1); y=     p(1)*xr+q(1)*expmm( xr); return; endif
if(x>=xs(n))then; xr=x-xs(n); y=m(n)+p(n)*xr+q(n)*expmm(-xr); return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle ! <- only consider intervals of positive width
   if(xs(ib)>=x)exit ! <- exit once finite interval straddling x is found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2   ! <- halfwidth of interval
shx =sinh (xh)
chmx=coshm(xh)
xcmsx=xcms(xh)
xr=x-xs(ia)-xh          ! <- x relative to interval midpoint
pa=(p(ib)+p(ia))*o2
pd=(p(ib)-p(ia))*o2/xh
qa=(q(ib)+q(ia))*o2
qd=(q(ib)-q(ia))*o2/shx
! a,b,c,d are analogous to the Taylor coefficients about the interval midpoint
c=qd
a=pa-c*chmx
d=(qa-pd)*xh/xcmsx
b=qa-d*chmx

t2=xr**2/2
shmx =sinhm (xr)
chmmx=coshmm(xr)
y=m(ia)+a*xr+b*t2+c*shmx+d*chmmx
end subroutine eval_itspline

!=============================================================================
subroutine eval_uspline(n,xs,p,q, x,y)!                         [eval_uspline]
!=============================================================================
! Assuming the 1st derivatives, q, are correctly given at the n nodes, xs,
! of the standardized untensioned spline, where p are the nodal values, 
! evaluate the (UNtensioned) spline function y at the location x.
! First find the nonvanishing interval in which x resides, then expand
! y using basis functions implied by the interval-end values of p and q
! using the interval midpoint as local origin when x is interior, or the
! single interval endpoint when it is exterior.
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y
!----------------------------------------------------------------------------
integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm
!============================================================================
if(x<=xs(1))then; xr=x-xs(1); y=p(1)+q(1)*xr; return; endif
if(x>=xs(n))then; xr=x-xs(n); y=p(n)+q(n)*xr; return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle ! <- only consider intervals of positive width
   if(xs(ib)>=x)exit ! <- exit once finite interval straddling x is found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2   ! <- halfwidth of interval
xr=x-xs(ia)-xh          ! <- x relative to interval midpoint
pm=(p(ia)+p(ib))*o2     ! average of end values
qm=(p(ib)-p(ia))/(xh*2) ! average gradient
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm ! Half the total excess q at interval ends
qdh=qbh-qah    ! Half the difference of q at interval ends
shh=xh;        chh=u1
sh =xr;        ch =u1
shm =xr**3/6;  chm =xr**2*o2
shhm=xh**3/6;  chhm=xh**2*o2
xcmsh=xh**3/3
qdh=qdh/shh; qxh=qxh/xcmsh ! <- rescale qdh, qxh
y=pm+xr*qm +qdh*(chm-chhm) +  qxh*(xh*shm-xr*shhm)

end subroutine eval_uspline
   
!=============================================================================
subroutine eval_usplined(n,xs,p,q, x,y,dydx)!                   [eval_uspline]
!=============================================================================
! Like eval_uspline, but also return the derivative dy/dx
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx
!----------------------------------------------------------------------------
integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm
!============================================================================
if(x<=xs(1))then; xr=x-xs(1); y=p(1)+q(1)*xr; dydx=q(1); return; endif
if(x>=xs(n))then; xr=x-xs(n); y=p(n)+q(n)*xr; dydx=q(n); return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle ! <- only consider intervals of positive width
   if(xs(ib)>=x)exit ! <- exit once finite interval straddling x is found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2    ! <- halfwidth of interval
xr=x-xs(ia)-xh           ! <- x relative to interval midpoint
pm=(p(ia)+p(ib))*o2      ! average of end values
qm=(p(ib)-p(ia))/(xh*2) ! average gradient
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm ! Half the total excess q at interval ends
qdh=qbh-qah    ! Half the difference of q at interval ends
shh=xh;        chh=u1
sh =xr;        ch =u1
shm =xr**3/6;  chm =xr**2*o2
shhm=xh**3/6;  chhm=xh**2*o2
xcmsh=xh**3/3
qdh=qdh/shh; qxh=qxh/xcmsh ! <- rescale qdh, qxh
y=pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx=qm+qdh*sh +qxh*(xh*chm-shhm)
end subroutine eval_usplined

!=============================================================================
subroutine eval_usplinedd(n,xs,p,q, x,y,dydx,ddydxx)!            [eval_uspline]
!=============================================================================
! Like eval_uspline, but also return the derivative dy/dx
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx,ddydxx
!----------------------------------------------------------------------------
integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm
!============================================================================
if(x<=xs(1))then; xr=x-xs(1); y=p(1)+q(1)*xr; dydx=q(1); return; endif
if(x>=xs(n))then; xr=x-xs(n); y=p(n)+q(n)*xr; dydx=q(n); return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle ! <- only consider intervals of positive width
   if(xs(ib)>=x)exit ! <- exit once finite interval straddling x is found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2    ! <- halfwidth of interval
xr=x-xs(ia)-xh           ! <- x relative to interval midpoint
pm=(p(ia)+p(ib))*o2      ! average of end values
qm=(p(ib)-p(ia))/(xh*2)  ! average gradient
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm ! Half the total excess q at interval ends
qdh=qbh-qah    ! Half the difference of q at interval ends
shh=xh;        chh=u1
sh =xr;        ch =u1
shm =xr**3/6;  chm =xr**2*o2
shhm=xh**3/6;  chhm=xh**2*o2
xcmsh=xh**3/3
qdh=qdh/shh; qxh=qxh/xcmsh ! <- rescale qdh, qxh
y=pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx=qm+qdh*sh +qxh*(xh*chm-shhm)
ddydxx=qdh +qxh*xh*sh
end subroutine eval_usplinedd

!=============================================================================
subroutine eval_usplineddd(n,xs,p,q, x,y,dydx,ddydxx,dddydxxx)! [eval_uspline]
!=============================================================================
! Like eval_uspline, but also return the derivative dy/dx
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx,ddydxx,dddydxxx
!----------------------------------------------------------------------------
integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm
!============================================================================
if(x<=xs(1))then; xr=x-xs(1); y=p(1)+q(1)*xr; dydx=q(1); return; endif
if(x>=xs(n))then; xr=x-xs(n); y=p(n)+q(n)*xr; dydx=q(n); return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle ! <- only consider intervals of positive width
   if(xs(ib)>=x)exit ! <- exit once finite interval straddling x is found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2    ! <- halfwidth of interval
xr=x-xs(ia)-xh           ! <- x relative to interval midpoint
pm=(p(ia)+p(ib))*o2      ! average of end values
qm=(p(ib)-p(ia))/(xh*2)  ! average gradient
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm ! Half the total excess q at interval ends
qdh=qbh-qah    ! Half the difference of q at interval ends
shh=xh;        chh=u1
sh =xr;        ch =u1
shm =xr**3/6;  chm =xr**2*o2
shhm=xh**3/6;  chhm=xh**2*o2
xcmsh=xh**3/3
qdh=qdh/shh; qxh=qxh/xcmsh ! <- rescale qdh, qxh
y=pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx=qm+qdh*sh +qxh*(xh*chm-shhm)
ddydxx=qdh +qxh*xh*sh
dddydxxx=qxh*xh
end subroutine eval_usplineddd

!=============================================================================
subroutine eval_iuspline(n,xs, p,q,m,  x,y)!                   [eval_iuspline]
!=============================================================================
! Evaluate the integrated untensioned spline at x, returning the value, y.
!=============================================================================
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q,m
real(dp),             intent(in ):: x
real(dp),             intent(out):: y
!-----------------------------------------------------------------------------
real(dp),parameter:: u3o2=3*o2
real(dp):: a,b,c,d,t2,t3,t4,xh,xr,pa,pd,qa,qd
integer :: ia,ib
!=============================================================================
if(x<=xs(1))then; xr=x-xs(1); y=p(1)*xr+q(1)*xr**2/2; return; endif
if(x>=xs(n))then; xr=x-xs(n); y=m(n)+p(n)*xr+q(n)*xr**2/2; return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle ! <- only consider intervals of positive width
   if(xs(ib)>=x)exit ! <- exit once finite interval straddling x is found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2   ! <- halfwidth of interval
xr=x-xs(ia)-xh          ! <- x relative to interval midpoint
t2=xh**2/2
t3=t2*xh/3
pa=(p(ib)+p(ia))*o2
pd=(p(ib)-p(ia))*o2/xh
qa=(q(ib)+q(ia))*o2
qd=(q(ib)-q(ia))*o2/xh
! a,b,c,d are the Taylor coefficients of the cubic about the interval midpoint:
c=qd
a=pa-c*t2
d=(qa-pd)*u3o2/t2
b=qa-d*t2
t2=xr**2/2
t3=t2*xr/3
t4=t3*xr/4
y=m(ia)+a*xr+b*t2+c*t3+d*t4
end subroutine eval_iuspline

!==============================================================================
subroutine best_tslalom(nh,mh,doru,ts,hs,halfgate,bigT, tp,hp, & ! [best_slalom]
   qbest,yabest,enbest,modebest,maxita,maxitb,maxit,maxrts,FF)
!==============================================================================
! Run through the different allowed routes between the slalom gates and
! select as the final solution the one whose spline has the smallest "energy".
!==============================================================================
integer,                 intent(in   ):: nh,mh,doru
real(dp),dimension(nh),  intent(in   ):: ts,hs
real(dp),                intent(in   ):: halfgate,bigT
real(dp),dimension(mh*2),intent(  out):: tp,hp
real(dp),dimension(mh*2),intent(  out):: qbest
real(dp),dimension(mh*2),intent(  out):: yabest
real(dp),                intent(  out):: enbest
integer,dimension(mh),   intent(  out):: modebest
integer,                 intent(inout):: maxita,maxitb,maxit,maxrts
logical,                 intent(  out):: FF
!-----------------------------------------------------------------------------
real(dp),dimension(mh*2)  :: q,ya
real(dp),dimension(2,mh)  :: tn
real(dp),dimension(2,2,mh):: hn
real(dp)                  :: en
integer, dimension(mh)    :: code,mode
integer, dimension(mh*2)  :: bend
integer                   :: i,k,m,kbest,route_count,ita,ittot
logical, dimension(mh*2)  :: off
logical                   :: flag
!==============================================================================
m=mh*2
call set_gates(nh,mh,doru,ts,hs,halfgate, tn,hn,code,FF)
if(FF)then
   write(41,*) 'In best_tslalom; failure flag was raised in call to set_gates'
   return
endif
call count_routes(mh,code,route_count,FF)
maxrts=max(maxrts,route_count)
if(FF)then
   write(41,*) 'In best_tslalom; failure flag was raised in call to count_routes'
   return
endif
call list_routes(mh,code)
enbest=hu
flag=T
do k=1,ihu
   call next_route(mh,code,mode,flag)
   if(flag)then; flag=F; exit; endif
   call set_posts(mh,mode,tn,hn,bend,tp,hp,off)
   call slalom_tspline(m,bend,tp,hp,off,bigT, q,ya,en,ita,maxitb,ittot,FF)
   if(FF)then
      write(41,*) &
           'In best_tslalom; failure flag was raised in call to slalom_tspline'
      return
   endif
   maxita=max(maxita,ita)
   maxit =max(maxit,ittot)
   if(en<enbest)then
      modebest=mode
      enbest  =en
      qbest   =q
      yabest  =ya
      kbest   =k
   endif
enddo
end subroutine best_tslalom
!==============================================================================
subroutine best_uslalom(nh,mh,doru,ts,hs,halfgate, tp,hp, & !     [best_slalom]
   qbest,yabest,enbest,modebest,maxita,maxitb,maxit,maxrts,FF)
!==============================================================================
! Like best_tslalom, except this treats the special limiting case where the
! spline tension vanishes
!==============================================================================
integer,                 intent(in   ):: nh,mh,doru
real(dp),dimension(nh),  intent(in   ):: ts,hs
real(dp),                intent(in   ):: halfgate
real(dp),dimension(mh*2),intent(  out):: tp,hp
real(dp),dimension(mh*2),intent(  out):: qbest
real(dp),dimension(mh*2),intent(  out):: yabest
real(dp),                intent(  out):: enbest
integer,dimension(mh),   intent(  out):: modebest
integer,                 intent(inout):: maxita,maxitb,maxit,maxrts
logical,                 intent(  out):: FF
!-----------------------------------------------------------------------------
real(dp),dimension(mh*2)  :: q,ya
real(dp),dimension(2,mh)  :: tn
real(dp),dimension(2,2,mh):: hn
real(dp)                  :: en
integer, dimension(mh)    :: code,mode
integer, dimension(mh*2)  :: bend
integer                   :: i,k,m,kbest,route_count,ita,ittot
logical, dimension(mh*2)  :: off
logical                   :: flag
!==============================================================================
m=mh*2
call set_gates(nh,mh,doru,ts,hs,halfgate, tn,hn,code,FF)
if(FF)then
   write(41,*) 'In best_uslalom; failure flag was raised in call to set_gates'
   return
endif
call count_routes(mh,code,route_count,FF)
maxrts=max(maxrts,route_count)
if(FF)then
   write(41,*) 'In best_uslalom; failure flag was raised in call to count_routes'
   return
endif
call list_routes(mh,code)
enbest=hu
flag=T
do k=1,ihu
   call next_route(mh,code,mode,flag)
   if(flag)then; flag=F; exit; endif
   call set_posts(mh,mode,tn,hn,bend,tp,hp,off)
   call slalom_uspline(m,bend,tp,hp,off, q,ya,en,ita,maxitb,ittot,FF)
   if(FF)then
      write(41,*) &
           'In best_uslalom; failure flag was raised in call to slalom_uspline'
      return
   endif
   maxita=max(maxita,ita)
   maxit =max(maxit,ittot)
   if(en<enbest)then
      modebest=mode
      enbest  =en
      qbest   =q
      yabest  =ya
      kbest   =k
   endif
enddo
end subroutine best_uslalom

!=============================================================================
subroutine count_gates(nh,ts,halfgate,mh)!                       [count_gates]
!=============================================================================
! Count the number of distinct "time gates" that can accommodate all the data
! from the given profile. This gate count is mh.
!=============================================================================
integer,                   intent(in ):: nh
real(dp),dimension(nh),    intent(in ):: ts
real(dp),                  intent(in ):: halfgate
integer,                   intent(out):: mh
!-----------------------------------------------------------------------------
real(dp):: tp
integer :: i
!=============================================================================
tp=-hu ! <- default "time at present"
mh=0
do i=1,nh
   if(ts(i)<=tp+geps)cycle
! A new nominal time of observation:
   mh=mh+1
   tp=ts(i)
enddo
end subroutine count_gates

!=============================================================================
subroutine set_gates(nh,mh,doru,ts,hs,halfgate, tn,hn,code,FF)!    [set_gates]
!=============================================================================
! Be sure to precede this routine by a call to "count_gates" to get a
! consistent tally of the number of time gates, mh.
! Set the locations of the "gateposts" and the routing codes of allowed
! trajectories that thread through them.
! Halfgate is half the (temporal) gate width (in seconds)
! The "inferior" gatepost is at t - halfgate, the "superior" at t + halfgate.
! The aggregated data lead to a tally of gates not exceeding the tally of obs.
! the gatepost times of the aggregated data are put into array tn(:,:) where
! tn(1,:) hold the inferior, and tn(2,:) the superior gatepost times.
! In general, it is not known a priori whether the trajectory will end up
! ascending or descending though a given gate, so both alternatives are
! accounted for, with hn(:,1,:) holding the height corresponding to tn(:,:)
! assuming descending passage; hn(:,2:) likewise assuming ascending passage.
! A running "attitude code" is maintained, atti between gates i-1 and i, and
! the previous attitude code, attim between gates i-2 and i-1, where applicable.
! This code is =2 when the later gate is wholly above the earlier, =1, when
! later is wholly below the earlier, and remains 0 when altitudes of 
! the consecutive gates overlap. When the consecutive attitude codes are
! both =2, we force the mode of passage through gate i-1 to be ascending (route
! code = 8) if its route code has not already been determined by one of the
! overriding contact conditions imposed by temporal contact between gate i-1
! and its predecessor. Likewise, if consecutive attitude codes are both =1,
! we force the mode of passage through gate i-1 to be descending (route code
! 4) unless overridden by a previous temporal contact condition.
! Temporal contact conditions between gates i-1 and i force one of three
! route codes at gate i: =2 when gate i is wholly above gate i-1; =3 when
! gate i is wholly below; =5 when the gate height ranges overlap. 
!
! The purpose of the route code is to specify the possible modes of passage
! (descending, ascending, or either of these alternatives) through gate i
! when the actual mode of passage through the preceeding gate is known. It is
! based on a 2-digit trinary code. The possible modes of passage through
! gate i are enumerated by the "option code" whose values are:
! 0 when passage may be either descending or ascending (indeterminate);
! 1 when passage is definitely descending;
! 2 when passage is definitely ascending.
! The "units" digit of the trinary expansion of the route code gives the
! option code for gate i when passage through gate i-1 is prescribed to be
! DESCENDING; the "threes" digit of the trinary expansion of the route code
! gives the option code for gate i when passage through gate i-1 is 
! prescribed to be ASCENDING. These possibilities for the option code for
! gate i are summarized in the table below.
!
!    Route        ;   DESCENDING at (i-1)   ;   ASCENDING at (i-1)
!   Code(i)       ;   ==> Option code(i)    ;   ==> Option Code(i)
!............................................................................
!      0                     0                        0
!      2                     2                        0
!      3                     0                        1
!      4                     1                        1
!      5                     2                        1
!      8                     2                        2
!............................................................................. 
!
! The first route code in a chain of gates, ie., code(1), is alway set
! to 0, so at the very least, two combinations of routes are always coded
! according as whether we choose to initialize the spline solution with
! descent through gate 1 or an ascent. If all the gates are temporally 
! separated, then then final gate's route_code also has this 0 value
! signifying an indeterminate mode of passage.
! 
! In the special case where mh=1 and the given hs data are not enough to
! decide whether this trajectory is descending or ascending, the tie-breaker
! code, doru ("down or up") forces the sense of the trajectory as follows:
! doru=1 ==> descending
! doru=2 ==> ascending
!=============================================================================
integer,                   intent(in ):: nh,mh,doru
real(dp),dimension(nh),    intent(in ):: ts,hs
real(dp),                  intent(in ):: halfgate
real(dp),dimension(2,  mh),intent(out):: tn
real(dp),dimension(2,2,mh),intent(out):: hn
integer, dimension(    mh),intent(out):: code
logical,                   intent(out):: FF
!-----------------------------------------------------------------------------
real(dp):: tp,hp
integer :: i,im,i2,i2m,imh,n,atti,attim,codeim
!=============================================================================
FF=F
n=nh*2
tp=-hu ! <- default "time at present"
imh=0
do i=1,nh
   i2=i*2
   i2m=i2-1
   hp=hs(i)
   if(ts(i)>tp+geps)then
! A new nominal time of observation:
      imh=imh+1
      tp=ts(i)
      tn(1,imh)=tp-halfgate
      tn(2,imh)=tp+halfgate
      hn(:,:,imh)=hp
   elseif(ts(i)<tp-geps) then
      FF=T
      write(41,*) 'In set_gates; data are not temporally monotonic'
      return
   else
! The same nominal time of observation:      
      hn(1,1,imh)=max(hn(1,1,imh),hp)
      hn(1,2,imh)=min(hn(1,2,imh),hp)
      hn(2,2,imh)=max(hn(2,2,imh),hp)
      hn(2,1,imh)=min(hn(2,1,imh),hp)
   endif
enddo
if(imh/=mh)stop 'In set_gates; inconsistent gate tallies, imh and mh'
! When consecutive gates' post times overlap, adjust their hn if the height
! ranges also overlap:
if(mh==1)then
   code(1)=4*doru
   return
endif
attim=0
code=0 ! <- Default code implies switchable independently for each i
codeim=9 ! <- arbitrary nonzero number
do i=2,mh
   atti=0 ! (default, until we learn anything more definite)
   im=i-1
   if(tn(1,i)<=tn(2,im)+geps)then
! No intermission separates these consecutive gates, im and i:
      if(hn(2,2,im)<=hn(1,2,i))then
         atti=2 ! <-ascending attitude at common time
         code(i)=2
         if(attim==2.and.(codeim==0.or.codeim==2))code(im)=8
      elseif(hn(2,1,im)>=hn(1,1,i))then
         atti=1 ! <-descending attitude at common time
         code(i)=3
         if(attim==1.and.(codeim==0.or.codeim==3))code(im)=4
      else
! Overlapping, attitude at common time neither ascending nor descending,
! but sense of passage through gates must alternate (code=5).
         code(i)=5
         if(hn(2,1,im)<=hn(1,2,i))then; hn(1,2,i)=hn(2,1,im)
                                  else; hn(2,1,im)=hn(1,2,i)
         endif
         if(hn(2,2,im)<=hn(1,1,i))then; hn(2,2,im)=hn(1,1,i)
                                  else; hn(1,1,i)=hn(2,2,im)
         endif
      endif
   else
! Gates im and i separated by an intermission:
      if(hn(2,2,im)<=hn(1,2,i))then
         atti=2 ! <-ascending attitude at intermission
         if(attim==2.and.(codeim==0.or.codeim==2))code(im)=8
      elseif(hn(2,1,im)>=hn(1,1,i))then
         atti=1 ! <-descending attitude at intermission
         if(attim==1.and.(codeim==0.or.codeim==3))code(im)=4
      endif
   endif
   attim=atti
   codeim=code(i)
enddo
end subroutine set_gates

!=============================================================================
subroutine set_posts(mh,mode,tn,hn, bend,tp,hp,off)!               [set_posts]
!=============================================================================
! Given a set of mh double-gates (both descending and ascending types) and
! the array of actual passage modes (i.e., the actual route threading
! the sequence of gates), set the array of actual gateposts coordinates,
! tp and hp, and the corresponding set of signs, bend, by which these
! gatepost constraints, when activatived, must alter the principal
! changed derivative of the optimal spline taking the prescribed route.
! Also, flag (using logical array, "off") those gateposts that, for this
! particular route, are redundant owing to existence of duplication of 
! consecutive pairs of (tp,hp) sometimes occurring when no intermission
! separates consecutive gates.
!=============================================================================
integer,                   intent(in ):: mh
integer, dimension(    mh),intent(in ):: mode
real(dp),dimension(2,  mh),intent(in ):: tn
real(dp),dimension(2,2,mh),intent(in ):: hn
integer, dimension(mh*2),  intent(out):: bend
real(dp),dimension(mh*2),  intent(out):: tp,hp
logical, dimension(mh*2),  intent(out):: off
!-----------------------------------------------------------------------------
real(dp):: tprev,hprev
integer :: i,i2,i2m,i2mm,im,modei
!=============================================================================
off=F
do i=1,mh
   im=i-1
   modei=mode(i)
   i2=i*2; i2m=i2-1; i2mm=i2-2
   tp(i2m)=tn(1,i)
   tp(i2 )=tn(2,i)
   hp(i2m)=hn(1,modei,i)
   hp(i2 )=hn(2,modei,i)
! Check whether gatepost duplications exist, or one dominates another at same t:
   if(i>1)then
      if(tprev+geps>=tp(i2m))then
         if(hprev==hp(i2m))off(i2m)=T
         if(mode(im)==2.and.modei==1)then
            if(hprev<=hp(i2m))then
               off(i2mm)=T
            else
               off(i2m)=T
            endif
         elseif(mode(im)==1.and.modei==2)then
            if(hprev<=hp(i2m))then
               off(i2m)=T
            else
               off(i2mm)=T
            endif
         endif
      endif
   endif
   bend(i2m)=modei*2-3 ! mode=1 ==> bend=-1; mode=2 ==> bend=+1 
   bend(i2 )=-bend(i2m)! mode=1 ==> bend=+1; mode=2 ==> bend=-1 
   tprev=tp(i2)
   hprev=hp(i2)
enddo
end subroutine set_posts

!=============================================================================
subroutine count_routes(n,code,count,FF)!                       [count_routes]
!=============================================================================
! Given the route code array, "code", list all the allowed combinations
! of passage modes (descending === 1; ascending === 2) through the sequence
! of slalom gates.
!=============================================================================
integer,             intent(in ):: n
integer,dimension(n),intent(in ):: code
integer,             intent(out):: count
logical,             intent(out):: FF
!-----------------------------------------------------------------------------
integer,dimension(n):: mode
logical             :: flag
!============================================================================
FF=F
flag=T
do count=0,ihu; call next_route(n,code,mode,flag); if(flag)return; enddo
FF=(count>ihu)
if(FF) write(41,*) 'In count_routes; number of routes exceeds allowance = ',ihu
end subroutine count_routes

!=============================================================================
subroutine list_routes(n,code)!                                  [list_routes]
!=============================================================================
! Given the route code array, "code", list all the allowed combinations
! of passage modes (descending === 1; ascending === 2) through the sequence
! of slalom gates.
!=============================================================================
integer,             intent(in ):: n
integer,dimension(n),intent(in ):: code
!-----------------------------------------------------------------------------
integer,dimension(n):: mode
integer             :: i
logical             :: flag
!============================================================================
write(41,'("List all route combinations of ",i4," allowed passage modes")'),n
flag=T
do i=1,ihu
   call next_route(n,code,mode,flag)
   if(flag)then
      write(41,'(" In list_routes; List of routes complete")'); flag=F; exit
   endif
   write(41,60)i,mode
enddo
if(i>ihu) write(41,'("This list is not necessarily complete")')
60 format(i5,3x,6(2x,5i2))
end subroutine list_routes

!=============================================================================
subroutine next_route(n,code,mode,flag)!                          [next_route]
!=============================================================================
! Given the combinatoric specification of sequentially-conditional
! allowable modes of passage through the n gates encoded in array
! codes, and generically given the present sequence, modes, (a series of
! 1's and 2's denoting respectively descents and ascents through the gates)
! return the next allowed combination defining the updated modes. If instead,
! the intent is to initialize the sequence of modes, input the flag to "true" 
! and the first route (array of modes) will be returned (and the flag lowered 
! to "false").
! If there is no "next" route, the sequence having been already exhausted,
! the flag is raised to "true" on output and the route encoded in array,
! modes, is not meaningful.
! When, at gate i, the preceding gate's mode is "modeim" ( = modes(i-1)) 
! and the present gate's given route code is code=codes(i), the options
! for choosing mode(i) are encoded in the options code, 
! option = options(code,
!=============================================================================
integer,             intent(in   ):: n
integer,dimension(n),intent(in   ):: code
integer,dimension(n),intent(inout):: mode
logical,             intent(inout):: flag
!-----------------------------------------------------------------------------
integer,dimension(0:8,2):: options ! <- evaluates the trinary digit of code
integer,dimension(0:2)  :: firstmode
integer                 :: i,im,j,modeim,modejm,option
data options/0,1,2,0,1,2,0,1,2, 0,0,0,1,1,1,2,2,2/
data firstmode/1,1,2/
!=============================================================================
modeim=1 ! <-arbitrarily set mode of previous gate passage to "descent"
if(flag)then
! Initialize the route sequence and reset the flag:
   do i=1,n
      option=options(code(i),modeim)
      mode(i)=firstmode(option)
      modeim=mode(i)
   enddo
   flag=F
   return
endif

! Use the present route (array of "mode" elements), and the route code, 
! to find the next allowed route, or return with the flag raised when 
! no more allowed routes are to be found:
do i=n,1,-1
   im=i-1
   if(i>1)then
      modeim=mode(im)
   else
      modeim=1
   endif
   option=options(code(i),modeim)
   if(option>0.or.mode(i)==2)cycle
   mode(i)=2
   modejm=mode(i)
   do j=i+1,n
      option=options(code(j),modejm)
      mode(j)=firstmode(option)
      modejm=mode(j)
   enddo
   return
enddo
flag=T
end subroutine next_route

!=============================================================================
subroutine slalom_tspline(n,bend,xn,yn,off,bigX,q, ya,en, &
     ita,maxitb,ittot,FF)                               !      [slalom_spline]
!=============================================================================
! Fit a tensioned spline, characteristic abscissa scale, bigX, between the
! "slalom gates" defined by successive pairs of abscissae, xn, and 
! corresponding ordinate values, yn. Even number n is the total number
! of inequality constraints, or twice the number of gates. There is no
! assumed conditional monotonicity for the gates, but the sense in which
! they are threaded is encoded in the array of signs (-1 or +1), "bend"
! which determines, when activated, the sense in which the gatepost constraint
! changes the principal non-continuous derivative (generally 3rd derivative)
! of the spline. Some gatepost inequality constraints are disabled, as flagged
! by logical array, "off", when two consecutive gateposts constraints are
! identical.
! Subject to the linear inequality constraints, we seek the tensioned
! spline with characteristic scale, bigX, whose energy is minimized.
! The energy of the tensioned spline in the infinitesimal segment [x,x+dx]
! is proportional to half*{ (dy/dx)**2 + (bigT**2)*(ddy/dxx)**2 }*dx.
! The problem is therefore of the type: minimize a quadratic functional
! subject to finitely many (n) linear inequality constraints.
!
! The problem is first standardized by rescaling xn (to xs=xn/bigT) so that
! the characteristic scale becomes unity. We start with a feasible spline
! fitted (equality constraints) to as many of the constraints with distinct
! xs as we can. We "A" iterate from one such feasible, conditionally minimum-
! energy solution to another with a different set of equality constraints
! via an "B" iteration" as follows. The "A" solution generally may have
! constraints at the gateposts that are "pushing" when they should be
! "pulling" (specifically, the sign of the discontinuity in the spline's
! third derivative is the opposite of what it should be at that point). Take
! the most egregious violation and simply switch "off" that constraint to
! initiate the "B" iteration. The energy-minimization "solution" lacking 
! the removed constraint is in general no longer in the feasible region of
! spline-space, but along the line segment to get there, is a conditionally-
! least energy feasible state, energy lower than at the start of this B-step,
! but prevented from going any further by running into a hitherto inactivative
! constraint (ie., bumping up against a new gatepost) which we therefore
! activate. From this point, we once again solve the spline solution, but now
! with this additional equality constraint, and perhaps we once again obtain
! a "solution" outside the feasible region; but we just repeat the remedy. Each
! step of the "B" iterations, where we are adding constraints one at a time,
! we are decreasing the energy of the spline, so the process must terminate.
! At termination of the B iteration, we go back to the A iteration, and 
! see whether there is any remaining constraint whose "jump" in 3rd derivative
! is of the wrong sign, indiciating the constraint to be redundant. But finally
! we must end with a solution that is both feasible and without redundant 
! constraints. Then the problem is solved. 
!
! In general, when the constraint of the final solution is not active, the
! value y of the spline differs from the yn there; it is therefore convenient
! to output what the actual y value of the spline is, which we do in the
! array, ya ("y actual").
!=============================================================================
integer,                intent(in   ):: n
integer, dimension(n),  intent(in   ):: bend
real(dp),dimension(n),  intent(in   ):: xn,yn
logical, dimension(n),  intent(in   ):: off
real(dp),               intent(in   ):: bigX
real(dp),dimension(n),  intent(  out):: q
real(dp),dimension(n),  intent(  out):: ya
real(dp),               intent(  out):: en
integer,                intent(  out):: ita,ittot
integer,                intent(inout):: maxitb
logical,                intent(  out):: FF
!-----------------------------------------------------------------------------
integer,parameter      :: nita=40,nitb=40
real(dp),dimension(n)  :: xs,jump,qt,yat
real(dp)               :: xp,sj,sjmin
integer                :: i,j,k,itb
logical,dimension(n)   :: on
!=============================================================================
FF=F
! For algebraic convenience, work in terms of rescaled times, xs, of 
! constraint
xs=xn/bigX

! Initialize the "A" iteration by fitting a feasible spline to as many
! "gateposts" as is possible with distinct xs. A constraint i is signified
! to be activated when logical array element, on(i), is true:
xp=-hu
do i=1,n
   if(off(i))then; on(i)=F; cycle; endif
   on(i)=(xs(i)>xp+geps); if(on(i))xp=xs(i)
enddo
ittot=1
call fit_gtspline(n,xs,yn,on,qt,jump,yat,en,FF)! <- Make the initial fit
if(FF)then
   write(41,*) 'In slalom_tspline; failure flag raised in call to fit_gtspline'
   write(41,*) 'at initialization of A loop'
   return
endif

! loop over steps of iteration "A" to check for jump-sign violations
do ita=1,nita
   q=qt   ! Copy solution vector q of nodal 1st-derivatives
   ya=yat ! Copy nodal intercepts

! Determine whether there exists a sign-violation in any active "jump"
! of the 3rd derviative and, if so, single out the worst violator, j:
   j=0
   k=0
   sjmin=0
   do i=1,n
      if(.not.on(i))cycle
      sj=-bend(i)*jump(i)
      if(sj<0)then
         j=i
         on(i)=F
      else
         k=k+1 ! <- new tally of constraints switched "on"
      endif
   enddo
   if(j==0)exit !<-Proper conditions for a solution are met
   if(k==0)on(j)=T ! <- must leave at least one constraint "on"

! Begin a new "B" iteration that adds as many new constraints as needed
! to keep the new conditional minimum energy spline in the feasible region:
   do itb=1,nitb
      call fit_gtspline(n,xs,yn,on,qt,jump,yat,en,FF)
      if(FF)then
         write(41,*) 'In slalom_tspline; failure flag raised in call to fit_gtspline'
         write(41,*) 'at B loop, iterations ita,itb = ',ita,itb
         return
      endif
       ittot=ittot+1 ! Increment the running total of calls to fit_tspline

! Determine whether this "solution" wanders outside any slalom gates at
! the unconstrained locations and identify and calibrate the worst violation.
! In this case, sjmin, ends up being the under-relaxation coefficient
! by which we need to multiply this new increment in order to just stay
! within the feasible region of spline space, and constraint j must be
! switched "on":
      j=0
      sjmin=u1
      do i=1,n
         if(on(i).or.off(i))cycle
         sj=bend(i)*(yn(i)-yat(i))
         if(sj<0)then
            sj=(yn(i)-ya(i))/(yat(i)-ya(i))
            if(sj<sjmin)then
               j=i
               sjmin=sj
            endif
         endif
      enddo
      if(j==0)exit !<- spline is feasible, exit B loop and adopt solution as A

! Back off to best feasible solution along this path, which modulates the
! change just made by an underrelaxation factor, sjmin, and activate 
! constraint j
      ya=ya+sjmin*(yat-ya)
      q=q+sjmin*(qt-q)
      on(j)=T
   enddo
   maxitb=max(maxitb,itb)
   if(itb>nitb) then
      FF=T
      write(41,*) 'In slalom_tspline; exceeding the allocation of B iterations'
      return
   end if
   q=qt
   ya=yat
enddo
if(ita>nita)then
   FF=T
   write(41,*) 'In slalom_tspline; exceeding the allocation of A iterations'
   return
endif
end subroutine slalom_tspline

!=============================================================================
subroutine slalom_uspline(n,bend,xn,yn,off,q, ya,en,                         &
     ita,maxitb,ittot,FF)                                !      [slalom_spline]
!=============================================================================
! Like slalom_tspline, except this treats the special case where the spline
! is untensioned, and therefore the characteristic scale in x become infinite,
! and the spline becomes piecewise cubic instead of involving hyperbolic
! (or exponential) function.
!=============================================================================
integer,                intent(in   ):: n
integer, dimension(n),  intent(in   ):: bend
real(dp),dimension(n),  intent(in   ):: xn,yn
logical, dimension(n),  intent(in   ):: off
real(dp),dimension(n),  intent(  out):: q
real(dp),dimension(n),  intent(  out):: ya
real(dp),               intent(  out):: en
integer,                intent(  out):: ita,ittot
integer,                intent(inout):: maxitb
logical,                intent(  out):: FF
!-----------------------------------------------------------------------------
integer,parameter      :: nita=40,nitb=40
real(dp),dimension(n)  :: jump,qt,yat
real(dp)               :: xp,sj,sjmin,ena
integer                :: i,j,k,itb
logical,dimension(n)   :: on
!=============================================================================
! Initialize the "A" iteration by fitting a feasible spline to as many
! "gateposts" as is possible with distinct xn. A constraint i is signified
! to be activated when logical array element, on(i), is true:
FF=F
xp=-hu
do i=1,n
   if(off(i))then
      on(i)=F
      cycle
   endif
   on(i)=(xn(i)>xp+geps)
   if(on(i))xp=xn(i)
enddo
ittot=1
ena=hu
call fit_guspline(n,xn,yn,on,qt,jump,yat,en,FF)! <- Make the initial fit
if(FF)then
   write(41,*) 'In slalom_uspline; failure flag raised in call to fit_guspline'
   write(41,*) 'at initialization of A loop'
   return
endif

! loop over steps of iteration "A" to check for jump-sign violations
do ita=1,nita
   q=qt   ! Copy solution vector q of nodal 1st-derivatives
   ya=yat ! Copy nodal intercepts

! Determine whether there exists a sign-violation in any active "jump"
! of the 3rd derviative and, if so, single out the worst violator, j:
   j=0
   k=0
   sjmin=0
   do i=1,n
      if(.not.on(i))cycle
      sj=-bend(i)*jump(i)
      if(sj<0)then
         j=i
         on(i)=F
      else
         k=k+1 ! <- new tally of constraints switched "on"
      endif
   enddo
   if(j==0)exit !<-Proper conditions for a solution are met
   if(en>=ena)exit
   if(k==0)on(j)=T ! <- must leave at least one constraint "on"
   ena=en

! Begin a new "B" iteration that adds as many new constraints as needed
! to keep the new conditional minimum energy spline in the feasible region:
   do itb=1,nitb
      call fit_guspline(n,xn,yn,on,qt,jump,yat,en,FF)
      if(FF)then
         write(41,*) 'In slalom_uspline; failure flag raised in call to fit_guspline'
         write(41,*) 'at B loop, iterations ita,itb = ',ita,itb
         return
      endif
      ittot=ittot+1 ! Increment the running total of calls to fit_uspline

! Determine whether this "solution" wanders outside any slalom gates at
! the unconstrained locations and identify and calibrate the worst violation.
! In this case, sjmin, ends up being the under-relaxation coefficient
! by which we need to multiply this new increment in order to just stay
! within the feasible region of spline space, and constraint j must be
! switched "on":
      j=0
      sjmin=u1
      do i=1,n
         if(on(i).or.off(i))cycle
         sj=bend(i)*(yn(i)-yat(i))
         if(sj<0)then
            sj=(yn(i)-ya(i))/(yat(i)-ya(i))
            if(sj<sjmin)then
               j=i
               sjmin=sj
            endif
         endif
      enddo
      if(j==0)exit !<- spline is feasible, exit B loop and adopt solution as A

! Back off to best feasible solution along this path, which modulates the
! change just made by an underrelaxation factor, sjmin, and activate 
! constraint j
      ya=ya+sjmin*(yat-ya)
      q=q+sjmin*(qt-q)
      on(j)=T
   enddo
   maxitb=max(maxitb,itb)
   if(itb>nitb) then
      FF=T
      write(41,*) 'In slalom_uspline; exceeding the allocation of B iterations'
      return
   end if
   q=qt
   ya=yat
enddo
if(ita>nita)then
   FF=T
   write(41,*) 'In slalom_uspline; exceeding the allocation of A iterations'
   return
endif
end subroutine slalom_uspline

!=============================================================================
subroutine convertd(n,halfgate,hdata,tdata,phof,                             &
     doru,idx,hs,ts,descending,FF)!                                 [convertd]
!=============================================================================
integer,              intent(in ):: n
real(dp),             intent(in ):: halfgate
real,    dimension(n),intent(in ):: hdata
real,    dimension(n),intent(in ):: tdata
integer, dimension(n),intent(in ):: phof
integer,              intent(out):: doru
integer, dimension(n),intent(out):: idx
real(dp),dimension(n),intent(out):: hs
real(dp),dimension(n),intent(out):: ts
logical,              intent(out):: descending
logical,              intent(out):: FF
!------------------------------------------------------------------------------
integer :: i,j,ii,upsign
real(dp):: s,gate
!=============================================================================
FF=F
if(size(hdata)/=n)stop 'In convertd; inconsistent dimensions of hdata'
if(size(tdata)/=n)stop 'In convertd; inconsistent dimensions of tdata'
if(size(hs)/=n)stop 'In convertd; inconsistent dimensions of hs'
if(size(ts)/=n)stop 'In convertd; inconsistent dimensions of ts'
hs=hdata
! convert to whole number of seconds rounded to the nearest gate=2*halfgate:
upsign=0
gate=halfgate*2
do i=1,n
   ts(i)=nint(tdata(i)*3600/gate)*gate
   if(phof(i)==5)upsign=1  ! Ascending flight
   if(phof(i)==6)upsign=-1 ! Descending flight
enddo
doru=0
if (upsign>0) then
   doru=2
else
   doru=1
endif
if(n==1)return
if(ts(1)>ts(n))then
   do i=1,n/2
      j=n+1-i
      s=ts(i); ts(i)=ts(j); ts(j)=s ! Swap ts
      s=hs(i); hs(i)=hs(j); hs(j)=s ! Swap hs
   enddo
endif
if(upsign==1)then
   descending=F
elseif(upsign==-1)then
   descending=T
else
   descending=(hs(n)<hs(1))
   if(descending)then; upsign=-1; write(41,'("mainly DESCENDING")')
                 else; upsign=1;  write(41,'("mainly ASCENDING")')
   endif
endif

! make sure the order is in time increasing order
do i=1,n
   idx(i)=i
end do
do i=2,n
   do ii=1,i-1
      if (ts(i)<ts(ii).or.(ts(i)==ts(ii).and.upsign*(hs(i)-hs(ii))<u0)) then  
         s=ts(i);  ts(i)=ts(ii);   ts(ii)=s   ! Swap ts
         s=hs(i);  hs(i)=hs(ii);   hs(ii)=s   ! Swap hs
         j=idx(i); idx(i)=idx(ii); idx(ii)=j  ! swap idx
      end if
   end do
end do
do i=2,n
   if(ts(i)<ts(i-1)) then
      write(41,*) 'In convertd; time sequence not monotonic', i, ts(i),ts(i-1)
      FF=T
      return
   end if
enddo
do i=2,n
   if(upsign*(hs(i)-hs(i-1))<u0)&
        write(41,*) 'In convertd; height sequence not monotonic'
enddo
end subroutine convertd

!=============================================================================
subroutine convertd_back(n,wdata,tdata,ws,ts,idx,descending)!  [convertd_back]
!=============================================================================
integer,              intent(in ):: n
integer, dimension(n),intent(in ):: idx
real(dp),dimension(n),intent(in ):: ws
real(dp),dimension(n),intent(in ):: ts
logical,              intent(in ):: descending
real,dimension(n),    intent(out):: wdata
real,dimension(n),    intent(out):: tdata
!------------------------------------------------------------------------------
integer          :: i,j,ii
real(dp)         :: s
real,dimension(n):: wn
real,dimension(n):: tn
!=============================================================================
if(size(wdata)/=n)stop 'In convertd; inconsistent dimensions of wdata'
if(size(tdata)/=n)stop 'In convertd; inconsistent dimensions of tdata'
if(size(ws)/=n)stop 'In convertd; inconsistent dimensions of ws'
if(size(ts)/=n)stop 'In convertd; inconsistent dimensions of ts'

do i = 1, n
   ii = idx(i); tn(ii) = ts(i);  wn(ii) = ws(i)
end do

wdata=wn; tdata=tn
if (descending.or.n==1) return
do i=1,n/2
   j=n+1-i
   s=tdata(i); tdata(i)=tdata(j); tdata(j)=s ! Swap tdata
   s=wdata(i); wdata(i)=wdata(j); wdata(j)=s ! Swap wdata
enddo
end subroutine convertd_back

end module pspl
