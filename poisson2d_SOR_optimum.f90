
!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Iterative Methods for solving 2D Poisson equation
!     d2u/dx2 + d2u/dy2 = f(x,y)
!     Drichlet b.c.
!     
!     Successive over relaxation (SOR)
!     Finding optimal omega
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) 
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Oct. 22, 2015
!-----------------------------------------------------------------------------!

program poisson2d
implicit none
integer::i,j,nitr,nx,ny,k,nom
real*8 ::dx,dy,x0,xL,y0,yL,omega,tol,dom,om1,om2,lmax,pi,optom
real*8,allocatable ::f(:,:),u(:,:),ue(:,:),x(:),y(:)

!Domain
x0 =-1.0d0 !left
xL = 1.0d0 !right

y0 =-1.0d0 !bottom
yL = 1.0d0 !up

!number of points
nx = 20
ny = 20

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)
dy = (yL-y0)/dfloat(ny)

!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

allocate(y(0:ny))
do j=0,ny
y(j) = y0 + dfloat(j)*dy
end do

!given source term and the exact solution on unit square domain
allocate(f(0:nx,0:ny))
allocate(ue(0:nx,0:ny))


do j=0,ny
do i=0,nx
f(i,j) =-2.0d0*(2.0d0-x(i)*x(i)-y(j)*y(j))
ue(i,j)= (x(i)*x(i)-1.0d0)*(y(j)*y(j)-1.0d0)
end do
end do

open(19,file='nitr.plt')
write(19,*) 'variables ="w","N-ITR"'

!Iterative schemes to solve Poisson Equation:
allocate(u(0:nx,0:ny))

!Tolerance
tol = 1.0d-4


!Compute with various omega and plot number of iterations required by the specified omega
nom = 50
om1=0.98
om2=1.98

dom = (om2-om1)/dfloat(nom)
do k=0,nom

  
	!initial guess which satisfies boundary condtions
	!Homogeneous Drichlet B.C.
	!we are not updating boundary points, they are zero
	do j=0,ny
	do i=0,nx
	u(i,j) = 0.0d0
	end do
	end do
  
  	omega = om1 + dfloat(k)*dom
       
	call SOR(nx,ny,dx,dy,f,ue,u,omega,tol,nitr) !Successive over relaxation (SOR)

	!write error
    write(19,*) omega, nitr
    write(*,*) omega, nitr
    
end do

close(19)


!compute optimum omega
pi = 4.0d0*datan(1.0d0)
!lmax = 1.0d0 - 0.25d0*((pi/dfloat(nx))**2 +(pi/dfloat(ny))**2 )
lmax=0.5d0*(dcos(pi/dfloat(nx))+dcos(pi/dfloat(ny)))
optom=2.0d0/(1.0d0 + dsqrt(1-lmax**2))


print*,"-------------------------"
print*,"optimum omega = ", optom
print*,"-------------------------"

open(200,file='field.plt')
write(200,*)'variables ="x","y","u","ue"'
write(200,*)'zone f=point i=',nx+1,',j=',ny+1
do j=0,ny
do i=0,nx
write(200,*) x(i),y(j),u(i,j),ue(i,j)
end do
end do
close(200)


end


!-----------------------------------------------------------------------------!
!Successive over relaxation (SOR)
!-----------------------------------------------------------------------------!
subroutine SOR(nx,ny,dx,dy,f,ue,u,omega,tol,nitr)  
implicit none
integer::nx,ny,i,j,nitr
real*8 ::dx,dy,tol,err,a,omega,up
real*8 ::f(0:nx,0:ny),ue(0:nx,0:ny),u(0:nx,0:ny),e(0:nx,0:ny)


a = -2.0d0/(dx*dx) - 2.0d0/(dy*dy)

nitr = 0

    !compute initial error
 	do j=0,ny
	do i=0,nx
	e(i,j) = dabs(u(i,j)-ue(i,j))
	end do
	end do   
    !max value of the error
    err = maxval(e)
    
do while(err.ge.tol)

	nitr = nitr + 1
  
    !update
	do j=1,ny-1
	do i=1,nx-1
	up = (1.0d0/a)*(f(i,j) &
                   - (u(i+1,j)+u(i-1,j))/(dx*dx) &
                   - (u(i,j+1)+u(i,j-1))/(dy*dy) )
    
	u(i,j) = u(i,j) + omega*(up-u(i,j))
	end do
	end do
    
    !compute error
 	do j=0,ny
	do i=0,nx
	e(i,j) = dabs(u(i,j)-ue(i,j))
	end do
	end do   
    !max value of the error
    err = maxval(e)        
end do

end 






