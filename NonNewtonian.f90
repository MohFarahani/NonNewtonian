Program ProjectionMethod
Implicit none 
Integer i,j,m,n,k
Real d,LperD,dx,dy,P0,V0,U0,dp,dtemp,Power
Real*8 du,MaxReError,MaxReU,MaxReV,Pr,re,dt,MaxReT,MaxReP,likeGrad2PX,likeGrad2PY,gradUstarX,gradVstarY,GradTX,GradTY,ubarT,VbarT,Grad2TX,Grad2TY
Real*8 Ux,Uy,Uxx,Uyy,Uxy,Uyx,Vx,Vy,Vxx,Vyy,Vxy,Vyx,VbarIJ,VbarIp1J,VbarIn1J,VIp5J,VIp5Jn1,UxIJp1,UxIJn1,UyIp1J,UyIn1J,VxIp5J,VxIp5Jn1,VyIp1Jn5,VyIJn5,UbarIJ,UbarIJp1,UbarIJn1,UIJp5,UIn1Jp5,UxIn5Jp1,UxIn5J,UyIJp5,UyIn1Jp5,VxIJp1,VxIJn1,VyIp1J,VyIn1J
Real*8 Pi5,Dpi5PerDY,Dpi5PerDx,DcPerDx,DdPerDy,DdPerDx,DfPerDy,P1,P2
Real,dimension(:,:),allocatable:: T,P,U,V,ReT,ReP,Ustar,Vstar,Told,Pold,Uold,Vold,ReU,ReV,Pcorrected,Tcorrected,Ucorrected,Vcorrected,x,y
!Print*,"Please Insert The Reynolds Number "
!Read*,re
!print*,"Please Insert The Pr Number "
!Read*,Pr
!Print*,"Please Insert The Ratio Of The Length To Plates Distance (L/D)"
!Read*,LperD
!Print*,"Please Insert The Number Of Horizontal Grid (m)"
!Read*,m
!Print*,"Please Insert The Number Of Vertival Grid (n) "
!Read*,n
!print*,"Please Insert The Initial Dimensionless Pressure(P0) "
!Read*,P0
!print*,"Please Insert The Initial Dimensionless Horizotal Velocity (U0) "
!Read*,U0
!print*,"Please Insert The Initial Dimensionless Vertical Velocity (V0) "
!Read*,V0
!print*,"Please Insert The Maximum Acceptable Difrence Between Valu for U"
!Read*,du
!print*,"Please Insert The Maximum Acceptable Difrence Between Valu for p "
!Read*,dp
!print*,"Please Insert The Maximum Acceptable Difrence Between Valu for T "
!Read*,dtemp
re=100
pr=0.7
Power=5
LperD=20
m=80
n=20
p0=0.1
u0=1
v0=0.0001
du=0.001
dp=0.00001
dtemp=0.000001
!k=1
!open(1,file='output.txt')
!open(3,file='outputP.txt')
d=1
dx=LperD/m    ! Horizontal Grid Distance
dy=d/n    ! Vertical Grid Distance

dt=(2*re*(dx**2)*(dy**2))/(dx**2+dy**2)
Print*,"Please Insert Time step (dt) that should be under this: ",dt
Read*,dt
print*, dx,dy,dt
Allocate (T(m+2,n+2),P(m+2,n+2),ReT(m+2,n+2),ReP(m+2,n+2),V(m+2,n+1),Vstar(m+2,n+1),U(m+1,n+2),Ustar(m+1,n+2),Told(m+2,n+2),Pold(m+2,n+2),Uold(m+1,n+2),Vold(m+2,n+1),ReU(m+1,n+2),ReV(m+2,n+1),Pcorrected(m+2,n+2),Tcorrected(m+2,n+2),Ucorrected(m+2,n+2),Vcorrected(m+2,n+2),x(m+2,n+2),y(m+2,n+2))

!------------------------------------------ Inital Value (Guess)
do i=2,m+1
 do j=2,n+1
  P(i,j)=P0
 end do
end do
!------ Wall (Neuman Condition)----
do i=2,m+1
p(i,1)=p(i,2)
p(i,n+2)=p(i,n+1)
end do
!--- inlet & outlet(Neuman Condition) ---
do j=1,n+2
 P(1,j)=p(2,j)
 P(m+2,j)=P(m+1,j)
end do
!---------- U --------
do i=2,m
  do j=2,n+1
  U(i,j)=u0
  end do
end do
!------- wall ----
do i=2,m
U(i,1)=-U(i,2)
U(i,n+2)=-U(i,n+1)
end do
!--- inlet & outlet ---
do j=1,n+2
U(1,j)=1
U(m+1,j)=U(m,j)
end do
!-------------- V --------
do i=2,m+1
  do j=2,n
  V(i,j)=v0
 end do 
end do
!---- wall ----
do i=2,m+1
V(i,1)=0
V(i,n+1)=0
end do
!---- Inlet & Outlet ----
do j=1,n+1
V(1,j)=-V(2,j)
V(m+2,j)=V(m+1,j)
end do
MaxReError=10000
!-------------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------- Loop -----------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------------
Do while (MaxReError>du)

  !------------------------------ Boundary Condition For  U   ------
  ! Botton & Top (With Wall Condition)
  do i=1,m+1
    U(i,1)=-U(i,2)
	U(i,n+2)=-U(i,n+1)
  end do
  ! Inlet & Outlet ( With Neumann Condition)
  do j=2,n+1
	 U(1,j)=1
	 U(m+1,j)=U(m,j)
  end do

    !---------------------------- Boundary Condition For  V   ------
  ! Botton & Top (With Wall Condition)
  do i=1,m+2
    V(i,1)=0
	V(i,n+1)=0
  end do
  ! Inlet & Outlet ( With Neumann Condition)
  do j=2,n
	 V(1,j)=-v(2,j)
	 V(m+2,j)=V(m+1,j)
  end do


  !------------Calculate The Vstar&Ustar In Between Two Plates------------------------------------------------Ustar & Vstar---------
     do i=2,m
      do j=2,n+1 
     
	  VbarIJ=(V(i+1,j)+V(i+1,j-1)+V(i,j-1)+V(i,j))/4         ! Vbar at point (i,j) according to U point
      VbarIp1J=(V(i+2,j)+V(i+2,j-1)+V(i+1,j)+V(i+1,j-1))/4   ! Vbar at point (i+1,j) according to U point
      VbarIn1J=(V(i-1,j)+V(i-1,j-1)+V(i,j)+V(i,j-1))/4       ! Vbar at point (i-1,j) according to U point
      VIp5J=(V(i,j)+V(i+1,j))/2                              ! V at point (i+1/2,j) according to V point
	  VIp5Jn1=(V(i,j-1)+V(i+1,j-1))/2                            ! V at point (i+1/2,j-1) according to V point
	  
	  UxIJp1=(u(i+1,j+1)-u(i-1,j+1))/(2*dx)                      ! du/dx at point (i,j+1) according to U point
	  UxIJn1=(u(i+1,j-1)-u(i-1,j-1))/(2*dx)                  ! du/dx at point (i,j-1) according to U point
      UyIp1J=(u(i+1,j+1)-u(i+1,j-1))/(2*dy)                  ! du/dy at point (i+1,j) according to U point
      UyIn1J=(u(i-1,j+1)-u(i-1,j-1))/(2*dy)                  ! du/dy at point (i-1,j) according to U point
      
	  VxIp5J=(V(i+1,j)-V(i,j))/dx                            ! dV/dx at point (i+1/2,j) according to V point
      VxIp5Jn1=(V(i+1,j-1)-V(i,j-1))/dx                      ! dV/dx at point (i+1/2,j-1) according to V point
      VyIp1Jn5=(V(i+1,j)-V(i+1,j-1))/dy                      ! dV/dy at point (i+1,j+1/2) according to V point
	  VyIJn5=(V(i,j)-V(i,j-1))/dy                            ! dV/dy at point (i,j+1/2) according to V point
      

	  Vx=(VbarIp1J-VbarIn1J)/(2*dx)                          ! dv/dx at point (i,j) according to U point
      Vy=(VIp5J-VIp5Jn1)/(dy)
      Vyy=(VIp5J-2*VbarIJ+VIp5Jn1)/((dy/2)**2)
      Vxx=(VbarIp1J-2*VbarIJ+VbarIn1J)/(dx**2)
	  Vxy=(VxIp5J-VxIp5Jn1)/dy
	  Vyx=(VyIp1Jn5-VyIJn5)/dx
	 
	  Uxy=(UxIJp1-UxIJn1)/(2*dy)
	  Uyx=(UyIp1J-UyIn1J)/(2*dx)	 
	  Ux=(U(i+1,j)-U(i-1,j))/(2*dx)
      Uy=(U(i,j+1)-U(i,j-1))/(2*dy)
	  Uxx=(U(i+1,j)-2*U(i,j)+U(i-1,j))/(dx**2)
	  Uyy=(U(i,j+1)-2*U(i,j)+U(i,j-1))/(dy**2)
      

	  Pi5=2*(Ux**2+Vy**2)+(Uy+Vx)**2
	  Dpi5PerDX=0.5*(4*(2*Ux*Uxx+2*Vy*Vyx)+4*(Uy+Vx)*(Uyx+Vxx))
      Dpi5PerDY=0.5*(4*(2*Ux*Uxy+2*Vy*Vyy)+4*(Uy+Vx)*(Uyy+Vxy))
      p1=0.5*(Power-1)*((Pi5)**((Power-1)/2-1))
	  p2=(Pi5)**((Power-1)/2)
	  DcPerDx=2*(P1*Dpi5PerDx*Ux+Uxx*P2)
      DdPerDY=P1*Dpi5PerDy*(Uy+Vx)+(Uyy+Vxy)*P2

      Ustar(i,j)=U(i,j)+dt*(DcPerDx+DdPerDy)/Re-dt*U(i,j)*Ux-dt*VbarIJ*Uy
     
	  end do
     end do

     do i=2,m+1
      do j=2,n
    
	  UbarIJ=(U(i,j)+U(i,j+1)+U(i-1,j)+U(i-1,j+1))/4                    ! Ubar at point (i,j) according to V point
      UbarIJp1=(U(i,j+2)+U(i,j+1)+U(i-1,j+2)+U(i-1,j+1))/4              ! Ubar at point (i,j+1) according to V point
      UbarIJn1=(U(i,j)+U(i,j-1)+U(i-1,j)+U(i-1,j-1))/4                  ! Ubar at point (i,j-1) according to V point
      UIJp5=(U(i,j+1)+U(i,j))/2                                         ! U at point (i,j+1/2) according to U point
	  UIn1Jp5=(U(i-1,j+1)+U(i-1,j))/2                                   ! U at point (i-1,j+1/2) according to U point

      UxIn5Jp1=(u(i,j+1)-u(i-1,j+1))/dx                                 ! du/dx at point (i-1/2,j+1) according to U point
	  UxIn5J=(u(i,j)-u(i-1,j))/dx                                       ! du/dx at point (i-1/2,j) according to U point
      UyIJp5=(u(i,j+1)-u(i,j))/dy                                       ! du/dy at point (i,j+1/2) according to U point
	  UyIn1Jp5=(u(i-1,j+1)-u(i-1,j))/dy                                 ! du/dy at point (i-1,j+1/2) according to U point
      
	  VxIJp1=(V(i+1,j+1)-V(i-1,j+1))/(2*dx)                                 ! dV/dx at point (i,j+1) according to V point
      VxIJn1=(V(i+1,j-1)-V(i-1,j-1))/(2*dx)                             ! dV/dx at point (i,j-1) according to V point
      VyIp1J=(V(i+1,j+1)-V(i+1,j-1))/(2*dy)                             ! dV/dx at point (i+1,j) according to V point
	  VyIn1J=(V(i-1,j+1)-V(i-1,j-1))/(2*dy)                             ! dV/dx at point (i-1,j) according to V point

      Ux=(UIJp5-UIn1Jp5)/dx
      Uy=(UbarIJp1-UbarIJn1)/(2*dy)
	  Uyy=(UbarIJp1-2*UbarIJ+UbarIJn1)/(dy**2)
	  Uxx=(UIJp5-2*UbarIJ+UIn1Jp5)/((dx/2)**2)
      Uxy=(UxIn5Jp1-UxIn5J)/dy
      Uyx=(UyIJp5-UyIn1Jp5)/dx


	  Vxy=(VxIJp1-VxIJn1)/(2*dy)
	  Vyx=(VyIp1J-VyIn1J)/(2*dx)
	  Vx=(V(i+1,j)-V(i-1,j))/(2*dx)
      Vy=(V(i,j+1)-V(i,j-1))/(2*dy)
	  Vxx=(V(i+1,j)-2*V(i,j)+V(i-1,j))/(dx**2)
	  Vyy=(V(i,j+1)-2*V(i,j)+V(i,j-1))/(dy**2)

	  Pi5=2*(Ux**2+Vy**2)+(Uy+Vx)**2
	  Dpi5PerDX=0.5*(4*(2*Ux*Uxx+2*Vy*Vyx)+4*(Uy+Vx)*(Uyx+Vxx))
      Dpi5PerDY=0.5*(4*(2*Ux*Uxy+2*Vy*Vyy)+4*(Uy+Vx)*(Uyy+Vxy))
      p1=0.5*(Power-1)*((Pi5)**((Power-1)/2-1))
	  p2=(Pi5)**((Power-1)/2)

	  DfPerDy=2*(P1*Dpi5PerDy*Vy+Vyy*P2)
      DdPerDx=(P1*Dpi5PerDx*(Uy+Vx)+(Uyx+Vxx)*P2)




      Vstar(i,j)=V(i,j)+dt*(DdPerDx+DfPerDy)/Re-dt*V(i,j)*Vy-dt*UbarIJ*Vx
      end do
     end do


!-------------Boundary Condition Ustar & Vstar -----------  
 
!----Wall---- 
    do i=2,m
	 Ustar(i,1)=-U(i,2)
	 Ustar(i,n+2)=-U(i,n+1)
	 end do
!--- Inlet & outlet
  do j=1,n+2
	 Ustar(1,j)=1
	 Ustar(m+1,j)=U(m,j)
     end do
!-------------- Wall ----
    do i=2,m+1
	 Vstar(i,1)=0
	 Vstar(i,n+1)=0
	 end do
!-------------- Inlet & Outlet --
	 do j=1,n+1
    	Vstar(1,j)=-V(2,j)
    	Vstar(m+2,j)=V(m+1,j) 
	end do


 !------------- Calculate The Pressure In Between Two Plates------------------------------------------------------PPPPPPPPPPPPPPPP------
MaxReP=100
 k=1
do While (k<20.or.MaxReP>dp)
Pold=P
!------Boundary Condition
!--- inlet & outlet
	do j=2,n+2
    P(1,j)=P(2,j)
	P(m+2,j)=P(m+1,j)
    end do
!---wall
    do i=1,m+2 
    P(i,1)=P(i,2)
	P(i,n+2)=P(i,n+1)
   end do
!--	----------------------
   do i=2,m+1
    do j=2,n+1
	gradUstarX=(Ustar(i,j)-Ustar(i-1,j))/(dx)
    gradVstarY=(Vstar(i,j)-Vstar(i,j-1))/(dy)
	likeGrad2PX=(P(i+1,j)+P(i-1,j))/(dx**2)
	likeGrad2PY=(P(i,j+1)+P(i,j-1))/(dy**2)
       P(i,j)=(-1/dt*(gradUstarX+gradVstarY)+likeGrad2PX+likeGrad2PY)/(2/dx**2+2/dy**2)
       ReP(i,j)=abs((P(i,j)-Pold(i,j))/(P(i,j)+0.0000001))
	  end do
    end do
MaxReP=maxval(Rep)
!Print*,'MaxReP=  ',MaxRep
k=k+1
   end do

 !-------------Boundary Condition For Pressure Inlet & Outlet & Wall -----------


 !---------------Calculate The Velocity In Between Two Plates --------------------------UUUUUUUUU&&&&&&&&&&&&&VVVVVVVVVVVVVVVVVVVVVVVVVVVVV
uold=u
do i=2,m
 do j=2,n+1
 U(i,j)=Ustar(i,j)-dt*(P(i+1,j)-P(i,j))/(dx)
 ReU(i,j)=abs((U(i,j)-Uold(i,j))/dt)
 end do 
end do
vold=v
do i=2,m+1
  do j=2,n
 V(i,j)=Vstar(i,j)-dt*(P(i,j+1)-P(i,j))/(dy)
 ReV(i,j)=abs((V(i,j)-Vold(i,j))/dt)
 end do
end do

 !----------- Boundary Condition For Velocity ---------
!---- wall
!do i=2,m
!u(i,1)=-u(i,2) 
!u(i,n+2)=-u(i,n+1)
! ReU(i,1)=abs((U(i,1)-Uold(i,1))/dt)
! ReU(i,n+2)=abs((U(i,n+2)-Uold(i,n+2))/dt)
!end do
!-- Inlet & Outlet
!do j=1,n+2
!u(1,j)=1 
!u(m+1,j)=u(m,j)
! ReU(1,j)=0
! ReU(m+1,j)=abs((U(m+1,j)-Uold(m+1,j))/dt)
!end do
!---- Wall
!do i=2,m+1
!V(i,1)=0
! v(i,n+1)=0
! ReV(i,1)=0
! ReV(i,n+1)=0
!end do
!--- Inlet & Outlet 
!do j=1,n+1
!V(1,j)=2*0-V(2,j)
!V(m+1,j)=V(m,j)
!ReV(1,j)=0
!ReV(m+1,j)=abs((V(m+1,j)-Vold(m+1,j))/dt)
!end do 
MaxReU= Maxval(ReU)
MaxReV= Maxval(ReV)
MaxReError=Max(MaxReU,MaxReV)

Print*,'MaxReError= ',MaxReError
End do
Print*,'End Loop'
!--------------------------------------------------------------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------- End The Loop --------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------------------- 

!--------------------------- Calculation T ------------------------

! T=(T-Tw)/(T0-Tw)
!----- Give initial Value
do i=2,m+1
 do j=2,n+1
 T(i,j)=2
 end do 
end do
!----wall
do i=2,m+2
T(i,1)=-T(i,2)
T(i,n+2)=-T(i,n+1)
end do
!----- inlet & outlet
do j=1,n+2
T(1,j)=2-T(2,j)
T(m+2,j)=T(m+1,j)
end do

!------------------
  k=1
   do while (k<20.or.MaxReT>dtemp)
   Told=T
do i=2,m+1
 do j=2,n+1 ! Relative Error T
 gradTX=(T(i,j)-T(i-1,j))/dx
 gradTY=(T(i,j)-T(i,j-1))/dy
 ubarT=(u(i,j)+u(i-1,j))/2
 VbarT=(v(i,j)+v(i,j-1))/2
 Grad2TX=(T(i+1,j)-2*T(i,j)+T(i-1,j))/(dx**2)
 Grad2TY=(T(i,j+1)-2*T(i,j)+T(i,j-1))/(dy**2)
    T(i,j)=Told(i,j)+dt*((Grad2TX+Grad2TY)/(re*pr)-(ubarT*GradTX+VbarT*GradTY))
    ReT(i,j)=abs((T(i,j)-Told(i,j))/dt)
	T(i,j)=T(i,j)
	T(m+2,j)=T(m+1,j)
	T(1,j)=2-T(2,j)
	T(i,1)=-T(i,2)
    T(i,n+2)=-T(i,n+1)
   end do
  end do
k=k+1
  MaxReT=maxval(ReT)
  Print*,MaxReT
end do
!-------- Boundary Condition
!------ wall
  do i=2,m+1
  T(i,1)=-T(i,2)
  T(i,n+2)=-T(i,n+1)
  end do
!----- inlet & Outlet
  do j=1,n+2
  T(m+2,j)=T(m+1,j)
  T(1,j)=2-T(2,j)
  end do

!------------------------------------------------------------------------------------------------------------------------
!------------------------------------- Transmission All Point To Pressure Point + TecPlot Format ------------------------
!------------------------------------------------------------------------------------------------------------------------

!---------------------------------------------pppppppppppppppppppppppppppp----------------
! ------------ P Corner -------------------------
Pcorrected(1,1)=0.25*(P(1,1)+P(1,2)+P(2,1)+P(2,2))
Pcorrected(1,n+2)=0.25*(P(1,n+2)+P(1,n+1)+P(2,n+1)+P(2,n+2))
Pcorrected(m+2,1)=0.25*(P(m+2,1)+P(m+1,1)+P(m+1,2)+P(m+2,2))
Pcorrected(m+2,n+2)=0.25*(P(m+2,n+1)+P(m+1,n+1)+P(m+1,n+2)+P(m+2,n+2))
!------------- P Inner (Not in Boundary) --------
do i=2,m+1
 do j=2,n+1
  Pcorrected(i,j)=P(i,j)
 end do
end do
!------------- P Top & Bottom Wall --------------
do i=2,m+1
 Pcorrected(i,1)=0.5*(P(i,1)+P(i,2))
 Pcorrected(i,n+2)=0.5*(P(i,n+2)+P(i,n+1))
end do
!------------- P Inlet & Outlet -----------------
do j=2,n+1
Pcorrected(1,j)=0.5*(P(1,j)+P(2,j))
Pcorrected(m+2,j)=0.5*(P(m+2,j)+P(m+1,j))
end do
!----------------------------------------------TTTTTTTTTTTTTTTTTTTTTTTTTTTT-------------------
! ------------ TCorner -------------------------
Tcorrected(1,1)=0.25*(T(1,1)+T(1,2)+T(2,1)+T(2,2))
Tcorrected(1,n+2)=0.25*(T(1,n+2)+T(1,n+1)+T(2,n+1)+T(2,n+2))
Tcorrected(m+2,1)=0.25*(T(m+2,1)+T(m+1,1)+T(m+1,2)+T(m+2,2))
Tcorrected(m+2,n+2)=0.25*(T(m+2,n+1)+T(m+1,n+1)+T(m+1,n+2)+T(m+2,n+2))
!------------- T Inner (Not in Boundary) --------
do i=2,m+1
 do j=2,n+1
  Tcorrected(i,j)=T(i,j)
 end do
end do
!------------- T Top & Bottom Wall --------------
do i=2,m+1
 Tcorrected(i,1)=0.5*(T(i,1)+T(i,2))
 Tcorrected(i,n+2)=0.5*(T(i,n+2)+T(i,n+1))
end do
!------------- T Inlet & Outlet -----------------
do j=2,n+1
Tcorrected(1,j)=0.5*(T(1,j)+T(2,j))
Tcorrected(m+2,j)=0.5*(T(m+2,j)+T(m+1,j))
end do
! ------------------------------------------------- U & V -------------------------

do i=2,m+1
  do j=2,n+1

    Ucorrected(i,j)=0.5*(U(i,j)+U(i-1,j))
    Vcorrected(i,j)=0.5*(V(i,j)+V(i,j-1))

  end do 
end do

! -------------------- U & V For Top & Bottom Wall -----

do i=2,m+1

!   Ucorrected(i,1)=0.25*(U(i,1)+U(i,2)+U(i-1,1)+U(i-1,2))
!   Ucorrected(i,n+2)=0.25*(U(i,n+1)+U(i,n+2)+U(i-1,n+1)+U(i-1,n+2))
!   Vcorrected(i,1)=V(i,1)
!   Vcorrected(i,n+2)=V(i,n+1)
Ucorrected(i,1)=0
Ucorrected(i,n+2)=0
Vcorrected(i,1)=0
Vcorrected(i,n+2)=0
end do

! -------------------- U & V For Inlet & Outlet -----

do j=2,n+1

  ! Vcorrected(1,j)=0.25*(V(1,j-1)+V(1,j)+V(2,j)+V(2,j-1))
  ! Vcorrected(m+2,j)=0.25*(V(m+2,j-1)+V(m+2,j)+V(m+1,j)+V(m+1,j-1))
   Vcorrected(1,j)=0
   Vcorrected(m+2,j)=0
   Ucorrected(1,j)=1
   Ucorrected(m+2,j)=UCorrected(m+1,j)
end do

!-------------- U & V Corner

   Ucorrected(1,1)=0 !0.5*(U(1,2)+U(1,1))
   Vcorrected(1,1)=0 !0.5*(V(2,1)+V(1,1))
  
   Ucorrected(1,n+2)=0 !0.5*(U(1,n+2)+U(1,n+1))
   Vcorrected(1,n+2)=0 !0.5*(V(2,n+1)+V(1,n+1))

   Ucorrected(m+2,1)=0 !0.5*(U(m+1,2)+U(m+1,1))
   Vcorrected(m+2,1)=0 !0.5*(V(m+2,1)+V(m+1,1))
  
   Ucorrected(m+2,n+2)=0 !0.5*(U(m+1,n+2)+U(m+1,n+1))
   Vcorrected(m+2,n+2)=0 !0.5*(V(m+2,n+1)+V(m+1,n+1))

!--------------------------------- X,Y   --------------
x(1,1)=0
y(1,1)=0

do j=1,n+2
x(1,j)=0
x(2,j)=dx/2
x(m+2,j)=LperD
end do
do i=1,m+2
y(i,2)=dy/2
y(i,1)=0
y(i,n+2)=1
end do 
do j=3,n+1
y(1,j)=y(1,j-1)+dy
y(2,j)=y(2,j-1)+dy
y(m+2,j)=y(m+2,j-1)+dy
end do
do i=3,m+1
x(i,1)=x(i-1,1)+dx
x(i,2)=x(i-1,2)+dx
x(i,n+2)=x(i-1,n+2)+dx
end do

do i=3,m+1
 do j=3,n+1
  x(i,j)=x(i-1,j)+dx
  y(i,j)=y(i,j-1)+dy
 end do
end do
 

! --------------------------------------------------  TecPlot ------------
open(2,file='output.dat')

write(2,'(a,/,a,I4,a,I4,a)') 'variables= "x","y","Ucorrected","Vcorrected","Pcorrected","Tcorrected"','zone I=' ,m+2,',j=',n+2,',F=point'


do j=1,n+2
 do i=1,m+2
 
 write(2,'(6(2x,e20.6))') x(i,j),y(i,j),Ucorrected(i,j),Vcorrected(i,j),Pcorrected(i,j),Tcorrected(i,j)

 end do
end do





end 


