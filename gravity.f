      PROGRAM Gravity
C*********************************************************
C  This code calculates the gravitational potential,
C  from which gravity vs Z can be calculated from
C  earth surface and up to 1000 meters.
C  MODEL: The Earth is assumed as a solid sphere with
C         3 mass densities:
C      
C  Rho = Rho0 ..........for r <= R-H
C      = Rho1 ..........for R-H < r < R and Theta < Alpha
C      = Rho2 ..........for R-H < r < R and Theta > Alpha
C
C  Theta is Azimuthal angle, R is radius of Earth,
C  H is the depth of Crust.
C*********************************************************
      PARAMETER (Alfa=10.0,High=1000.0)
      PARAMETER (N=1001)
      REAL Z(N),V(N),grav(N),VV0,VV1,VV1p,VV2
      REAL PI,G,M,R,Rho0,Rho1,Rho2,Rho1av,Dz,zz,Alpha,H
C
      PI=4.0*ATAN(1.0)
      Alpha=PI*Alfa/180.0
      G=6.67408E-11
      M=5.972E24
      R=6.378E6
      Rho1av=2600.0
      Rho0=5545.0
      Rho1=2750.0
      Rho2=(PI/(PI-Alpha))*(Rho1av-Rho1*Alpha/PI)
      print*,'Please Give me H ='
      read*,H
c---------------------------------
c  Measure g just above the ground 
c  from Z=0 to Z= 1000 meters.
c
      Dz=High/FLOAT(N-1)
      DO i=1,N
       Z(i)=R+Dz*FLOAT(i-1)
      ENDDO
c-----------------------------------------
c  Calculation of gravitational potential
      DO i=1,N
       zz=Z(i)
       VV0=V0(zz,Rho0,G,R,H)
       VV1=V1(zz,Rho1,G,R,H,Alpha)
       VV1p=V1(zz,Rho2,G,R,H,Alpha)
       VV2=V2(zz,Rho1,Rho2,G,R,H,Alpha,VV1p)
       V(i)=VV0+VV1+VV2
      ENDDO
c-----------------------------------------
c  Calculation of gravity
C
      OPEN(UNIT=10,FILE='g.dat')
      DO i=1,N
       WRITE(10,1)Z(i)-R,V(i)-V(1)
      ENDDO
1     FORMAT(F10.3,3x,E15.3)
      print*,'Gravity = ',(V(1001)-V(1))/1000
C
      STOP
      END


      REAL FUNCTION V0(zz,Rho0,G,R,H)
C******************************************************
C  Calculation of contribution: V0
C******************************************************
      REAL zz,Rho0,fact1,fact2
      REAL G,R,H,PI,Rm,Mm,gm,Vm
      PI=4.0*ATAN(1.0)
      IF(zz.GE.(R-H))THEN
        fact1=G*Rho0*4.0*PI/(3.0*zz)
        fact2=(R-H)**3
        V0=-fact1*fact2
      ELSE
        Rm=(R-H)
        Mm=Rho0*4.0*PI*Rm**3/3.0
        gm=G*Mm/(Rm**2)
        Vm=-G*Mm/Rm
        V0=Vm+(gm/(2.0*Rm))*(zz**2-Rm**2)
      ENDIF
      RETURN
      END


      REAL FUNCTION V1(zz,Rho1,G,R,H,Alpha)
C******************************************************
C  Calculation of contribution: V1
C******************************************************
      REAL zz,Rho1,F1,F2,F3,Term1,Term2,Term3
      REAL T1A,T1B,T2A,T2B
      REAL G,R,H,Alpha,PI
      PI=4.0*ATAN(1.0)
      F1=2.0*PI*G*Rho1/zz
      F2=-2.0*PI*G*Rho1/(3.0*zz)
      F3=-2.0*PI*G*Rho1*COS(Alpha)
      T1A=zz*(R**2)/2.0 -(R**3)/3.0
      T1B=-zz*((R-H)**2)/2.0 + ((R-H)**3)/3.0
      Term1=F1*(T1A+T1B)
      T2A=(R**2+zz**2-2.0*zz*R*COS(Alpha))**(3.0/2.0)
      T2B=((R-H)**2+zz**2-2.0*zz*(R-H)*COS(Alpha))**(3.0/2.0)
      Term2=F2*(T2A-T2B)
C  To do the integration:
      Term3=F3*Simson(R,H,zz,Alpha)
      V1=Term1+Term2+Term3
      RETURN
      END

      REAL FUNCTION SIMSON(R,H,zz,Alpha)
C*****************************************
C Do Integration using Simpson Rule.
C*****************************************
      PARAMETER (Nx=1001)
      REAL R,H,zz,Alpha,Dx,rr,fact
      REAL Y(Nx)
      Dx=H/FLOAT(Nx-1)
      DO i=1,Nx
       rr=(R-H)+Dx*FLOAT(i-1)
       Y(i)=SQRT(rr**2+zz**2-2.0*rr*zz*COS(Alpha))
      ENDDO
      fact=Dx/3.0
      SIMSON=0.0
      DO i=1,N
       IP=MOD(i,2)
       IF(i.EQ.1.OR.i.EQ.N)THEN
         SIMSON=SIMSON+Y(i)*fact
       ELSEIF(IP.EQ.0)THEN
         SIMSON=SIMSON+4.0*Y(i)*fact
       ELSE
         SIMSON=SIMSON+2.0*Y(i)*fact
       ENDIF
      ENDDO
      RETURN
      END

   
      REAL FUNCTION V2(zz,Rho1,Rho2,G,R,H,Alpha,VV1p)
C******************************************************
C  Calculation of contribution: V2
C******************************************************
      REAL zz,Rho1,Rho2,Fact
      REAL G,R,H,Alpha,PI,VV1p
      PI=4.0*ATAN(1.0)
      Fact=-4.0*PI*G*Rho2/(3.0*zz)
      V2=Fact*(R**3-(R-H)**3)-VV1p
      RETURN
      END
