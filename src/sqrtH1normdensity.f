c-----------------------------------------------------------------------
c calculate square root of the integrand of the H_AdS^(1,0) norm (defined in arxiv:1110.6794v2)
c where we substitute the coordinates of arxiv:1110.6794v2 by uncompactified quasi-spherical Kerr-Schild coords.
c We write the (Euclidean) volume measure dV=R^2sin(theta)dRdthetadphi as 
c R^2/rho^2dR_drho rho^2 sin(theta)drhodthetadphi=R^2/rho^2 dR_drho dxdydz.
c In the definition of sqrth1normdensity_n we only include the factors R^2/rho^2 dR_drho. Therefore, in order
c to obtain the H_AdS^(1,0) norm for quasi-spherical coordinates, 
c sqrth1normdensity_n needs to be squared and integrated in dxdydz to obtain the H_AdS^(1,0) norm.
c This can be done (up to a numerical factor coming from the integral over the entire space), for example, 
c by using the L^2-norm calculation in DV and then squaring the result.
c-----------------------------------------------------------------------
        subroutine sqrth1normdensity(
     &                  sqrth1normdensity_n,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &                  phys_bdy,ghost_width,
     &                  ief_bh_r0,a_rot)

        implicit none

        logical is_nan
        logical calc_der,calc_adv_quant
        data calc_der/.false./
        data calc_adv_quant/.false./
        real*8  ief_bh_r0,a_rot,M0,M0_min
        integer Nx,Ny,Nz
        integer i,j,k
        integer phys_bdy(6),ghost_width(6)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,ct,L
        real*8 lambda4
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 sqrth1normdensity_np1(Nx,Ny,Nz)
        real*8 sqrth1normdensity_n(Nx,Ny,Nz)
        real*8 sqrth1normdensity_nm1(Nx,Ny,Nz)
        real*8 phi10,h1normdensity0
        real*8 df_drho,df_dtheta,df_dphi
        real*8 dphi1_drho_n,dphi1_dtheta_n,dphi1_dphi_n
        real*8 derphi_dRad

        real*8 rblhor,Radhor,rhohor


        integer is,ie,js,je,ks,ke

        integer i1,j1,k1,a,b,c,d,e,f,g,h,p,q,r
        integer ic,jc,kc
        real*8 efe_ires(4,4)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy,dz
        real*8 x0,y0,z0
        real*8 rho0,theta0,phi0   
        real*8 Rad0
        real*8 drho_dRad,dRad_drho

        real*8 dtdt,dtdrho,dtdtheta,dtdphi
        real*8 dxdt,dxdrho,dxdtheta,dxdphi
        real*8 dydt,dydrho,dydtheta,dydphi
        real*8 dzdt,dzdrho,dzdtheta,dzdphi

        real*8 dxcar_dxqssph(4,4)

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,y,theta,phi)
        !--------------------------------------------------------------
        real*8 gkerrads_ll(4,4),gkerrads_uu(4,4)
        real*8 gkerrads_ll_x(4,4,4),gkerrads_uu_x(4,4,4)
        real*8 gkerrads_ll_xx(4,4,4,4)
        real*8 gammakerrads_ull(4,4,4)
        real*8 gammakerrads_ull_x(4,4,4,4)
        real*8 Hkerrads_l(4)
        real*8 phi1kerrads, phi1kerrads_x(4)
        real*8 detg0_kerrads_qssph0
        real*8 gkerrads_ll_qssph(4,4),gkerrads_uu_qssph(4,4)



        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je,ks,ke/0,0,0,0,0,0,0,0,0/
        data ic,jc,kc/0,0,0/
        data i1,j1,k1,a,b,c,d,e,p,q,r/0,0,0,0,0,0,0,0,0,0,0/

        data dx,dy,dz/0.0,0.0,0.0/
        data x0,y0,rho0/0.0,0.0,0.0/    

        data gkerrads_ll,gkerrads_uu/16*0.0,16*0.0/
        data gkerrads_ll_x,gkerrads_uu_x/64*0.0,64*0.0/
        data gkerrads_ll_xx/256*0.0/
        data gammakerrads_ull/64*0.0/
        data gammakerrads_ull_x/256*0.0/
        data Hkerrads_l/4*0.0/
        data phi1kerrads_x/4*0.0/
        data dxcar_dxqssph/16*0.0/
        data gkerrads_ll_qssph/16*0.0/
        data gkerrads_uu_qssph/16*0.0/



!----------------------------------------------------------------------


      ! Black hole mass
        M0=ief_bh_r0/2
      ! Minimum black hole mass. For M0 below this value, there is a naked singularity
        M0_min=((2*(1 + a_rot**2/L**2) + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2))*Sqrt(-1 + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2) - a_rot**2/L**2))/(3.*Sqrt(6.))

        if (a_rot.ge.L) then
         write (*,*) "ERROR in choice of Kerr-AdS initial parameters: 
     &   the rotation parameter a must be smaller than the AdS radius L"
          write (*,*) "a_rot,L=",a_rot,L
          stop
        end if

        if (M0.le.M0_min) then
          write (*,*) "ERROR in choice of Kerr-AdS initial parameters: 
     &      the black hole mass M0=2*r0 must be larger or equal 
     &      than the M0_min value"
          write (*,*) "M0,M0_min=",M0,M0_min
          stop
        end if

      !event horizon radius in Boyer-Lindquist coordinates rotating at the boundary (non-spherical coordinates)
        rblhor=(Sqrt(-2*a_rot**2 - 2*L**2 + 
     -      (a_rot**4 + 14*a_rot**2*L**2 + L**4)/
     - (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 54*L**4*M0**2 +
     -          Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -   4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 + L**6 + 
     -                54*L**4*M0**2)**2)/2.)**0.3333333333333333 + 
     -   (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 + L**6 
     -  + 54*L**4*M0**2 + 
     -         Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -            4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 
     -                54*L**4*M0**2)**2)/2.)**0.3333333333333333) + 
     -    Sqrt(-4*a_rot**2 - 4*L**2 - 
     -      (a_rot**4 + 14*a_rot**2*L**2 + L**4)/
     -   (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 + L**6 
     -  + 54*L**4*M0**2 +
     -          Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -             4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 
     -                 54*L**4*M0**2)**2)/2.)**0.3333333333333333 - 
     -   (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 + L**6 
     -  + 54*L**4*M0**2 + 
     -         Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -            4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 
     -                54*L**4*M0**2)**2)/2.)**0.3333333333333333 + 
     -      (12*Sqrt(3.)*L**2*M0)/
     -       Sqrt(-2*a_rot**2 - 2*L**2 + 
     -         (a_rot**4 + 14*a_rot**2*L**2 + L**4)/
     -          (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 + L**6 + 
     -             54*L**4*M0**2 + 
     -             Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -                4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 
     -                    54*L**4*M0**2)**2)/2.)**0.3333333333333333 + 
     -    (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  +L**6+54*L**4*M0**2+
     -            Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -               4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 
     -                   54*L**4*M0**2)**2)/2.)**0.3333333333333333)))/
     -  (2.*Sqrt(3.))


        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)


        ! (MAIN LOOP) loop through spacetime points x(i),y(j)
        do i=is,ie
          do j=js,je
           do k=ks,ke

              x0=x(i)
              y0=y(j)
              z0=z(k)
              rho0=sqrt(x0**2+y0**2+z0**2)
              theta0=acos(x0/rho0)
              if (z0.lt.0) then
               phi0=atan2(z0,y0)+2*PI
              else
               phi0=atan2(z0,y0)
              end if
              !uncompactified quasi-sherical Kerr-Schild radial coordinate
              Rad0=2*rho0/(1-rho0**2)
              drho_dRad=((1-rho0**2)**2)/(2*(1+rho0**2))
              dRad_drho=(-1 + rho0)**(-2) 
     &           + (1 + rho0)**(-2)

            !position of event horizon in uncompactified spherical coordinates (non-rotating at the boundary)
            Radhor=(Sqrt(2.)*L*
     -       Sqrt(rblhor**2*(a_rot**2 + rblhor**2)))/
     -       Sqrt(2*L**2*rblhor**2 + a_rot**2*(L**2 - rblhor**2) 
     -       + a_rot**2*(L**2 + rblhor**2)*Cos(2*theta0))

            !position of event horizon in compactified spherical coordinates: rho=sqrt(x**2+y**2+z**2) where x,y,z are the code coordinates
            rhohor=(-1 + Sqrt(1 + Radhor**2))/
     -            Radhor

            if ((rho0.gt.rhohor).and.
     &         (chr(i,j,k).ne.ex)) then

                ! set phi1 value
                phi10=phi1_n(i,j,k)

                call kerrads_derivs_kerrschildcoords(
     &          gkerrads_ll,gkerrads_uu,gkerrads_ll_x,
     &          gkerrads_uu_x,gkerrads_ll_xx,
     &          Hkerrads_l,
     &          gammakerrads_ull,
     &          phi1kerrads,
     &          phi1kerrads_x,
     &          x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k,
     &          ief_bh_r0,a_rot,
     &          calc_der,calc_adv_quant)

!define transformation matrix between Cartesian coordinates to compactified (quasi-)spherical coordinates, 
!        !e.g. dxcar_dxqssph(3,2)=dy/drho

             dtdt=1
             dtdrho=0
             dtdtheta=0
             dtdphi=0
             dxdt=0
             dxdrho=cos(theta0)
             dxdtheta=-rho0*sin(theta0)
             dxdphi=0
             dydt=0
             dydrho=sin(theta0)*cos(phi0)
             dydtheta=rho0*cos(theta0)*cos(phi0)
             dydphi=-rho0*sin(theta0)*sin(phi0)
             dzdt=0
             dzdrho=sin(theta0)*sin(phi0)
             dzdtheta=rho0*cos(theta0)*sin(phi0)
             dzdphi=rho0*sin(theta0)*cos(phi0)

             dxcar_dxqssph(1,1)=dtdt
             dxcar_dxqssph(1,2)=dtdrho
             dxcar_dxqssph(1,3)=dtdtheta
             dxcar_dxqssph(1,4)=dtdphi
             dxcar_dxqssph(2,1)=dxdt
             dxcar_dxqssph(2,2)=dxdrho
             dxcar_dxqssph(2,3)=dxdtheta
             dxcar_dxqssph(2,4)=dxdphi
             dxcar_dxqssph(3,1)=dydt
             dxcar_dxqssph(3,2)=dydrho
             dxcar_dxqssph(3,3)=dydtheta
             dxcar_dxqssph(3,4)=dydphi
             dxcar_dxqssph(4,1)=dzdt
             dxcar_dxqssph(4,2)=dzdrho
             dxcar_dxqssph(4,3)=dzdtheta
             dxcar_dxqssph(4,4)=dzdphi

             do a=1,4
               do b=1,4
                gkerrads_ll_qssph(a,b)=0.0d0
                do c=1,4
                 do d=1,4
                  gkerrads_ll_qssph(a,b)=gkerrads_ll_qssph(a,b)
     &                        +dxcar_dxqssph(c,a)*dxcar_dxqssph(d,b)
     &                          *gkerrads_ll(c,d)
                 end do
                end do
               end do
             end do

                call calc_g0uu(gkerrads_ll_qssph(1,1),
     &              gkerrads_ll_qssph(1,2),
     &              gkerrads_ll_qssph(1,3),
     &              gkerrads_ll_qssph(1,4),
     &              gkerrads_ll_qssph(2,2),
     &              gkerrads_ll_qssph(2,3),
     &              gkerrads_ll_qssph(2,4),
     &              gkerrads_ll_qssph(3,3),
     &              gkerrads_ll_qssph(3,4),
     &              gkerrads_ll_qssph(4,4),
     &              gkerrads_uu_qssph(1,1),
     &              gkerrads_uu_qssph(1,2),
     &              gkerrads_uu_qssph(1,3),
     &              gkerrads_uu_qssph(1,4),
     &              gkerrads_uu_qssph(2,2),
     &              gkerrads_uu_qssph(2,3),
     &              gkerrads_uu_qssph(2,4),
     &              gkerrads_uu_qssph(3,3),
     &              gkerrads_uu_qssph(3,4),
     &              gkerrads_uu_qssph(4,4),
     &              detg0_kerrads_qssph0)


        do a=1,3
          do b=a+1,4
            gkerrads_uu_qssph(b,a)=gkerrads_uu_qssph(a,b) 
          end do
        end do

             dphi1_drho_n  =
     &           df_drho(phi1_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
             dphi1_dtheta_n  =
     &           df_dtheta(phi1_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
             dphi1_dphi_n  =
     &           df_dphi(phi1_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)


            !phi=(1-rho^2)^2*phi1
            derphi_dRad=
     &       drho_dRad*(-4*rho0*(1-rho0**2)*phi10
     &       +(1-rho0**2)**2*dphi1_drho_n)


            h1normdensity0=
     &        (Rad0**2*derphi_dRad**2
     &        +(1-rho0**2)**4*
     &    (
     &      gkerrads_uu_qssph(3,3)
     &         *
     &          dphi1_dtheta_n**2
     &        +2*gkerrads_uu_qssph(3,4)
     &          *
     &          dphi1_dtheta_n*dphi1_dphi_n
     &    +gkerrads_uu_qssph(4,4)
     &      *dphi1_dphi_n**2
     &          )
     &        +((1-rho0**2)**2*phi10)**2)
     &          *Rad0**2/rho0**2
     &          *dRad_drho

!TEST logarithm dependence
        h1normdensity0=log(ct)


                sqrth1normdensity_n(i,j,k)=
     &             sqrt(h1normdensity0)

                 !the phi coordinate is not defined at y=z=0, i.e., theta=0,PI 
                !(not surprising since spherical coordinates are not valid everywhere on the sphere), 
                !hence we set sqrth1normdensity_n(i,j,k) to 0 at these points

                if ((abs(y0).lt.10.0d0**(-10)).and.
     &          (abs(z0).lt.10.0d0**(-10))) then
                    sqrth1normdensity_n(i,j,k)=0.0d0
                end if


            else !i.e., the point is either excised or inside the Kerr-AdS analytic horizon
               sqrth1normdensity_n(i,j,k)=0.0d0
            end if

           end do
          end do
        end do


        return
        end