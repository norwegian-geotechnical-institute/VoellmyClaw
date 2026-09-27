subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)

    use geoclaw_module, only: grav, spherical_distance, rho
    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad

    use geoclaw_module, only: coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: manning_coefficient
    use geoclaw_module, only: manning_break, num_manning
    use geoclaw_module, only: rad2deg

    use slidefriction_module

    implicit none

    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt

    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j, nman
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), coeff
    real(kind=8) :: xm, xc, xp, ym, yc, yp, dx_meters, dy_meters
    real(kind=8) :: u, v, hu0, hv0, u_norm, g, one_minus_r
    real(kind=8) :: tau, wind_speed, theta, phi, psi, P_gradient(2), S(2)
    real(kind=8) :: Ddt, sloc(2)

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    ! Physics
    ! Nominal density of water
    !real(kind=8), parameter :: rho = 1025.d0

    !real(kind=8) :: u_norm, g

    ! Coriolis source term
    ! TODO: May want to remove the internal calls to coriolis as this could
    !       lead to slow downs.
    if (coriolis_forcing) then
        do j=1,my
            y = ylower + (j - 0.5d0) * dy
            fdt = coriolis(y) * dt ! Calculate f dependent on coordinate system

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)

            do i=1,mx
                q(2,i,j) = q(2, i, j) * a(1,1) + q(3, i, j) * a(1,2)
                q(3,i,j) = q(2, i, j) * a(2,1) + q(3, i, j) * a(2,2)
            enddo
        enddo
    endif
    ! End of coriolis source term
    one_minus_r = (1. - rho(1)/rho_s)

    ! Coulomb and Voellmy friction
    if (Coulomb_fric.or.Voellmy_fric) then
        do j=1,my
            do i=1,mx
                if (q(1,i,j) < depth_tolerance) then
                    q(2:3,i,j) = 0.d0
                else
                    u = q(2,i,j)/q(1,i,j)
                    v = q(3,i,j)/q(1,i,j)
                    u_norm = sqrt(u**2 + v**2)
                    if (aux(2,i,j) == 1) then
                        ! water cell. Adjust gravity
                        g = one_minus_r*grav
                    else
                        g = grav
                    endif
                    if (u_norm>1d-10) then

                        if (Coulomb_fric) then
                            gamma = dt*mu*g/u_norm
                            dgamma = 1.d0 + gamma*cos(atan((aux(1,i+1,j)-aux(1,i-1,j))/2.d0/dx))
                            q(2,i,j) = q(2,i,j)/dgamma

                            dgamma = 1.d0 + gamma*cos(atan((aux(1,i,j+1)-aux(1,i,j-1))/2.d0/dy))
                            q(3,i,j) = q(3,i,j)/dgamma
                        endif

                        if (Voellmy_fric) then
                            gamma  = g*u_norm/xi/q(1,i,j)
                            dgamma = 1.d0 + dt*gamma
                            q(2,i,j) = q(2,i,j) / dgamma
                            q(3,i,j) = q(3,i,j) / dgamma
                        endif

                    endif
                endif
            enddo
        enddo
    endif
    ! End of Coulomb and Voellmy friction

end subroutine src2
