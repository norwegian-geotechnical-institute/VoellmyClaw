subroutine setprob()

    use qinit_module, only: qinit_style
    use slidefriction_module, only: Coulomb_fric, mu, Voellmy_fric, xi, rho_s, water_level

    implicit none

    integer :: iunit
    character(len=25) fname

    iunit = 7
    fname = 'setprob.data'

    call opendatafile(iunit, fname)

    read(7,*) Coulomb_fric
    read(7,*) mu
    read(7,*) Voellmy_fric
    read(7,*) xi
    read(7,*) rho_s
    read(7,*) water_level
    read(7,*) qinit_style

end subroutine setprob
