! Scattering

! Ti, Si parameters
! Einzel & Wolfle JLTP32 (1987) page 34
      subroutine he3_s0s1t0t1(P, S0,S1,T0,T1)
        implicit none
        include 'he3.fh'
        real*8 P
        real*8 f0s, f0a, f1s, f1a
        real*8 A0s, A0a, A1s, A1a
        real*8 S0,S1,T0,T1, W, Wa
        f0s = he3_f0s(P)
        f0a = he3_f0a(P)
        f1s = he3_f1s(P)
        f1a = he3_f1a(P)

        A0s = f0s/(1D0+f0s)
        A0a = f0a/(1D0+f0a)
        A1s = f1s/(1D0+f1s/3D0)
        A1a = f1a/(1D0+f1a/3D0)

        S0 = A0s - 3D0*A0a
        S1 = A1s - 3D0*A1a
        T0 = A0s + A0a
        T1 = A1s + A1a
      end

! Crossection averages over Abrikosov angles
! Einzel & Wolfle JLTP32 (1987) page 34
      function he3_crsect_w(P, S0,S1,T0,T1)
        implicit none
        include 'he3.fh'
        real*8 P, S0,S1,T0,T1, he3_crsect_w
        he3_crsect_w = const_pi/2D0 *
     .    (S0**2 - 2D0/3D0*S0*S1 + 7D0/15D0*S1**2
     .     + 1.5D0*T0**2 - T0*T1 + 7D0/10D0*T1**2)
      end
      function he3_crsect_wi(P, S0,S1,T0,T1)
        implicit none
        include 'he3.fh'
        real*8 P, S0,S1,T0,T1, he3_crsect_wi
        he3_crsect_wi = const_pi/2D0 *
     .    ((S0**2)/3D0 + 2D0/15D0*S0*S1 - 29D0/105D0*S1**2
     .     + 2D0/3D0*S0*T0 - 2D0/5D0*S0*T1
     .     + 5D0/3D0*(25D0-36D0*dlog(2D0))*T0**2
     .     + (84D0-120D0*dlog(2D0))*T0*T1
     .     + 5D0/21D0*(173D0-252D0*dlog(2D0))*T1**2)
      end
      function he3_crsect_wd(P, S0,S1,T0,T1)
        implicit none
        include 'he3.fh'
        real*8 P, S0,S1,T0,T1, he3_crsect_wd
        he3_crsect_wd = const_pi/2D0 *
     .    (7D0/15D0*S0**2 - 18D0/35D0*S0*S1 + 107D0/315D0*S1**2
     .     + 8D0/15D0*S0*T0 - 8D0/105D0*(S0*T1 + S1*T0)
     .     + 8D0/63D0*S1*T1 + 29D0/30D0*T0**2
     .     - 19D0/35D0*T0*T1 + 33D0/70D0*T1**2)
      end

! Scattering factors
! Einzel & Wolfle JLTP32 (1987) f.74
! code from Samuli
      function He3_scatt_l1a(P)
        implicit none
        include 'he3.fh'
        real*8 P
        real*8 S0,S1,T0,T1, W, Wl
        real*8 he3_crsect_w
        call he3_s0s1t0t1(P, S0,S1,T0,T1)
        W = he3_crsect_w(P, S0,S1,T0,T1)
        Wl = const_pi * 1D0/420D0 *
     .    (- 70D0*S0**2 - 54D0*S1**2 + 175D0*T0**2
     .     + 28D0*S0*(7D0*S1 + 10D0*T0 - 6D0*T1)
     .     - 42D0*T0*T1 + 71D0*T1**2
     .     + 8D0*S1*(-21D0*T0 + 19D0*T1))
        He3_scatt_l1a = Wl / W
      end

      function He3_scatt_g0(P)
        implicit none
        include 'he3.fh'
        real*8 P
        real*8 S0,S1,T0,T1, W, Wi
        real*8 he3_crsect_w, he3_crsect_wi
        call he3_s0s1t0t1(P, S0,S1,T0,T1)
        W  = he3_crsect_w(P, S0,S1,T0,T1)
        Wi = he3_crsect_wi(P, S0,S1,T0,T1)
        He3_scatt_g0 = Wi / W
      end

      function He3_scatt_d0(P)
        implicit none
        include 'he3.fh'
        real*8 P
        real*8 S0,S1,T0,T1, W, Wd
        real*8 he3_crsect_w, he3_crsect_wd
        call he3_s0s1t0t1(P, S0,S1,T0,T1)
        W  = he3_crsect_w(P, S0,S1,T0,T1)
        Wd = he3_crsect_wd(P, S0,S1,T0,T1)
        He3_scatt_d0 = Wd / W
      end

! Normal state quasiparticle lifetime at the Fermi level
! \tau_N(0,T)
! Einzel JLTP32 (1978) p.28,34
! Einzel JLTP84 (1991) f.4
      function He3_tau_n0(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 p,ttc
        real*8 S0,S1,T0,T1, W, Wa
        real*8 he3_crsect_w
        call he3_s0s1t0t1(p, S0,S1,T0,T1)
        W = he3_crsect_w(p, S0,S1,T0,T1)
        He3_tau_n0 = 32 *
     .    he3_tfeff(P)*const_hbar/const_pi**2
     .    / W / const_kb / (1D-3*ttc*he3_tc(P))**2
      end

! Thermal average of normal state quasiparticle lifetime
! \bar\tau_N(T)
! Einzel JLTP84 (1991) f.5
      function He3_tau_n_av(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 p,ttc
        He3_tau_n_av = 0.75D0 * He3_tau_n0(ttc, p)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Collision integral for Bogoliubov quasiparticles

! Collision integral in Einzel approximation
! Einzel, Wolfle, Hirschfeld, JLTP80 (1990), Appendix, p.66
      function he3_coll_int(xi,ttc, gap, g0, d0)
        implicit none
        include 'he3.fh'
        real*8 xi, ttc, gap, x
        real*8 J0,J1,J2,J3, K0,K1,K2,K3, I0,I1,I2,I3
        real*8 a0,a1,a2,a3, b0,b1,b2, g0,d0

        x = gap/ttc
        a0 = -0.5768D0
        a1 =  0.2694D0
        a2 =  0.29D0
        a3 = -0.08D0
        J0 = 3D0/4D0/const_pi**0.5D0
     .       * (1D0+2D0*x)**1.5D0
     .       / (dexp(x) + a0 + a1*x + a2*x**2 + a3*x**3)
        J1 = x**2 / (2D0*const_pi)**0.5D0
     .       * (0.5D0 + x)**1.5D0
     .       / (1.3D0 + x**2) * dexp(-x)
        J2 = x**2 / (2D0*const_pi)**0.5D0
     .       / (0.5D0 + x)**0.5D0
     .       / (dexp(x)+2.3D0)
        J3 = 3D0*x**4 / (2D0*const_pi)**0.5D0
     .       / (0.5D0 + x)**0.5D0
     .       / (0.2D0 + x**2)
     .       / (dexp(x)+1D0+0.1D0*x**2)
        b0 =  3.4296D0
        b1 = -3.2148D0
        b2 =  2.375D0
        K0 = 9D0/8D0/(2D0*const_pi)**0.5D0
     .       / (0.5D0 + x)**0.5D0
     .       / (dexp(x) + b0 + b1*x + b2*x**2)
        K1 = - 5D0/8D0/(2D0*const_pi)**0.5D0
     .       * x**2 / (1D0+x)**0.5D0
     .       * dexp(-x) / (const_pi**2 + x**2)
        K2 = 3D0/8D0/(2D0*const_pi)**0.5D0
     .       * x**2 / (1D0+x)**0.5D0
     .       * dexp(-x) / (127D0/150D0 * const_pi**2 + x**2)
        K3 = - 15D0/8D0/(2D0*const_pi)**0.5D0
     .       * x**4/(1D0+x)**0.5D0
     .       * dexp(-x) / (const_pi**2 + x**2)**2
!       Possible error in K3: in the paper it contains (1+x**2)
!       but maybe it is (1+x), because Ki sould be ~exp(x)/x at low temp.
!       Anyway change looks negligable.
        I0 = J0 + K0 * (xi/ttc)**2
        I1 = J1 + K1 * (xi/ttc)**2
        I2 = J2 + K2 * (xi/ttc)**2
        I3 = J3 + K3 * (xi/ttc)**2
        he3_coll_int = ( I0 - g0*(I1+I2) + d0*I3 )
      end

! Collision integral for low temp (good for < 0.7Tc)
! Einzel, JLTP84 (1991), p.345
      function he3_coll_int_lt(xi,ttc, gap, g0, d0)
        implicit none
        include 'he3.fh'
        real*8 xi, ttc, gap, x, g0,d0,w0
        x=xi/dsqrt(2*ttc*gap)
        w0 = (1D0 - 2D0/3D0*g0 + d0)
        he3_coll_int_lt =
     .    3D0/2D0/const_pi * gap/ttc * he3_yosida(ttc,gap, 0D0)
     .    * (w0 + ttc/gap*(0.75D0*(1D0 + x**2) * w0
     .                      - (1D0+2D0*x**2)*(g0/3D0+d0) ))
      end

! Collision integral for high temp (good above 0.95 Tc)
! Einzel, JLTP84 (1991), p.345
      function he3_coll_int_ht(xi,ttc, gap, g0, d0)
        implicit none
        include 'he3.fh'
        real*8 xi, ttc, gap, x, g0,d0,w0
        he3_coll_int_ht = 1D0 + (xi**2 + gap**2)/(ttc*const_pi)**2
      end

! Quasiparticle lifetime at fermi level
! Einzel JLTP84 (1991) p.344
      function he3_tau0(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, g0, d0, tn
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_tau0 = NaN
          return
        endif
        g0  = he3_scatt_g0(p)
        d0  = he3_scatt_d0(p)
        gap = he3_trivgap(ttc, p)
        tn  = he3_tau_n0(ttc, p)
        he3_tau0 = tn / he3_coll_int(0D0, ttc, gap, g0, d0);
      end

! Integrand for tau_av calculations
! Integration is similar to Y0 calculation in he3_gap.f
      function he3_tau_av_int(x,ttc, gap, g0, d0)
        implicit none
        real*8 x,ttc, gap, g0, d0, xi, ek, C
        real*8 he3_tau_av_int, he3_coll_int
        C=3D0 ! see tests/plot_tauav_int.m
        xi = datanh(x)*C
        ek=dsqrt(xi**2 + gap**2)
        he3_tau_av_int = he3_coll_int(xi,ttc, gap, g0, d0)
     .   / (dcosh(ek/(2D0*ttc)))**2
     .   * C/(1D0-x**2)
      end

! Averaged quasiparticle lifetime
! Einzel JLTP84 (1991) p.345
      function he3_tau_av(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, sum, g0, d0, Y0, tn
        real*8 dx, xp, xm
        real*8 he3_tau_av_int
        integer i, maxi
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_tau_av=NaN
          return
        endif
        g0  = he3_scatt_g0(p)
        d0  = he3_scatt_d0(p)
        gap = he3_trivgap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        tn  = he3_tau_n0(ttc, p)
        sum = 0D0
        maxi=1000
        dx=1D0/dble(maxi)
        ! intergation from 0 to 1 using Gaussian quadrature
        do i=1,maxi 
          xp = dx * (dble(i) - 0.5D0 + 0.5D0/dsqrt(3D0))
          xm = dx * (dble(i) - 0.5D0 - 0.5D0/dsqrt(3D0))
          sum = sum
     .       + he3_tau_av_int(xp, ttc, gap, g0, d0) * dx/2D0
     .       + he3_tau_av_int(xm, ttc, gap, g0, d0) * dx/2D0
        enddo
        he3_tau_av = Y0 * tn / (2D0*sum)
      end

! Integrand for he3_fpath calculation
! Integration is similar to Y0 calculation in he3_gap.f
      function he3_fpath_int(x,ttc, gap, g0, d0)
        implicit none
        real*8 x,ttc, gap, g0, d0, xi, ek, I, C
        real*8 he3_fpath_int, he3_coll_int
        C=2D0 ! see tests/plot_tauav_int.m
        xi = datanh(x)*C
        Ek=dsqrt(xi**2 + gap**2)
        I = he3_coll_int(xi,ttc, gap, g0, d0)
        he3_fpath_int =
     .   1/I**2  ! (t/tN)^2
     .   * (xi/Ek)**2
     .   / (1D0 + dexp(Ek/ttc)) ! Fermi function
     .   * C/(1D0-x**2)
      end

! Mean free path of Bogoliubov quasiparticles
! Einzel JLTP32 (1978) f.84
      function he3_fpath(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, sum, g0, d0, Y0, tn, vf, Ife
        real*8 dx, xp, xm
        real*8 he3_fpath_int
        integer i, maxi
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_fpath=NaN
          return
        endif
        g0  = he3_scatt_g0(p)
        d0  = he3_scatt_d0(p)
        gap = he3_trivgap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        tn  = he3_tau_n0(ttc, p)
        vf  = he3_vf(p)
        sum = 0D0
        maxi=100
        dx=1D0/dble(maxi)
        ! intergation from 0 to 1 using Gaussian quadrature
        do i=1,maxi 
          xp = dx * (dble(i) - 0.5D0 + 0.5D0/dsqrt(3D0))
          xm = dx * (dble(i) - 0.5D0 - 0.5D0/dsqrt(3D0))
          sum = sum
     .       + he3_fpath_int(xp, ttc, gap, g0, d0) * dx/2D0
     .       + he3_fpath_int(xm, ttc, gap, g0, d0) * dx/2D0
        enddo
        Ife = ttc * dlog(1+dexp(-gap/ttc)) ! integral of fermi funvtion
        he3_fpath = vf * tn * dsqrt(sum/Ife)
      end

! Spin diffusion perpendicular transport time, s
! Einzel JLTP84 (1991) f.90,96
      function he3_tau_dperp(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, l1a, gap, y0, y2
        l1a = he3_scatt_l1a(p)
        gap = he3_trivgap(ttc, p)
        y0  = he3_yosida(ttc, gap, 0D0)
        y2  = he3_yosida(ttc, gap, 2D0)
        he3_tau_dperp = he3_tau_av(ttc,p)
     .   /(1D0 - l1a*(4D0*y0 + y2)/5D0/y0)
      end

! Spin diffusion parallel transport time, s
! Einzel JLTP84 (1991) f.90,96
      function he3_tau_dpar(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, l1a, gap, y0, y2
        l1a = he3_scatt_l1a(p)
        gap = he3_trivgap(ttc, p)
        y0  = he3_yosida(ttc, gap, 0D0)
        y2  = he3_yosida(ttc, gap, 2D0)
        he3_tau_dperp = he3_tau_av(ttc,p)
     .   /(1D0 - l1a*(2D0*y0 + 3D0*y2)/5D0/y0)
      end


! Integrand for spin diffusion calculation
! Integration is similar to Y0 calculation in he3_gap.f
! 2D integration needed: x [0:1] and th angle [0:pi]
! Einzel JLTP84 (191) f.108
! o0 -- Larmor freq (rad/s)
! oe -- Exchange freq (rad/s)
! td -- Spin diffusion perpendicular transport time, s
      function he3_sdiff_int(x, th, ttc, gap, o0, oe, td)
! tau_dp w0 w_exch Vf suspar
        implicit none
        real*8 x,th, ttc, gap, o0, oe, td
        real*8 C, xi, Ek, kz, kp, phi, u
        complex*16 t,s, Sp2
        complex*16 he3_sdiff_int
        C=2D0
        xi = datanh(x)*C
        Ek=dsqrt(xi**2 + gap**2)
        phi = (dcosh(Ek/(2D0*ttc)))**(-2)
        u = xi/Ek

        kz=dsin(th)
        kp=dcos(th)

        t = td / dcmplx(1D0, o0 * td)
        s = (o0 + oe) * t

        Sp2 = (u*(1-kz**2)+kz**2)**2


        he3_sdiff_int =  kz**3 * t
     .    * (3D0/8D0*(u-1D0)**2 * kp**4
     .         + (u-1D0)*kp**2 + 1 - (0,1)*s*Sp2)
     .    / (1 + s**2*Sp2)
     .    * phi * C/(1D0-x**2)
      end

! Integrand for spin diffusion calculation - 2
! Integrate he3_sdiff_int by th angle [0:pi],
! keep x for future integration
      function he3_sdiff_int2(x, ttc, gap, o0, oe, td)
        implicit none
        include 'he3.fh'
        real*8 x, ttc, gap, o0, oe, td
        real*8 dt, tp, tm
        complex*16 he3_sdiff_int, he3_sdiff_int2, sum
        integer i, maxi

        sum = (0D0, 0D0)
        maxi = 100
        dt = const_pi/dble(maxi)
        ! intergation th from 0 to pi using Gaussian quadrature
        do i=1,maxi
          tp = dt * (dble(i) - 0.5D0 + 0.5D0/dsqrt(3D0))
          tm = dt * (dble(i) - 0.5D0 - 0.5D0/dsqrt(3D0))
          sum = sum
     .       + he3_sdiff_int(x, tp, ttc,gap,o0,oe,td) * dt/2D0
     .       + he3_sdiff_int(x, tm, ttc,gap,o0,oe,td) * dt/2D0
        enddo
        he3_sdiff_int2 = sum
      end

! Spin diffusion coefficient D_perp, cm/s
! Integrate he3_sdiff_int2 by th angle [0:pi],
      function he3_sdiff(ttc, p, nu0)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        real*8 gap, Vf, suspar, Y0, f0a
        real*8 o0, oe, td
        real*8 dx, xp, xm
        complex*16 he3_sdiff_int2, sum
        integer i, maxi
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_sdiff=NaN
          return
        endif
        gap = he3_trivgap(ttc, p)
        Vf  = he3_vf(p)
        f0a = he3_f0a(p)
        Y0  = he3_yosida(ttc, gap, 0D0);
        suspar = (2D0 + Y0) / (3D0 + f0a*(2D0 + Y0))
        o0  = nu0/2/const_pi
        oe  = -f0a*o0*suspar
        td  = he3_tau_dperp(ttc, p)

        sum = (0D0, 0D0)
        maxi=100
        dx=1D0/dble(maxi)
        ! intergation x from 0 to 1 using Gaussian quadrature
        do i=1,maxi
          xp = dx * (dble(i) - 0.5D0 + 0.5D0/dsqrt(3D0))
          xm = dx * (dble(i) - 0.5D0 - 0.5D0/dsqrt(3D0))
          sum = sum
     .       + he3_sdiff_int2(xp, ttc, gap, o0, oe, td) * dx/2D0
     .       + he3_sdiff_int2(xm, ttc, gap, o0, oe, td) * dx/2D0
        enddo

        he3_sdiff = dreal(sum) / 4D0 * Vf**2 / suspar / ttc
      end
