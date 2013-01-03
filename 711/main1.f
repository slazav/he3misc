c********************************************************************
c BTN: Sample problem 1 (simple usage of BTNEZ, CHKDER)
c main program
c last changed: 10/10/90
c********************************************************************
c
      program          main
      implicit         double precision (a-h,o-z)
      parameter        (n = 100, ndmx = 10,
     *                 lw = (3+7*ndmx)*(n+ndmx))
      double precision x(n), g(n), w(lw)
      external         sfun
c
c specify initial guess of parameters
c
      call xstart (x, n)
c
c check derivative values at initial point
c
      call chkez (n, x, f, g, w, lw, sfun, iflag, errmax, imax)
c
c solve optimization problem
c
      call btnez (n, x, f, g, w, lw, sfun, iflag)
c
c Print results (using output unit 10)
c
      if (iflag .eq. 999) write (*,*) ' DISASTER'
      iunit = 10
      write (iunit,810)
      if (iflag .ne.   0) then
          write (    *,*) ' Error code = ', iflag
          write (iunit,*) ' Error code = ', iflag
      end if
      ig   = idamax (n, g, 1)
      write (iunit,820) f
      write (iunit,830) abs(g(ig))
      write (iunit,840)
      write (iunit,850) (x(i), i = 1,n)
c
      stop
810   format (' Results of optimization problem', /)
820   format (' Final function value = ', 1pd24.16)
830   format (' Norm of gradient     = ', 1pd12.4)
840   format (' Parameters')
850   format (4(1pd16.8))
      end
c
c----------------------------------------------------------------------
c initial guess
c----------------------------------------------------------------------
c
      subroutine xstart (x, n)
      double precision    x(n)
c
c specify initial values of variables
c
      do 10 i = 1,n
          x(i) = 0.d0
10    continue
c
      return
      end
c
c----------------------------------------------------------------------
c function evaluation
c----------------------------------------------------------------------
c
      subroutine sfun (n ,x, f, g)
c
c evaluate nonlinear function and gradient
c
      implicit  double precision (a-h,o-z)
      double precision    x(n), f, g(n)
C$    PRIVATE
c
      f = 0.d0
      do 30 i = 1,n
          b = i
          if (i .gt. n/2) b = -b
          f = f + 5.d-1*i*x(i)*x(i) - b*x(i)
          g(i) = i*x(i) - b
30    continue
c
      return
      end
