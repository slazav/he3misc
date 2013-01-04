c********************************************************************
c BTN: Sample problem 4 (using BTN to solve a constrained problem)
c
c A logarithmic barrier method is used to minimize a nonlinear
c function subject to simple bounds on the variables ("box" constraints)
c
c main node program
c last changed: 10/10/90
c********************************************************************
c
      program          node
      implicit         double precision (a-h,o-z)
      parameter        (n = 100, ndmx = 16, liw = 3*ndmx,
     *                 lw = 14*n + ndmx*(4*ndmx + 9))
      double precision x(n), g(n), w(lw), lb(100), ub(100), mu
      integer          iw(liw), pid
      logical          idzero
      common /bnd/     lb, ub
      common /bar/     mu
      external         sfun, maxstp
c
c set up node information
c
      id     = mynode ()
      pid    = 0
      idzero = id .eq. 0
      iu = 9
      open (iu, file = 'node.out', status = 'unknown')
c
c Initialize nonlinear parameters x and barrier parameter mu
c
      call xstart (x, n)
      mu    = 1.d0
c
c set up parameters for derivative checker
c
      call chkpar (n, kmax, msglvl, iunit, istart, iend, istep)
c
c check derivative values at initial point
c
      call chkder (n, x, f, g, w, lw, iw, sfun, iflag,
     *             ndmx, kmax, errmax, imax,
     *             msglvl, iunit, istart, iend, istep)
      if (idzero) write (iu,800) errmax, imax
      if (iflag .ne. 0) stop
c
c Set customizing parameters for BTN (the user-supplied routine
c maxstp will be used to bound the step length in the line search)
c
      call btnpar (nodemx, kmax, maxit, msglvl, pid, iprec, nlmax,
     *          initv, tolq, iunit, rnktol, maxfun, accrcy, stepmx,
     *          eta, ireset, indstp)
      indstp = 1
c
c solve optimization problem
c
      scale = 10.d0
10    if (idzero) write (*,860) mu
      call btn (n, x, f, g, w, lw, iw, sfun, iflag,
     *          nodemx, kmax, maxit, msglvl, pid, iprec, nlmax,
     *          initv, tolq, iunit, rnktol, maxfun, accrcy,
     *          stepmx, eta, ireset, indstp, maxstp)
      if (iflag .eq. 999) then
          if (idzero) write (*,*) ' Fatal error in BTN'
          stop
      end if
c
c update barrier parameter
c
      mu = mu / scale
      if (mu .gt. 1.d-6) go to 10
c
c Print results (node 0 only)
c
      if (.not. idzero) stop
      write (iu,810)
      if (iflag .ne.   0) then
          write ( *,*) ' Error code = ', iflag
          write (iu,*) ' Error code = ', iflag
      end if
      ig   = idamax (n, g, 1)
      write (iu,820) f
      write (iu,830) abs(g(ig))
      write (iu,840)
      write (iu,850) (x(i), i = 1,n)
      close (iu)
      stop
800   format (' CHKDER: Max. error in gradient = ', 1pd12.4,
     *        ' at index ', i5)
810   format (//, ' Results of optimization problem', /)
820   format (' Final function value = ', 1pd24.16)
830   format (' Norm of gradient     = ', 1pd12.4)
840   format (' Parameters')
850   format (4(1pd16.8))
860   format (//,' mu =', 1pd12.4)
      end
c
c----------------------------------------------------------------------
c initial guess
c----------------------------------------------------------------------
c
      subroutine xstart (x, n)
      double precision x(n)
      double precision lb(100), ub(100), mu
      common /bnd/     lb, ub
      common /bar/     mu
c
c specify initial values of variables and their bounds
c
      do 10 i = 1,n
          x(i)  = 0.5d0
          lb(i) = 0.d0
          ub(i) = 1.d0
10    continue
c
      return
      end
c
c----------------------------------------------------------------------
c function evaluation
c----------------------------------------------------------------------
c
      subroutine sfun (n ,x, f, g, active, k)
c
c evaluate nonlinear function and gradient
c
      implicit         double precision (a-h,o-z)
      double precision x(n), f, g(n)
      double precision lb(100), ub(100), mu
      common /bnd/     lb, ub
      common /bar/     mu
      logical          active
c
      if (.not. active) return
c
      f = 0.d0
      do 30 i = 1,n
          b = 1.d0
          if (i .gt. n/2) b = -b
          d1 = x(i) - lb(i)
          d2 = ub(i) - x(i)
          if (d1 .lt. ep) write (*,800) d1
          if (d2 .lt. ep) write (*,800) d2
          f = f + b * x(i)
     *       - mu * (log (d1)) - mu * (log (d2))
          g(i) =  b - mu / d1 + mu / d2
30    continue
c
      return
800   format (' Warning: negative slack value ', d16.8, /,
     *        ' Parameter mu may be too small')
      end
c********************************************************************
c maxstp
c computes maximum feasible step in a given direction d
c under box constraints
c********************************************************************
c
      subroutine maxstp (n, x, d, stepmx)
c
c Computes maximum stepsize moving from the (feasible) point x
c in direction d, under box constraints
c
c Parameters:
c n       --> integer, number of variables
c x       --> double precision, size n, current vector of parameters
c d       --> double precision, size n, search direction vector
c stepmx <--  double precision, maximum step allowed
c
      implicit         double precision (a-h, o-z)
      double precision d(*), x(*), lb(100), ub(100)
      common /bnd/     lb, ub
c
c  set up
c
      alpha   = 1.d8
      t       = alpha
      ep      = 1.d-8
c
c find distance to nearest bound
c (this is done on all processors, but could be done in parallel)
c
      do 10 i = 1,n
          if (d(i) .gt. ep) then
              t = (ub(i) - x(i)) / d(i)
           else if (d(i) .lt. -ep) then
              t = (x(i) - lb(i)) /(-d(i))
           endif
           if (t .lt. alpha) alpha = t
10    continue
      stepmx = alpha*0.99999
c
      return
      end
