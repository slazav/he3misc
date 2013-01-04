c********************************************************************
c BTN: Sample problem 3 (complex usage of BTN, CHKDER)
c main node program
c last changed: 10/5/90
c********************************************************************
c
      program          node
      implicit         double precision (a-h,o-z)
      parameter        (nn = 100, ndmx = 16,
     *                 liw = 3*ndmx,
     *                 lw = 14*nn + ndmx*(4*ndmx + 9),
     *                 nmsgi = 8, lmsgi = 4*nmsgi, idpsze = 8,
     *                 nmsgr = 2, lmsgr = idpsze*nmsgr)
      double precision x(nn), g(nn), w(lw), msgr(nmsgr)
      integer          iw(liw), msgi(nmsgi), hid, pid
      logical          idzero
      external         sfun
      common /cube/    kmax, pid
c
c set up node information
c
      id     = mynode ()
      hid    = myhost ()
      idzero = id .eq. 0
c
c Receive problem specification from host processor
c
      call crecv (101, msgi, lmsgi)
c
c set default parameters for chkder and btn
c
      n      = msgi(1)
      call chkpar (n, kmax, msglvl, iunit, istart, iend, istep)
      call btnpar (nodemx, kmax, maxit, msglvl, pid, iprec, nlmax,
     *    initv, tolq, iunit, rnktol, maxfun, accrcy, stepmx, eta,
     *    ireset, indstp)
c
c modify selected parameters
c
      ichk   = msgi(2)
      istart = msgi(3)
      iend   = msgi(4)
      istep  = msgi(5)
      msglvl = msgi(6)
      iunit  = msgi(7)
      kmax   = msgi(8)
c
      call crecv (102, msgr, lmsgr)
      rnktol = msgr(1)
      accrcy = msgr(2)
c
c check derivative values (if desired)
c
      call xstart (x, n)
      if (ichk .ne. 0) then
          call chkder (n, x, f, g, w, lw, iw, sfun, iflag,
     *                 nodemx, kmax, errmax, imax,
     *                 msglvl, iunit, istart, iend, istep)
          if (idzero) then
              open  (iunit, file = 'node.out', status = 'unknown')
              write (iunit,800) errmax, imax
              close (iunit)
          end if
          stop
      end if
c
c solve optimization problem
c
      call btn (n, x, f, g, w, lw, iw, sfun, iflag,
     *          nodemx, kmax, maxit, msglvl, pid, iprec, nlmax,
     *          initv, tolq, iunit, rnktol, maxfun, accrcy,
     *          stepmx, eta, ireset, indstp, sfun)
      if (.not. idzero) stop
c
c Send results to host processor (node 0 only)
c
      call csend (103, iflag, 4, hid, pid)
      call csend (104, f, idpsze, hid, pid)
      call csend (105, x, n*idpsze, hid, pid)
      call csend (106, g, n*idpsze, hid, pid)
c
      stop
800   format (' CHKDER: Max. error in gradient = ', 1pd12.4,
     *        ' at index ', i5)
      end
c
c----------------------------------------------------------------------
c initial guess
c----------------------------------------------------------------------
c
      subroutine xstart (x, n)
      double precision x(n)
c
c specify initial values of the variables
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
c This example illustrates how the function evaluation can be
c spread over several processors (with blocksize less than cubesize).
c
c This example assumes that the cube size (cubesz) is exactly divisible
c by kmax (the maximum block size), and that n (the number of variables)
c is exactly divisible by m (the number of processors working on a
c function evaluation).
c
c The block size and the cube size can be varied.
c
c This example uses message passing.  The message types used are
c 200, 201, ..., 200+(cubesize), 401, 402, ..., 400+(cubesize)
c----------------------------------------------------------------------
c
      subroutine sfun (n ,x, f, g, active, kblock)
c
c evaluate the nonlinear function and its gradient [in parallel]
c
      parameter        (idpsze = 8)
      implicit         double precision (a-h,o-z)
      double precision x(*), g(*)
      integer          cubesz, pid
      logical          active
      common /cube/    kmax, pid
c
c get cube information, and set up indexes for this processor
c    kblock = current block size
c    id     = # of this processor
c    cubesz = # of processors in hypercube
c    itotal = # of processors working together to compute f(x)
c    icount = counter for the processors working to compute f(x)
c             [icount = 1, 2, ..., itotal]
c    nt     = # of variables handled by this processor
c    i1     = subscript of first variable handled by this processor
c    i2     = subscript of  last variable handled by this processor
c    id0    = # of processor that accumulates f(x) and g(x)
c
      id     = mynode ()
      cubesz = numnodes ()
      itotal = cubesz / kmax
      icount = id/kmax + 1
      nt     = n / itotal
      i1     = nt*(icount-1) + 1
      i2     = nt*icount
      id0    = mod (id,kmax)
c
c determine if a function evaluation is necessary
c
      if (id0 .ge. kblock) return
c
c send portion of the x vector to the other nodes computing f(x)
c
      if (id+1 .le. kmax) then
          icount = 1
          do 20 id1 = id+kmax,cubesz-1,kmax
              icount = icount + 1
              j1     = nt*(icount-1) + 1
              call csend (200,  x(j1), nt*idpsze, id1, pid)
20        continue
          icount = 1
      else
          call crecv (200,  x(i1), nt*idpsze)
      end if
c
c compute portion of function value and gradient
c
      f = 0.d0
      do 30 i = i1,i2
          b = i
          if (i .gt. n/2) b = -b
          f = f + 5.d-1*i*x(i)*x(i) - b*x(i)
          g(i) = i*x(i) - b
30    continue
c
c send partial results to main processor
c
      if (id+1 .le. kmax) then
          icount = 1
          do 40 id1 = id+kmax,cubesz-1,kmax
              icount = icount + 1
              i1 = nt*(icount-1) + 1
              call crecv (201+id1,    f1,    idpsze)
              call crecv (401+id1, g(i1), nt*idpsze)
              f = f + f1
40        continue
      else
          call csend (201+id,     f,    idpsze, id0, pid)
          call csend (401+id, g(i1), nt*idpsze, id0, pid)
      end if
c
      return
      end
