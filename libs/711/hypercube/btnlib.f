c********************************************************************
c BTN: parallel software for unconstrained optimization
c last changed: 09/30/91
c********************************************************************
c
      subroutine btnez (n, x, f, g, w, lw, iw, sfun, iflag)
c
c parallel block truncated-Newton routine for unconstrained optimization
c
c Parameters
c n        -> integer, number of variables
c x       <-> double precision, size n, nonlinear variables (initial guess
c                 and solution)
c f       <-  double precision, nonlinear function value (output)
c g       <-  double precision, size n, gradient (output)
c w        -  double precision, size lw, work space; currently must be of
c                 length at least 10n + 9k + 4k^2 + 2(n/k + 1),
c                 where k = # of processors on hypercube
c lw       -> integer, length of w
c iw       -  integer, size 3k, array for storing index information for
c                 work arrays (k = # of processors on hypercube)
c sfun     -  subroutine to evaluation nonlinear function:
c                 subroutine sfun (n, x, f, g, active, kblock)
c             Routine sfun must be declared EXTERNAL in the calling program.
c             The parameters n, x, f, and g are as above; x and n must
c             not be modified by routine sfun; sfun must return the
c             nonlinear function value in f and the gradient vector in g.
c             The logical variable active indicates if the processor is
c             working (active = .true.) or idle (active = .false.).  When
c             BTN is used via BTNEZ, then routine sfun can begin with the line
c                     if (.not. active) return
c             since in this case idle processors are ignored. [For more
c             detailed information, see routine BTN.] Parameter kblock is the
c             current block size (this parameter can be ignored when using
c             BTNEZ).
c iflag   <-  integer, error code:
c                 0 => normal return: successful termination
c                 1 => linesearch failed to find lower point
c                 2 => search direction points uphill
c                 3 => function appears to be unbounded below
c                 4 => more than maxfun evaluations of f(x) in linesearch
c                 5 => not converged after nlmax iterations
c                 6 => n <= 0
c               900 => insufficient work space provided
c               999 => disaster (see file btnout.iunit for message)
c             See the user manual for more information about the meaning
c             of these error codes.
c
      implicit         double precision (a-h,o-z)
      double precision x(*), g(*), w(*)
      integer          iw(*), pid
      external         sfun
c
c set parameters for BTN
c
      call btnpar (nodemx, kmax, maxit, msglvl, pid, iprec, nlmax,
     *   initv, tolq, iunit, rnktol, maxfun, accrcy, stepmx, eta,
     *   ireset, indstp)
c
c solve optimization problem
c
      call btn (n, x, f, g, w, lw, iw, sfun, iflag,
     *          nodemx, kmax, maxit, msglvl, pid, iprec, nlmax,
     *          initv, tolq, iunit, rnktol, maxfun, accrcy,
     *          stepmx, eta, ireset, indstp, sfun)
c
      return
      end
c
c
      subroutine btnpar (nodemx, kmax, maxit, msglvl, pid, iprec,
     *    nlmax, initv, tolq, iunit, rnktol, maxfun, accrcy, stepmx,
     *    eta, ireset, indstp)
      implicit double precision (a-h,o-z)
      integer  pid
c
c Sets parameters for BTN
c
c Parameters:
c
c nodemx   integer, size of hypercube
c kmax     integer, block size
c          default:  kmax = numnodes()
c maxit    integer, maximum number of inner iterations allowed
c          default:  maxit = 20
c msglvl   integer,  amount of printing desired
c          (<-1: none,
c            -1: warning messages only
c             0: one line per outer iteration).
c          default:  msglvl = 0
c pid      integer, process i.d. number (used for message passing).
c          Unless multiple processes are being run simultaneously, pid
c          can be set to 0.
c          default:  pid = 0
c iprec    integer, type of preconditioner
c          (0 = none, 1 = BFGS, 2 = approx. BFGS)
c          default:  iprec = 2.
c nlmax    integer, maximum number of outer iterations allowed
c          default:  nlmax = 500
c initv    integer, type of initialization
c           (1 = random,2 = rhs + random, 3 = limited memory)
c          default:  initv = 3
c tolq     double precision, tolerance for quadratic-based test
c          default:  tolq =  5.d-1
c iunit    integer, output unit # for printing from this processor
c          (if this unit has not already been opened, messages go to
c          file btnout.<iunit>)
c rnktol   double precision, tolerance for rank test in inner iteration
c          default:  rnktol = 1.d-9]
c maxfun   integer, maximum number of function evaluations allowed in
c          the linesearch.
c          default:  maxfun = 2*nlmax
c accrcy   double precision, user-chosen convergence tolerance (stop the
c          algorithm if the infinity-norm of the gradient is  <=
c          accrcy*(1.d0 + |f(x)|).
c          default:  accrcy = 0.d0
c stepmx   double precision, maximum step allowed in line search
c          default:  stepmx = 1.d3
c eta      double precision, accuracy parameter for line search (must be
c          strictly between 0 and 1, and it is recommended that
c          it be chosen between .1 and .9)
c          default:  eta = 0.2
c ireset   integer, reset preconditioner after ireset iterations.
c          If ireset <= 0, then no resetting is done.
c          default:  ireset = 0 (no resetting)
c indstp   integer, indicates if the user wishes to control the
c          maximum step size in the line search with a user provided
c          subroutine MAXSTP.  If indstp=0, BTN automatically determines
c          the maximum step.  If indstp<>0, then subroutine MAXSTP must
c          be provided.
c          default:  indstp = 0 (automatic setting of step size)
c
      nodemx  = numnodes()
      kmax    = numnodes()
      maxit   = 20
      msglvl  = 0
      pid     = 0
      iprec   = 2
      nlmax   = 500
      initv   = 3
      tolq    = 5.0d-1
      iunit   = 10
      rnktol  = 1.0d-9
      maxfun  = 2*nlmax
      accrcy  = 0.d0
      stepmx  = 1.d3
      eta     = 0.2d0
      ireset  = 0
      indstp  = 0
c
      return
      end
c
c
      subroutine btn (n, x, f, g, w, lw, rowinf, sfun, iflag,
     *            nodemx, kmax, maxit, msglvl, pid, iprec, nlmax,
     *            initv, tolq, iunit, rnktol, maxfun, accrcy,
     *            stepmx, eta, ireset, indstp, maxstp)
c
c parallel block truncated-Newton routine
c
c Parameters
c n        -> integer, number of variables
c x        -> double precision, size n, nonlinear variables (initial guess
c                  and solution)
c f       <-  double precision, nonlinear function value (output)
c g       <-  double precision, size n, gradient (output)
c w        -  double precision, size lw, work space; currently must be of
c                 length at least 10n + 9kmax + 4kmax^2 + 2(n/kmax + 1);
c lw       -> integer, length of w
c rowinf   -  integer, size nodemx*3, array for storing index
c                 information for work arrays
c sfun     -  subroutine to evaluation nonlinear function:
c                 subroutine sfun (n, x, f, g, active, kblock)
c             Routine sfun must be declared EXTERNAL in the calling program.
c             The parameters n, x, f, and g are as above; x and n must
c             not be modified by routine sfun; sfun must return the
c             nonlinear function value in f and the gradient vector in g.
c             The logical variable active indicates if the processor is
c             working (active = .true.) or idle (active = .false.).  Under
c             normal circumstances (i.e., if each function evaluation is
c             performed using only one processor) then routine sfun can
c             begin with the line
c                     if (.not. active) return
c             since in this case idle processors are ignored.
c
c             If each function evaluation is being performed in parallel by
c             several processors, then active, kblock (the current block size),
c             and kmax can be used to determine what each processor should
c             do.  In this case, processors 0, 1, ..., (kblock-1) will call
c             sfun with active=.true., and all other processors will use
c             active=.false.  An idle processor # (id) will call sfun
c             at the same time as processor # (mod(id,kmax)).  Only the
c             active processor will be passed the vector x, and only
c             this processor need return the values of f and g.  An example
c             illustrating this idea is provided in the user manual.
c iflag   <-  integer, error code:
c                 0 => ok
c                 1 => linesearch failed to find lower point
c                 2 => search direction points uphill
c                 3 => function appears to be unbounded below
c                 4 => more than maxfun evaluations of f(x)
c                 5 => not converged after nlmax iterations
c                 6 => n <= 0
c               900 => insufficient work space provided
c               999 => disaster (see file btnout.iunit for message)
c            See the user manual for more information about the meaning
c            of these error codes.
c nodemx   -> integer, size of hypercube
c kmax     -> integer, block size [usually equal to nodemx]
c maxit    -> integer, maximum number of inner iterations allowed
c                 [typical value: 20]
c msglvl   -> integer, amount of printing desired (<0 = none;
c                 0 = 1 line per outer iteration)
c pid      -> integer, process id number [typical value: 0]
c iprec    -> integer, type of preconditioner (0 = none, 1 = BFGS,
c                 2 = approx. BFGS) [it is strongly suggested that
c                 iprec = 2 be used]
c nlmax    -> integer, maximum number of outer iterations allowed
c                 [typical value: 500]
c initv    -> integer, type of initialization (1 = random,
c                 (2 = rhs + random, 3 = limited memory) [it is
c                 suggested that initv = 3 be used]
c tolq     -> double precision, tolerance for quadratic-based test
c                 [typical value: 5.d-1]
c iunit    -> integer, output unit # for printing from this processor
c                 (if this unit has not already been opened, output
c                 will go to file btnout.<iunit>); iunit must be
c                 less than 1000
c rnktol   -> double precision, tolerance for rank test in inner iteration
c                 [typical value: 1.d-9]
c maxfun   -> integer, maximum number of function evaluations allowed in
c                 the linesearch [typical value: 2*nlmax]
c accrcy   -> double precision, user-chosen convergence tolerance (stop the
c                 algorithm if the infinity-norm of the gradient is  <=
c                 accrcy*(1.d0 + |f(x)|).  Under normal circumstances,
c                 the user should set accrcy=0.0 and let the stringent
c                 convergence test built-in to BTN be used.  If the function
c                 f(x) or its gradient cannot be obtained accurately, then
c                 the user may wish to override this test by setting
c                 accrcy = 1.d-6 (say), or perhaps some slightly larger value.
c                 [typical value: 0.d0]
c stepmx   -> double precision, maximum step allowed in line search
c                 [typical value: 1.d3]
c eta      -> double precision, accuracy parameter for line search (must be
c                 strictly between 0 and 1, and it is recommended that
c                 it be chosen between .1 and .9) [typical value: 0.2]
c ireset   -> integer, reset preconditioner after ireset iterations; if
c                 no resetting is desired, set ireset=0; normally it
c                 will not be necessary to reset the preconditioner based
c                 on the iteration count, but it can occasionally improve
c                 performance on certain specially structured problems; see
c                 the user manual for further guidance.
c                 [recommended value: 0]
c indstp   -> integer, indicates if the user wishes to select the maximum
c                 step in the line search manually at every iteration.
c                 See sample main program 4 for an example.
c                 If indstp=0, then BTN sets the maximum step automatically.
c                 If indstp<>0, then the user-supplied routine MAXSTP is
c                 called (see below).
c                 [typical value: 0]
c maxstp   -  (optional) subroutine to set the maximum step for the line
c                 search.  It must have the following calling sequence:
c                     subroutine maxstp (n, x, d, stepmx)
c                 where n and x are as for subroutine sfun, d is the current
c                 search direction (the new vector of parameters will be of
c                 the form x+a*d for some scalar a), and stepmx is the largest
c                 allowable step (computed by the user within subroutine
c                 maxstp).  The values of n, x, and d must not be modified by
c                 the user, and stepmx must be defined by the user.  Subroutine
c                 maxstp must be declared EXTERNAL in the calling program.
c                 NOTE: most applications will not require this feature, and
c                 should set indstp=0 when calling BTN (see above); when
c                 indstp=0, BTN will not attempt to call maxstp, and no
c                 subroutine need be provided.  In this case, the call to BTN
c                 can have the form:  call btn (..., indstp, sfun), where sfun
c                 is the name of the function evaluation routine described
c                 above.  See the sample main programs for further examples.
c                 [typical value: routine maxstp not provided]
c
      implicit         double precision (a-h,o-z)
      double precision x(*), g(*), w(*)
      integer          rowinf(nodemx,3), pid, idpsze, type0,
     *                 flgbcg, cubesz
      logical          active, msg, idzero, inuse
      character        outfle*10
      external         sfun
c
c get cube information
c
      cubesz   = numnodes()
      id       = mynode ()
c
c check for inappropriate parameter settings
c
      if (maxit .lt. 0) then
          if (id .eq. 0 .and. msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: maxit < 0'
              write (*,*) '       maxit = ', maxit
              write (*,*) '       Resetting maxit = 20'
          end if
          maxit = 20
      end if
      if (iprec .lt. 0 .or. iprec .gt. 2) then
          if (id .eq. 0 .and. msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: iprec out of range'
              write (*,*) '       iprec = ', iprec
              write (*,*) '       Resetting iprec = 2'
          end if
          iprec = 2
      end if
      if (initv .lt. 1 .or. initv .gt. 3) then
          if (id .eq. 0 .and. msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: initv out of range'
              write (*,*) '       initv = ', initv
              write (*,*) '       Resetting initv = 3'
          end if
          initv = 3
      end if
      if (tolq .le. 0.d0 .or. tolq .ge. 1.d0) then
          if (id .eq. 0 .and. msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: tolq out of range'
              write (*,*) '       tolq = ', tolq
              write (*,*) '       Resetting tolq = 0.5'
          end if
          tolq = 5.d-1
      end if
      if (rnktol .lt. 0.d0) then
          if (id .eq. 0 .and. msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: rnktol out of range'
              write (*,*) '       rnktol = ', rnktol
              write (*,*) '       Resetting rnktol = 1.d-9'
          end if
          rnktol = 1.d-9
      end if
      if (eta .le. 0.d0 .or. eta .ge. 1.d0) then
          if (id .eq. 0 .and. msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: eta out of range'
              write (*,*) '       eta = ', eta
              write (*,*) '       Resetting eta = 0.2'
          end if
          eta = 2.d-1
      end if
c
c set up
c
      inquire (iunit, opened = inuse)
      if (.not. inuse) then
          if (iunit .gt. 99) then
              write (outfle,830) iunit
          else if (iunit .gt. 9) then
              write (outfle,840) iunit
          else
              write (outfle,850) iunit
          end if
          open (iunit, file = outfle, status = 'unknown')
      end if
      if (n .le. 0) then
          iflag = 6
          go to 500
      end if
      mmax     = n/cubesz + 1
      maxit1   = maxit
      idpsze   = 8
      idzero   = id .eq. 0
c
c compare n with blocksize
c
      if (n .lt. kmax) then
          if (id .eq. 0 .and. msglvl .ge. -1) then
              write (*,*) ' BTN - Warning: blocksize too small'
              write (*,*) '       kmax = ', kmax
              write (*,*) '       n    = ', n
              write (*,*) '       Setting kmax = n'
          end if
          kmax = n
      end if
      itsave   = 0
      type0    = 11002
      kxk      = kmax * kmax
      k        = kmax
      nk       = n + k
      nf       = 1
c
c subdivide the work array
c
      lp      = 1
      lpre    = lp     + n
      ivinit  = lpre   + mmax
      iv0     = ivinit + nk
      iv1     = iv0    + nk
      ir      = iv1    + nk
      iar     = ir     + nk
      ip      = iar    + nk
      ialpha  = ip     + nk
      ibeta   = ialpha + kxk
      il1     = ibeta  + kxk
      il2     = il1    + kxk
      id1     = il2    + kxk
      id2     = id1    + kmax
      is      = id2    + kmax
      iwk     = is     + kmax
      iwk1    = iwk    + n
      iwk2    = iwk1   + n
      lprenw  = iwk2   + n
      lwtest  = lprenw + mmax - 1
      if (lw .lt. lwtest) then
          call disast ('BTN', 'INSUFFICIENT STORAGE', iunit, iflag)
          write (iunit,*) ' # of locations required = ', lwtest
          iflag = 900
          go to 500
      end if
c
      msg     = msglvl .ge.  0
      ncg     = 0
      iter    = 0
      iflag   = 0
      flgbcg  = 0
      intflg  = 0
c
c determine the columns under control of each node
c    m defines the number of rows handled by each processor
c    myrows is the starting location of rows handled.
c    rowinf(i,1) and rowinf(i,2) define myrows and m for processor i
c
      call getmr (n, cubesz, rowinf, nodemx, 1, n, 1)
      myrows = rowinf(id+1,1)
      m      = rowinf(id+1,2)
      mm     = m
c
c initial printouts
c
      if (id .lt. kmax) then
          active = .true.
      else
          active = .false.
      end if
      call sfun (n, x, f, g, active, kmax)
      fsave = abs(f)
      if (kmax .lt. cubesz) then
          itest = cubesz - kmax - 1
          if (id .le. itest) then
              idtmp = id + kmax
              call csend (type0  , g, n*idpsze, idtmp, pid)
              call csend (type0+1, f,   idpsze, idtmp, pid)
          else if (id .ge. kmax) then
              call crecv (type0  , g, n*idpsze)
              call crecv (type0+1, f,   idpsze)
          end if
          type0 = type0 + 2
      end if
      gnrm = pdnrmi (m, g(myrows), 1,
     *               cubesz, id, type0, rowinf(1,3), pid)
      if (idzero .and. msg) then
          write (*,800)
          write (*,810) iter, nf, ncg, f, gnrm
      end if
c
c set tolerances for convergence tests
c
      call settol (tol0,tol1,tol2,tol3,tol4,id)
c
c convergence test at initial point
c
      ftest = 1.d0 + dabs(f)
      if (gnrm .lt. tol0*ftest) go to 500
      type0   = 10003
c
c initialize preconditioner to the identity
c
      if (iprec .gt. 0) call dfill (m, 1.d0, w(lprenw), 1)
c
c  main loop
c
      do 10 iter = 1,nlmax
c
c  set the preconditioner
c
          itsave = itsave + 1
          if (iprec .eq. 0) then
              call dfill (m, 1.d0, w(lpre), 1)
          else
              call dcopy (m, w(lprenw), 1, w(lpre), 1)
          end if
c
c  get search direction
c  direction resides on w(lp) ... w(lp+n-1)
c
          call precg (w(iv1), m, m, myrows, g, kmax, iter,
     *        initv, flgbcg)
          call bcg (n, w(lp), g, gnrm, k, kmax, flgbcg,
     *        maxit1, x, nbl, w(lpre), w(lprenw), pid,
     *        tolq, mm, kmax, w(iv0),  w(iv1),  w(ir),
     *        w(iar), w(ip), w(ialpha),  w(ibeta),  w(il1),
     *        w(il2), w(id1),  w(id2),  w(is), w(iwk), w(iwk1),
     *        w(iwk2), rowinf, nodemx, myrows, m, iprec,
     *        rowinf(1,3), sfun, type0, iunit, rnktol)
          k = kmax
          if (iflag .eq. 999) go to 500
          call dneg (m, w(lp+myrows-1), 1)
          dg = pddot (m, w(lp+myrows-1), 1, g(myrows), 1,
     *        cubesz, id, type0, rowinf(1,3), pid)
          ncg = ncg + nbl
          if (initv .eq. 3 .and.
     *        flgbcg .gt. 0 .and. intflg .eq. flgbcg)
     *        initv = 2
          call postcg (w(iv1), m, n, m, myrows, k, iter, initv,
     *        w(ivinit), m, w(lp), g, gnrm, flgbcg, id, cubesz,
     *        rowinf(1,3), type0, pid)
          if (flgbcg .ge. 0) intflg = flgbcg
c
c collect the search direction on all nodes (even those not
c involved in computing the search direction)
c
          type0 = 10001
          call xgcol (w(lp+myrows-1), m*idpsze, w(lp), n*idpsze, ncnt,
     *            cubesz, m, id, type0, n, rowinf(1,3), pid)
c
c parallel line search
c
          oldf = f
          if (indstp .ne. 0) call maxstp (n, x, w(lp), stepmx)
          call psrch (iflag, n, x, f, g, w(lp), w(iwk), w(iwk1),
     *               w(iwk2), sfun, nf, eta, stepmx, id, pid, kmax,
     *               cubesz, rowinf(1,3), type0, dg, alpha)
          if (iflag .lt. 0) then
              maxit1 = 1 + 3/kmax
          else
              maxit1 = maxit
          end if
          if (iflag .lt. 0) iflag = 0
          gnrm = pdnrmi (m, g(myrows), 1,
     *              cubesz, id, type0, rowinf(1,3), pid)
          if (idzero .and. msg) then
              write (*,810) iter, nf, ncg, f, gnrm
          end if
          if (iflag .ne. 0) go to 500
c
c test for convergence
c
          diff  = abs(oldf - f)
          ftest = 1.d0 + dabs(f)
          xnrm = pdnrmi (m, x, 1,
     *              cubesz, id, type0, rowinf(1,3), pid)
          xtest = 1.d0 + xnrm
          pnrm = pdnrmi (m, w(lp), 1,
     *              cubesz, id, type0, rowinf(1,3), pid)
c
          if (alpha*pnrm  .lt.   tol1*xtest
     *        .and. diff  .lt.   tol2*ftest
     *        .and. gnrm  .lt.   tol3*ftest) go to 500
          if (      gnrm  .lt.   tol4*ftest) go to 500
          if (      gnrm  .lt. accrcy*ftest) go to 500
          if (nf .gt. maxfun) then
              iflag = 4
              go to 500
          end if
          ftest = abs(f)
          if (ftest .le. fsave*1.d-5 .or. itsave .eq. ireset) then
              if (iprec .gt. 0) call dfill (n, 1.d0, w(lprenw), 1)
              fsave = ftest
              itsave = 0
          end if
10    continue
      iflag = 5
c
500   if (idzero .and. msg) write (*,820) iflag
      if (.not. inuse) close (iunit)
      return
800   format (/, 3x, 'it', 3x, 'nf', 2x, 'ncg', 11x, 'f', 14x, '|g|')
810   format (' ', i4, 1x, i4, 1x, i4, 2x, 1pd16.8, 2x, 1pd12.4)
820   format (' BTN terminating with error code = ', i4)
830   format ('btnout.',   i3)
840   format ('btnout.0',  i2)
850   format ('btnout.00', i1)
      end
c
c
      subroutine settol (tol0,tol1,tol2,tol3,tol4,id)
      implicit double precision (a-h,o-z)
c
c set tolerances for convergence tests
c
      eps    = d1mach(3)
      if (1.d0 + 2.d0*eps .eq. 1.d0) then
          if (id .eq. 0) then
              write (*,*) ' ### BTN: ERROR ###'
              write (*,*) '     Machine epsilon too small'
              write (*,*) '     Value = ', eps
              write (*,*) '     Modify routine D1MACH'
              write (*,*) '     Terminating execution'
          end if
          stop
      end if
      rteps  = sqrt(eps)
      rtol   = 1.d1*rteps
      rtolsq = rtol*rtol
      peps   = eps**0.6666d0
c
      tol0 = 1.d-2  * rteps
      tol1 = rtol   + rteps
      tol2 = rtolsq + eps
      tol3 = sqrt(peps)
      tol4 = 1.d-2*rtol
c
      return
      end
c
c
      subroutine bcg (n, x, b, bnrm, k, kmax, iflag, maxit,
     *                xnl, ncg, PC, PCnew, pid,
     *                tolq, mm, kk, V0, V1, R, AR, U,
     *                alpha, beta, L1, L2, D1, D2, s,
     *                wk, wk1, wk2, rowprc, nodemx, myrows, m,
     *                iprec, nodenm, sfun, type0, iunit, rnktol)
c
c Parallel block conjugate-gradient iteration (to compute a
c search direction for a truncated-Newton optimization algorithm)
c
c PARAMETERS
c n          -> integer, dimension of problem
c x         <-  double precision, size n, search direction
c b          -> double precision, size n, right-hand side vector
c bnrm       -> double precision, infinity-norm of right-hand side
c k         <-> integer, current block size
c kmax       -> integer, maximum block size
c iflag     <-  integer, flag:
c                       0 => okay
c                      -1 => steepest-descent direction
c                       i => linear depend. in col. i at iteration 1
c                     999 => disaster (message to unit iunit)
c maxit      -> integer, maximum number of inner iterations allowed
c                     (if maxit<0, then nonlinearity test is also used)
c xnl        -> double precision, size n, vector of nonlinear parameters
c ncg       <-  integer, counter for number of inner iterations
c PC        <-> double precision, size m, current and new preconditioner
c PCnew      -  double precision, size m, work array (new PC)
c pid        -> integer, process i.d. number
c tolq       -> double precision, tolerance for quadratic test
c mm         -> integer, leading dimension of m*k matrices
c kk         -> integer, leading dimension of k*k matrices
c V0         -  double precision, size m*k, work array (old Lanczos vectors)
c V1         -> double precision, size m*k, work array (current Lanczos vectors)
c R          -  double precision, size m*k, work array (preconditioned V1)
c AR         -  double precision, size m*k, work array (Hessian times R)
c U          -  double precision, size m*k, work array (inner search directions)
c alpha      -  double precision, size k*k, work array (diagonal block of V'AV)
c beta       -  double precision, size k*k, work array (off-diag block of V'AV)
c L1         -  double precision, size k*k, work array (Cholesky factor of V'AV)
c L2         -  double precision, size k*k, work array (Cholesky factor of V'AV)
c D1         -  double precision, size   k, work array (Cholesky factor of V'AV)
c D2         -  double precision, size   k, work array (Cholesky factor of V'AV)
c s          -  double precision, size   k, work array (inner step sizes)
c wk         -  double precision, size   n, work array
c wk1        -  double precision, size   n, work array
c wk2        -  double precision, size   n, work array
c rowprc     -> integer, size kmax*2, information on processor allocation
c nodemx     -> integer, maximum number of processors available
c myrows     -> integer, index of first row on this processor
c m          -> integer, number of rows on this processor
c iprec      -> integer, choice of preconditioner:
c                       0 => no preconditioner
c                       1 => LMQN diagonal PC
c                       2 => approximate LMQN diagonal PC
c nodenm     -  integer, size kmax, array used by lower-level routines
c sfun       -> external, subroutine sfun (n,x,f,g,a,k), evaluate f(x)
c type0     <-> integer, type number for message passing
c iunit      -> integer, output unit # for printing
c               (see file btnout.iunit for output)
c rnktol     -> double precision, tolerance for rank test in inner iteration
c
      implicit         double precision (a-h,o-z)
      parameter        (idpsze = 8)
      double precision b(n), xnl(n), x(n), PC(m), PCnew(m),
     *                 V0(mm,kk), V1(mm,kk), R(mm,kk), AR(mm,kk),
     *                 U(mm,kk), alpha(kk,kk), beta(kk,kk),
     *                 L1(kk,kk), L2(kk,kk), D1(kk), D2(kk), s(kk),
     *                 wk(n), wk1(n), wk2(n)
      integer          nodenm(*), rowprc(nodemx,2),
     *                 type0, kmax, pid, cubesz
      logical          active, indef
      external         sfun
c
c set up
c
      cubesz = numnodes ()
      id     = mynode ()
      kold   = k
      iflag  = 0
      info   = 0
      iend   = 0
      indef  = .false.
      xAx    = 0.d0
      nblock = n/k
      if (nblock*k .lt. n) nblock = nblock + 1
c
c**********************************************************************
c Initialization
c**********************************************************************
c
      call dfill (m, 0.d0, x(myrows), 1)
      call dfill (k, 0.d0, s, 1)
      if (bnrm .eq. 0.d0) then
          call disast ('BCG', 'bnrm = 0', iunit, iflag)
          go to 500
      end if
      call dcopy (m, b(myrows), 1, V1, 1)
      call dscal (m, 1.d0/bnrm, V1, 1)
c
c**********************************************************************
c Main loop
c**********************************************************************
c
      maxnbl = nblock
      if (maxit  .eq.     0) maxit  = 1
      if (maxnbl .gt. maxit) maxnbl = maxit
c
      do 100 nbl = 1, maxnbl
          ncg  = nbl
          kold = k
c
c**********************************************************************
c Get QR of V1
c**********************************************************************
c
          call msolve (V1, mm, m, k, R, mm, PC)
          call pgetch (V1, mm, R, mm, beta, kk, m, cubesz, kold, id,
     *           L1, kk, k, type0, nodenm, iunit, iflag, rnktol, pid)
          if (iflag .eq. 999) return
          if (nbl .eq. 1) then
              if (k .eq. 0) then
                  call disast ('BCG', 'K=0 AT ITER=1', iunit, iflag)
                  iflag = -1
                  call dcopy  (m, b(myrows), 1, x(myrows), 1)
                  return
              else
                  if (k. lt. kold) iflag = k+1
                  kold = k
              endif
          else
              if (k .eq. 0) go to 500
          endif
c
c Compute V1 x (Beta')^(-1)
c
          call vltinv (V1, mm, m, k, V1, mm, beta, kk, iunit, iflag)
          if (iflag .eq. 999) return
          call msolve (V1, mm, m, k, R, mm, PC)
          if (nbl .eq. 1) then
              s(1) = pddot (m, b(myrows), 1, R(1,1), 1,
     *            cubesz, id, type0, nodenm, pid)
              s(1) = s(1) / bnrm
          end if
c
c**********************************************************************
c Compute AR, alpha
c**********************************************************************
c
          call getar(R, mm, k, AR, mm, cubesz, m, myrows, rowprc,
     *           nodemx, wk, wk1, wk2, type0, id, pid, n, xnl, b,
     *           sfun)
          call AtrBA (R, mm, k, AR, mm, alpha, kk, L1, kk,
     *           id, m, cubesz, type0, nodenm, pid)
c
c update diagonal preconditioner (if desired)
c
          if (iprec .eq. 1) then
              call frmpc (PCnew, alpha, R, AR, mm, m, kk, k,
     *            id, iunit, cubesz, wk, wk1, nodenm, type0, pid)
          endif
          if (iprec .eq. 2) then
              call frmpc1 (PCnew, alpha, R, AR, mm, m, kk, k,
     *            id, iunit, cubesz, wk, wk1, nodenm, type0, pid)
          endif
c
c**********************************************************************
c Cholesky for L1
c set D1 = D2 (kold * kold), then solve for L1 (kold * k)
c L2 * D1 * L1 = beta, where L2, D1 are of dimension kold
c**********************************************************************
c
          if (nbl .gt. 1) then
              call dcopy (kold, D2, 1, D1, 1)
              do 30 i = 1,k
                  call lsol (L2(i,i), kk, kold-i+1, beta(i,i),
     *                  L1(i,i), iunit, iflag)
                  if (iflag .eq. 999) return
                  call dvdiv (kold-i+1,L1(i,i),1,D1(i),1,L1(i,i),1)
30            continue
          end if
c
c**********************************************************************
c Cholesky for alpha
c Cholesky factorization of block tridiagonal matrix.
c In the following L2 will be the factor of the diagonal block.
c All processors compute the Cholesky simultaneously.
c**********************************************************************
c
          call matcpy (alpha, kk, L2, kk, k, k)
          if (nbl .gt. 1) then
              do 50 i = 1,k
              do 50 j = 1,i
                  do 40 l = i,kold
                      L2(i,j) = L2(i,j) - L1(l,i) * D1(l) * L1(l,j)
40                continue
                  if (j .lt. i) L2(j,i) = L2 (i,j)
50            continue
          end if
          call dpofa2 (L2, kk, k, info, id, iunit, iflag, rnktol)
          if (iflag .eq. 999) return
          if (info .ne. 0) then
              indef = .true.
              k = abs(info) - 1
              if (nbl .eq. 1 .and. k .eq. 0) then
                  call dcopy (m, b(myrows), 1, x(myrows), 1)
                  iflag = -1
                  return
              endif
              if (nbl .gt. 1 .and. k .eq. 0) go to 500
          endif
c
c get L2 and D2
c
          call getld2 (L2, kk, k, D2, iunit, iflag)
          if (iflag .eq. 999) return
c
c**********************************************************************
c Search direction U
c**********************************************************************
c compute  V1 - U L1
c     V1 (n * k), U (n * kold), L1 (kold * k) lower triangular
c     V1 - U L1 will be stored in R
c     U*L1 is stored in U (=U)
c
          if (nbl .gt. 1) then
              call timesl (U, mm, m, kold, L1, kk, k)
              call aminsb (R, mm, U, mm, R, mm, m, k)
          end if
c
c solve for U:
c    U * L2' = V1
c    V1 is n * k
c    L2 is k * k lower triangular
c
          call vltinv (R, mm, m, k, U, mm, L2, kk, iunit, iflag)
          if (iflag .eq. 999) return
c
c**********************************************************************
c Step size
c**********************************************************************
c
          if (nbl .eq. 1) then
              call lsol  (L2, kk, k, s, s, iunit, iflag)
              if (iflag .eq. 999) return
              call dvdiv (k, s, 1, D2, 1, s, 1)
          else
              call dvmul  (kold, D1, 1, s, 1, s, 1)
              call premlt (s, mm, kold, 1, L1, kk, k, wk2, n)
              call lsol   (L2, kk, k, wk2, s, iunit, iflag)
              if (iflag .eq. 999) return
              call dvdiv  (k, s, 1, D2, 1, s, 1)
              call dscal  (k, -1.d0, s, 1)
          end if
c
c**********************************************************************
c Update x
c**********************************************************************
c
          call maxpy (x(myrows), n, m, 1, 1.d0, U, mm, k, s, kk)
c
c**********************************************************************
c compute termination criterion using quadratic test
c**********************************************************************
c
          if (indef) go to 500
          if (nbl .eq. maxnbl) go to 500
          if (nbl .gt. 1 .and. k .eq. 0) go to 500
c
c quadratic test
c
          call dvmul  (k, D2, 1, s, 1, wk2, 1)
          xAx  = xAx + ddot(k, s, 1, wk2, 1)
          bx   = pddot (m, b(myrows), 1, x(myrows), 1, cubesz, id,
     *                type0, nodenm, pid)
          qnew = xAx * 0.5d0 - bx
          if (nbl .gt. 1) then
              diff = qnew - qold
              if (qnew .eq. 0.d0) then
                  call disast ('BCG','Q(P) = 0',iunit, iflag)
                  go to 500
              end if
              if ((nbl * diff / qnew) .lt. tolq) iend = 1
          end if
          qold = qnew
c
          if (iend .eq. 1) go to 500
c
c**********************************************************************
c Update for next iteration:  form new V
c**********************************************************************
c AR:= AR - V0 * beta; V0 (n * kold), beta (kold * k)
c
          if (nbl .gt. 1) then
              call timesl (V0, mm, m, kold, beta, kk, k)
              call aminsb (AR, mm, V0, mm, AR, mm, m, k)
          end if
c
c V1:= AR - V1 * alpha; update V0
c
          call matcpy (V1, mm,  V0, mm,  m, k)
          call atimsb (V0, mm, m, k, alpha, kk, k, V1, mm)
          call aminsb (AR, mm, V1, mm, V1, mm, m, k)
c
100   continue
c
500   call dscal (m, bnrm, x(myrows), 1)
      return
      end
c
c
      subroutine pgetch (V, ldv, Av, ldav, Beta, ldb, m, kmax, k,
     *             id, buf, ldbf, kr, type0, nodenm, iunit,
     *             iflag, rnktol, pid)
c
c Cholesky factorization of V'A V, where
c
c PARAMETERS
c V         -> double precision, size m*k, input matrix
c ldv       -> integer, leading dimension of V
c Av        -> double precision, size m*k, input matrix A*V
c ldav      -> integer, leading dimension of AV
c Beta     <-  double precision, size k*k, Cholesky factor
c ldb       -> integer, leading dimension of Beta
c m         -> integer, number of rows on this processor
c kmax      -> integer, number of active processors
c k         -> integer, block size (current)
c id        -> integer, i.d. number of this processor
c buf       -  double precision, size k*k, work matrix
c ldbf      -> integer, leading dimension of buf
c kr       <-  integer, rank of Beta (new block size)
c type0    <-> integer, message type for message passing
c nodenm    -  integer, size p, array for lower-level routines
c iunit     -> integer, output unit # for printing
c iflag    <-  integer, flag (=999 in case of disaster)
c rnktol    -> double precision, tolerance for rank test in inner iteration
c pid       -> integer, process i.d. number
c
      implicit    double precision (a-h,o-z)
      parameter   (izero = 0)
      dimension   V(ldv,*), Av(ldav,*), Beta(ldb,*), buf(ldbf,*)
      integer     nodenm(*), kmax, type0, pid
c
c form V'AV
c
      call AtrBA (V, ldv, k, Av, ldav, Beta, ldb, buf, ldbf,
     *            id, m, kmax, type0, nodenm, pid)
c
c find Cholesky factorization of V'AV (simultaneously on all nodes)
c and determine new block size (i.e., rank of Beta)
c
      call dpofa2 (Beta, ldb, k, info, id, iunit, iflag, rnktol)
      if (iflag .eq. 999) return
      if (info .ne. 0) then
          kr = abs(info) - 1
      else if (info .eq. 0) then
          kr = k
      endif
c
      return
      end
c
c
      subroutine AtrBA (A, lda, k, BA, ldb, C, ldc, buf, ldbf,
     *                   id, m, kmax, type0, nodenm, pid)
c
c computes C = A'* BA, where BA is B*A with B symmetric.
c A (n x k), BA (n * k) and C (k * k)
c This node controls m rows of A and the respective m rows of BA
c The result is explicitly symmetrized, and both halves are computed.
c
c PARAMETERS
c A         -> double precision, size m*k, input matrix
c lda       -> integer, leading dimension of A
c k         -> integer, block size
c BA        -> double precision, size m*k, input matrix B*A
c ldb       -> integer, leading dimension of BA
c C        <-  double precision, size k*k, result matrix
c ldc       -> integer, leading dimension of C
c buf       -  double precision, size k*k, work matrix
c ldbf      -> integer, leading dimension of buf
c id        -> integer, i.d. number of this processor
c m         -> integer, number of rows on this processor
c kmax      -> integer, number of active processors
c type0    <-> integer, type number for message passing
c nodenm    -  integer, size kmax, array used by lower-level routines
c pid       -> integer, process i.d. number
c
      implicit         double precision (a-h,o-z)
      parameter        (idpsze = 8, izero = 0)
      double precision A(lda,*), BA(ldb,*), C(ldc,*)
      double precision buf(ldbf,*)
      integer          nodenm(*), kmax, type0, pid
c
      msgl = ldc * k
c
c form piece of A'BA
c
      do 10 i = 1, k
      do 10 j = 1, k
          C(i,j) = ddot(m, A(1,i), 1, BA(1,j), 1)
10    continue
c
c symmetrize result
c
      do 15 i = 1,k
      do 15 j = 1,i-1
           C(j,i) = (C(i,j)+C(j,i))*5.d-1
           C(i,j) = C(j,i)
15    continue
c
c  get sum of C from all nodes
c
      call xgdsum (C, msgl, buf, kmax, id, type0, nodenm, pid)
c
      return
      end
c
c
      subroutine getar (R, ldr, k, AR, ldar, kmax, m, myrows, rowprc,
     *                  nodemx, wk, wk1, wk2, type0, id, pid, n,
     *                  xnl, g, sfun)
c
c getar: forms A*R, where R is a block matrix, by finite differencing
c
c PARAMETERS
c
c R             -> double precision, size m*k, original block matrix
c ldr           -> integer, leading dimension of R
c k             -> integer, block size
c AR           <-  double precision, size m*k, result: A*R (A = Hessian)
c ldar          -> integer, leading dimension of AR
c kmax          -> integer, # of active processors
c m             -> integer, # of rows in R and AR
c myrows        -> integer, initial row of block on this processor
c rowprc        -> integer, size kmax*2, blocking information
c nodemx        -> integer, leading dimension of rowprc
c wk, wk1, wk2  -  double precision, size n, work arrays
c type0         -> integer, type # for sending messages
c id            -> integer, index of this processor
c pid           -> integer, process id
c n             -> integer, # of variables
c xnl           -> double precision, current estimate of nonlinear variables
c g             -> double precision, current gradient
c
      implicit         double precision (a-h,o-z)
      parameter        (idpsze = 8)
      double precision R(ldr,*), AR(ldar,*), wk(*), wk1(*),
     *                 wk2(*), g(*), xnl(*)
      integer          rowprc(nodemx,2), kmax, pid, type0, type1
      logical          active
      external         sfun
c
c transpose R (sending one sub-column at a time)
c
      type1 = type0 + kmax
      do 10 j = 1,k
          if (j .ne. id+1) then
              call csend (type0+id, R(1,j), m*idpsze, j-1, pid)
          else
              call dcopy (m, R(1,j),1, wk(myrows),1)
          endif
10    continue
c
c receive one column of R in wk, and form one column of AR in wk1
c
      if (id .lt. k) then
          do 20 j = 1,kmax
          if (j .ne. id+1)
     *         call crecv (type0+j-1, wk(rowprc(j,1)),
     *                  rowprc(j,2)*idpsze)
20        continue
          active = .true.
          call atimes  (wk, wk1, n, xnl, g, wk2, sfun, active, k)
c
c 'untranspose' the result wk1, sending sub-columns to rest of cube
c
          do 30 j = 1,kmax
              if (j. ne. id+1) then
                  call csend (type1+id, wk1(rowprc(j,1)),
     *                    rowprc(j,2)*idpsze, j-1, pid)
              else
                  call dcopy(m, wk1(myrows), 1, AR(1,j), 1)
              endif
30        continue
      else
          active = .false.
          call atimes  (wk, wk1, n, xnl, g, wk2, sfun, active, k)
      endif
c
c store result in AR
c
      do 40 j = 1,k
          if (j .ne. id+1) call crecv(type1+j-1, AR(1,j), m*idpsze)
40    continue
      type0 = type0 + 2*kmax
c
      return
      end
c
c
      subroutine atimes (v, av, n, x, g, wk, sfun, active, k)
      implicit         double precision  (a-h,o-z)
      double precision v(*), av(*), x(*), g(*), wk(*)
      logical          active
c
c compute matrix/vector product via finite-differencing
c
c PARAMETERS
c v       -> double precision, size n, input for matrix/vector product
c av      <- double precision, size n, result of matrix/vector product (Av)
c n       -> integer, dimension of problem
c x       -> double precision, size n, current nonlinear parameter vector
c g       -> double precision, size n, current nonlinear gradient
c wk      -  double precision, size n, work array
c sfun    -> subroutine sfun (n,x,f,g,a,k) to evaluate nonlinear function
c active  -> logical, true if this processor must evaluate f(x)
c k       -> integer, the current block size
c
c REQUIRES
c d1mach -  function to specify machine constants
c BLAS   -  basic linear algebra subroutines
c
      if (.not. active) then
          call sfun (n, wk, fw, av, active, k)
          return
      end if
      eps = d1mach(3)
      h   = 1.d1*sqrt(eps)
      rh  = 1.d0 / h
      call dcopy (n, x, 1, wk, 1)
      call daxpy (n, h, v, 1, wk, 1)
      call sfun  (n, wk, fw, av, active, k)
      call daxpy (n, -1.d0, g, 1, av, 1)
      call dscal (n, rh, av, 1)
c
      return
      end
c
c
      subroutine frmpc (D, alpha, R, AR, mm, m, kk, k,
     *            id, iunit, kmax, wk, wk1,
     *            nodenm, type0, pid)
c
c Update the diagonal preconditioner, based on BFGS formula
c
c PARAMETERS
c D        <-> double precision, size m, diagonal preconditioner
c alpha     -> double precision, size k*k, matrix from routine BCG
c R         -> double precision, size m*k, matrix from routine BCG
c AR        -> double precision, size m*k, matrix from routine BCG
c mm        -> integer, leading dimension of R and AR
c m         -> integer, number of rows handled by this processor
c kk        -> integer, leading dimension of alpha
c k         -> integer, current block size
c id        -> integer, i.d. number of this processor
c iunit     -> integer, unit number for output
c kmax      -> integer, number of processors in use
c wk        -  double precision, size 1, work variable
c wk1       -  double precision, size 1, dummy argument (to match
c              with FORMM1)
c nodenm    -  integer, size kmax, work array for lower routines
c type0    <-> integer, records current message type
c pid       -> integer, process i.d. number
c
      implicit         double precision (a-h,o-z)
      double precision D(*), R(mm,*), AR(mm,*), alpha(kk,*),
     *                 wk(*), wk1(*)
      integer          nodenm(*), type0, pid
c
c update D for each vector in the block
c
      do 30 ind = 1,k
c
c form v'Dv for this node, and then produce global result
c
          vdv = 0.d0
          do 10 i = 1,m
              vdv = vdv + r(i,ind)*D(i)*r(i,ind)
10        continue
          call xgdsum (vdv, 1, wk, kmax, id, type0, nodenm, pid)
c
c update D
c
          vgv = alpha(ind,ind)
          if (abs(vdv) .le. 1.d-12) go to 30
          if (abs(vgv) .le. 1.d-12) go to 30
          do 20 i = 1,m
              t1 = D(i)*R(i,ind)
              t2 = AR(i,ind)
              D(i) = D(i) - t1*t1/vdv + t2*t2/vgv
              if (D(i) .le. 1.d-6) D(i) = 1.d0
20        continue
c
30    continue
c
      return
      end
c
c
      subroutine frmpc1 (D, alpha, R, AR, mm, m, kk, k,
     *            id, iunit, kmax, wk, wk1,
     *            nodenm, type0, pid)
c
c Update the diagonal preconditioner, BFGS (approximate) formula
c
c PARAMETERS
c D        <-> double precision, size m, diagonal preconditioner
c alpha     -> double precision, size k*k, matrix from routine BCG
c R         -> double precision, size m*k, matrix from routine BCG
c AR        -> double precision, size m*k, matrix from routine BCG
c mm        -> integer, leading dimension of R and AR
c m         -> integer, number of rows handled by this processor
c kk        -> integer, leading dimension of alpha
c k         -> integer, current block size
c id        -> integer, i.d. number of this processor
c iunit     -> integer, unit number for output
c kmax      -> integer, number of processors in use
c wk        -  double precision, size k, work array
c wk1       -  double precision, size k, work array
c nodenm    -  integer, size kmax, work array for lower routines
c type0    <-> integer, records current message type
c pid       -> integer, process i.d. number
c
      implicit         double precision (a-h,o-z)
      double precision D(*), R(mm,*), AR(mm,*), alpha(kk,*),
     *                 wk(*), wk1(*)
      integer          nodenm(*), type0, pid
c
c form v'Dv for this node, and then produce global result
c
      do 20 ind = 1,k
          vdv = 0.d0
          do 10 i = 1,m
              vdv = vdv + r(i,ind)*D(i)*r(i,ind)
10        continue
          wk(ind) = vdv
20    continue
      call xgdsum (wk, k, wk1, kmax, id, type0, nodenm, pid)
c
c update D
c
      do 40  ind = 1,k
          vgv = alpha(ind,ind)
          vdv = wk(ind)
          if (abs(vdv) .le. 1.d-12) go to 40
          if (abs(vgv) .le. 1.d-12) go to 40
          do 30 i = 1,m
              t1 = D(i)*R(i,ind)
              t2 = AR(i,ind)
              D(i) = D(i) - t1*t1/vdv + t2*t2/vgv
              if (D(i) .le. 1.d-6) D(i) = 1.d0
30        continue
c
40    continue
c
      return
      end
c
c
      subroutine precg (V, ldv, m, myrows, b, kmax, iter, initv, iflag)
c
c initialize the matrix V
c      initv = 1  col 1 = b, all others random
c      initv = 2  col 1 = b, col 2 = previous direc, all others random
c      initv = 3  col 1 = b, alternate previous direc and previous grad
c
c PARAMETERS
c V      <-> double precision, size m*kmax, initial matrix for block CG
c            iteration
c ldv     -> integer, leading dimension of V
c m       -> integer, # of rows on this processor
c myrows  -> integer, index of first row on this processor
c b       -> double precision, size n [size of problem], right-hand side
c            vector
c kmax    -> integer, # of active processors
c iter    -> integer, outer iteration number
c initv   -> integer, type of initialization desired
c iflag   -> integer, indicator flag (tests if initv=3 is failing)
c
      implicit         double precision (a-h,o-z)
      double precision V(ldv,*), b(*)
c
c Store b as first column of V
c
      call dcopy (m, b(myrows), 1, V, 1)
c
c INITV=1,2,3:  Dynamic (and random) initialization
c
      krand = 2
      if (initv .eq. 2) then
          if (iter .gt. 1) krand = 3
      else if (initv .eq. 3) then
          if (iflag .eq. -1) then
              krand = 3
              go to 20
          endif
          krand = 2*iter
          if (krand .gt. kmax) go to 40
      endif
20    do 30 j = krand,kmax
          call rancol (V, ldv, j, m, myrows)
30    continue
c
40    return
      end
c
c
      subroutine postcg (V, ldv, n, m, myrows, k, iter, initv,
     *                  Vinit, ldvi, x, g, gnrm, iflag, id, kmax,
     *                  nodenm, type0, pid)
c
c initialize the matrix V
c    initv = 2  col 1 = b, col 2 = previous direc(x), all others random
c    initv = 3  col 1 = b, alternate previous direcs(x) and grads(g)
c If directions are linearly dependent (at the first inner iteration)
c the offending column is replaced by a random column.  If this happens
c for the same column two outer iterations in a row, subroutine btn
c switches to initialization strategy 3.
c
c V        <-  double precision, size m*k, new initialization matrix
c ldv       -> integer, leading dimension of Vf
c n         -> integer, size of problem
c m         -> integer, number of rows on this processor
c myrows    -> integer, index of first row on this processor
c k         -> integer, block size
c iter      -> integer, outer iteration number
c initv     -> integer, choice of initialization scheme
c Vinit    <-> double precision, size m*k, storage of initialization vectors
c ldvi      -> integer, leading dimension of Vinit
c x         -> double precision, size n, search-direction vector
c g         -> double precision, size n, current gradient
c gnrm      -> double precision, infinity-norm of g
c iflag     -> integer, flag from BCG (indicates linear dependence)
c id        -> integer, i.d. number of this processor
c kmax      -> integer, number of active processors
c nodenm    -  integer, size kmax, array used by lower-level routines
c type0    <-> integer, type number for message passing
c pid       -> integer, process i.d. number
c
      implicit         double precision (a-h,o-z)
      double precision V(ldv,*), Vinit(ldvi,*), x(*), g(*)
      integer          nodenm(*), type0, pid
c
      if (k .eq. 1) return
      if (initv .eq. 1) return
c
c store search direction in column 2
c
      xnrm = pdnrmi (m, x(myrows), 1, kmax,
     *            id, type0, nodenm, pid)
      call dcopy (m, x(myrows), 1, V(1,2), 1)
      call dscal (m, 1.d0/xnrm, V(1,2), 1)
      if (initv .eq.  2) return
      if (iflag .eq. -1) return
      if (k     .eq.  2) return
c
c store current gradient in column 3
c
      call dcopy (m, g(myrows), 1, V(1,3), 1)
      call dscal (m, 1.d0/gnrm, V(1,3), 1)
      if (k .eq. 3) return
      if (iflag .gt. 0)
     *     call rancol (Vinit, ldvi, iflag, m, myrows)
c
c update remaining columns
c
      do 10 j = k, 4, -1
          call dcopy (m, Vinit(1,j-2), 1, Vinit(1,j), 1)
          call dcopy (m, Vinit(1,j  ), 1,     V(1,j), 1)
10    continue
c
      call dcopy (m, V(1,2), 1, Vinit(1,2), 1)
      call dcopy (m, V(1,3), 1, Vinit(1,3), 1)
c
      return
      end
c
c
      subroutine psrch (iflag, n, x, f, g, d, x1, g1, gopt,
     *            sfun, nf, eta, stepmx, id, pid, kmax, cubesz,
     *            nodenm, type0, dg, alpha)
c
c Parallel line search
c      Using Armijo convergence test based on function values only
c      See Luenberger (2nd edition), p. 212
c
c Parameters:
c iflag  <--  integer, error code:
c                 0 => okay
c                -1 => okay, but alpha <> 1
c                 1 => no acceptable point found
c                 2 => d is a direction of ascent
c                 3 => function may be unbounded below
c n       --> integer, number of variables
c x      <--> double precision, size n, current and new vector of parameters
c f      <--> double precision, current and new function value
c g      <--> double precision, size n, current and new gradient vector
c d       --> double precision, size n, search direction vector
c x1      --  double precision, size n, work array to store temporary x
c g1      --  double precision, size n, work array to store temporary g
c gopt    --  double precision, size n, work array to store temporary g
c sfun    --> subroutine to evaluate f(x) (call sfun (n,x,f,g,a,k))
c nf     <--> integer, total number of function/gradient evaluations
c eta     --> double precision, parameter for line search
c stepmx  --> double precision, maximum step allowed
c id      --> integer, id number of processor
c pid     --> integer, process id number
c kmax    --> integer, number of processors available to evaluate f(x)
c cubesz  --> integer, number of processors
c alpha  <--  double precision, final step size
c
      implicit         double precision (a-h, o-z)
      double precision d(*), g(*), g1(*), gopt(*),
     *                 x(*), x1(*), ww(3), wx(3)
      integer          nodenm(*), kmax, pid, type0, cubesz
      logical          active, aup, adown
      external         fcomp
c
c  set up
c
      id      = mynode ()
      aup     = .true.
      adown   = .true.
      idpsze  = 8
      nidpsze = n * idpsze
      alfopt  = 0.d0
      fopt    = f
      itmax   = 2 + 30/kmax
c
      iflag   = 0
      if (dg .ge. 0.d0) then
          iflag = 2
          return
      endif
      icount = 0
c
c  get maximum step and initialize alpha
c
      if (kmax .eq. 1) then
          alpha = 1.d0
          if (stepmx .lt. alpha) alpha = stepmx
          amax = alpha
          amin = amax
      else
          amax = kmax
          if (amax .gt. stepmx) amax = stepmx
          if (amax .gt. 1.d0) then
              amin = 1.d0 / dfloat(kmax)
              if (kmax .eq. 2) amin = 1.d0
          else
              amin = amax / dfloat(kmax)
          end if
          call setalf (amax, amin, alpha, kmax, id, 1)
      endif
c
c  test function values at initial points
c
20    nf = nf + 1
      if (id .lt. kmax) then
          active = .true.
          call dsvtvp (n, alpha, d, 1, x, 1, x1, 1)
          call sfun (n, x1, f1, g1, active, kmax)
      else
          active = .false.
          call sfun (n, x1, f1, g1, active, kmax)
          f1 = f + 1.d0
      end if
      ftest = f + eta * alpha * dg
c
c determine minimum function value using the subroutine gopf
c
      fmin = f
      f2   = f1
      if (f2 .gt. ftest) f2 = fmin + 1.d0
      if (cubesz .eq. 1) then
          fmin = f2
          imin = 1
      else
          wx(1) = f2
          wx(2) = alpha
          wx(3) = id + 1
          nw = 3 * idpsze
          call gopf (wx, nw, ww, fcomp)
          fmin  = wx(1)
          alpha = wx(2)
          imin  = wx(3)
      endif
      if (fmin .ge. f) then
          imin  = 0
          fmin  = f
          alpha = 0.d0
      end if
      icount = icount + 1
c
c Test for success at this iteration and failure at last iteration.
c In this case, use the previous optimal step and exit.
c
      if (fmin .ge. fopt .and. alfopt .ne. 0.d0) then
          imin  = kmax
          alpha = alfopt
          fmin  = fopt
          if (id .eq. kmax-1) call dcopy (n, gopt, 1, g1, 1)
          go to 30
      end if
c
c Test for too many iterations
c
      if (icount .ge. itmax .and. imin .eq. 0) then
          iflag = 1
          return
      endif
      if (icount .ge. itmax .and. imin .eq. kmax) then
          iflag = 3
      endif
c
c Test to see if step must be reduced
c
      if (imin .eq. 0 .and. adown) then
          aup  = .false.
          amax = 0.9d0 * amin
          amin = amax / dfloat (kmax)
          if (kmax .eq. 1) then
              amax = amin / 2.d0
              amin = amax
          end if
          call setalf (amax, amin, alpha, kmax, id, 2)
          go to 20
      endif
c
c Test to see if step must be increased
c
      if (kmax .gt. 1 .and. imin .eq. kmax
     *          .and. amax .lt. stepmx .and. aup) then
          adown = .false.
          if (fmin .lt. fopt) then
              alfopt = alpha
              fopt   = fmin
              if (id .eq. kmax-1) call dcopy (n, g1, 1, gopt, 1)
          end if
          amin = amax * 1.1d0
          amax = amin * kmax
          if (amax .gt. stepmx) amax = stepmx
          call setalf (amax, amin, alpha, kmax, id, 1)
          go to 20
      endif
c
c form the new x, and send g to all the nodes
c
30    call dsvtvp (n, alpha, d, 1, x, 1, x, 1)
      if (alpha .ne. 1.d0) iflag = -1
      f = fmin
      if (imin .eq. id + 1) then
          call dcopy (n, g1, 1, g, 1)
          do 60 j = 1, cubesz-1
              nodenm(j) = j-1
              if (j .gt. id) nodenm(j) = nodenm(j) + 1
60        continue
          call gsendx (type0, g, nidpsze, nodenm, cubesz-1)
      else
          call crecv (type0, g, nidpsze)
      endif
      type0 = type0 + 1
c
      return
      end
c
c
      subroutine setalf (amax, amin, alpha, cubesz, id, ind)
c
c determine a set of steplengths alpha for the line search
c
c Parameters
c amax    --> maximum value of alpha
c amin   <--  minimum value of alpha
c alpha  <--  value of alpha for this processor
c cubesz  --> number of processors
c id      --> id number of this processor
c ind     --> indicator: =1 (normal) =2 (shrinking step)
c
      implicit   double precision (a-h, o-z)
      integer    id, cubesz
c
c Set up
c
      if (id .ge. cubesz) then
          alpha = 0.0d0
          return
      end if
      np1 = cubesz - 1
      if (amin .ge. amax) amin = amax / dfloat(cubesz)
      if (ind .eq. 2) go to 100
c
      if (cubesz .gt. 2) then
          if (amin .lt. 1.d0 .and. amax .gt. 1.d0) then
              imid = cubesz/2
              if (id .lt. imid) then
                  alpha = amin + id * (1.d0-amin)/dfloat(imid)
              else if (id .gt. imid) then
                  alpha = 1.d0 +
     *                  (id-imid) * (amax-1.d0)/dfloat(np1-imid)
              else if (id .eq. imid) then
                  alpha = 1.d0
              end if
          else
              alpha = amin + id * (amax-amin)/dfloat(np1)
          end if
      else
          if (id .eq. 0) alpha = amin
          if (id .eq. 1) alpha = amax
      end if
      return
c
c Shrinking of step (more rapid reduction of step to zero)
c
100   if (cubesz .ge. 2) then
          alpha = amax / 2.d0**(cubesz-id-1)
          amin  = amax / 2.d0**(cubesz-1)
      else
          alpha = amax
          amin  = amax
      end if
      return
c
      end
c
c
      subroutine fcomp (wx, ww)
c
c Comparison routine used by gopf above.
c Test for smallest function value, and set related parameters
c
      double precision wx(*), ww(*)
c
      if (ww(1) .lt. wx(1)) then
          wx(1) = ww(1)
          wx(2) = ww(2)
          wx(3) = ww(3)
      else if (ww(1) .eq. wx(1) .and. ww(3) .lt. wx(3)) then
          wx(1) = ww(1)
          wx(2) = ww(2)
          wx(3) = ww(3)
      end if
c
      return
      end
c
c
      subroutine aminsb (A, lda, B, ldb, C, ldc, m, n)
      double precision A(lda,*), B(ldb,*), C(ldc,*)
c
c computes C = A - B, for m*n matrices A, B, and C
c [C can overwrite either A or B]
c
      do 10 j = 1, n
          call dvsub (m, A(1,j), 1, B(1,j), 1, C(1,j), 1)
10    continue
c
      return
      end
c
c
      subroutine atimsb (A, lda, n, m, B, ldb, k, C, ldc)
      implicit         double precision (a-h, o-z)
      double precision A(lda,*), B(ldb,*), C(ldc,*)
c
c forms C = A x B,  for A(n*m), B(m*k), and C(n*k)
c
      do 10 j = 1,k
      do 10 i = 1,n
          C(i,j) = ddot (m,A(i,1),lda,B(1,j),1)
10    continue
c
      return
      end
c
c
      subroutine disast (routne, msg, iunit, iflag)
c
c print fatal error message to unit iunit, then return
c
c PARAMETERS
c routne -> character*(*), name of routine where error was detected
c msg    -> character*(*), error message
c iunit  -> integer, unit for output
c iflag <-  integer, flag=999 upon return from disaster
c
      character*(*) routne, msg
c
      write (iunit,800) routne, msg
      iflag = 999
c
      return
800   format (' ********************',         /,
     *        ' ERROR, ERROR, ERROR',          /,
     *        ' Fatal error in routine ', a10, /,
     *        ' ', a40,                        /,
     *        ' Terminating',                  /,
     *        ' ********************')
      end
c
c
      subroutine dpofa2 (A, lda, n, info, id, iunit, iflag, rnktol)
      implicit         double precision (a-h,o-z)
      double precision A(lda,*),  t, s, tol, mnorm
      integer          info, n
c
c Cholesky factor of a double precision positive definite matrix.
c this is a modification of the LINPACK routine DPOFA, designed to
c produce the lower triangular matrix a column at a time.
c
c info:   integer = 0 for normal return
c                 = k signals an error condition - the leading minor
c                     of order k is not positive definite.
c
      tol = rnktol * mnorm (A, lda, n, n)
c
c  get LDL decomposition from a Cholesky factorization
c
      do 30 i = 1,n
          s = ddot (i-1, A(i,1), lda, A(i,1), lda)
          s = A(i,i) - s
c
          info = i
          if (s .lt. -tol) return
          if (s .le.  tol) then
              info = - info
              return
          endif
          A(i,i) = sqrt(s)
c
          do 20 k = i+1,n
              t = a(k,i)  - ddot (i-1, A(i,1), lda, A(k,1), lda)
              if (A(i,i) .eq. 0.d0) then
                  call disast ('DPOFA2', 'A(i,i) = 0', iunit, iflag)
                  return
              end if
              t = t / A(i,i)
              A(k,i)  = t
20        continue
30    continue
      info = 0
c
      return
      end
c
c
      subroutine getld2 (L, ldl, k, D, iunit, iflag)
      implicit         double precision (a-h, o-z)
      double precision L(ldl,*), D(*)
c
c  get LDL decomposition from a Cholesky factorization
c
      do 10 i = 1,k
          D(i) = L(i,i)
          if (D(i) .eq. 0.d0) then
              call disast ('GETLD2', 'D(i) = 0', iunit, iflag)
              return
          end if
          do 10 j = 1,k
              if (j .gt. i) L(i,j) = 0.d0
              if (j .le. i) L(i,j) = L(i,j) / D(j)
10    continue
      call dvmul (k, D, 1, D, 1, D, 1)
c
      return
      end
c
c
      subroutine getmr (n, p, rowprc, nodemx, istart, iend, istep)
      integer rowprc(nodemx,2), p
c
c determine the rows handled by each processor
c
c PARAMETERS
c n        -> # of variables
c p        -> # of processors
c rowprc  <-  array to store indexing information
c nodemx   -> declared dimension of rowprc
c istart   -> beginning index
c iend     -> ending index
c istep    -> stride for index
c
      nn = (iend - istart + istep) / istep
c     m = n / p
      m = nn / p
c     rowprc(1,1) = 1
      rowprc(1,1) = istart
      rowprc(1,2) = m
      mnp = mod(nn,p)
c
      do 10 i = 2,p
c         rowprc(i,1) = rowprc(i-1,1) + rowprc(i-1,2)
          rowprc(i,1) = rowprc(i-1,1) + rowprc(i-1,2)*istep
          if (mnp .ge. i-1) then
              rowprc(i,2) = m + 1
          else
              rowprc(i,2) = m
          endif
10    continue
c
      return
      end
c
c
      subroutine lsol (l, ldl, k, y, x, iunit, iflag)
      implicit         double precision (a-h, o-z)
      double precision l(ldl,*), y(*), x(*)
c
c solve Lx = y where L is lower triangular
c the vectors x and y may be the same (overwrite solution on rhs)
c
      do 10 i = 1,k
          x(i) = y(i) - ddot (i-1, l(i,1), ldl, x, 1)
          if (l(i,i) .eq. 0.d0) then
              call disast ('LSOL', 'l(i,i) = 0', iunit, iflag)
              return
          end if
          x(i) = x(i) / l(i,i)
10    continue
c
      return
      end
c
c
      subroutine matcpy (A, lda, B, ldb, m, n)
      double precision A(lda,*), B(ldb,*)
c
c copies matrix A  (m x n) onto matrix B
c
      do 10 j = 1,n
          call dcopy (m, A(1,j), 1, B(1,j), 1)
10    continue
c
      return
      end
c
c
      subroutine maxpy (y, ldy, n, k, alpha, a, lda, m, x, ldx)
      double precision y(ldy,*), a(lda,*), x(ldx,*), alpha, t
c
c form y = y + alpha * Ax
c      y      n x k
c      A      n x m
c      x      m x k
c      alpha  scalar
c
      do 10 i = 1,n
      do 10 l = 1,m
          t = alpha*A(i,l)
          call daxpy (k,t,x(l,1),ldx,y(i,1),ldy)
10    continue
c
      return
      end
c
c
      double precision function mnorm (A, lda, m, n)
      implicit double precision (a-h,o-z)
      double precision A(lda,*), colnrm, zero
      data   zero /0.d0/
c
c compute the 1-norm of the m x n matrix A
c
      mnorm = zero
      do 10 j = 1,n
          colnrm = dasum (m, A(1,j), 1)
          if (mnorm .lt. colnrm) mnorm = colnrm
10    continue
c
      return
      end
c
c
      subroutine msolve (A, lda, m, n, Am, ldam, Rm)
      double precision A(lda,*), Am(ldam,*), Rm(*)
c
c  preconditions the m*n matrix A, using the vector Rm.
c  resulting matrix is in Am [A and Am can be the same]
c
      do 10 j = 1,n
          call dvdiv (m, A(1,j), 1, Rm, 1, Am(1,j), 1)
10    continue
c
      return
      end
c
c
      double precision function pddot (m, x, lx, y, ly, p, id,
     *                  type0, nodenm, pid)
c
c Parallel inner product routine
c Each node forms the inner product of its part of x and y in parallel
c Node 0 accumulates the inner products and sends result to other nodes.
c
c PARAMETERS
c m         -> integer, size of arrays on each processor
c x         -> double precision, size m, input array 1
c lx        -> integer, stride for x (as in DDOT)
c y         -> double precision, size m, input array 2
c ly        -> integer, stride for y (as in DDOT)
c p         -> integer, number of active processors
c id        -> integer, i.d. number of this processor
c type0    <-> integer, message number for communication
c nodenm    -  integer, size p, array for lower-level routines
c pid       -> integer, process i.d. number
c
      implicit         double precision (a-h,o-z)
      double precision x(*), y(*)
      integer          nodenm(*), p, type0, pid
c
      prod  = ddot (m, x, lx, y, ly)
      call xgdsum (prod, 1, sval, p, id, type0, nodenm, pid)
      pddot = prod
c
      return
      end
c
c
      double precision function pdnrmi (m, x, lx, p, id,
     *                  type0, nodenm, pid)
c
c Parallel infinity-norm routine
c Each node computes the norm of its portion of x in parallel
c Node 0 accumulates the norms and sends result to other nodes.
c
c PARAMETERS
c m         -> integer, size of arrays on each processor
c x         -> double precision, size m, input array
c lx        -> integer, stride for x (as in DNRM2)
c p         -> integer, number of active processors
c id        -> integer, i.d. number of this processor
c type0    <-> integer, message number for communication
c nodenm    -  integer, size p, array for lower-level routines
c pid       -> integer, process i.d. number
c
      implicit         double precision (a-h,o-z)
      double precision x(*)
      integer          nodenm(*), p, type0, pid
c
      if (m .gt. 0) then
          ix   = idamax (m, x, lx)
          xnrm = abs(x(ix))
      else
          xnrm = 0.d0
      end if
      call xgdmax (xnrm, 1, sval, p, id, type0, nodenm, pid)
      pdnrmi = xnrm
c
      return
      end
c
c
      subroutine premlt (A, lda, k1, n,  L, ldl, k2, B,ldb)
      implicit         double precision (a-h,o-z)
      double precision A(lda,*), L(ldl,*), B(ldb,*)
c
c forms B =  L' x A
c for A(k1*n) and L(k1*k2) lower triangular (with k1 .ge. k2)
c
      do 10 j = 1,n
      do 10 i = 1,k2
          B(i,j) = ddot (k1-i+1,A(i,j),1,L(i,i),1)
10    continue
c
      return
      end
c
c
      subroutine rancol (V, ldv, j, m, myrows)
c
c  generate a random vector for the j-th column of the matrix V
c
c PARAMETERS
c V      <-  double precision, size m*j, initialization matrix (see
c            PREBCG)
c ldv     -> integer, leading dimension of V
c j       -> integer, column to be randomized
c m       -> integer, # of rows on this processor
c myrows  -> integer, index of first row on this processor
c
      implicit         double precision (a-h,o-z)
      double precision V(ldv,*)
      integer          click, seed, rand
c
c setup
c
      click = 52343
      icol  = j
      seed  = mod(2*(icol)*click+1,65536)
      rand  = seed
c
c increment seed
c
      do 15 i = 1, myrows-1
         rand = mod(3125*rand,65536)
15    continue
c
c generate random column
c
      do 20 i = 1, m
          rand = mod(3125*rand,65536)
          V(i,j) = (rand - 32768.0d0)/16384.0d0
20    continue
c
      return
      end
c
c
      subroutine timesl (A, lda, n, k1, L, ldl, k2)
      implicit         double precision (a-h,o-z)
      double precision A(lda,*), L(ldl,*), t
c
c performs A = A x L
c for A(n*k1) and L(k1*k2) lower triangular (with k1 .ge. k2)
c
      do 20 j = 1,k2
      do 20 i = 1,n
          t      = ddot(k1-j+1,A(i,j),lda,L(j,j),1)
          A(i,j) = t
20    continue
c
      return
      end
c
c
      subroutine vltinv (V, ldv, m, k, P, ldp, L, ldl, iunit, iflag)
      implicit         double precision (a-h,o-z)
      double precision L(ldl,*),  V(ldv,*), P(ldp,*), t
c
c Solves P L' = V, for L(k*k) lower triangular, P(m*k) and V(m*k)
c V may be overwritten
c
      do 20 i = 1,m
          do 10 j = 1,k
            t = V(i,j) - ddot (j-1,P(i,1),ldp,L(j,1),ldl)
            if (L(j,j) .eq. 0.d0) then
                call disast ('VLTINV', 'L(j,j) = 0', iunit, iflag)
                return
            end if
            P(i,j) = t / L(j,j)
10        continue
20    continue
c
      return
      end
c
c
      subroutine xgcol (x, lx, y, ly, ncnt, p, m, id, type0, n,
     *                  nodenm, pid)
c
c If p=cubesz, then call GCOL.
c Otherwise, collect a vector x in the vector y on node 0, and then
c     send y to all nodes.
c
c PARAMETERS
c x         -> double precision, size lx, input vector (one on each
c              processor)
c lx        -> integer, length of x (in bytes)
c y         -> double precision, size ly, output vector
c ly        -> integer, length of y (in bytes)
c ncnt      -  integer, dummy used by GCOL
c p         -> integer, number of active processors
c m         -> integer, dimension of x
c id        -> integer, i.d. number of this processor
c type0    <-> integer, type number for messages
c n         -> integer, dimension of y
c nodenm    -  integer, size p, array for global send
c pid       -> integer, process i.d. number
c
      implicit         double precision (a-h,o-z)
      double precision x(*), y(*)
      integer          nodenm(*), type0, p, pid, cubesz
      parameter        (izero = 0, idpsze=8)
c
      cubesz = numnodes()
      if (p .eq. cubesz) then
          call gcol (x, lx, y, ly, ncnt)
          return
      end if
c
      if (id .eq. izero) then
          call dcopy (m, x, 1, y, 1)
          mnp = mod(n,p)
          loc = m + 1
          do 20 iproc = 1, p-1
              msgl = m
              if (iproc .le. mnp) msgl = m + 1
              call crecv (type0+iproc, y(loc), msgl*idpsze)
              loc = loc + msgl
20        continue
          do 30 i = 1,cubesz-1
              nodenm(i) = i
30        continue
          ndpsize = n * idpsze
          call gsendx (type0, y, ndpsize, nodenm, cubesz-1)
      else if (id .lt. p) then
          call csend (type0+id, x, lx, izero, pid)
          call crecv (type0,  y, ly)
      else
          call crecv (type0,  y, ly)
      endif
      type0 = type0 + p
c
      return
      end
c
c
      subroutine xgdmax (x, n, work, p, id, type0, nodenm, pid)
c
c If p=cubesz, call GDHIGH.
c Otherwise, maximize the components of the vectors x on the nodes, and
c    store the result on all the nodes.
c
c PARAMETERS
c x        <-> double precision, size n, vector to be maxed, as well
c              as the result
c n         -> integer, dimension of x
c work      -  double precision, size n, work array
c p         -> integer, number of active processors
c id        -> integer, i.d. number of processor
c type0    <-> integer, type number for messages
c nodenm    -  integer, size p, array for global send
c pid       -> integer, process i.d. number
c
      implicit         double precision (a-h, o-z)
      double precision x(*), work(*)
      integer          nodenm(*), p, type0, pid, cubesz
      parameter        (izero=0, idpsze=8)
c
      cubesz = numnodes()
      if (p .eq. cubesz) then
          call gdhigh (x, n, work)
          return
      end if
c
      ndpsize = n * idpsze
      if (id .eq. izero) then
          do 10 iproc = 1, p-1
              call crecv (type0+iproc, work, ndpsize)
              call dvmax (n, x, 1, work, 1, x, 1)
10        continue
          do 20 i = 1,p-1
              nodenm(i) = i
20        continue
          call gsendx (type0, x, ndpsize, nodenm, p-1)
      else
          call csend (type0+id, x, ndpsize, izero, pid)
          call crecv (type0,  x, ndpsize)
      endif
      type0 = type0 + p
c
      return
      end
c
c
      subroutine xgdsum (x, n, work, p, id, type0, nodenm, pid)
c
c If p=cubesz, call GDSUM.
c Otherwise, sum the components of the vectors x on the nodes, and
c    store the result on all the nodes.
c
c PARAMETERS
c x        <-> double precision, size n, vector to be summed and
c              the result
c n         -> integer, dimension of x
c work      -  double precision, size n, work array
c p         -> integer, number of active processors
c id        -> integer, i.d. number of processor
c type0    <-> integer, type number for messages
c nodenm    -  integer, size p, array for global send
c pid       -> integer, process i.d. number
c
      implicit         double precision (a-h, o-z)
      double precision x(*), work(*)
      integer          nodenm(*), p, type0, pid, cubesz
      parameter        (izero=0, idpsze=8)
c
      cubesz = numnodes()
      if (p .eq. cubesz) then
          call gdsum (x, n, work)
          return
      end if
c
      ndpsize = n * idpsze
      if (id .eq. izero) then
          do 10 iproc = 1, p-1
              call crecv (type0+iproc, work, ndpsize)
              call dvadd (n, x, 1, work, 1, x, 1)
10        continue
          do 20 i = 1,p-1
              nodenm(i) = i
20        continue
          call gsendx (type0, x, ndpsize, nodenm, p-1)
      else
          call csend (type0+id, x, ndpsize, izero, pid)
          call crecv (type0,  x, ndpsize)
      endif
      type0 = type0 + p
c
      return
      end
