c********************************************************************
c BTN: CHKDER, parallel derivative checker for nonlinear functions
c last changed: 09/27/91
c********************************************************************
c
      subroutine chkez (n, x, f, g, w, lw, iw, sfun, iflag, errmax,
     *                  imax)
      implicit         double precision (a-h,o-z)
      double precision x(*), g(*), w(*)
      integer          iw(*)
      external         sfun
c
c Check derivatives of objective function via finite-differencing.
c Easy to use version
c
c PARAMETERS
c n       -> integer, dimension of problem
c x       -> double precision, size n, point where derivatives are checked
c f      <-  double precision, value of function at the point x (from sfun)
c g      <-  double precision, size n, gradient at the point x (from sfun)
c w       -  double precision, size n, work array
c lw      -> integer, declared length of array w [not used]
c iw      -> integer, size 2k, integer work array (k = # of processors
c            on hypercube)
c sfun    -> subroutine sfun (n,x,f,g) to evaluate f(x), g(x);
c            see subroutine BTNEZ for details. Routine SFUN must
c            be declared EXTERNAL in the calling program.
c iflag  <-  integer, error code:   0 => gradient values appear accurate
c                                 100 => gradient error > .01
c errmax <-  double precision, size of largest error in the gradient
c imax   <-  integer, index where errmax was obtained
c
c set default values of customizing parameters
c
      call chkpar (n, kmax, msglvl, iunit, istart, iend,
     *      istep)
      nodemx  = numnodes()
c
c check selected gradient values
c
      call chkder (n, x, f, g, w, lw, iw, sfun, iflag,
     *             nodemx, kmax, errmax, imax,
     *             msglvl, iunit, istart, iend, istep)
      return
      end
c
c
      subroutine chkpar (n, kmax, msglvl, iunit, istart, iend,
     *                   istep)
c
c
c  sets parameters for  derivative checker
c
c n         -> integer, number of variables
c kmax     <-  integer, block size
c              default:   kmax = numnodes()
c msglvl   <-  integer,  amount of printing desired,
c              (<-1 => none
c                -1 => warning messages only
c                 0 => one line with maximum error
c                 1 => individual components are printed.
c                      Processor i will print to iunit + i
c                      (if unit is not already opened, then
c                      output will be sent to file chkout.<iunit+i>)
c              default:   msglvl = 0
c iunit    <-  see discussion for msglvl
c              default:   iunit = 10
c istart   <-  integer, index of the first partial derivative
c              to be checked
c              default:  istart = 1
c iend     <-  integer, index of the last partial derivative to be
c              checked
c              default:  iend = istep*kmax
c istep    <-  integer, gradient is tested every istep components.
c              default: istep = n/kmax
c              Warning: if n is large or if gradient values are
c                       expensive, checking the gradient values can
c                       be very time consuming.  The default values
c                       result in each processor checking exactly
c                       one gradient component.
c
      kmax    = numnodes()
      msglvl  = 0
      iunit   = 10
      istart  = 1
      istep   = n/kmax
      if (istep .le. 0) istep = 1
      iend    = istart + istep*kmax - 1
      if (iend .gt. n) iend = n
c
      return
      end
c
c
      subroutine chkder (n, x, f, g, w, lw, rowinf, sfun, iflag,
     *                   nproc, kmax, errmax, imax,
     *                   msglvl, iunit, istart, iend, istep)
c
c Check derivatives of objective function via finite-differencing.
c
c PARAMETERS
c n       -> integer, dimension of problem
c x       -> double precision, size n, point where derivatives are checked
c f      <-  double precision, value of function at the point x (from sfun)
c g      <-  double precision, size n, gradient at the point x (from sfun)
c w       -  double precision, size n, work array
c lw      -> integer, declared length of array w [not used]
c rowinf <-  integer, size (nproc,2), information on processor allocation
c sfun    -> subroutine sfun (n,x,f,g,a,k) to evaluate f(x), g(x);
c            see subroutine BTNEZ for details. [Note: if more than one
c            processor is being used to compute each function value,
c            then an 'inactive' processor # (id) will call sfun at the
c            same time as processor # (mod(id,kmax)).] Routine sfun must
c            be declared EXTERNAL in the calling program.
c iflag  <-  integer, error code:   0 => gradient values appear accurate
c                                 100 => gradient error > .01
c nproc   -> integer, maximum number of processors allowed (also the
c            leading dimension of rowinf); need not equal actual number
c            of processors
c kmax    -> integer, block size [normally equal to nproc, unless more
c            than one processor is being used to compute each function
c            value; see subroutine BTN for further details]
c errmax <-  double precision, size of largest error in the gradient
c imax   <-  integer, index where errmax was obtained
c msglvl  -> integer, specify amount of printing desired:
c                   <0 => none
c                  <-1 => none
c                   -1 => warning messages only
c                    1 => individual components to unit (iunit+id)
c                         (if unit is not already opened, then output
c                         will be sent to file chkout.<iunit+id>)
c            where id is the processor #
c iunit   -> integer, see discussion of parameter msglvl
c istart  -> integer, test gradient starting at component istart
c iend    -> integer, test gradient starting at component iend
c istep   -> integer, test gradient at every istep component
c                    e.g., if istep=2 then test every 2nd component
c
c REQUIRES
c d1mach  -  function to specify machine constants
c
      implicit         double precision (a-h,o-z)
      double precision x(*), g(*), w(*), ww(2), wx(2), f
      integer          rowinf(nproc,2), cubesz
      logical          active, msg0, msg1, inuse
      character        outfle*10
      external         fcomp1
c
c Setup: obtain cube and processor information
c
      iflag  = 0
      id     = mynode()
      cubesz = numnodes()
      idpsze = 8
      eps = d1mach(3)
      if (1.d0 + 2.d0*eps .eq. 1.d0) then
          if (id .eq. 0) then
              write (*,*) ' ### CHKDER: ERROR ###'
              write (*,*) '     Machine epsilon too small'
              write (*,*) '     Value = ', eps
              write (*,*) '     Modify routine D1MACH'
              write (*,*) '     Terminating execution'
          end if
          stop
      end if
      h      = sqrt(eps)
      errmax = 0.d0
      imax   = 0
      msg0   = msglvl .ge. 0
      msg1   = msglvl .ge. 1
      iu     = iunit + id
      inquire (iu, opened = inuse)
      if (msg1 .and. .not. inuse) then
          if (id .gt. 99) then
              write (outfle,830) id
          else if (id .gt. 9) then
              write (outfle,840) id
          else
              write (outfle,850) id
          end if
          open (iu, file = outfle, status = 'unknown')
      end if
      if (msg1) write (iu,800)
c
c compare n with blocksize
c
      if (n .lt. kmax) then
          if (id .eq. 0 .and. msglvl .ge. -1) then
              write (*,*) ' CHKDER - Warning: blocksize too small'
              write (*,*) '          kmax = ', kmax
              write (*,*) '          n    = ', n
              write (*,*) '          Setting kmax = n'
          end if
          kmax = n
      end if
      if (id .lt. kmax) then
          active = .true.
      else
          active = .false.
          id0 = mod (id,kmax)
      end if
c
c get row information for this processor
c
      call getmr (n, kmax, rowinf, nproc, istart, iend, istep)
      if (active) then
          i1 = rowinf(id+1,1)
          if (id+1 .lt. kmax) then
              i2 = rowinf(id+2,1)-1
          else
              i2 = iend
          end if
c
c evaluate function and gradient at x
c
          call sfun (n, x, f, g, active, kmax)
c
c compute finite-difference estimate of gradient
c
          do 10 i = i1, i2, istep
              xi = x(i)
              x(i) = x(i) + h
              call sfun (n, x, fx, w, active, kmax)
              gi = (fx - f) / h
              erri = abs(gi - g(i)) / (1.d0 + abs(g(i)))
              if (msg1) write (iu,810) i, g(i), gi, erri
              if (erri .gt. errmax) then
                  errmax = erri
                  imax = i
              end if
              x(i) = xi
10        continue
      else
          i1 = rowinf(id0+1,1)
          if (id0+1 .lt. kmax) then
              i2 = rowinf(id0+2,1)-1
          else
              i2 = iend
          end if
          call sfun (n, x, f, g, active, kmax)
          do 15 i = i1, i2, istep
              call sfun (n, x, fx, w, active, kmax)
15        continue
      end if
c
c determine maximum error (and its index) using the subroutine GOPF
c with subroutine FCOMP1 used to make the comparison
c
      if (cubesz .gt. 1) then
          wx(1) = errmax
          wx(2) = imax
          nw = 2 * idpsze
          call gopf (wx, nw, ww, fcomp1)
          errmax  = wx(1)
          imax    = wx(2)
      endif
      if (errmax .gt. 1.d-2) iflag = 100
      if (msg0 .and. id .eq. 0) write (*,820) errmax, imax
      if (msg1) write (iu,820) errmax, imax
c
      if (.not. inuse) close (iu)
      return
800   format (' CHKDER:  Testing Derivatives', //,
     *        '   i', 9x, 'g(i)', 14x, 'gi', 13x, 'error')
810   format (' ', i3, 2x, d16.8, 2x, d16.8, 2x, d12.4)
820   format (/, ' CHKDER:  Max. error in gradient = ', d12.4, /
     *           '          observed at component ', i6)
830   format ('chkout.',   i3)
840   format ('chkout.0',  i2)
850   format ('chkout.00', i1)
      end
c
c
      subroutine fcomp1 (wx, ww)
c
c Comparison subroutine used by GOPF to find maximum error
c See subroutine CHKDER
c
      double precision wx(*), ww(*)
c
      if (ww(1) .gt. wx(1)) then
          wx(1) = ww(1)
          wx(2) = ww(2)
      end if
c
      return
      end
