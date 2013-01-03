c********************************************************************
c BTN: Sample problem 3 (complex usage of BTN, CHKDER)
c main host program
c last changed: 10/5/90
c********************************************************************
c
c problem declarations
c
      program          host
      parameter        (nmax = 100, idpsze = 8,
     *                 nmsgi = 8, lmsgi = 4*nmsgi,
     *                 nmsgr = 2, lmsgr = idpsze*nmsgr)
      implicit         double precision (a-h,o-z)
      integer          msgi(nmsgi), pid
      double precision x(nmax), g(nmax), msgr(nmsgr)
c
c set up hypercube
c
      pid = 0
      call setpid   (pid)
      call killcube (-1, pid)
      call load     ('node', -1, pid)
c
c read data from file
c
      open (7, file = 'host3.d', status = 'old')
      read (7,800) n
      read (7,800) ichk
      read (7,800) istart
      read (7,800) iend
      read (7,800) istep
      read (7,800) msglvl
      read (7,800) iunit
      read (7,800) kblock
c
      read (7,810) rnktol
      read (7,810) accrcy
c
c Send problem specification to the cube
c
      msgi(1)  = n
      msgi(2)  = ichk
      msgi(3)  = istart
      msgi(4)  = iend
      msgi(5)  = istep
      msgi(6)  = msglvl
      msgi(7)  = iunit
      msgi(8)  = kblock
      call csend (101, msgi, lmsgi, -1, pid)
c
      msgr(1)  = rnktol
      msgr(2)  = accrcy
      call csend (102, msgr, lmsgr, -1, pid)
c
      if (ichk .ne. 0) stop
c
c receive solution
c
      call crecv (103, iflag, 4)
      call crecv (104, f, idpsze)
      call crecv (105, x, n*idpsze)
      call crecv (106, g, n*idpsze)
c
c write results to file
c
      open  (iunit, file = 'host.out', status = 'unknown')
      write (iunit,820)
      write (iunit,830) iflag
      write (iunit,840) f
      write (iunit,850)
      write (iunit,860) (x(i), i = 1,n)
      write (iunit,870)
      write (iunit,860) (g(i), i = 1,n)
      close (iunit)
c
      stop
800   format (10x, i5)
810   format (10x, d16.4)
820   format (' Results of BTN', /)
830   format (' Error code = ', i4)
840   format (' Function value = ', 1pd24.16)
850   format (/, ' x = ')
860   format (3(1pd16.8))
870   format (/, ' g = ')
      end
