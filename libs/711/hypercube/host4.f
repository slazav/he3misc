c********************************************************************
c BTN: Sample problem 4 (using BTN to solve a constrained problem)
c main host program
c last changed: 7/30/90
c********************************************************************
c
c problem declarations
c
      program    host
      integer    pid
c
c set up hypercube
c
      pid = 0
      call setpid   (pid)
      call killcube (-1, pid)
      call load     ('node', -1, pid)
c
      stop
      end
