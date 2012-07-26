program texture
  implicit none

  common /texture_pars/ XiH, Lambda, KappaH
  double precision XiH, Lambda, KappaH
  data XiH    /1D0/    &
       Lambda /5D0/ &
       KappaH /0.1D0/

  integer nr,nz, ir,iz
  parameter (nr=20, nz=20)

  double precision U(nr,nz,2), HS, HF
  double precision rmin,rmax,zmin,zmax
  data rmin /0D0/ &
       rmax /1D0/ &
       zmin /0D0/ &
       zmax /1D0/ &
       HS/1D-1/   &
       HF/1D-5/

  do ir=1,nr
    do iz=1,nz
      U(ir,iz,1) = 0D0;
      U(ir,iz,2) = 0D0;
    enddo
  enddo

  call VAR2D(rmin,rmax,zmin,zmax,U,nr,nz,2, HF)

  open (20, FILE='data.txt')
  do ir=1,nr
    do iz=1,nz
      write(20,*) ir, iz, U(ir,iz,1), U(ir,iz,2)
    enddo
    write(20,*)
  enddo
  close(20)

end
