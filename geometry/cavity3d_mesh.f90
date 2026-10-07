program cavity3d_mesh
!
!  Build high order meshes of the axisymmetric cavity (the 3d analogue of
!  the 2d elliptical cavity in
!  helm-high-freq-conv/paper_examples/cavity_example.f) and write them out
!  as vtk plots and, for the triangular meshes, go3 files.
!
!  See cavity3d_geom.f90 for the parametrization and the meshing strategy.
!
!  Resolution is wavelength driven: for each wavenumber zk in zks and each
!  order in norders_run, the target patch arclength is
!  h = (norder+1)*lambda/ppw.
!
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)

  real *8, allocatable :: srcvals(:,:), srccoefs(:,:), qwts(:)
  integer *8, allocatable :: norders(:), ixyzs(:), iptype(:)
  real *8, allocatable :: tchse(:)
  integer *8, allocatable :: nrts(:,:)

  real *8 pars(3)
  real *8 zks(4)
  integer *8 norders_run(2), iptypes_run(2)

  character *300 fname, dir_out
  character *12 ptag

  external funcurve_cavity3d

  call prini(6,13)

  done = 1.0d0
  pi = atan(done)*4

  dir_out = 'meshes/'

!
!  generating curve parameters: defaults reproduce the 2d cavity
!
  np = 3
  call cavity3d_defpars(np, pars)

!
!  meshing controls
!
  ppw = 10.0d0
  dthmax = pi/6
  ntbase = 8
!
!  iort = -1 gives outward normals on this geometry (exterior scattering)
!
  iort = -1

!
!  The quadratic skeleton plot (_msh.vtk) is small and is written for
!  every mesh.  The full high order surface plot and the normals plot are
!  tens to hundreds of megabytes at the finer resolutions, so they are
!  only written for the coarsest wavenumber, which is all a visual check
!  of the lip, the caps and the normal orientation needs.
!
  ifplot_full = 1

  zks(1) = 5.0d0
  zks(2) = 10.0d0
  zks(3) = 20.0d0
  zks(4) = 40.0d0
  nzk = 4

  norders_run(1) = 4
  norders_run(2) = 8
  nord = 2

  iptypes_run(1) = 1
  iptypes_run(2) = 11
  nipt = 2

!
!  reference quantities for the default parameters, from high order
!  quadrature on the generating curve
!
  arefl = 6.56986269594066d0
  aref = 29.53229570305616d0
  vref = 5.88149777703912d0

  print *, ""
  print *, "axisymmetric cavity: generating curve parameters"
  print *, "  a   =", pars(1)
  print *, "  b   =", pars(2)
  print *, "  zsc =", pars(3)
  print *, ""
  call cavity3d_report_curve(np, pars)
  print *, ""

  do iik = 1,nzk
    zk = zks(iik)
    do iio = 1,nord
      norder = norders_run(iio)

      call cavity3d_h_from_k(zk, ppw, norder, h)

      do iip = 1,nipt
        iptype0 = iptypes_run(iip)

        call get_cavity3d_npat_mem(np, pars, h, dthmax, ntbase, iort, &
          norder, iptype0, nch2d, npatches, npts)

        allocate(srcvals(12,npts), srccoefs(9,npts), qwts(npts))
        allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

        call get_cavity3d_npat(np, pars, h, dthmax, ntbase, iort, &
          norder, iptype0, npatches, npts, norders, ixyzs, iptype, &
          srccoefs, srcvals)

!
!  surface area and enclosed volume, the latter by the divergence
!  theorem  V = (1/3) \int x . n dA
!
        call get_qwts(npatches, norders, ixyzs, iptype, npts, srcvals, &
          qwts)

        area = 0
        vol = 0
        do i = 1,npts
          area = area + qwts(i)
          vol = vol + (srcvals(1,i)*srcvals(10,i) + &
            srcvals(2,i)*srcvals(11,i) + &
            srcvals(3,i)*srcvals(12,i))*qwts(i)
        enddo
        vol = vol/3

        if (iptype0.eq.1) then
          ptag = 'tri'
        else
          ptag = 'quad'
        endif

        write(6,'(a)') ' '
        write(6,'(a,f6.2,a,i2,a,a,a,f9.6)') ' zk=', zk, '  norder=', &
          norder, '  ', trim(ptag), '   h=', h
        write(6,'(a,i6,a,i8,a,i10)') '   nch2d=', nch2d, &
          '  npatches=', npatches, '  npts=', npts
        write(6,'(a,e12.5,a,e11.4)') '   area  =', area, &
          '   rel err =', abs(area-aref)/aref
        write(6,'(a,e12.5,a,e11.4)') '   volume=', vol, &
          '   rel err =', abs(abs(vol)-vref)/vref
        if (vol.lt.0) then
          write(6,'(a)') '   WARNING: normals point inward'
        endif

!
!  vtk plots.  The quadratic skeleton mesh is always written; the full
!  high order surface and the normals only for the coarsest wavenumber.
!
        dlam = 2*pi/zk

        write(fname,'(a,a,a,i2.2,a,i2.2,a)') trim(dir_out), &
          'cavity3d_', trim(ptag), int(zk,8), '_p', norder, '_msh.vtk'
        call surf_quadratic_msh_vtk_plot(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, trim(fname), 'a')

        if (ifplot_full.eq.1 .and. iik.eq.1) then
          write(fname,'(a,a,a,i2.2,a,i2.2,a)') trim(dir_out), &
            'cavity3d_', trim(ptag), int(zk,8), '_p', norder, '.vtk'
          call plot_surface_info_all(dlam, npatches, norders, ixyzs, &
            iptype, npts, srccoefs, srcvals, trim(fname), 'a')

          write(fname,'(a,a,a,i2.2,a,i2.2,a)') trim(dir_out), &
            'cavity3d_', trim(ptag), int(zk,8), '_p', norder, '_nrm.vtk'
          call surf_vtk_plot_vec(npatches, norders, ixyzs, iptype, &
            npts, srccoefs, srcvals, srcvals(10:12,1:npts), &
            trim(fname), 'a')
        endif

!
!  go3 files.  The go3 reader (in_go3.f90) hardwires iptype = 1, so only
!  the triangular meshes can round trip through that format.
!
        if (iptype0.eq.1) then
          write(fname,'(a,a,i2.2,a,i2.2,a)') trim(dir_out), &
            'cavity3d_tri', int(zk,8), '_p', norder, '.go3'
          call write_go3(trim(fname), norder, npatches, npts, srcvals)
        endif

        deallocate(srcvals, srccoefs, qwts, norders, ixyzs, iptype)
      enddo
    enddo
  enddo

  print *, ""
  print *, "meshes written to ", trim(dir_out)
  print *, ""

  stop
end program cavity3d_mesh
!
!
!
!
subroutine cavity3d_report_curve(np, pars)
!
!  Print the characteristic dimensions of the generating curve, plus its
!  arclength, surface area of revolution and enclosed volume computed by
!  composite Gauss-Legendre quadrature.  Useful as a sanity check when
!  the curve parameters are changed from the defaults.
!
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  integer *8, intent(in) :: np
  real *8, intent(in) :: pars(np)

  real *8, allocatable :: ts(:), ws(:), umat(:,:), vmat(:,:)

  done = 1.0d0
  pi = atan(done)*4

  k = 32
  npan = 400
  itype = 1
  allocate(ts(k), ws(k), umat(k,k), vmat(k,k))
  call legeexps(itype, k, ts, umat, vmat, ws)

  rl = 0
  area = 0
  vol = 0
  rmax = 0
  rinmax = 0
  do ipan = 1,npan
    slo = pi*(ipan - done)/npan
    shi = pi*ipan/npan
    sm = (shi + slo)/2
    sh = (shi - slo)/2
    do j = 1,k
      ss = sm + sh*ts(j)
      w = ws(j)*sh
      call funcurve_cavity3d(ss, np, pars, r, z, drds, dzds, d2rds2, &
        d2zds2)
      spd = sqrt(drds*drds + dzds*dzds)
      rl = rl + w*spd
      area = area + w*2*pi*r*spd
      vol = vol + w*pi*r*r*dzds
      if (ss.lt.pi/2 .and. r.gt.rmax) rmax = r
      if (ss.gt.pi/2 .and. r.gt.rinmax) rinmax = r
    enddo
  enddo

  call funcurve_cavity3d(0.0d0, np, pars, r0, z0, d1, d2, d3, d4)
  call funcurve_cavity3d(pi/2, np, pars, rl2, zl2, d1, d2, d3, d4)
  call funcurve_cavity3d(pi, np, pars, r1, z1, d1, d2, d3, d4)

  write(6,'(a,e11.4,a,f10.6)') '  outer pole      r=', r0, '  z=', z0
  write(6,'(a,e11.4,a,f10.6)') '  inner apex      r=', r1, '  z=', z1
  write(6,'(a,f11.6,a,f10.6)') '  aperture lip    r=', rl2, '  z=', zl2
  write(6,'(a,f11.6)')         '  outer rho_max    =', rmax
  write(6,'(a,f11.6)')         '  inner rho_max    =', rinmax
  write(6,'(a,f14.9)')         '  arclength        =', rl
  write(6,'(a,f14.9)')         '  surface area     =', area
  write(6,'(a,f14.9)')         '  enclosed volume  =', abs(vol)

  return
end subroutine cavity3d_report_curve
