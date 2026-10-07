!
!  Analytic axisymmetric ("3d") version of the elliptical cavity used in
!  helm-high-freq-conv/paper_examples/cavity_example.f
!
!  ---------------------------------------------------------------------
!  Geometry
!  ---------------------------------------------------------------------
!
!  The 2d cavity in cavity_example.f is built by load_cavity_zpars, which
!  samples a generating curve in the half plane rho >= 0, reflects it
!  through rho -> -rho, and then represents the resulting closed curve by
!  a Fourier series (funcurv_zfft).  The Fourier step is only a smoothing
!  / resampling device; the curve being sampled is fully analytic.  With
!  the smoothing width h = a/sqrt(2) that load_cavity_zpars uses, the
!  Gaussian convolutions collapse and the generating curve is
!
!     G(x)    = x*erf(x/a) + (a/sqrt(pi))*exp(-(x/a)^2)      (smoothed |x|)
!     theta(s)= (b-a) + 2*(1 - (b-a)/pi)*G(s - pi/2)
!     r0(s)   = 1 - a*erf((s - pi/2)/a)
!
!     rho(s)  = r0(s)*sin(theta(s))
!     z(s)    = zsc*r0(s)*cos(theta(s))            s in [0, pi]
!
!  with a = 0.2, b = pi/12, zsc = 1.2 reproducing cavity_example.f.
!
!  The derivatives are elementary because G'(x) = erf(x/a):
!
!     theta' = c*erf(x/a),                 c = 2*(1 - (b-a)/pi)
!     theta''= c*(2/(a*sqrt(pi)))*exp(-(x/a)^2)
!     r0'    = -(2/sqrt(pi))*exp(-(x/a)^2)
!     r0''   = (4*x/(a^2*sqrt(pi)))*exp(-(x/a)^2)
!
!  The key structural fact that makes a 3d version possible is that the
!  ramp coefficient telescopes exactly,
!
!     theta(0) = (b-a) + c*(pi/2) = (b-a) + (pi - (b-a)) = pi = theta(pi),
!
!  (up to the O(exp(-(pi/(2a))^2)) Gaussian tail, i.e. ~1d-27 for a=0.2)
!  so rho(0) = rho(pi) = 0 and the generating curve both starts and ends
!  ON the axis of rotation.  Revolving it about the z axis therefore gives
!  a closed, smooth, analytic surface rather than something with a seam:
!
!      (rho(s)*cos(t), rho(s)*sin(t), z(s)),   s in [0,pi], t in [0,2pi)
!
!  which is a hollow spheroidal shell with a circular aperture at the top
!  -- the axisymmetric open cavity / Helmholtz resonator analogue of the
!  2d geometry.  For the default parameters:
!
!      outer shell      rho_max = 1.2000 at z ~ 0, south pole z = -1.4400
!      aperture lip     rho     = 0.2793 at z =  1.1523  (diameter 0.5585)
!      inner cavity     rho_max = 0.7996, apex at z = -0.9600
!      shell thickness at the bottom                      0.4800
!      generating curve arclength                          6.5698627
!      (2*6.5698627 = 13.1397 is the "13.14" in the 2d nch formula)
!      surface area                                       29.5322957
!      enclosed volume                                     5.8814978
!
!  Both axis points sit at the bottom (z = -1.44 outer, z = -0.96 inner),
!  so the revolved surface carries two polar caps there.
!
!  ---------------------------------------------------------------------
!  Meshing
!  ---------------------------------------------------------------------
!
!  The mesh is produced by fmm3dbie's get_axissym_fcurve_npat, which takes
!  the generating curve as a function handle and handles the two polar
!  caps with a circle mesh (circ_mesh.f90) so there is no 1/rho
!  coordinate singularity on the axis.
!
!  Resolution is wavelength driven.  For a wavenumber zk and a target of
!  ppw points per wavelength at order norder,
!
!      h = (norder+1)*lambda/ppw,        lambda = 2*pi/|zk|
!
!  is the target patch arclength (each patch carries norder+1 nodes per
!  direction).  The chunking of [0,pi] equidistributes the monitor
!  function
!
!      M(s) = speed(s)*( 1/h + kappa(s)/dthmax )
!
!  so that a chunk is at most ~h in arclength AND turns by at most
!  ~dthmax.  The curvature term matters here: the aperture lip has
!  kappa_max = 7.43 (radius of curvature 0.1345) at s = 1.519, so at
!  modest frequencies the lip, not the wavelength, sets the local mesh
!  size.
!
!  The azimuthal count is graded per chunk, nt ~ 2*pi*rho/h rounded to
!  ntbase*2^k, so patches stay roughly isotropic instead of becoming
!  azimuthally thin near the axis (a constant azimuthal count costs about
!  2x the patches for the same resolution).  ntbase is a multiple of 8 so
!  that the cap rim, which carries 4*nrts(2,cap) uniformly spaced angular
!  intervals, lines up with the adjacent band chunk.
!
!  ---------------------------------------------------------------------
!  Routines
!  ---------------------------------------------------------------------
!    cavity3d_defpars        - default (a, b, zsc) reproducing the 2d cavity
!    funcurve_cavity3d       - the generating curve, fmm3dbie fcurve interface
!    cavity3d_h_from_k       - target patch arclength from (zk, ppw, norder)
!    get_cavity3d_tchse_mem  - number of 2d chunks for a given h, dthmax
!    get_cavity3d_tchse      - the graded chunk breakpoints in s
!    get_cavity3d_nrts       - per-chunk patch counts, and the cap parameters
!    get_cavity3d_npat_mem   - npatches, npts  (convenience wrapper)
!    get_cavity3d_npat       - the discretized surface (convenience wrapper)
!
!  get_cavity3d_npat_mem / get_cavity3d_npat are the user facing pair and
!  follow the usual fmm3dbie two stage (size, then fill) convention.
!

      subroutine cavity3d_defpars(np, pars)
!
!  Default generating curve parameters, chosen to reproduce the 2d cavity
!  cross section of cavity_example.f exactly.
!
!  Input arguments:
!    - np: integer *8
!        length of pars, must be at least 3
!
!  Output arguments:
!    - pars: real *8 (np)
!        * pars(1) = a,   radial step / corner smoothing width  (0.2)
!        * pars(2) = b,   aperture angle parameter              (pi/12)
!        * pars(3) = zsc, axial elongation of the shell         (1.2)
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: np
      real *8, intent(out) :: pars(np)

      done = 1.0d0
      pi = atan(done)*4

      pars(1) = 0.2d0
      pars(2) = pi/12.0d0
      pars(3) = 1.2d0

      return
      end
!
!
!
!
      subroutine funcurve_cavity3d(s, np, pars, r, z, drds, dzds, &
        d2rds2, d2zds2)
!
!  Generating curve of the axisymmetric cavity, in the fmm3dbie fcurve
!  calling convention expected by get_axissym_fcurve_npat.
!
!  Input arguments:
!    - s: real *8
!        curve parameter, s in [0, pi].  s = 0 is the outer (south) pole,
!        s = pi/2 the aperture lip, s = pi the inner cavity apex.  Both
!        endpoints lie on the axis of rotation.
!    - np: integer *8
!        number of parameters, must be at least 3
!    - pars: real *8 (np)
!        see cavity3d_defpars
!
!  Output arguments:
!    - r, z: real *8
!        cylindrical radius and height of the point
!    - drds, dzds: real *8
!        first derivatives with respect to s
!    - d2rds2, d2zds2: real *8
!        second derivatives with respect to s
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: np
      real *8, intent(in) :: s, pars(np)
      real *8, intent(out) :: r, z, drds, dzds, d2rds2, d2zds2

      intrinsic erf

      done = 1.0d0
      pi = atan(done)*4
      rsqpi = done/sqrt(pi)

      a = pars(1)
      b = pars(2)
      zsc = pars(3)

      x = s - pi/2
      ee = erf(x/a)
      gg = exp(-(x/a)**2)

!
!  smoothed corner G and its derivatives; G' = erf(x/a) exactly
!
      cth = 2*(done - (b-a)/pi)

      gval = x*ee + a*rsqpi*gg
      th = (b-a) + cth*gval
      thp = cth*ee
      thpp = cth*2*rsqpi/a*gg

!
!  radial profile r0 and its derivatives
!
      r0 = done - a*ee
      r0p = -2*rsqpi*gg
      r0pp = 4*x/(a*a)*rsqpi*gg

      st = sin(th)
      ct = cos(th)

      r = r0*st
      drds = r0p*st + r0*ct*thp
      d2rds2 = r0pp*st + 2*r0p*ct*thp + r0*(ct*thpp - st*thp*thp)

      z = zsc*r0*ct
      dzds = zsc*(r0p*ct - r0*st*thp)
      d2zds2 = zsc*(r0pp*ct - 2*r0p*st*thp - r0*st*thpp - &
        r0*ct*thp*thp)

      return
      end
!
!
!
!
      subroutine cavity3d_h_from_k(zk, ppw, norder, h)
!
!  Target patch arclength for ppw points per wavelength at order norder.
!
!  Each patch carries norder+1 nodes in each parameter direction, so a
!  patch of arclength h resolves (norder+1)*lambda/h points per
!  wavelength; inverting gives h = (norder+1)*lambda/ppw.
!
!  Input arguments:
!    - zk: real *8
!        wavenumber (magnitude is used)
!    - ppw: real *8
!        desired points per wavelength
!    - norder: integer *8
!        order of discretization on each patch
!
!  Output arguments:
!    - h: real *8
!        target patch arclength
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: norder
      real *8, intent(in) :: zk, ppw
      real *8, intent(out) :: h

      done = 1.0d0
      pi = atan(done)*4

      rlam = 2*pi/abs(zk)
      h = (norder + done)*rlam/ppw

      return
      end
!
!
!
!
      subroutine cavity3d_monitor_tab(np, pars, h, dthmax, nst, stab, &
        rtab, ctab)
!
!  Tabulate the mesh density monitor function
!
!      M(s) = speed(s)*( 1/h + kappa(s)/dthmax )
!
!  and its running integral on a uniform grid of nst points covering
!  [0, pi].  Equidistributing the running integral gives chunks that are
!  at most ~h long and turn by at most ~dthmax.
!
!  This is an internal helper for get_cavity3d_tchse / _mem.
!
!  Input arguments:
!    - np: integer *8
!        number of curve parameters
!    - pars: real *8 (np)
!        curve parameters, see cavity3d_defpars
!    - h: real *8
!        target patch arclength
!    - dthmax: real *8
!        maximum turning angle per chunk, in radians
!    - nst: integer *8
!        number of tabulation points
!
!  Output arguments:
!    - stab: real *8 (nst)
!        uniform grid on [0, pi]
!    - rtab: real *8 (nst)
!        cylindrical radius r(s) at the grid points
!    - ctab: real *8 (nst)
!        running trapezoidal integral of M, ctab(1) = 0
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: np, nst
      real *8, intent(in) :: pars(np), h, dthmax
      real *8, intent(out) :: stab(nst), rtab(nst), ctab(nst)

      real *8, allocatable :: mtab(:)

      done = 1.0d0
      pi = atan(done)*4

      allocate(mtab(nst))

      ds = pi/(nst - done)

      do i = 1,nst
        stab(i) = (i - done)*ds
        call funcurve_cavity3d(stab(i), np, pars, r, z, drds, dzds, &
          d2rds2, d2zds2)
        rtab(i) = r
        spd = sqrt(drds*drds + dzds*dzds)
        rkap = abs(drds*d2zds2 - dzds*d2rds2)/spd**3
        mtab(i) = spd*(done/h + rkap/dthmax)
      enddo

      ctab(1) = 0
      do i = 2,nst
        ctab(i) = ctab(i-1) + (mtab(i) + mtab(i-1))*ds/2
      enddo

      return
      end
!
!
!
!
      subroutine get_cavity3d_tchse_mem(np, pars, h, dthmax, nch2d)
!
!  Number of chunks in the generating curve discretization for a given
!  target patch arclength and turning angle.  The first and last chunks
!  become the two polar caps.
!
!  Input arguments:
!    - np: integer *8
!        number of curve parameters
!    - pars: real *8 (np)
!        curve parameters, see cavity3d_defpars
!    - h: real *8
!        target patch arclength, e.g. from cavity3d_h_from_k
!    - dthmax: real *8
!        maximum turning angle per chunk, in radians.  pi/6 is a
!        reasonable default
!
!  Output arguments:
!    - nch2d: integer *8
!        number of chunks, at least 5
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: np
      real *8, intent(in) :: pars(np), h, dthmax
      integer *8, intent(out) :: nch2d

      real *8, allocatable :: stab(:), rtab(:), ctab(:)
      integer *8 nst

      nst = 20001
      allocate(stab(nst), rtab(nst), ctab(nst))

      call cavity3d_monitor_tab(np, pars, h, dthmax, nst, stab, rtab, &
        ctab)

      nch2d = ceiling(ctab(nst))
      if (nch2d.lt.5) nch2d = 5

      return
      end
!
!
!
!
      subroutine get_cavity3d_tchse(np, pars, h, dthmax, nch2d, tchse)
!
!  Graded chunk breakpoints in s for the generating curve, obtained by
!  equidistributing the monitor function of cavity3d_monitor_tab.  Chunk 1
!  and chunk nch2d are the polar caps, as expected by
!  get_axissym_fcurve_npat.
!
!  Input arguments:
!    - np: integer *8
!        number of curve parameters
!    - pars: real *8 (np)
!        curve parameters, see cavity3d_defpars
!    - h: real *8
!        target patch arclength
!    - dthmax: real *8
!        maximum turning angle per chunk, in radians
!    - nch2d: integer *8
!        number of chunks, from get_cavity3d_tchse_mem
!
!  Output arguments:
!    - tchse: real *8 (nch2d+1)
!        chunk breakpoints, tchse(1) = 0, tchse(nch2d+1) = pi
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: np, nch2d
      real *8, intent(in) :: pars(np), h, dthmax
      real *8, intent(out) :: tchse(nch2d+1)

      real *8, allocatable :: stab(:), rtab(:), ctab(:)
      integer *8 nst, i, j

      done = 1.0d0
      pi = atan(done)*4

      nst = 20001
      allocate(stab(nst), rtab(nst), ctab(nst))

      call cavity3d_monitor_tab(np, pars, h, dthmax, nst, stab, rtab, &
        ctab)

!
!  invert the (monotone) running integral at nch2d equal increments
!
      tchse(1) = 0
      tchse(nch2d+1) = pi

      j = 1
      do i = 2,nch2d
        ctarg = ctab(nst)*(i - done)/(nch2d + done)
        do while (j.lt.nst-1 .and. ctab(j+1).lt.ctarg)
          j = j + 1
        enddo
        dc = ctab(j+1) - ctab(j)
        if (dc.gt.0) then
          tchse(i) = stab(j) + (stab(j+1) - stab(j))*(ctarg - ctab(j))/dc
        else
          tchse(i) = stab(j)
        endif
      enddo

      return
      end
!
!
!
!
      subroutine get_cavity3d_nrts(np, pars, nch2d, tchse, h, ntbase, &
        nrts, nmid, rmid)
!
!  Per chunk patch counts for the body of revolution, plus the polar cap
!  parameters, in the form get_axissym_fcurve_npat expects.
!
!  Band chunks (2 .. nch2d-1) get one patch in the s direction, since the
!  chunking in s already targets h, and a graded azimuthal count
!
!      nrts(2,i) = ntbase*2^k  ~  2*pi*max(r)/h on the chunk
!
!  so that patches stay roughly isotropic instead of becoming azimuthally
!  thin as r -> 0.  The cap rim carries 4*nrts(2,cap) uniformly spaced
!  angular intervals, so the caps are given a quarter of the azimuthal
!  count of the band chunk they adjoin; with ntbase a multiple of 8 the
!  cap rim and the band then line up.
!
!  Input arguments:
!    - np: integer *8
!        number of curve parameters
!    - pars: real *8 (np)
!        curve parameters, see cavity3d_defpars
!    - nch2d: integer *8
!        number of chunks
!    - tchse: real *8 (nch2d+1)
!        chunk breakpoints, from get_cavity3d_tchse
!    - h: real *8
!        target patch arclength
!    - ntbase: integer *8
!        smallest azimuthal patch count, should be a positive multiple
!        of 8.  8 is a reasonable default
!
!  Output arguments:
!    - nrts: integer *8 (2,nch2d)
!        nrts(1,i) patches in the s direction, nrts(2,i) in the azimuthal
!        direction.  For the caps, nrts(1,i) is the number of radial
!        intervals and nrts(2,i) the angular intervals per quadrant
!    - nmid: integer *8
!        number of intervals across the middle square of the cap mesh
!    - rmid: real *8
!        radius of the middle square of the cap mesh
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: np, nch2d, ntbase
      real *8, intent(in) :: pars(np), tchse(nch2d+1), h
      integer *8, intent(out) :: nrts(2,nch2d), nmid
      real *8, intent(out) :: rmid

      integer *8 i, j, nsamp, nt, ntb
      real *8 rbar

      done = 1.0d0
      pi = atan(done)*4

      ntb = ntbase
      if (ntb.lt.8) ntb = 8

      nsamp = 33

      do i = 1,nch2d
!
!  largest cylindrical radius on this chunk
!
        rbar = 0
        do j = 1,nsamp
          ss = tchse(i) + (tchse(i+1) - tchse(i))*(j - done)/(nsamp - done)
          call funcurve_cavity3d(ss, np, pars, r, z, drds, dzds, &
            d2rds2, d2zds2)
          if (r.gt.rbar) rbar = r
        enddo

!
!  nt ~ 2*pi*rbar/h, rounded to ntb*2^k with k >= 0
!
        nt = ntb
        do while (nt.lt.2*pi*rbar/h)
          nt = nt*2
        enddo
        if (nt.gt.ntb) then
          if (abs(log(nt/(2*pi*rbar/h))).gt.abs(log(2*pi*rbar/h/(nt/2)))) then
            nt = nt/2
          endif
        endif

        nrts(1,i) = 1
        nrts(2,i) = nt
      enddo

!
!  polar caps: one radial interval, a quarter of the adjacent azimuthal
!  count, and a 2x2 middle square of comparable size
!
      nrts(1,1) = 1
      nrts(2,1) = max(nrts(2,2)/4, 2_8)

      nrts(1,nch2d) = 1
      nrts(2,nch2d) = max(nrts(2,nch2d-1)/4, 2_8)

      nmid = 2
      rmid = 0.5d0

      return
      end
!
!
!
!
      subroutine get_cavity3d_npat_mem(np, pars, h, dthmax, ntbase, &
        iort, norder, iptype0, nch2d, npatches, npts)
!
!  Number of patches and points in a high order mesh of the axisymmetric
!  cavity.  Convenience wrapper that builds the graded chunking and patch
!  counts internally and defers to get_axissym_fcurve_npat_mem.
!
!  Input arguments:
!    - np: integer *8
!        number of curve parameters, at least 3
!    - pars: real *8 (np)
!        curve parameters, see cavity3d_defpars
!    - h: real *8
!        target patch arclength, e.g. from cavity3d_h_from_k
!    - dthmax: real *8
!        maximum turning angle per chunk, in radians.  pi/6 is a
!        reasonable default
!    - ntbase: integer *8
!        smallest azimuthal patch count, a positive multiple of 8
!    - iort: integer *8
!        orientation flag.  iort = -1 gives outward normals on this
!        geometry (see the note in get_cavity3d_npat)
!    - norder: integer *8
!        order of discretization on each patch
!    - iptype0: integer *8
!        * iptype0 = 1,  triangular patches, Vioreanu-Rokhlin nodes
!        * iptype0 = 11, quadrangular patches, Gauss-Legendre nodes
!        * iptype0 = 12, quadrangular patches, Chebyshev nodes
!
!  Output arguments:
!    - nch2d: integer *8
!        number of chunks used in the generating curve
!    - npatches: integer *8
!        number of patches
!    - npts: integer *8
!        number of discretization points on the surface
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: np, ntbase, iort, norder, iptype0
      real *8, intent(in) :: pars(np), h, dthmax
      integer *8, intent(out) :: nch2d, npatches, npts

      real *8, allocatable :: tchse(:)
      integer *8, allocatable :: nrts(:,:)
      integer *8 nmid
      real *8 rmid

      external funcurve_cavity3d

      call get_cavity3d_tchse_mem(np, pars, h, dthmax, nch2d)

      allocate(tchse(nch2d+1), nrts(2,nch2d))

      call get_cavity3d_tchse(np, pars, h, dthmax, nch2d, tchse)
      call get_cavity3d_nrts(np, pars, nch2d, tchse, h, ntbase, nrts, &
        nmid, rmid)

      call get_axissym_fcurve_npat_mem(nch2d, tchse, &
        funcurve_cavity3d, np, pars, nrts, rmid, nmid, iort, norder, &
        iptype0, npatches, npts)

      return
      end
!
!
!
!
      subroutine get_cavity3d_npat(np, pars, h, dthmax, ntbase, iort, &
        norder, iptype0, npatches, npts, norders, ixyzs, iptype, &
        srccoefs, srcvals)
!
!  High order mesh of the axisymmetric cavity.  Convenience wrapper that
!  builds the graded chunking and patch counts internally and defers to
!  get_axissym_fcurve_npat.
!
!  Note on orientation: with u = s and v = t the surface normal is
!  X_s x X_t, which points INTO the solid for this generating curve, so
!  iort = -1 is the choice that gives outward normals (the convention for
!  exterior scattering).
!
!  Input arguments:
!    - np, pars, h, dthmax, ntbase, iort, norder, iptype0:
!        as in get_cavity3d_npat_mem
!    - npatches: integer *8
!        number of patches, from get_cavity3d_npat_mem
!    - npts: integer *8
!        number of discretization points, from get_cavity3d_npat_mem
!
!  Output arguments:
!    - norders: integer *8 (npatches)
!        order of discretization on each patch, = norder
!    - ixyzs: integer *8 (npatches+1)
!        starting location in srccoefs, srcvals for each patch
!    - iptype: integer *8 (npatches)
!        patch type, = iptype0
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du, dxyz/dv
!    - srcvals: real *8 (12,npts)
!        xyz, dxyz/du, dxyz/dv and the normal at the discretization nodes
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: np, ntbase, iort, norder, iptype0
      integer *8, intent(in) :: npatches, npts
      real *8, intent(in) :: pars(np), h, dthmax
      integer *8, intent(out) :: norders(npatches), iptype(npatches)
      integer *8, intent(out) :: ixyzs(npatches+1)
      real *8, intent(out) :: srccoefs(9,npts), srcvals(12,npts)

      real *8, allocatable :: tchse(:)
      integer *8, allocatable :: nrts(:,:)
      integer *8 nch2d, nmid
      real *8 rmid

      external funcurve_cavity3d

      call get_cavity3d_tchse_mem(np, pars, h, dthmax, nch2d)

      allocate(tchse(nch2d+1), nrts(2,nch2d))

      call get_cavity3d_tchse(np, pars, h, dthmax, nch2d, tchse)
      call get_cavity3d_nrts(np, pars, nch2d, tchse, h, ntbase, nrts, &
        nmid, rmid)

      call get_axissym_fcurve_npat(nch2d, tchse, funcurve_cavity3d, &
        np, pars, nrts, rmid, nmid, iort, norder, iptype0, npatches, &
        npts, norders, ixyzs, iptype, srccoefs, srcvals)

      return
      end
