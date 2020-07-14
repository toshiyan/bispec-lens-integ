!///////////////////////////////////////////////////////////////////////!
! Lensing Bispectrum
!///////////////////////////////////////////////////////////////////////!

module bispec
  use cosmofunc, only: cosmoparams
  use bstool,    only: bispecfunc, bispec_lens_pb_kernel, bispec_lens_pb_init, bispec_lens_lss_kernel, bispec_lens_lss_init, zinterp
  implicit none

contains


subroutine bispec_lens_init(lmax,zs,k,pk0,fitmodel,pkmodel,zn,zmin,zmaxcut,zspacetype,kn,Omegam,H0,b1d,b2d,b3d,pb_wp,pb_ck)
  implicit none
  character(8), intent(in) :: fitmodel, pkmodel
  integer, intent(in) :: lmax, zn, kn, zspacetype
  double precision, intent(in) :: zmin, zmaxcut, Omegam, H0
  double precision, intent(in), dimension(3) :: zs
  double precision, intent(in), dimension(1:kn) :: k, pk0
  double precision, intent(out), dimension(14,zn) :: b1d
  double precision, intent(out), dimension(6,zn,lmax) :: b2d
  double precision, intent(out), dimension(3,zn,lmax) :: b3d
  double precision, intent(out), dimension(3,zn,lmax) :: pb_ck
  double precision, intent(out), dimension(3,3,zn,lmax) :: pb_wp
  !opt4py :: fitmodel = 'RT'
  !opt4py :: pkmodel = 'T12'
  !opt4py :: zn = 100
  !opt4py :: zmin = 0.0001
  !opt4py :: zmaxcut = 1100
  !opt4py :: zspacetype = 1
  !opt4py :: Omegam = 0.279
  !opt4py :: H0 = 70.
  !opt4py :: kn = 0
  !add2py :: kn = len(k)
  !internal
  double precision :: z(zn), dz(zn), zmax
  type(cosmoparams) :: cp
  type(bispecfunc)  :: b

  ! cosmological parameters (to compute background quantities) matched to RT fullsky sims
  cp%Om = Omegam
  cp%H0 = H0
  cp%w0 = -1d0
  cp%wa = 0d0
  cp%nu = 0d0
  cp%h  = cp%H0/100d0
  cp%Ov = 1d0 - cp%om !flat
  cp%ns = 0.9645d0

  ! z points
  zmax = minval(zs)  ! delta_m at only z < min(zs) are significantly correlated
  if (zmaxcut<zmax)  zmax = zmaxcut  ! integrate z upto zmaxcut
  call zinterp(zmin,zmax,zn,zspacetype,z,dz)

  ! precompute quantities relevant to bispectrum
  call bispec_lens_lss_init(cp,b,z,dz,zs,k*cp%h,pk0/cp%h**3,(/1,lmax/),fitmodel,pkmodel) !correction for h/Mpc to /Mpc
  call bispec_lens_pb_init(cp,b%kl,b%pl,z,dz,zs,(/1,lmax/),pb_wp,pb_ck)

  b1d = 0d0
  b2d = 0d0
  b3d = 0d0

  b1d(1,:) = b%z
  b1d(2,:) = b%D
  b1d(3,:) = b%sz
  b1d(4,:) = b%knl
  b1d(5,:) = b%zker
  
  b2d(1,:,:) = b%kl
  b2d(2,:,:) = b%plL
  b2d(3,:,:) = b%pl
  b2d(4,:,:) = b%q

  select case(fitmodel)
  case ('SC','GM','3B')
    b3d(:,:,:) = b%abc
  case ('RT')
    b1d(6,:) = b%gammanz
    b1d(7,:) = b%anz
    b1d(8,:) = b%bnz
    b1d(9,:) = b%cnz
    b1d(10,:) = b%alphanz1
    b1d(11,:) = b%alphanz2
    b1d(12,:) = b%betanz
    b1d(13,:) = b%dnq
    b1d(14,:) = b%n
    b2d(5,:,:) = b%Ik
    b2d(6,:,:) = b%PE
  end select
 
end subroutine bispec_lens_init


subroutine bispec_lens(l1,l2,l3,b1d,b2d,b3d,pb_wp,pb_ck,fitmodel,pb,Omegam,H0,lmax,zn,bl)
!*  Compute lensing bispectrum analytically
!*
!*  Args
!*    :l1/l2/l3 (int) : multipoles of the bispectrum
!*    :b1d/b2d/b3d (double) : functions
!*    :pb_wp/pb_ck (double) : functions for post-Bron
!*
!*  Args(optional):
!*    :fitmodel (str) : fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
!*
!*  Returns:
!*    :bl (double) : lensing bispectrum
!*
  implicit none
  !I/O
  character(8), intent(in) :: fitmodel
  integer, intent(in) :: l1, l2, l3, zn, lmax
  logical, intent(in) :: pb
  double precision, intent(in) :: Omegam, H0
  double precision, intent(in), dimension(14,zn) :: b1d
  double precision, intent(in), dimension(6,zn,lmax) :: b2d
  double precision, intent(in), dimension(3,zn,lmax) :: b3d
  double precision, intent(in), dimension(3,zn,lmax) :: pb_ck
  double precision, intent(in), dimension(3,3,zn,lmax) :: pb_wp
  double precision, intent(out) :: bl
  !opt4py :: fitmodel = 'RT'
  !opt4py :: pb = True
  !opt4py :: Omegam = 0.279
  !opt4py :: H0 = 70.
  !opt4py :: zn = 0
  !opt4py :: lmax = 0
  !add2py :: zn = len(b1d[0,:])
  !add2py :: lmax = len(b2d[0,0,:])
  !internal
  double precision :: bl_lss, bl_pb
  type(cosmoparams) :: cp
  type(bispecfunc)  :: b

  bl = 0d0
  if (mod(l1+l2+l3,2)/=0) return

  ! cosmological parameters (to compute background quantities)
  cp%Om = Omegam
  cp%H0 = H0
  cp%w0 = -1d0
  cp%wa = 0d0
  cp%nu = 0d0
  cp%h  = cp%H0/100d0
  cp%Ov = 1d0 - cp%om !flat
  cp%ns = 0.9645d0

  ! pass bispectrum related arrays
  allocate(b%zker(zn),b%z(zn),b%sz(zn),b%D(zn),b%knl(zn),b%q(zn,lmax),b%kl(zn,lmax),b%plL(zn,lmax),b%pl(zn,lmax))
  b%z = b1d(1,:) 
  b%D = b1d(2,:) 
  b%sz = b1d(3,:) 
  b%knl = b1d(4,:) 
  b%zker = b1d(5,:) 
  b%kl = b2d(1,:,:) 
  b%plL = b2d(2,:,:) 
  b%pl = b2d(3,:,:) 
  b%q = b2d(4,:,:) 
  select case(fitmodel)
  case ('SC','GM','3B')
    allocate(b%abc(3,zn,lmax))
    b%abc = b3d(:,:,:)
  case ('RT')
    allocate(b%dnq(zn),b%PE(zn,lmax),b%Ik(zn,lmax),b%n(zn),b%gammanz(zn),b%anz(zn),b%bnz(zn),b%cnz(zn),b%alphanz1(zn),b%alphanz2(zn),b%betanz(zn))
    b%gammanz = b1d(6,:)
    b%anz = b1d(7,:)
    b%bnz = b1d(8,:)
    b%cnz = b1d(9,:) 
    b%alphanz1 = b1d(10,:) 
    b%alphanz2 = b1d(11,:) 
    b%betanz = b1d(12,:) 
    b%dnq = b1d(13,:) 
    b%n = b1d(14,:) 
    b%Ik = b2d(5,:,:) 
    b%PE = b2d(6,:,:) 
  end select
 
  ! bispectrum
  bl_lss = 0d0
  bl_pb  = 0d0
  call bispec_lens_lss_kernel(cp,b,l1,l2,l3,bl_lss,fitmodel)
  if (pb)  call bispec_lens_pb_kernel(l1,l2,l3,pb_wp,pb_ck,bl_pb)

  bl = bl_lss + bl_pb

end subroutine bispec_lens


end module bispec


