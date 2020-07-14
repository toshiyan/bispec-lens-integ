import lib

def bispec_lens_init(lmax,zs,k,pk0,fitmodel='RT',pkmodel='T12',zn=100,zmin=0.0001,zmaxcut=1100,zspacetype=1,kn=0,Omegam=0.279,H0=70.):
  """
  Usage:
    :b1d,b2d,b3d,pb_ck,pb_wp = .bispec.bispec_lens_init(lmax,zs,k,pk0,fitmodel,pkmodel,zn,zmin,zmaxcut,zspacetype,kn,Omegam,H0):
  """
  kn = len(k)
  return lib.bispec.bispec_lens_init(lmax,zs,k,pk0,fitmodel,pkmodel,zn,zmin,zmaxcut,zspacetype,kn,Omegam,H0)

def bispec_lens(l1,l2,l3,b1d,b2d,b3d,pb_wp,pb_ck,fitmodel='RT',pb=True,Omegam=0.279,H0=70.,lmax=0,zn=0):
  """
  Compute lensing bispectrum analytically

  Args
    :l1/l2/l3 (*int*): multipoles of the bispectrum
    :b1d/b2d/b3d (*double*): functions
    :pb_wp/pb_ck (*double*): functions for post-Bron

  Args(optional):
    :fitmodel (*str*): fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)

  Returns:
    :bl (*double*): lensing bispectrum

  Usage:
    :bl = .bispec.bispec_lens(l1,l2,l3,b1d,b2d,b3d,pb_wp,pb_ck,fitmodel,pb,Omegam,H0,lmax,zn):
  """
  zn = len(b1d[0,:])
  lmax = len(b2d[0,0,:])
  return lib.bispec.bispec_lens(l1,l2,l3,b1d,b2d,b3d,pb_wp,pb_ck,fitmodel,pb,Omegam,H0,lmax,zn)

