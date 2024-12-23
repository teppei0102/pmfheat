module parameters
  implicit none
  ! control parameter:
  integer,parameter::input_switch=0
  ! 0 means w/o magnetic fields, 1 means with input data.
  ! 2 means WIMP with UCMH.
  ! 3 means shock heating with UCMH.
  double precision,parameter::B_1Mpc=0.50
  double precision,parameter::n_B=-2.9
  double precision,parameter::sp_G=19.47
  ! B_1Mpc is in the unit of nano Gauss.
  ! n_B magnetic spectral index, -3 is scale-free.
  ! sp_G Gamma function with an argument of (n_B+3)/2 = [(n_B+1)/2]!
  double precision,parameter::B_cut=1.00  ! [nG], only used for TS06.
  character(9),parameter::index="05nG/n-29"
  ! character(14),parameter::dirname="k_cut/n-29/B05"
  character(8),parameter::dirname="./output"
  ! character(13),parameter::dirname="01nG/n-2/5Mpc"
  integer,parameter::ngrid=1
  integer,parameter::dgrid=1
  ! number of grid in x-direction calculated by each core;(dgrid*n_cpu=ngrid)
  integer,parameter::diss_switch=1
  ! 0:amplitude decreasing, 1:cut-off increasing
  integer,parameter::heat_switch=3
  ! 0:no magnetic heating, 1:ambipolar diffusion, 2:decaying turbulence, 3:both of AD and DT
  ! 4:heating model in Schleicher, Banerjee & Klessen 2008 (arXiv:0807.3802)
  ! 5:heating model with AD in Sethi & Subramanian 2005 (arXiv:astro-ph/0405413)
  ! 6:heating model with AD in Tashiro & Sugiyama 2006 (arXiv:astro-ph/0607169)
  integer,parameter::cool_switch=0
  ! 0 means absence of cooling, 1 means presence.
  double precision,parameter::z_abs=17.0
  ! the position of the absorption (17.0 for EDGES)

  double precision,parameter::pi=3.141592653589793238
  ! physical parameters
  double precision,parameter::nano=1.0d-9        !
  double precision,parameter::erg=6.2415d11      ! value of 1 erg (eV)
  double precision,parameter::Mpc=3.0857d24      ! value of 1 Mpc (cm)
  double precision,parameter::c=2.99792d10       ! speed of light (cm/sec)
  double precision,parameter::c2_inv=1.0/(c*c)   !
  double precision,parameter::kB=8.617d-5        ! Boltzmann constant (eV/K)
  double precision,parameter::h=4.135667d-15     ! Planck constant (eV sec)
  double precision,parameter::h2_inv=1.0/(h*h)   !
  double precision,parameter::G=6.674d-8         ! gravitational constant (cm^3/g/sec^2)
  double precision,parameter::mc2_e=511.0d3      ! static energy of electron (eV)
  double precision,parameter::mc2_p=938.272d6    ! static energy of proton (eV)
  double precision,parameter::E_1s=13.6          ! binding energy of 1s hydrogen (eV)
  double precision,parameter::E_2s=E_1s/4.0      ! binding energy of 2s hydrogen (eV)
  double precision,parameter::SB=4.72d-3         ! a=4\sigma/c, radiation energy density constant (eV/cm^3/K^4)
  double precision,parameter::sigma=6.6524d-25   ! Thomson scattering cross section (cm^2)
  ! cosmological parameters
  double precision,parameter::T_cmb_now=2.7255
  double precision,parameter::z_rec=1089.0
  double precision,parameter::T_rec=(1.0+z_rec)*T_cmb_now
  double precision,parameter::time_rec=1.192d+13 ! age of the universe at the recombination (sec) (= 378,000 (yr))
  double precision,parameter::x_rec=0.5
  double precision,parameter::a_rec=1.0/(1.0+z_rec)
  ! data from Planck 2018
  ! double precision,parameter::reduced_h=0.500 ! default value in RECFAST code
  double precision,parameter::reduced_h=0.6727
  double precision,parameter::omega_k=0.000
  double precision,parameter::omega_m=0.120/reduced_h/reduced_h   !
  double precision,parameter::omega_b=0.02237/reduced_h/reduced_h !
  double precision,parameter::omega_l=1.0-omega_b-omega_m
  double precision,parameter::h0=100.0*reduced_h*0.3d-19 ! Hubble constant (1/sec)

  double precision::a, la, a_old, la_old, d_com, eta
  double precision::Tg1, Tg2, Tg3, Tg4
  double precision::nH, w

  ! simulation parameters
  integer(kind(8))::plan
  ! integer,parameter::FFTW_FORWARD=-1
  ! integer,parameter::FFTW_BACKWARD=1
  ! integer,parameter::FFTW_ESTIMATE=0

  integer::loop
  integer,parameter::nloop=1    ! number of realizations
  character(2)::loop_number
  integer::step=0
  double precision,parameter::boxsize=1.0        ! fixed for input!!! boxsize (Mpc)
  double precision,parameter::B_amp=1.0d0        ! fixed for input!!! amplitude of comoving B-field in Mpc (nano G)
  double precision,parameter::cfl=0.05           ! default is 0.01 (2017.01.17)
  integer::xi,yi,zi
  double precision::xgrid,ygrid,zgrid
  double precision::x1,x2,x3,x4
  double precision::fdata1=1.0
  double precision::dfdata1=0.0
  double precision::dfdata2=0.0
  double precision,dimension(ngrid,ngrid,ngrid)::data_in1, data_in2, Tg, dTg, x, dx
  ! minoda added 2023.10.11
  double precision::Ts
  double precision,dimension(ngrid,ngrid,ngrid)::fb, dfb ! fb is the real baryon density fluctuation, and dfb is time-differential fb.
  ! minoda added on 25th Feb. 2019
  double precision,dimension(ngrid,ngrid,ngrid)::fm,fm_out
  ! fm is the matter density contrast with lower bound, fm_out is without bound.
  double precision::E_ad,E_dt,E_ad_0,E_dt_0
  double precision,dimension(ngrid,ngrid,ngrid)::Tg_one_one,Tg_two_one,Tg_two_two
  double precision,dimension(ngrid,ngrid,ngrid)::dTg_one_one,dTg_two_one,dTg_two_two

  ! initial condition
  ! double precision,parameter::data_in1_init=8.0e2
  double precision,parameter::z_start=1000.00  ! default is 1000.0
  double precision,parameter::a_start=1.0/(1.0+z_start)
  double precision,parameter::T_start=T_rec*(a_rec/a_start)
  double precision,parameter::x_start=0.9999     ! default is 0.9999
  double precision,parameter::x_test=0.10900
  ! double precision,parameter::f0b_start=a_start
  double precision,parameter::fb_start=0.0
  ! minoda added on 25th Feb. 2019
  double precision,parameter::fm_start=0.0
  ! file input/output parameter
  integer,parameter::foutput=20000
  integer,parameter::noutput=100               ! how many times you write data from z=z_start to z=0
  integer::ioutput=0
  integer::i
  character(2)::output_number
  character(3)::cpu_number
  logical,dimension(noutput)::print_check=.true.
  double precision,dimension(noutput)::print_time
  !
  double precision,parameter::dla_start=1.0d-6 ! default is 1.0d-6 (2017.01.17)
  double precision::dla,dla_old,dla_Tg,dla_x
  !  double precision::dla_new
  double precision,parameter::z_tc=950.0       ! when tight coupling is violated.
  double precision,parameter::z_end=10.0
  double precision,parameter::a_end=1.0/(1.0+z_end)
  double precision,parameter::a_now=1.0

  double precision::T_photon                   ! tempolariliy using since 6/13 if input=0
  double precision::lambda_cut                 ! in the unit of Mpc.
  double precision::m_B                        ! =2(n_B+3)/(n_B+5)

end module parameters
