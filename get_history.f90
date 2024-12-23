program get_history
  use parameters
  implicit none
  include "mpif.h"
  integer ierr,nprocs,myrank
  integer,parameter::checkstep=1

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,nprocs,ierr)
  call mpi_comm_rank(mpi_comm_world,myrank,ierr)

  if(input_switch.eq.0)then
     print*,""
     print*,"   //////////////////////////////////////////////////////"
     print*,"   ///  B_1Mpc =",B_1Mpc,"nG   ///"
     print*,"   ///  n_B =",n_B,"   ///"
     print*,"   ///  heat switch =",heat_switch," ///"
     print*,"   //////////////////////////////////////////////////////"
     print*,""
  else if(input_switch.eq.1)then
     if(myrank==0)print*,""
     if(myrank==0)print*,"   //////////////////////////////////////////////////////"
     if(myrank==0)print*,"   ///  amplitude of B  ",B_amp,"nG   ///"
     if(myrank==0)print*,"   ///  power index of B  ",index,"                    ///"
     if(myrank==0)print*,"   ///  boxsize         ",boxsize,"Mpc  ///"
     if(myrank==0)print*,"   ///  number of grid  ",ngrid,"                 ///"
     if(myrank==0)print*,"   //////////////////////////////////////////////////////"
     if(myrank==0)print*,""
  else if(input_switch.eq.2.or.input_switch.eq.3)then
     print*,""
     print*,"   //////////////////////////////////////////////////////"
     print*,"   ///   Thermal history with WIMP and UCMH model     ///"
     print*,"   //////////////////////////////////////////////////////"
     print*,""
  else
     print*,"!! --- Error --- Error --- Error --- !!"
     print*,"Please select appropriate input_switch"
  end if

  do loop=1,nloop
     write(loop_number,'(I2.2)')loop
     ! Read input file
     if(input_switch.eq.0)then
        if(heat_switch.eq.0)then
           call system("mkdir -p "//dirname//"/no")
           open(1,file=dirname//"/no/check.dat",status="replace")
           open(2,file=dirname//"/no/gamma_of_z.dat",status="replace")
           open(3,file=dirname//"/no/B_1Mpc.dat",status="replace")
           open(4,file=dirname//"/no/x_ion.dat",status="replace")
        else if(heat_switch.eq.1)then
           call system("mkdir -p "//dirname//"/ad")
           open(1,file=dirname//"/ad/check.dat",status="replace")
           open(2,file=dirname//"/ad/gamma_of_z.dat",status="replace")
           open(3,file=dirname//"/ad/B_1Mpc.dat",status="replace")
           open(4,file=dirname//"/ad/x_ion.dat",status="replace")
        else if(heat_switch.eq.2)then
           call system("mkdir -p "//dirname//"/dt")
           open(1,file=dirname//"/dt/check.dat",status="replace")
           open(2,file=dirname//"/dt/gamma_of_z.dat",status="replace")
           open(3,file=dirname//"/dt/B_1Mpc.dat",status="replace")
           open(4,file=dirname//"/dt/x_ion.dat",status="replace")
        else if(heat_switch.eq.3)then
           if(diss_switch.eq.0)then
              call system("mkdir -p "//dirname//"/both/k_all")
              open(1,file=dirname//"/both/k_all/check.dat",status="replace")
              open(2,file=dirname//"/both/k_all/gamma_of_z.dat",status="replace")
              open(3,file=dirname//"/both/k_all/B_1Mpc.dat",status="replace")
              open(4,file=dirname//"/both/k_all/x_ion.dat",status="replace")
           else if(diss_switch.eq.1)then
              call system("mkdir -p "//dirname//"/both/k_cut")
              open(1,file=dirname//"/both/k_cut/check.dat",status="replace")
              write(1,*)"####    z,       Tgas,       TCMB,       Tspin"
              open(2,file=dirname//"/both/k_cut/gamma_of_z.dat",status="replace")
              open(3,file=dirname//"/both/k_cut/B_1Mpc.dat",status="replace")
              open(4,file=dirname//"/both/k_cut/x_ion.dat",status="replace")
           else
              call system("mkdir -p "//dirname//"/both")
              open(1,file=dirname//"/both/check.dat",status="replace")
              open(2,file=dirname//"/both/gamma_of_z.dat",status="replace")
              open(3,file=dirname//"/both/B_1Mpc.dat",status="replace")
              open(4,file=dirname//"/both/x_ion.dat",status="replace")
           end if
        else if(heat_switch.ge.4.and.heat_switch.le.6)then
           call system("mkdir -p "//dirname)
           open(1,file=dirname//"/check.dat",status="replace")
           open(2,file=dirname//"/gamma_of_z.dat",status="replace")
           open(3,file=dirname//"/B_1Mpc.dat",status="replace")
           open(4,file=dirname//"/x_ion.dat",status="replace")
        else
           call system("mkdir -p "//dirname)
           open(1,file=dirname//"/check.dat",status="replace")
           open(2,file=dirname//"/gamma_of_z.dat",status="replace")
           open(3,file=dirname//"/B_1Mpc.dat",status="replace")
           open(4,file=dirname//"/x_ion.dat",status="replace")
        end if
        ! This part calculates the heating rate with B_n & n_B.
        ! Minoda added for zero-field calculation
        if(B_1Mpc.le.0.0.or.B_cut.le.0.0)then
           print*,"B_1Mpc =",B_1Mpc
           print*,"B_cut =",B_cut
           print*,"This is not positive."
           lambda_cut = 0.0
           m_B = 0.0
           E_ad_0 = 0.0
           E_dt_0 = 0.0
        else
           ! Minoda added: according to Fedeli & Moscardini, 2012
           if(heat_switch.eq.4)then
              lambda_cut = 2.0*pi/234.0*B_1Mpc
              m_B = 1.5
           else if(heat_switch.eq.6)then
              lambda_cut = 2.0*pi/52.0*B_cut
              print*,"lambda_c",lambda_cut
              lambda_cut = c*B_cut*nano/sqrt(4.0*pi*rho_baryon(a_now)/erg*a_rec)
              lambda_cut = lambda_cut/hubble(a_rec)/a_rec/Mpc
              print*,"lambda_c",lambda_cut
              lambda_cut = B_cut*nano/hubble(a_now)/sqrt(4.0*pi*omega_m*omega_b*1.0d-29)/Mpc
              print*,"lambda_c",lambda_cut
           else
              lambda_cut = (1.319e-3*B_1Mpc*B_1Mpc)**(1.0/(5.0+n_B))
              m_B = 2.0*(n_B+3.0)/(n_B+5.0)
           end if
           print*,"lambda_c =",lambda_cut," Mpc"
           ! data_in1(1,1,1)=(n_B+3.0)/(n_B+5.0)*(B_1Mpc**4.0)*4.0*pi*pi
           ! data_in1(1,1,1)=data_in1(1,1,1)/(lambda_cut**(2.0*(n_B+4.0)))
           ! print*,data_in1(1,1,1)
           data_in1(:,:,:)=0.0
           ! data_in1(:,:,:)=data_in1_init
           ! data_in1(:,:,:)=0.05*(1.05**loop)
           data_in2(:,:,:)=0.0

           if(heat_switch.eq.4)then
              E_ad_0 = B_1Mpc*B_1Mpc*B_1Mpc*B_1Mpc/lambda_cut/lambda_cut
              E_dt_0 = B_1Mpc*B_1Mpc/(8.0*pi)
           else if(heat_switch.eq.6)then
              E_ad_0 = B_cut*B_cut*B_cut*B_cut/lambda_cut/lambda_cut*(n_B+3.0)/(n_B+5.0)
              E_dt_0 = 0.0
           else
              E_ad_0 = 4.0*(2.0*pi)**(2.0*(n_B+4.0))*B_1Mpc*B_1Mpc*B_1Mpc*B_1Mpc
              E_ad_0 = E_ad_0/(n_B+3.0)/(n_B+5.0)/sp_G/sp_G/(lambda_cut**(2.0*(n_B+4.0)))
              E_dt_0 = (2.0*pi)**(n_B+2.0)*B_1Mpc*B_1Mpc/2.0/(n_B+3.0)/sp_G/(lambda_cut**(n_B+3.0))
           end if
        end if
        print*,"|rotBxB|^2 =",E_ad_0,"nG^4/Mpc^2"
        print*,"|B|^2/8pi =",E_dt_0,"nG^2"
        ! write(1,*)"### heating term Gamma = ",data_in1(1,1,1),"###"
     else if(input_switch.eq.1)then
        call system("mkdir -p "//dirname//"/output-"//loop_number)
        write(cpu_number,"(I3.3)")myrank
        ! open(1000+myrank,file="output/"//dirname//"/output-"//loop_number//"/delta_Tb_glob.dat"//cpu_number,status="replace")
        open(100,file="input/B_fields/"//index//"/absL2_"//loop_number//".dat",status="old")
        open(200,file="input/B_fields/"//index//"/divL_"//loop_number//".dat",status="old")
        do xi=1,ngrid
           do yi=1,ngrid
              do zi=1,ngrid
                 read(100,*)xgrid,ygrid,zgrid,data_in1(zi,yi,xi)
                 read(200,*)xgrid,ygrid,zgrid,data_in2(zi,yi,xi)
                 ! data_in1(zi,yi,xi)=0.001*data_in1(zi,yi,xi)
                 ! data_in2(zi,yi,xi)=0.001*data_in2(zi,yi,xi)
                 ! data_in2(zi,yi,xi)=0.0
              end do
              read(100,*)
              read(200,*)
           end do
        end do
        close(100)
        close(200)
     else
        call system("mkdir -p "//dirname)
        open(1,file=dirname//"/T_gas.dat",status="replace")
        open(4,file=dirname//"/x_ion.dat",status="replace")
     end if

     ! Set initial parameters
     step=0
     ioutput=0
     a=a_start
     la=log(a)
     do i=1,noutput
        print_check(i)=.true.
        print_time(i)=dble(i)*la/dble(noutput)
     end do

     Tg(:,:,:)=T_start
     x(:,:,:)=x_start
     fb(:,:,:)=fb_start
     dfb(:,:,:)=fb_start
     fm(:,:,:)=fm_start
     fm_out(:,:,:)=fm_start
     if(input_switch.eq.0)then
        Tg_one_one(:,:,:)=T_start
        Tg_two_one(:,:,:)=T_start
        Tg_two_two(:,:,:)=T_start
        E_ad=E_ad_0
        E_dt=E_dt_0
     end if

     if(myrank==0)print*,"      Start calculation loop = ",loop_number
     if(myrank==0)print*,"       step                  z                      dla"

     !modified at 16th (Tue) Oct. 2018
     if(input_switch.eq.1.and.ioutput==0)then
        call system("mkdir -p "//dirname//"/output-"//loop_number) ! modified at 20:15 5th(Thu) Jan 2017.
     end if

     do
        ! if(mod(step,10000).eq.0)print*,step,z(la),dla
        ! make output files
        if(input_switch.eq.1)then
           do i=noutput,1,-1
              if(print_check(i).and.la>print_time(i))then
                 call make_output(ioutput,la,Tg,x,fb,fm,fm_out)
                 print_check(i)=.false.
              end if
           end do
        else
           ! minoda changed 13th. June, 2018
           T_photon=T_rec*a_rec/exp(la)
           ! if(mod(step,1000)==0.and.z(la)<990.0)write(1,*)z(la),Tg(1,1,1),T_photon
           do i=noutput,1,-1
              if(print_check(i).and.la>print_time(i))then
                 ! write(1,*)z(la),Tg(1,1,1),T_photon
                 write(1,"(f15.5,f15.5,f15.5,f15.5)")z(la),Tg(1,1,1),T_photon,Ts
                 write(2,*)z(la),fdata1
                 write(3,*)z(la),B_1Mpc*sqrt(fdata1)
                 write(4,*)z(la),x(1,1,1)
                 print_check(i)=.false.
                 print*,step,z(la),dla
              end if
           end do
           if(z(la)<z_abs+0.01.and.z(la)>z_abs-0.01)then
              if(Tg(1,1,1).ge.T_photon)print*,z(la),Tg(1,1,1)-T_photon,"emission"
              if(Tg(1,1,1).le.T_photon)print*,z(la),Tg(1,1,1)-T_photon,"absorption"
           end if
        end if

        ! Calculate appropriate time step
        if(input_switch.eq.0.or.input_switch.eq.2)then
           dla=1.0e-4
           ! print*,"minoda (1)"
           ! dla=dla_start
           ! ! do
           ! !    dla=dla*0.5
           ! !    if(mod(step,1000).eq.0)print*,dla
           ! !    if(abs(dT_gas(la,Tg(1,1,1),x(1,1,1),fb_start,fb_start))<Tg(1,1,1))exit
           ! ! end do
           ! minoda modified 22nd Nov. 2018
           do
              dla=dla*0.5
              call calc_dTg(la,x,Tg,dTg_one_one)
              ! print*,"minoda (2)"
              Tg_one_one(1,1,1)=Tg(1,1,1)+dTg_one_one(1,1,1)
              if(Tg_one_one(1,1,1)>0.0)exit
           end do
           ! if(dTg_one_one(1,1,1)>0.01*Tg(1,1,1))print*,"dla=",dla
           ! print*,"minoda (3)"
           call calc_dTg(la,x,Tg,dTg_two_one)
           Tg_two_one(1,1,1)=Tg(1,1,1)+dTg_two_one(1,1,1)
           call calc_dTg(la+dla,x,Tg_two_one,dTg_two_two)
           Tg_two_two(1,1,1)=Tg_two_one(1,1,1)+dTg_two_two(1,1,1)
           do
              ! print*,"Tg_one_one",Tg_one_one(1,1,1)
              ! print*,"Tg_two_two",Tg_two_two(1,1,1)
              ! print*,"log(a_start)/dla",log(a_start)/dla
              ! print*,"Tg(1,1,1)",Tg(1,1,1)
              ! if(step.eq.0)print*,"T_gas 0",Tg(1,1,1)
              ! if(step.eq.0)print*,"T_gas error",Tg_one_one(1,1,1)-Tg_two_two(1,1,1)
              ! if(step.eq.0)print*,"predicted step number",log(a_start)/dla
              ! if(step.eq.0)print*,"predicted total error",abs((Tg_one_one(1,1,1)-Tg_two_two(1,1,1))*log(a_start)/dla)
              if(abs((Tg_one_one(1,1,1)-Tg_two_two(1,1,1))*log(a_start)/dla)<0.05*Tg(1,1,1))exit
              Tg_one_one(1,1,1)=Tg_two_one(1,1,1)
              dla=0.5*dla
              call calc_dTg(la,x,Tg,dTg_two_one)
              Tg_two_one(1,1,1)=Tg(1,1,1)+dTg_two_one(1,1,1)
              call calc_dTg(la+dla,x,Tg_two_one,dTg_two_two)
              Tg_two_two(1,1,1)=Tg_two_one(1,1,1)+dTg_two_two(1,1,1)
           end do
           ! if(dla<1.0e-6)print*,step,z(la),dla
        else
           if(step.le.10)then
              dla=dla_start
           else
              dla_old=1.0d-2
              !           if(step<10000.and.mod(step,100).eq.0)print*,"Now step =",step,", so write down dla for debug"
              do xi=dgrid*myrank+1,dgrid*myrank+dgrid
                 do yi=1,ngrid
                    do zi=1,ngrid
                       x1=dx_ion(la,x(zi,yi,xi),Tg(zi,yi,xi),fb(zi,yi,xi))
                       x2=dx_ion(la+0.5*dla,x(zi,yi,xi)+0.5*dla*x1,Tg(zi,yi,xi),fb(zi,yi,xi))
                       x3=dx_ion(la+0.5*dla,x(zi,yi,xi)+0.5*dla*x2,Tg(zi,yi,xi),fb(zi,yi,xi))
                       x4=dx_ion(la+dla,x(zi,yi,xi)+dla*x3,Tg(zi,yi,xi),fb(zi,yi,xi))
                       dla_x=x(zi,yi,xi)*6.0/abs(x1+2.0*x2+2.0*x3+x4)
                       if(dx_ion(la,x(zi,yi,xi),Tg(zi,yi,xi),max(-0.9,fb(zi,yi,xi)))==0.0)then
                          dla_x = dla_old
                       else
                          dla_x = abs(x(zi,yi,xi)/dx_ion(la,x(zi,yi,xi),Tg(zi,yi,xi),max(-0.9,fb(zi,yi,xi))))
                       end if
                       if(dT_gas(la,Tg(zi,yi,xi),x(zi,yi,xi),fb(zi,yi,xi),max(-0.9,dfb(zi,yi,xi)))==0.0)then
                          dla_Tg = dla_old
                       else
                          dla_Tg = abs(Tg(zi,yi,xi)/dT_gas(la,Tg(zi,yi,xi),x(zi,yi,xi),max(-0.9,fb(zi,yi,xi)),dfb(zi,yi,xi)))
                       end if
                       dla_old = min(dla_old,min(dla_x,dla_Tg)*cfl)
                    end do
                 end do
              end do
              call mpi_allreduce(dla_old,dla,1,mpi_double_precision,mpi_min,mpi_comm_world,ierr)
           end if
        end if

        ! write time step for check
        if(myrank==0.and.mod(step,100000)==0)then
        end if

        ! Calculate differential values
        call calc_delta(la,fb,dfb,fm,fm_out)
        do xi=dgrid*myrank+1,dgrid*myrank+dgrid
           do yi=1,ngrid
              do zi=1,ngrid
!                 if(myrank==0.and.step>checkstep)print*,"start runge-kutta" ! for check
                 x1=dx_ion(la,x(zi,yi,xi),Tg(zi,yi,xi),max(-0.9,fb(zi,yi,xi)))
                 x2=dx_ion(la+0.5*dla,x(zi,yi,xi)+0.5*dla*x1,Tg(zi,yi,xi),max(-0.9,fb(zi,yi,xi)))
                 x3=dx_ion(la+0.5*dla,x(zi,yi,xi)+0.5*dla*x2,Tg(zi,yi,xi),max(-0.9,fb(zi,yi,xi)))
                 x4=dx_ion(la+dla,x(zi,yi,xi)+dla*x3,Tg(zi,yi,xi),max(-0.9,fb(zi,yi,xi)))
                 dx(zi,yi,xi)=(x1+2.0*x2+2.0*x3+x4)*dla/6.0
!                 print*,"this is check",fb(1,1,1),dfb(1,1,1)
!                 if(myrank==0.and.step>checkstep)print*,"finish runge-kutta",dx(zi,yi,xi)
              end do
           end do
        end do
        call calc_dTg(la,x,Tg,dTg)

        ! Update
        ! print*,step,z(la),dla
        step=step+1
        la=la+dla
        ! minoda comment 10.16 for bug-fix
        do xi=dgrid*myrank+1,dgrid*myrank+dgrid
           do yi=1,ngrid
              do zi=1,ngrid
                 x(zi,yi,xi)=x(zi,yi,xi)+dx(zi,yi,xi)
                 if(x(zi,yi,xi).ge.x_start)x(zi,yi,xi)=x_start
                 Tg(zi,yi,xi)=Tg(zi,yi,xi)+dTg(zi,yi,xi)
              end do
           end do
        end do

        ! minoda added 1st. Oct. 2018
        ! reduce the magnetic energy due to the ambipolar diffusion
        if(input_switch.eq.0.and.B_1Mpc>0.0.and.B_cut>0.0)then
           call dissipate_energy(la,x(1,1,1),Tg(1,1,1),E_ad,E_dt)
        end if
        ! minoda added 11th Oct. 2023
        ! calculate the 21cm differential brightness temperature
        call calc_deltaTb(la,Tg(1,1,1),T_photon,Ts,nH,x(1,1,1))
        if(la>=log(a_end))exit
     end do

     if(input_switch.eq.0)then
        ! minoda changed 13th. June, 2018
        ! write(1,*),z(la),Tg(1,1,1),x(1,1,1)
        T_photon=T_rec*a_rec/exp(la)
        ! write(1,*),z(la),Tg(1,1,1),T_photon
        write(1,"(f15.5,f15.5,f15.5,f15.5)")z(la),Tg(1,1,1),T_photon,Ts
        write(2,*)z(la),fdata1
        write(3,*)z(la),B_1Mpc*sqrt(fdata1)
        write(4,*)z(la),x(1,1,1)
        close(1)
        close(2)
        close(3)
        close(4)
     else
        call make_output(ioutput,la,Tg,x,fb,fm,fm_out)
        ! close(1000+myrank)
     end if
  end do


  call mpi_finalize(ierr)

contains

  double precision function z(la)
    implicit none
    double precision::la
    z=1.0/exp(la)-1.0
  end function z

  double precision function hubble(a)
    implicit none
    double precision::a,a_inv
    a_inv=1.0/a
    hubble=h0*sqrt(omega_m*a_inv*a_inv*a_inv+omega_l)
  end function hubble

  double precision function rho_baryon(a)
    implicit none
    double precision::a,a_inv,rho_matter
    a_inv=1.0/a
    rho_matter=1.0540d4*reduced_h*reduced_h*a_inv*a_inv*a_inv
    rho_baryon=rho_matter*(omega_b/omega_m) !!! note: this value has (eV/cm^3) dimension !!!
  end function rho_baryon

  double precision function rho_cdm(a)
    implicit none
    double precision::a,a_inv,rho_matter
    a_inv=1.0/a
    rho_matter=1.0540d4*reduced_h*reduced_h*a_inv*a_inv*a_inv
    rho_cdm=rho_matter*(1.0-omega_b/omega_m) !!! note: this value has (eV/cm^3) dimension !!!
  end function rho_cdm

  double precision function chi(Temp)
    implicit none
    double precision::Temp
    ! print*,"called chi(Temp), Temp=",Temp
    if(Temp>0)then
       chi = 1.9e14*(Temp**0.375)   ! drag coefficient (cm^3/g/sec)
    else
       chi = 3.5e13
    end if
    ! chi = 3.5e13
    ! chi = 1.9e14*(Temp**0.375)   ! drag coefficient (cm^3/g/sec)
  end function chi


  subroutine calc_delta(la,fb,dfb,fm,fm_out)
    implicit none
    double precision::la,a,eta,delta,ddelta,source,delta_matter
    double precision,parameter::divL0=nano*nano/Mpc/Mpc
    double precision,dimension(ngrid,ngrid,ngrid)::fb,dfb,fm,fm_out
    ! compute local baryon density in real space
    a=exp(la)
    if(a<a_rec)then
       fb(:,:,:)=0.0
       dfb(:,:,:)=0.0
       fm(:,:,:)=0.0
       fm_out(:,:,:)=0.0
    else
       eta=a/a_rec
       delta=(3.0*eta+2.0/sqrt(eta*eta*eta)-15.0*log(eta))*omega_b/omega_m
       delta=delta + 15.0*log(eta) + 30.0*(1.0-omega_b/omega_m)/sqrt(eta) - (30.0-25.0*omega_b/omega_m)
       delta=delta*2.0/15.0/hubble(a)/hubble(a)
       ddelta=(eta-1.0/sqrt(eta*eta*eta)-5.0)*omega_b/omega_m + 5.0*(1.0-(1.0-omega_b/omega_m)/sqrt(eta))
       ddelta=ddelta*2.0/5.0/hubble(a)
       ! minoda modified at 28th. Feb. 2019
       delta_matter=(2.0/sqrt(eta*eta*eta)+3.0*eta-5.0)*omega_b/omega_m*2.0/15.0/hubble(a)/hubble(a)
!       delta=delta*2.0/15.0/h0/h0
       do xi=dgrid*myrank+1,dgrid*myrank+dgrid
          do yi=1,ngrid
             do zi=1,ngrid
                source=data_in2(zi,yi,xi)*divL0/(4.0*pi*rho_baryon(a_now)*c2_inv/erg*a*a*a)
!                source=data_in2(zi,yi,xi)*divL0/(4.0*pi*rho_baryon(a_now)*c2_inv/erg)
!                print*,rho_baryon(a_now)*c2_inv/erg,"g/cc"
!                print*,h0*h0,hubble(a)*hubble(a)*a*a*a
                fb(zi,yi,xi)=delta*source
                dfb(zi,yi,xi)=ddelta*source
                !minoda added on 25th. Feb. 2019
                fm(zi,yi,xi)=delta_matter*source
                fm_out(zi,yi,xi)=delta_matter*source
                !!!! condition delta > -1 !!!!
                if(fb(zi,yi,xi)<-0.9)then
                   fb(zi,yi,xi)=-0.9
                   dfb(zi,yi,xi)=0.0
                   fm(zi,yi,xi)=-0.9
!                   print*,fb(zi,yi,xi),dfb(zi,yi,xi)
                end if
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             end do
          end do
       end do
    end if
  end subroutine calc_delta


  ! subroutine calc_dx(la,x,Tg,dx)
  !   implicit none
  !   double precision::la
  !   double precision,dimension(ngrid,ngrid,ngrid)::x, Tg,  dx
  !   double precision::x1, x2, x3, x4
  ! end subroutine calc_dx

  subroutine calc_dTg(la,x,Tg,dTg)
    implicit none
    double precision::la
    double precision,dimension(ngrid,ngrid,ngrid)::x, Tg, dTg
    double precision::Tg1, Tg2, Tg3, Tg4
    do xi=dgrid*myrank+1,dgrid*myrank+dgrid
       do yi=1,ngrid
          do zi=1,ngrid
             Tg1=dT_gas(la,Tg(zi,yi,xi),x(zi,yi,xi),max(-0.9,fb(zi,yi,xi)),dfb(zi,yi,xi))
             Tg2=dT_gas(la+0.5*dla,Tg(zi,yi,xi)+0.5*dla*Tg1,x(zi,yi,xi),max(-0.9,fb(zi,yi,xi)),dfb(zi,yi,xi))
             Tg3=dT_gas(la+0.5*dla,Tg(zi,yi,xi)+0.5*dla*Tg2,x(zi,yi,xi),max(-0.9,fb(zi,yi,xi)),dfb(zi,yi,xi))
             ! if(Tg3<0.0)then
             !    dTg(zi,yi,xi)=-2*Tg(zi,yi,xi)
             !    exit
             ! end if
             ! print*,"minoda (1-1)"
             Tg4=dT_gas(la+dla,Tg(zi,yi,xi)+dla*Tg3,x(zi,yi,xi),max(-0.9,fb(zi,yi,xi)),dfb(zi,yi,xi))
             ! print*,"minoda (1-2)"
             dTg(zi,yi,xi)=(Tg1+2.0*Tg2+2.0*Tg3+Tg4)*dla/6.0
             ! print*,"minoda (1-3)"
          end do
       end do
    end do
  end subroutine calc_dTg

  subroutine make_output(ioutput,la,Tg,x,fb,fm,fm_out)
    implicit none
    integer::ioutput
    double precision::la
    double precision::a,a_inv,T_photon
    double precision,dimension(ngrid,ngrid,ngrid)::Tg,x,fb,dTb,fm,fm_out
    double precision::dTb_g=0.0
    double precision::kSZ_fac
    ioutput=ioutput+1
    write(cpu_number,"(I3.3)")myrank
    write(output_number,"(I2.2)")ioutput
    if(myrank==0)print*,"       step                  z           ioutput"
    if(myrank==0)print*,step,z(la),ioutput
    if(myrank==0)print*,"       step                  z                      dla"
    call system("mkdir -p "//dirname//"/time"//output_number)
    open(200000+100*myrank+ioutput,file=dirname//"/time"//output_number//"/delta_Tb.dat"//cpu_number,status="replace")
    open(300000+100*myrank+ioutput,file=dirname//"/time"//output_number//"/pressure_e.dat"//cpu_number,status="replace")
    open(400000+100*myrank+ioutput,file=dirname//"/time"//output_number//"/ion_rate.dat"//cpu_number,status="replace")
    open(500000+100*myrank+ioutput,file=dirname//"/time"//output_number//"/diff_temp.dat"//cpu_number,status="replace")
    open(700000+100*myrank+ioutput,file=dirname//"/time"//output_number//"/number_baryon.dat"//cpu_number,status="replace")
    open(800000+100*myrank+ioutput,file=dirname//"/time"//output_number//"/delta_baryon.dat"//cpu_number,status="replace")
    ! minoda added 5th. Feb. 2019
    open(900000+100*myrank+ioutput,file=dirname//"/time"//output_number//"/Doppler.dat"//cpu_number,status="replace")
    open(1000000+100*myrank+ioutput,file=dirname//"/time"//output_number//"/Doppler_v.dat"//cpu_number,status="replace")
    open(1100000+100*myrank+ioutput,file=dirname//"/time"//output_number//"/delta_matter.dat"//cpu_number,status="replace")
    open(1200000+100*myrank+ioutput,file=dirname//"/time"//output_number//"/delta_matter_out.dat"//cpu_number,status="replace")
    ! calculate comoving distance (in the unit of Mpc)
    la_old=la
    dla_old=0.0001*la_old
    d_com=0.0
    do
       a_old=exp(la_old)
       la_old=la_old-dla_old
       d_com=d_com-c*dla_old/a_old/hubble(a_old)
       if(la_old.ge.0.0)exit
    end do
    d_com=d_com/Mpc             ! modified at 18:08, 1st(Sun), Jan.
    ! calculate average of y-parameter/(ion_rate * T_e * n_b)
    a=exp(la)
    w=sigma*kB/mc2_e*Mpc     ! updated at 17:10, 2nd(Thu), Mar. 2017
    ! minoda added for test 2019.08.31
    ! write(300000+100*myrank+ioutput,*)la,d_com,w
    write(300000+100*myrank+ioutput,*)la,d_com,sigma*Mpc*1.0d-4*rho_baryon(a)/mc2_p
    ! minoda added 5th. Feb. 2019
    ! Factor for the Doppler effect
    write(400000+100*myrank+ioutput,*)la,d_com,sigma*Mpc*hubble(a)*Mpc*rho_baryon(a)/mc2_p
    write(900000+100*myrank+ioutput,*)la,d_com,sigma*Mpc*hubble(a)*Mpc/c
    write(1000000+100*myrank+ioutput,*)la,d_com,sigma*Mpc*1.0d-4*rho_baryon(a)/mc2_p
    ! minoda changed 30th. Aug. 2019
    ! hubble has unit of (1/cm), therefore multiplied by Mpc to convert (cm/sec/Mpc)
    ! minoda added 12th. Sept. 2019
    kSZ_fac=(a*hubble(a)*Mpc/c)

    a_inv=1.0/a
    T_photon=T_rec*a_rec*a_inv
    dTb_g=0.0
    do xi=dgrid*myrank+1,dgrid*myrank+dgrid
       do yi=1,ngrid
          do zi=1,ngrid
             nH=rho_baryon(a)/mc2_p*(1.0+fb(zi,yi,xi)) ! added at 17:10, 2nd(Thu), Mar. 2017
!             nH=rho_baryon(a)/mc2_p ! added at 17:10, 2nd(Thu), Mar. 2017
             dTb(zi,yi,xi)=28.0*sqrt(0.1*a_inv)*(1.0-x(zi,yi,xi))*(1.0+fb(zi,yi,xi))*(Tg(zi,yi,xi)-T_photon)/Tg(zi,yi,xi)
             dTb_g = dTb_g+Tg(zi,yi,xi)/ngrid/ngrid/dgrid
             ! dTb_g = dTb_g+dTb(zi,yi,xi)/ngrid/ngrid/dgrid
             write(200000+100*myrank+ioutput,*)xi,yi,zi,dTb(zi,yi,xi)
             ! added at 15:10, 19th(Thu), June 2018
             ! minoda added for test 2019.08.31
             write(300000+100*myrank+ioutput,*)xi,yi,zi,kB*(Tg(zi,yi,xi)-T_photon)/mc2_e ! updated at 17:10, 2nd(Thu), Mar. 2017
             ! write(300000+100*myrank+ioutput,*)xi,yi,zi,x(zi,yi,xi)*(Tg(zi,yi,xi)-T_photon)*nH ! updated at 17:10, 2nd(Thu), Mar. 2017
             write(400000+100*myrank+ioutput,*)xi,yi,zi,x(zi,yi,xi)
             write(500000+100*myrank+ioutput,*)xi,yi,zi,Tg(zi,yi,xi),T_photon
!             write(50000+100*myrank+ioutput,*)xi,yi,zi,Tg(zi,yi,xi)-T_photon
             write(700000+100*myrank+ioutput,*)xi,yi,zi,nH
             write(800000+100*myrank+ioutput,*)xi,yi,zi,fb(zi,yi,xi)
             ! minoda added 5th. Feb. 2019
             write(900000+100*myrank+ioutput,*)xi,yi,zi,x(zi,yi,xi)*nH*fb(zi,yi,xi)
             write(1000000+100*myrank+ioutput,*)xi,yi,zi,3.5*kSZ_fac*kSZ_fac*fb(zi,yi,xi)
             write(1100000+100*myrank+ioutput,*)xi,yi,zi,fm(zi,yi,xi)
             write(1200000+100*myrank+ioutput,*)xi,yi,zi,fm_out(zi,yi,xi)
          end do
       end do
       write(300000+100*myrank+ioutput,*)""
       write(400000+100*myrank+ioutput,*)""
       write(500000+100*myrank+ioutput,*)""
       write(700000+100*myrank+ioutput,*)""
       write(800000+100*myrank+ioutput,*)""
       write(900000+100*myrank+ioutput,*)""
       write(1000000+100*myrank+ioutput,*)""
       write(1100000+100*myrank+ioutput,*)""
       write(1200000+100*myrank+ioutput,*)""
    end do
    ! write(1000+myrank,*)a_inv,dTb_g
    close(300000+100*myrank+ioutput)
    close(400000+100*myrank+ioutput)
    close(500000+100*myrank+ioutput)
    close(700000+100*myrank+ioutput)
    close(800000+100*myrank+ioutput)
    close(900000+100*myrank+ioutput)
    close(1000000+100*myrank+ioutput)
    close(1100000+100*myrank+ioutput)
    close(1200000+100*myrank+ioutput)
  end subroutine make_output

  subroutine calc_deltaTb(la,T_mat,T_CMB,T_spin,nH,x_ion)
    implicit none
    double precision::la,T_mat,T_CMB,T_spin,nH,x_ion
    double precision::KHH,KeH,KpH,lnt
    double precision,parameter::A10 = 2.85d-15
    double precision::xcoll
    a=exp(la)
    nH=rho_baryon(a)/mc2_p
    if(T_mat<1.0)then
       KHH = 1.380d-13
       KeH = 0.239d-9
       KpH = 0.4028d-9
    else
       lnt=dLog(T_mat)
       KHH=-4.34289510d1+2.16894228*lnt-6.73666541*lnt*lnt
       KHH=KHH+7.38234512 *lnt*lnt*lnt
       KHH=KHH-3.46789178*lnt*lnt*lnt*lnt
       KHH=KHH+8.87827433d-1 *lnt*lnt*lnt*lnt*lnt
       KHH=KHH-1.34182754d-1*lnt*lnt*lnt*lnt*lnt*lnt
       KHH=KHH+1.19968634d-2 *lnt*lnt*lnt*lnt*lnt*lnt*lnt
       KHH=KHH-5.88215431d-4*lnt*lnt*lnt*lnt*lnt*lnt*lnt*lnt
       KHH=KHH+1.22111242d-5*lnt*lnt*lnt*lnt*lnt*lnt*lnt*lnt*lnt
       KHH=dexp(KHH)*1.0d6

       KeH=-3.59759221d1+5.17874468d-1*lnt+1.48896712d-3*lnt*lnt
       KeH=KeH+3.66801839d-4*lnt*lnt*lnt
       KeH=KeH-3.86512231d-4*lnt*lnt*lnt*lnt
       KeH=dexp(KeH)*1.0d6

       KpH=-3.54503981d1+2.86461006d-1*lnt-1.22528505d-1*lnt*lnt
       KpH=KpH-6.79646651d-2 *lnt*lnt*lnt
       KpH=KpH-3.93059229d-2*lnt*lnt*lnt*lnt
       KpH=KpH-6.92090993d-3 *lnt*lnt*lnt*lnt*lnt
       KpH=KpH+5.30493146d-4*lnt*lnt*lnt*lnt*lnt*lnt
       KpH=KpH-1.52387862d-5 *lnt*lnt*lnt*lnt*lnt*lnt*lnt
       KpH=dexp(KpH)*1.0d6
    end if

    if(T_mat>10000.0)then
       KHH = 7.870e-10
    end if
    if(T_mat>20000.0)then
       KpH = 2.201d-9
    end if
    if(T_mat>100000.0)then
       KeH = 3.72254d-9
    end if

    xcoll = 0.0628/T_CMB/A10*(nH*(1.0-x_ion)*KHH+nH*x_ion*KeH+nH*x_ion*KpH)
    T_spin=(1.0+xcoll)*T_CMB*T_mat/(T_mat+xcoll*T_CMB)
  end subroutine calc_deltaTb

  subroutine dissipate_energy(la,x_ion,T_mat,E_ambi,E_decay)
    implicit none
    double precision::la,x_ion,T_mat,E_ambi,E_decay
    double precision::a_inv, rho_b_now, B_sq, time_dt
    ! double precision::la,data_input1,x_ion,T_mat,E_ambi,E_decay
    a=exp(la)
    a_inv=1.0/a
    ! minoda changed on 18.10.17
    rho_b_now = (rho_baryon(a_now)/erg)*c2_inv
    B_sq = (nano*B_amp*a_inv*a_inv)*(nano*B_amp*a_inv*a_inv)
    ! minoda changed on 18.10.26
    ! time_dt = lambda_cut*Mpc*sqrt(rho_baryon(a)/erg*c2_inv/2.0/E_dt_0/B_sq)
    ! dfdata1=4.4*(2.0*pi)**(n_B+4.0)/(n_B+5.0)/sp_G
    ! ! dfdata1=dfdata1/Gamma_function((n_B+3.0)/2.0)/Gamma_function((n_B+3.0)/2.0)
    ! dfdata1=dfdata1*(1.0-x_ion)/x_ion/(T_mat**0.375)
    ! minoda changed
    ! if(heat_switch>3.and.a<a_rec)then
    if(heat_switch>3.and.a<a_rec)then
       dfdata1=0.0
       dfdata2=0.0
    else
       ! minoda wrote (2019/6/7): modified increasing cut-off scale for time-scale of decaying turbulence
       ! time_dt = a*lambda_cut*Mpc*sqrt(rho_baryon(a)/erg*c2_inv/2.0/E_dt/B_sq)
       time_dt = a*lambda_cut*fdata1*Mpc*sqrt(rho_baryon(a)/erg*c2_inv/2.0/E_dt/B_sq)
       dfdata1=(1.0-x_ion)/x_ion/(16.0*pi*pi*chi(T_mat)*rho_b_now*rho_b_now)
       dfdata2=(log(1.0+time_dt/time_rec)**m_B)/((log(1.0+time_dt/time_rec)+1.5*log(a/a_rec))**(1.0+m_B))
       if(diss_switch.eq.0)then
          dfdata1=dfdata1/hubble(a)*E_ad_0/E_dt_0*(nano*nano/Mpc/Mpc)*fdata1*fdata1*dla
          dfdata2=dfdata2*1.5*m_B*fdata1*dla
          ! minoda modified 2019/05/21
          ! dfdata2=dfdata2*1.5*m_B*a*a*a*a*fdata1*dla
       else if(diss_switch.eq.1)then
          ! minoda put a bug (nB+6 => nB+5)
          dfdata1=dfdata1/hubble(a)*E_ad_0/E_dt_0*(nano*nano/Mpc/Mpc)*(fdata1**(n_B+6.0))*dla/(n_B+3.0)
          dfdata2=dfdata2*1.5*m_B*fdata1*dla/(n_B+3.0)
          ! minoda wrote to test the importance of (n+3) term
          ! dfdata1=dfdata1/hubble(a)*E_ad_0/E_dt_0*(nano*nano/Mpc/Mpc)*(fdata1**(n_B+6.0))*dla
          ! dfdata2=dfdata2*1.5*m_B*fdata1*dla
       end if
    end if
    ! dfdata2=2.0*(n_B+3.0)/(n_B+5.0)/log(3.8e9/a_start*(a**1.5))*fdata1*dla
    ! print*,"fdata1",fdata1
    ! print*,"dfdata1",dfdata1
    ! print*,"dfdata2",dfdata2
    if(heat_switch.eq.0)then
       fdata1=fdata1
    else if(heat_switch.eq.1)then
       fdata1=fdata1-dfdata1
    else if(heat_switch.eq.2)then
       fdata1=fdata1-dfdata2
    else if(heat_switch.eq.3.or.heat_switch.eq.4)then
       fdata1=fdata1-dfdata1-dfdata2
    else if(heat_switch.eq.5.or.heat_switch.eq.6)then
       fdata1=fdata1
    end if
    ! data_input1=data_in1_init
       if(diss_switch.eq.0)then
          E_ambi = E_ad_0*fdata1*fdata1
          E_decay = E_dt_0*fdata1
       else if(diss_switch.eq.1)then
          E_ambi = E_ad_0*(fdata1**(2.0*n_B+8.0))
          E_decay = E_dt_0*(fdata1**(n_B+3.0))
       end if
    ! ! if(mod(step,200)==0.and.z(la)>950.0)write(2,*)z(la),sqrt(fdata1)
    ! ! if(mod(step,20)==0.and.z(la)<950.0.and.z(la)>900.0)write(2,*)z(la),sqrt(fdata1)
    ! ! if(mod(step,2)==0.and.z(la)<900.0.and.z(la)>500.0)write(2,*)z(la),sqrt(fdata1)
    ! ! if(mod(step,1)==0.and.z(la)<500.0.and.z(la)>10.0)write(2,*)z(la),sqrt(fdata1)
    ! if(mod(step,200)==0.and.z(la)>950.0)write(2,*)z(la),data_input1,x_ion
    ! if(mod(step,20)==0.and.z(la)<950.0.and.z(la)>900.0)write(2,*)z(la),data_input1,x_ion
    ! if(mod(step,10)==0.and.z(la)<900.0)write(2,*)z(la),data_input1,x_ion
    ! ! if(mod(step,1)==0.and.z(la)<500.0.and.z(la)>100.0)write(2,*)z(la),data_input1,x_ion

    ! if(mod(step,200)==0.and.z(la)>950.0)write(2,*)z(la),fdata1
    ! if(mod(step,20)==0.and.z(la)<950.0.and.z(la)>900.0)write(2,*)z(la),fdata1
    ! if(mod(step,10)==0.and.z(la)<900.0)write(2,*)z(la),fdata1
  end subroutine dissipate_energy

  double precision function dT_gas(la,T_gas,x_ion,delta_baryon,delta_delta)
    implicit none
    double precision::delta_baryon,delta_delta
    double precision::la, a_inv, T_gas, x_ion, T_photon, rho_photon, rho_b2
    double precision::dT_gas1, dT_gas2, dT_gas3, dT_gas4
    double precision::eps
    double precision::rotBxB,gamma
    double precision::Theta_ff,Eta_rec,Phi,Zeta
    double precision::B_sq
    double precision::time_dt
    double precision,parameter::chi_SS = 3.5e13
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! if(delta_baryon<-0.9)then
    !    delta_baryon=-0.9
    !    delta_delta=0.0
    ! end if
    ! print*,"minoda (1-1-1)"
    if(x_ion.ge.x_start)x_ion=x_start
    a=exp(la)
    a_inv=1.0/a
    T_photon=T_rec*a_rec*a_inv
    nH=rho_baryon(a)/mc2_p*(1.0+delta_baryon)
    rho_photon=SB*T_photon*T_photon*T_photon*T_photon
    rotBxB = (nano*B_amp*a_inv*a_inv)*(nano*B_amp*a_inv*a_inv)/(Mpc*boxsize*a)
    rho_b2 = (rho_baryon(a)/erg)*c2_inv*(1.0+delta_baryon)*(rho_baryon(a)/erg)*c2_inv*(1.0+delta_baryon)
    ! minoda wrote for check 4th. Oct.
    ! gamma = data_in1_init*rotBxB*rotBxB/(16.0*pi*pi*chi(T_gas)*rho_b2)*(1.0-x_ion)/x_ion
    ! minoda changed on 18.10.04
    if(input_switch.eq.0)then
       B_sq = (nano*B_amp*a_inv*a_inv)*(nano*B_amp*a_inv*a_inv)
       if(B_1Mpc.le.0.0.or.B_cut.le.0.0)then
          time_dt = 0.0
          gamma = 0.0
       else
          if(E_dt.eq.0.0)then
             time_dt = 0.0
          else
             ! minoda wrote (2019/6/7): modified increasing cut-off scale for time-scale of decaying turbulence
             ! time_dt = a*lambda_cut*Mpc*sqrt(rho_baryon(a)/erg*c2_inv/2.0/E_dt/B_sq)
             time_dt = a*lambda_cut*fdata1*Mpc*sqrt(rho_baryon(a)/erg*c2_inv/2.0/E_dt/B_sq)
          end if

          if(heat_switch>3.and.a<a_rec)then
             gamma=0.0
          else if(heat_switch.eq.0)then
             gamma = 0.0
          else if(heat_switch.eq.1)then
             gamma = E_ad*rotBxB*rotBxB/(16.0*pi*pi*chi(T_gas)*rho_b2)*(1.0-x_ion)/x_ion
             ! minoda changed 18.10.17
          else if(heat_switch.eq.2)then
             gamma = ((log(1.0+time_dt/time_rec))**m_B)
             gamma = gamma/((log(1.0+time_dt/time_rec)+1.5*log(a/a_rec))**(1.0+m_B))
             gamma = gamma*1.5*m_B*hubble(a)*E_dt*B_sq
             ! gamma = 3.0*(n_B+3.0)/(n_B+5.0)/log(3.8e9*a*sqrt(1.0+z_rec))*hubble(a)*E_dt*B_sq
          else if(heat_switch.eq.3.or.heat_switch.eq.4)then
             gamma = (log(1.0+time_dt/time_rec))**m_B
             gamma = gamma/((log(1.0+time_dt/time_rec)+1.5*log(a/a_rec))**(1.0+m_B))
             gamma = gamma*1.5*m_B*hubble(a)*E_dt*B_sq
             gamma = gamma + E_ad*rotBxB*rotBxB/(16.0*pi*pi*chi(T_gas)*rho_b2)*(1.0-x_ion)/x_ion
          else if(heat_switch.eq.5)then
             gamma = E_ad*rotBxB*rotBxB/(16.0*pi*pi*chi_SS*rho_b2)*(1.0-x_ion)/x_ion
          else if(heat_switch.eq.6)then
             gamma = 7.0*E_ad*rotBxB*rotBxB/(768.0*pi*pi*chi_SS*rho_b2*x_ion)
          else
             gamma = 0.0
          end if
       end if
    else
       gamma = data_in1(zi,yi,xi)*rotBxB*rotBxB/(16.0*pi*pi*chi(T_gas)*rho_b2)*(1.0-x_ion)/x_ion
    end if

    ! print*,"minoda (1-1-2)"
    eps= (x_ion/(1.0+x_ion))*rho_photon*((8.0*sigma*c)/(3.0*mc2_e))
    dT_gas1 = -2.0*T_gas
    ! modified 5th May 2017
!    dT_gas2 = (T_photon-T_gas)*eps
    ! minoda wrote for check 1015
    ! print*,"z",1.0/a
    ! print*,"delta_delta/(1.0+delta_baryon)*T_gas",delta_delta/(1.0+delta_baryon)*T_gas
    ! print*,"eps",eps
    ! print*,"(T_photon-T_gas)*eps",(T_photon-T_gas)*eps
    ! print*,"hubble",hubble(a)
    dT_gas2 = delta_delta/(1.0+delta_baryon)*T_gas + (T_photon-T_gas)*eps

    ! if(x(zi,yi,xi).ge.x_start)then
    !    dT_gas3 = 0.0
    ! else
    !    dT_gas3 = gamma/(1.5*(kB/erg)*nH)
    ! end if
    dT_gas3 = gamma/(1.5*(kB/erg)*nH)

    ! wrote 05:47 14th Nov 2018
    ! X-ray heating term
    if(heat_switch.eq.4)then
       dT_gas3 = dT_gas3 + 0.0
       ! dT_gas3 = dT_gas3 + 0.5d3*hubble(a)*(0.1*a_inv)
    end if

    ! print*,"minoda (1-1-3)"
    if(cool_switch.eq.0.or.T_gas.le.0.0)then
       dT_gas4 = 0.0
    else
       Theta_ff = 1.42d-27*sqrt(T_gas)
       Eta_rec = 6.5d-27*sqrt(T_gas)*(T_gas*1.0d-3)**(-0.2)/(1.0+(T_gas*1.0d-6)**0.7)
       Phi = 7.5d-19/(1.0+sqrt(T_gas*1.0d-5))*exp(-1.18d5/T_gas)
       Zeta = 1.27d-21*sqrt(T_gas)/(1.0+sqrt(1.0d-5*T_gas))*exp(-1.58d5/T_gas)
       dT_gas4 = (Theta_ff+Eta_rec)*x_ion + (Phi+Zeta)*(1.0-x_ion)
       dT_gas4 = -2.0/3.0*x_ion*nH/(kB/erg)*dT_gas4
    end if
    ! print*,"minoda (1-1-4)"
!    if(mod(step,100).eq.0)print*,dT_gas1,dT_gas2/hubble(a),dT_gas3/hubble(a)
    ! if(z(la)>z_tc.and.eps*sqrt(a*a*a).gt.0.5d3*h0)then ! tight coupling approximation
    !    dT_gas1 = 3.0+1.5*sqrt(omega_m/(omega_m+omega_l*a*a*a))
    !    dT_gas2 = dx_ion(la,x_ion,T_gas,delta_baryon)/(x_ion*(1.0+x_ion))
    !    dT_gas = -(dT_gas1+dT_gas2)/eps*T_photon*hubble(a)-T_photon+dT_gas3/hubble(a)
    ! else
    !    dT_gas1 = -2.0*T_gas
    !    dT_gas2 = (T_photon-T_gas)*eps
    !    dT_gas = dT_gas1+(dT_gas2+dT_gas3)/hubble(a)
    ! end if
    ! print*,"dT_gas1",dT_gas1
    ! print*,"dT_gas2",dT_gas2
    ! print*,"dT_gas3",dT_gas3
    ! print*,"dT_gas4",dT_gas4
    ! print*,"sum 2-4",dT_gas2+dT_gas3+dT_gas4
    ! print*,"1/hubble",1.0/hubble(a)
    ! dT_gas3=0.0
    dT_gas = dT_gas1+(dT_gas2+dT_gas3+dT_gas4)/hubble(a)
    ! print*,"minoda (1-1-5)"
  end function dT_gas

  double precision function dx_ion(la,x_ion,T_gas,delta_baryon)
    implicit none
    double precision::delta_baryon
    double precision::la, a_inv, x_ion, T_gas, T_photon!, rho_photon
    double precision::dx_ion1, dx_ion2, dx_ion3
    double precision::K, U, alpha_e, beta_e, gamma_e, Con
    double precision,parameter::L21=8.22458 ! Goldman 1989
    ! if(delta_baryon<-0.9)then
    !    delta_baryon=-0.9
    ! end if
    ! if(x_ion.ge.x_start)x_ion=x_start ! comment at 10th Feb. 2018
    a=exp(la)
    a_inv=1.0/a
    T_photon=T_rec*a_rec*a_inv
    nH=rho_baryon(a)/mc2_p*(1.0+delta_baryon)
    K = 1.21567d-5*1.21567d-5*1.21567d-5/(8.0*pi*hubble(a))
    U = E_1s/(kB*T_gas)
    ! ! minoda comment:
    ! ! Now check the dependence on the case B recombination rate.
    ! ! Peebles 1993 (used in TS 2006)
    ! alpha_e = 2.84d-13/sqrt(T_gas*1.0d-4)
    ! ! ! Peebles 1993 (used in SS 2004)
    ! ! alpha_e = 2.6d-13*((T_gas*1.0d-4)**(-0.8))
    ! ! Hummer 1994 (used in RECFAST)
    ! alpha_e = 1.14d-13*(4.309*(T_gas*1.0d-4)**(-0.6166))/(1.0+0.6703*(T_gas*1.0d-4)**0.5300)
    ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    alpha_e = 1.14d-13*(4.309*(T_gas*1.0d-4)**(-0.6166))/(1.0+0.6703*(T_gas*1.0d-4)**0.5300)
    beta_e = alpha_e * (2.0*pi*mc2_e*kB*T_photon*h2_inv*c2_inv)**1.5 * exp(-E_2s/(kB*T_photon))
    gamma_e = 0.291d-7 * U**0.39 / (0.232+U) * exp(-U) ! Voronov 1997
    ! if(x_ion.gt.0.985)then
    !    Con = 1.0
    ! else
    !    Con = (1.0+K*L21*nH*(1.0-x_ion)) / (1.0+K*(L21+beta_e)*nH*(1.0-x_ion))
    ! end if
    Con = (1.0+K*L21*nH*(1.0-x_ion)) / (1.0+K*(L21+beta_e)*nH*(1.0-x_ion))
    ! if(T_gas.lt.1.0d5)then
    !    dx_ion1 = beta_e*(1.0-x_ion)*exp(-0.75*E_1s/(kB*T_gas)) ! modified at 22nd March 2017
    ! else
    !    dx_ion1 = 0.0
    ! end if

    ! minoda commentted on 2018.11.28 !!!!!!!!!!!!!!!!
    ! Equation in the original RECFAST code
    ! dx_ion1 = beta_e*(1.0-x_ion)*exp(-0.75*E_1s/(kB*T_gas))
    ! J. Chluba's modification, written 2017.03.22
    dx_ion1 = beta_e*(1.0-x_ion)*exp(-0.75*E_1s/(kB*T_photon))
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dx_ion2 = alpha_e*nH*x_ion*x_ion
    ! dx_ion3 = gamma_e*nH*(1.0-x_ion)*x_ion
    dx_ion3 = 0.0

    ! for check
!    if(myrank.le.0.and.mod(step,100).eq.0)then
! !        print*,step
! !        print*,"x=",x_ion
! !        print*,"T=",T_gas
! ! !       print*,"beta_e",beta_e
!        print*,"dx1=",dx_ion1*Con
!        print*,"dx2=",dx_ion2*Con
!        print*,"dx3=",dx_ion3
!        !        print*,"hubble=",hubble(a)
! !        print*,"dx",((dx_ion1-dx_ion2)*Con + dx_ion3)/hubble(a)
!    end if
    dx_ion = ((dx_ion1-dx_ion2)*Con + dx_ion3)/hubble(a)
  end function dx_ion


end program get_history
