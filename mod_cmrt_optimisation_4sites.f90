Module mod_cmrt
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: mass_h=1.007825d0,kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=4184.d0
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! system paramters
integer nquant
real*8,allocatable:: H_diab(:,:),E_exc(:),lambda_diab(:,:)
real*8,allocatable::lambda_exc(:,:,:,:),c_tr(:,:)
real*8,allocatable::RR(:,:)
complex*16,allocatable::gg(:,:,:,:),gdot(:,:,:,:),g2dot(:,:,:,:)
complex*16,allocatable::g_diab(:,:),gdot_diab(:,:),g2dot_diab(:,:)
complex*16,allocatable::sigma(:,:),sigma_diab(:,:)
real*8 p2collapse,p3collapse,H_2,H_3,H_4,Ham_12,Ham_13,Ham_23

real*8 gamma,temperature,lambda
integer nsteps
real*8 dt,tot_time,curr_time

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer nold
integer,allocatable :: seed(:)
real*8,allocatable:: work(:)
complex*16,allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i,size_seed,seed2(2)

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar

  open(10,file="cmrt.inp")
  read(10,*) nquant
  read(10,*) tot_time
  read(10,*) dt
  read(10,*) seed2
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  !---------------------------------------------------------- 

  allocate(H_diab(nquant,nquant),E_exc(nquant),lambda_diab(nquant,nquant))
  allocate(lambda_exc(nquant,nquant,nquant,nquant),c_tr(nquant,nquant))
  allocate(gg(nquant,nquant,nquant,nquant),gdot(nquant,nquant,nquant,nquant),g2dot(nquant,nquant,nquant,nquant))
  allocate(g_diab(nquant,nquant),gdot_diab(nquant,nquant),g2dot_diab(nquant,nquant))
  allocate(sigma(nquant,nquant),RR(nquant,nquant),sigma_diab(nquant,nquant))

  nsteps=nint(tot_time/dt)


  !-----------------------------------------------------------------  
  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)
  enddo
  call random_seed(put=seed)
  !-----------------------------------------------------------------  

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i,j,k
  real::max_value,max_index,acceptance_ratio
  real*8 r5,x5 
  real,dimension(1500)::P_4
  real,dimension(1500)::robustness
  real,dimension(1500)::Hamil_2,Hamil_3,Hamil_4,Hamil_12,Hamil_13,Hamil_23
  real::ita=0.0005d+12
  real::ita1=0.0006d+12
  open(1,file="rob_15.out")
  open(102,file="max_rob_15.out")
  open(103,file="final_pop_15.out")
  open(104,file="final_coup_15.out")
  open(105,file="final_popmatrix_15.out")
  open(106,file="final_coupmatrix_15.out")
  print*, nsteps
!Implementing Metropolis
  do j=1,1500
     call setup_parameters  
     call init_cond
     open(100,file="sigma_15.out")
     do i=1,nsteps
        call write_sigma
        call evolve
        P_4(i)=real(sigma_diab(4,4)) 
     enddo 
!Calculating robustness value
     robustness(j)=P_4(nsteps)*((exp(-(p2collapse*ita)))*(exp(-(p3collapse*ita1))))
     acceptance_ratio=robustness(j)/robustness(j-1)
     call gaussian_random_number(r5)
     x5=r5
     if(acceptance_ratio>x5) then
       write(1,*) robustness(j)
     end if
     close(100)
     Hamil_2(j)=H_2/wave_to_J
     Hamil_3(j)=H_3/wave_to_J
     Hamil_4(j)=H_4/wave_to_J
     Hamil_12(j)=Ham_12/wave_to_J
     Hamil_13(j)=Ham_13/wave_to_J
     Hamil_23(j)=Ham_23/wave_to_J  
  end do
  max_value=robustness(1)
  max_index=1
  do k=1,1500
     if(robustness(k)>max_value)then
       max_value=robustness(k)
       max_index=k
     end if
  end do
  write(105,*) Hamil_2(max_index),Hamil_3(max_index),Hamil_4(max_index)
  write(106,*) Hamil_12(max_index),Hamil_13(max_index),Hamil_23(max_index)
  write(102,*) max_value,max_index,P_4(nsteps)
  close(1)
  close(102)
  close(103)
  close(104)
  close(105)
  close(106) 
end subroutine main

subroutine setup_parameters
  implicit none
  integer i
  real*8 en(nquant),vect(nquant,nquant)
  real*8 r1,x1,r2,x2,r3,x3,r4,x4,r5,x5,r6,x6,r7,x7,r8,x8
  real*8 H_12,H_13,H_23
  gamma=1/50.d-15
  temperature=77.d0
  lambda=35.d0*wave_to_J/pi
  

  H_diab(1,1)= 0.d0*wave_to_J
  H_diab(2,2)= 120.d0*wave_to_J
  H_diab(3,3)= -310.d0*wave_to_J
  H_diab(4,4)= -510.d0*wave_to_J
  H_diab(1,2)= 30.9d0*wave_to_J ;  H_diab(2,1)= H_diab(1,2)
  H_diab(1,3)= 0.d0*wave_to_J    ;  H_diab(3,1)= H_diab(1,3)
  H_diab(1,4)= 0.d0*wave_to_J    ;  H_diab(4,1)= H_diab(1,4)
  H_diab(2,3)= -87.7d0*wave_to_J  ;  H_diab(3,2)= H_diab(2,3)
  H_diab(2,4)= 0.d0*wave_to_J    ;  H_diab(4,2)= H_diab(2,4)
  H_diab(3,4)= 2.927d0*wave_to_J ;  H_diab(4,3)= H_diab(3,4)

  lambda_diab=0.d0
  do i=1,nquant
    lambda_diab(i,i)=lambda*pi
  enddo
!Changing parameter values
  call gaussian_random_number(r1)
  x1=45*r1
  H_diab(2,2)=H_diab(2,2)+x1*wave_to_J
  H_2=H_diab(2,2)
  call gaussian_random_number(r5)
  x5=70*r5
  H_diab(3,3)=H_diab(3,3)+x5*wave_to_J
  H_diab(4,4)=H_diab(3,3)-200*wave_to_J
  H_3=H_diab(3,3)
  H_4=H_diab(4,4)
  call gaussian_random_number(r2)
  x2=25*r2
  H_12=H_diab(1,2)+x2*wave_to_J*5
  if(H_12/wave_to_J<100 .and. H_12/wave_to_J>-100) then
    H_diab(1,2)=H_12
    H_diab(2,1)=H_12
  end if
  Ham_12=H_diab(1,2)
  call gaussian_random_number(r3)
  x3=25*r3
  H_23=H_diab(2,3)+x3*wave_to_J*5
  if(H_23/wave_to_J<100 .and. H_23/wave_to_J>-100) then
    H_diab(2,3)=H_23
    H_diab(3,2)=H_23
  end if
  Ham_23=H_diab(2,3)
  call gaussian_random_number(r4)
  x4=25*r4
  H_13=H_diab(1,3)+x4*wave_to_J*5
  if(H_13/wave_to_J<20 .and. H_13/wave_to_J>-20) then
    H_diab(1,3)=H_13
    H_diab(3,1)=H_13
  end if
  Ham_13=H_diab(1,3)  

  write(103,*) H_diab(1,1)/wave_to_J,H_diab(2,2)/wave_to_J,H_diab(3,3)/wave_to_J,H_diab(4,4)/wave_to_J
  write(104,*) H_diab(1,2)/wave_to_J,H_diab(2,3)/wave_to_J,H_diab(1,3)/wave_to_J

  call diag(H_diab,nquant,en,vect,nquant)
  c_tr = vect
  E_exc = en

end subroutine setup_parameters
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none

  sigma_diab=0.d0
  sigma_diab(1,1)=1.d0

  sigma=matmul(transpose(c_tr),matmul(sigma_diab,c_tr))

  curr_time=0.d0
  RR=0.d0

end subroutine init_cond
!-----------------------------------------------------------------  

subroutine evolve
  implicit none
  complex*16 sigma_dot(nquant,nquant)
  
  call calcualte_gg
  call calculate_sigmadot(sigma_dot)

!Calculating the probability of collapse of site 2 and site 3 
  p2collapse=p2collapse+(real(sigma_diab(2,2))*dt)
  p3collapse=p3collapse+(real(sigma_diab(3,3))*dt)

  sigma=sigma+sigma_dot*dt
  curr_time=curr_time+dt

end subroutine evolve
!-----------------------------------------------------------------  


subroutine calcualte_gg
  implicit none
  integer i,j,k,l,n,m

  call calculate_g_diab

  do i=1,nquant
    do j=1,nquant
      do k=1,nquant
        do l=1,nquant
          gg(i,j,k,l)=0.d0
          gdot(i,j,k,l)=0.d0
          g2dot(i,j,k,l)=0.d0
          lambda_Exc(i,j,k,l)=0.d0
          do n=1,nquant
            do m=1,nquant
              gg(i,j,k,l)=gg(i,j,k,l)+c_tr(n,i)*c_tr(n,j)*c_tr(m,k)*c_tr(m,l)*g_diab(n,m)
              lambda_exc(i,j,k,l)=lambda_exc(i,j,k,l)+c_tr(n,i)*c_tr(n,j)*c_tr(m,k)*c_tr(m,l)*lambda_diab(n,m)
              gdot(i,j,k,l)=gdot(i,j,k,l)+c_tr(n,i)*c_tr(n,j)*c_tr(m,k)*c_tr(m,l)*gdot_diab(n,m)
              g2dot(i,j,k,l)=g2dot(i,j,k,l)+c_tr(n,i)*c_tr(n,j)*c_tr(m,k)*c_tr(m,l)*g2dot_diab(n,m)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine calcualte_gg
!-----------------------------------------------------------------  

subroutine calculate_g_diab
  implicit none
  integer i,j,k
  integer,parameter::nw=5000
  real*8 dw,w_max
  real*8 omg(nw),wt(nw),spec(nquant,nquant,nw)

  w_max=25*gamma
  dw=w_max/real(nw)
  do i=1,nw
    omg(i)=i*dw
    do j=1,nquant
      do k=1,nquant
        spec(j,k,i)=spectral(j,k,omg(i))
!write(20,*) omg(i)/(1*pi*clight),spec(1,1,i)/wave_to_J
      enddo
    enddo
  enddo
!  stop
  wt=omg*curr_time

  do i=1,nquant
    do j=1,nquant
      g_diab(i,j)=dw*sum(spec(i,j,:)/omg**2*(1/tanh(hbar*omg/(2*kb*temperature))*(1-cos(wt))+iota*(sin(wt)-wt)))
      g_diab(i,j)=g_diab(i,j)/(hbar)
      gdot_diab(i,j)=dw*sum(spec(i,j,:)/omg**2*(1/tanh(hbar*omg/(2*kb*temperature))*(sin(wt)*omg)+iota*(cos(wt)*omg-omg)))
      gdot_diab(i,j)=gdot_diab(i,j)/(hbar)
      g2dot_diab(i,j)=dw*sum(spec(i,j,:)/omg**2*(1/tanh(hbar*omg/(2*kb*temperature))*(cos(wt)*omg*omg)-iota*(sin(wt)*omg*omg)))
      g2dot_diab(i,j)=g2dot_diab(i,j)/(hbar)
    enddo
  enddo



end subroutine calculate_g_diab
!-----------------------------------------------------------------  

real*8 function spectral(i,j,w)
  implicit none
  integer,intent(in)::i,j
  real*8,intent(in)::w

  if(i==j) then
    spectral=2*lambda * gamma*w/(w**2+gamma**2)
  else
    !spectral=2*lambda * gamma*w/(w**2+gamma**2)
    spectral=0.d0
  endif


end function spectral
!-----------------------------------------------------------------  

subroutine calculate_sigmadot(sigma_dot)
  implicit none
  complex*16,intent(out):: sigma_dot(nquant,nquant)
  complex*16 RR_pd(nquant,nquant)
  integer i,j,f

  call calcualte_RR
  call calcualte_RR_pd(RR_pd)

  do i=1,nquant
    do j=1,nquant
      sigma_dot(i,j)=-iota*(E_exc(i)-E_exc(j))/hbar * sigma(i,j)
      if(i==j) then
        do f=1,nquant
          sigma_dot(i,j)=sigma_dot(i,j)+RR(i,f)*sigma(f,f)-RR(f,i)*sigma(i,i)
        enddo
      endif
      sigma_dot(i,j)=sigma_dot(i,j)-RR_pd(i,j)*sigma(i,j)
      if(i.ne.j) then
        do f=1,nquant
          sigma_dot(i,j)=sigma_dot(i,j)-0.5*(RR(f,i)+RR(f,j))*sigma(i,j)
        enddo
      endif
    enddo
  enddo


end subroutine calculate_sigmadot
!-----------------------------------------------------------------  

subroutine calcualte_RR
  implicit none
  complex*16 AA(nquant),FF(nquant),XX(nquant,nquant),fac
  integer i,j
  real*8 tt

  tt=curr_time

  do i=1,nquant
    AA(i) = exp(-iota*E_exc(i)*tt/hbar-gg(i,i,i,i))
    FF(i) = exp(-iota*(E_exc(i)-2*lambda_exc(i,i,i,i))*tt/hbar-conjg(gg(i,i,i,i)))
  enddo

  do i=1,nquant
    do j=1,nquant
      XX(i,j)=exp(2*(gg(i,i,j,j)+iota*lambda_exc(i,i,j,j)*tt/hbar))
      fac=(gdot(j,i,i,i)-gdot(j,i,j,j)-2*iota*lambda_exc(j,i,j,j)/hbar)
      fac=fac*(gdot(i,j,i,i)-gdot(i,j,j,j)-2*iota*lambda_exc(i,j,j,j)/hbar)
      XX(i,j)=XX(i,j)*(g2dot(j,i,i,j)-fac)
    enddo
  enddo

  !RR=0.d0

  do i=1,nquant
    do j=1,nquant
      if(i.ne.j)RR(i,j)=RR(i,j)+2*dt*real(conjg(FF(j))*AA(i)*XX(i,j))
    enddo
  enddo

end subroutine calcualte_RR
!-----------------------------------------------------------------  

subroutine calcualte_RR_pd(RR_pd)
  implicit none
  complex*16,intent(out)::RR_pd(nquant,nquant)
  integer i,j

  RR_pd=0.d0
  do i=1,nquant
    do j=1,nquant
      if(i.ne.j) then
        RR_pd(i,j) = real(gdot(i,i,i,i)+gdot(j,j,j,j)-2*gdot(i,i,j,j))
        RR_pd(i,j)=RR_pd(i,j)+iota*aimag(gdot(i,i,i,i)-gdot(j,j,j,j))
      endif
    enddo
  enddo

end subroutine calcualte_RR_pd
!-----------------------------------------------------------------  

subroutine write_sigma
  implicit none

  sigma_diab=matmul(c_tr,matmul(sigma,transpose(c_tr)))
  write(100,'(5es)') curr_time*1.d15,real(sigma_diab(1,1)),real(sigma_diab(2,2)),real(sigma_diab(3,3)),real(sigma_diab(4,4))
end subroutine write_sigma
!----------------------------------------------------------------

!----------------------------------------------------------------



 
!----------------------------------------------------------------- 
subroutine compute_KE_matrix_dvr(KE_matrix,ngrid,delq,mass)
  !! computes KE matrix in DVR basis
  !! Appendix A of JCP 96, 1982 (1991)
  implicit none
  integer,intent(in) :: ngrid       !! size of DVR grid
  real*8,intent(inout) :: KE_matrix(ngrid,ngrid)
  real*8,intent(in) :: delq         !! step size in position of DVR basis
  real*8,intent(in) :: mass         !! mass
  integer i,j
  real*8 pi,hbar

  pi=dacos(-1.d0);hbar=1.05457266D-34

  KE_matrix=hbar**2/(2*mass*delq**2)
  do i=1,ngrid
    do j=1,ngrid
      KE_matrix(i,j)=KE_matrix(i,j)*(-1.d0)**(i-j)
      if(i==j) then
        KE_matrix(i,j)=KE_matrix(i,j)*pi**2/3.d0
      else
        KE_matrix(i,j)=KE_matrix(i,j)*2.d0/real(i-j)**2
      endif
    enddo
  enddo
end subroutine compute_KE_matrix_dvr
!---------------------------------------------------------- 

subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!---------------------------------------------------------- 
subroutine gaussian_random_number(rnd)
  !! generates gaussian distribution with center 0, sigma 1
  !! q0+sig*rnd gives center=q0, sigma=sig
  implicit none
  real*8,intent(out)::rnd
  real*8 rnd1,rnd2,pi

  pi=dacos(-1.d0)

  call random_number(rnd1)
  call random_number(rnd2)
  rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!----------------------------------------------------------
End Module mod_cmrt
