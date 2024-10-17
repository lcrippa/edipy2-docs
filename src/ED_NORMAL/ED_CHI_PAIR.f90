MODULE ED_CHI_PAIR
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,reg,txtfy
  USE SF_LINALG,  only: inv,eigh
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_NORMAL
  USE ED_AUX_FUNX

  implicit none
  private


  public :: build_chi_pair_normal

  integer                          :: istate,iorb,jorb,ispin,jspin
  integer                          :: isector,jsector,ksector
  real(8),allocatable              :: vvinit(:),vvinit_tmp(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: ialfa
  integer                          :: jalfa
  integer                          :: ipos,jpos
  integer                          :: i,j,k
  real(8)                          :: sgn,norm2
  real(8),dimension(:),allocatable :: state_dvec
  real(8)                          :: state_e

contains


  !+------------------------------------------------------------------+
  !                            PAIR
  !PURPOSE  : Evaluate the pair susceptibility \Chi_pair for a 
  ! \chi_ab = <Delta*_a(\tau)Delta_b(0)>
  !+------------------------------------------------------------------+
  subroutine build_chi_pair_normal()
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG build_Chi_spin_normal: build pair-Chi"
#endif
    write(LOGfile,"(A)")"Get impurity pair Chi:"
    do iorb=1,Norb
       write(LOGfile,"(A)")"Get Chi_pair_l"//reg(txtfy(iorb))
       if(MPIMASTER)call start_timer(unit=LOGfile)
       call lanc_ed_build_pairChi_diag(iorb)
       if(MPIMASTER)call stop_timer
#ifdef _DEBUG
       if(ed_verbose>1)write(Logfile,"(A)")""
#endif
    enddo

    if(Norb>1)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             write(LOGfile,"(A)")"Get Chi_pair_mix_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
             if(MPIMASTER)call start_timer(unit=LOGfile)
             call lanc_ed_build_pairChi_mix(iorb,jorb)
             if(MPIMASTER)call stop_timer
#ifdef _DEBUG
             if(ed_verbose>1)write(Logfile,"(A)")""
#endif
          end do
       end do
       !
       !
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             select case(ed_diag_type)
             case default
                pairChi_w(iorb,jorb,:)   = 0.5d0*(pairChi_w(iorb,jorb,:) - pairChi_w(iorb,iorb,:) - pairChi_w(jorb,jorb,:))
                pairChi_tau(iorb,jorb,:) = 0.5d0*(pairChi_tau(iorb,jorb,:) - pairChi_tau(iorb,iorb,:) - pairChi_tau(jorb,jorb,:))
                pairChi_iv(iorb,jorb,:)  = 0.5d0*(pairChi_iv(iorb,jorb,:) - pairChi_iv(iorb,iorb,:) - pairChi_iv(jorb,jorb,:))
                !
             case ("full")
                !
             end select
             !
             pairChi_w(jorb,iorb,:)   = pairChi_w(iorb,jorb,:)
             pairChi_tau(jorb,iorb,:) = pairChi_tau(iorb,jorb,:)
             pairChi_iv(jorb,iorb,:)  = pairChi_iv(iorb,jorb,:)
          enddo
       enddo
    endif
    !

  end subroutine build_chi_pair_normal




  ! \chi_aa = <Delta*_a(\tau)Delta_a(0)>
  !         = <[C^+_a(\tau)C^+_a(\tau)][C_a(0)C_a(0)]>
  subroutine lanc_ed_build_pairChi_diag(iorb)
    integer                     :: iorb
    type(sector)                :: sectorI,sectorJ,sectorK
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG lanc_ed_build_pairChi diag: Lanczos build pair Chi l"//str(iorb)
#endif
    !
    if(ed_total_ud)then
       ialfa = 1
       ipos  = iorb
    else
       ialfa = iorb
       ipos  = 1
    endif
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          call es_return_dvector(MpiComm,state_list,istate,state_dvec)
       else
          call es_return_dvector(state_list,istate,state_dvec)
       endif
#else
       call es_return_dvector(state_list,istate,state_dvec)
#endif
       !
       ksector = getCsector(ialfa,2,isector)
       jsector = getCsector(ialfa,1,ksector)
       if(jsector/=0.AND.ksector/=0)then
          !
          if(MpiMaster)then
             call build_sector(isector,sectorI)
             call build_sector(ksector,sectorK)
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A30,I6,20I4)")&
                  'Apply Cup*Cdw',isector,sectorI%Nups,sectorI%Ndws

             allocate(vvinit_tmp(sectorK%Dim)) ;  vvinit_tmp=0d0
             allocate(vvinit(sectorJ%Dim))     ;  vvinit=0d0
             !C_dw|gs>=|tmp>
             do i=1,sectorI%Dim
                call apply_op_C(i,k,sgn,ipos,ialfa,2,sectorI,sectorK)
                if(sgn==0d0.OR.k==0)cycle
                vvinit_tmp(k) = sgn*state_dvec(i)
             enddo
             !C_up|tmp>=C_up[C_dw|gs>]
             do k=1,sectorK%Dim
                call apply_op_C(k,j,sgn,ipos,ialfa,1,sectorK,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*vvinit_tmp(k)
             enddo
             deallocate(vvinit_tmp)
             call delete_sector(sectorI)
             call delete_sector(sectorK)
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_pairChi(norm2,state_e,alfa_,beta_,iorb,iorb)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       endif
       if(allocated(state_dvec))deallocate(state_dvec)
       !
    enddo
    return
  end subroutine lanc_ed_build_pairChi_diag







  ! \chi_ab = <Delta*_a(\tau)Delta_b(0)>
  !         = <[C^+_a(\tau)C^+_a(\tau)][C_b(0)C_b(0)]>
  !from aux: <[C^+_a C^+_a + C^+_b C^+_b][C_a C_a + C_b C_b]>  
  subroutine lanc_ed_build_pairChi_mix(iorb,jorb)
    integer                     :: iorb,jorb
    type(sector)                :: sectorI,sectorJ,sectorK
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG lanc_ed_build_pairChi mix: Lanczos build pair Chi l"//str(iorb)//",m"//str(jorb)
#endif
    !    
    if(ed_total_ud)then
       ialfa = 1
       jalfa = 1
       ipos  = iorb
       jpos  = jorb
    else
       write(LOGfile,"(A)")"ED_CHI_PAIR warning: can not evaluate \Chi_pair_ab with ed_total_ud=F"
       return
    endif
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          call es_return_dvector(MpiComm,state_list,istate,state_dvec)
       else
          call es_return_dvector(state_list,istate,state_dvec)
       endif
#else
       call es_return_dvector(state_list,istate,state_dvec)
#endif
       !
       !
       ! --> Apply [C_b C_b + C_a C_a]|state>
       ksector = getCsector(ialfa,2,isector)
       jsector = getCsector(jalfa,1,ksector)
       if(jsector/=0.AND.ksector/=0)then
          if(MpiMaster)then
             call build_sector(isector,sectorI)
             call build_sector(ksector,sectorK)
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A30,I6,20I4)")&
                  'Apply C_bu*C_bd+C_au*C_ad',isector,sectorI%Nups,sectorI%Ndws
             allocate(vvinit_tmp(sectorK%Dim)) ;  vvinit_tmp=0d0
             allocate(vvinit(sectorJ%Dim))     ;  vvinit=0d0
             !
             !Apply C_a,up*C_a,dw:
             do i=1,sectorI%Dim
                call apply_op_C(i,k,sgn,ipos,ialfa,2,sectorI,sectorK)
                if(sgn==0d0.OR.k==0)cycle
                vvinit_tmp(k) = sgn*state_dvec(i)
             enddo
             do k=1,sectorK%Dim
                call apply_op_C(k,j,sgn,ipos,ialfa,1,sectorK,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*vvinit_tmp(k)
             enddo
             !
             !+ C_b,up*C_b,dw
             vvinit_tmp=0d0
             do i=1,sectorI%Dim
                call apply_op_C(i,k,sgn,jpos,jalfa,2,sectorI,sectorK)
                if(sgn==0d0.OR.k==0)cycle
                vvinit_tmp(k) = sgn*state_dvec(i)
             enddo
             do k=1,sectorK%Dim
                call apply_op_C(k,j,sgn,jpos,jalfa,1,sectorK,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*vvinit_tmp(k)
             enddo
             deallocate(vvinit_tmp)
             call delete_sector(sectorI)
             call delete_sector(sectorK)
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_pairChi(norm2,state_e,alfa_,beta_,iorb,jorb)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       endif
       if(allocated(state_dvec))deallocate(state_dvec)
       !
    enddo
    return
  end subroutine lanc_ed_build_pairChi_mix







  subroutine add_to_lanczos_pairChi(vnorm2,Ei,alanc,blanc,iorb,jorb)
    integer                                    :: iorb,jorb
    real(8)                                    :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
    real(8)                                    :: Ei,Ej,Egs,de
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_pairChi: add-up to GF"
#endif
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm2/zeta_function 
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_pairChi: LApack tridiagonalization"
#endif
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       ! the correct behavior for beta*dE << 1 is recovered only by assuming that v_n is still finite
       ! beta*dE << v_n for v_n--> 0 slower. First limit beta*dE--> 0 and only then v_n -->0.
       ! This ensures that the correct null contribution is obtained.
       ! So we impose that: if (beta*dE is larger than a small qty) we sum up the contribution, else
       ! we do not include the contribution (because we are in the situation described above).
       ! For the real-axis case this problem is circumvented by the usual i*0+ = xi*eps
       if(beta*dE > 1d-3)pairChi_iv(iorb,jorb,0)=pairChi_iv(iorb,jorb,0) + peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          pairChi_iv(iorb,jorb,i)=pairChi_iv(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=0,Ltau
          pairChi_tau(iorb,jorb,i)=pairChi_tau(iorb,jorb,i) + exp(-tau(i)*dE)*peso
       enddo
       do i=1,Lreal
          pairChi_w(iorb,jorb,i)=pairChi_w(iorb,jorb,i) - &
               peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
    !
  end subroutine add_to_lanczos_pairChi



END MODULE ED_CHI_PAIR
























