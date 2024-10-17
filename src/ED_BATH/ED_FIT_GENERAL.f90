MODULE ED_FIT_GENERAL
  USE ED_FIT_COMMON
  USE SF_SPIN, only: pauli_sigma_z
  USE SF_LINALG, only: kron, diag, trace

  implicit none
  private


  public :: chi2_fitgf_general
  public :: chi2_fitgf_general_superc

  complex(8),dimension(:,:,:,:,:),allocatable       :: FGmatrix


contains


  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for GENERAL BATH
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_general(fg,bath_)
    complex(8),dimension(:,:,:,:,:)             :: fg ![Nspin][Nspin][Norb][Norb][Lmats]
    logical,dimension(Nspin,Nspin,Norb,Norb)    :: Hmask
    real(8),dimension(:),intent(inout)          :: bath_
    real(8),dimension(:),allocatable            :: array_bath
    integer                                     :: i,j,iorb,jorb,ispin,jspin,io,jo,ibath
    integer                                     :: iter,stride,counter,Asize
    real(8)                                     :: chi
    logical                                     :: check
    type(effective_bath)                        :: dmft_bath
    character(len=256)                          :: suffix
    integer                                     :: unit
    complex(8),dimension(:,:,:,:,:),allocatable :: fgand ![Nspin][][Norb][][Ldelta]  
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG chi2_fitgf_general: Fit"
#endif
    !
    if(size(fg,1)/=Nspin)stop "chi2_fitgf_general error: size[fg,1]!=Nspin"
    if(size(fg,2)/=Nspin)stop "chi2_fitgf_general error: size[fg,2]!=Nspin"
    if(size(fg,3)/=Norb)stop "chi2_fitgf_general error: size[fg,3]!=Norb"
    if(size(fg,4)/=Norb)stop "chi2_fitgf_general error: size[fg,4]!=Norb"
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_general error: wrong bath dimensions"
    !
    if(cg_pow/=2.AND.cg_norm=="frobenius")then
       print *, "WARNING: CG_POW must be 2 for a meaningful definition of the Frobenius norm."
       print *, "         we'll still let you go ahead with the desired input, but please be "
       print *, "         be aware that CG_POW is not doing what you would expect for a chi^q"
    endif
    !
    call allocate_dmft_bath(dmft_bath)
    call set_dmft_bath(bath_,dmft_bath)
    allocate(array_bath(size(bath_)-1))
    Nlambdas  =bath_(1)
    array_bath=bath_(2:)
    !
    if(ed_all_g)then
      Hmask=.true.
    else
      Hmask=Hgeneral_mask(wdiag=.true.,uplo=.false.)
    endif 
    !
    totNso=count(Hmask)
    !
    allocate(getIspin(totNso),getJspin(totNso))
    allocate(getIorb(totNso) ,getJorb(totNso))
    counter=0
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                if (Hmask(ispin,jspin,iorb,jorb))then
                   counter=counter+1
                   getIspin(counter) = ispin
                   getIorb(counter)  = iorb
                   getJspin(counter) = jspin
                   getJorb(counter)  = jorb
                endif
             enddo
          enddo
       enddo
    enddo
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,5))Ldelta=size(fg,5)
    !
    allocate(FGmatrix(Nspin,Nspin,Norb,Norb,Ldelta))
    allocate(Gdelta(totNso,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
    !
    select case(cg_weight)
    case default
       Wdelta=1d0
    case(2)
       Wdelta=1d0*arange(1,Ldelta)
    case(3)
       Wdelta=Xdelta
    end select
    !
    write(LOGfile,*)"Fitted functions",totNso
    do i=1,totNso
       Gdelta(i,1:Ldelta) = fg(getIspin(i),getJspin(i),getIorb(i),getJorb(i),1:Ldelta)
    enddo
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG chi2_fitgf_general: cg_method:"//str(cg_method)//&
         ", cg_grad:"//str(cg_grad)//&
         ", cg_scheme:"//str(cg_scheme)//&
         ", cg_norm:"//str(cg_norm)
#endif
    !
    FGmatrix=fg
    !
    select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
    case default
       if(cg_grad==0)then
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,chi2_weiss_general,grad_chi2_weiss_general,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,chi2_delta_general,grad_chi2_delta_general,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_general error: cg_scheme != [weiss,delta]"
          end select
       else
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,chi2_weiss_general,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,chi2_delta_general,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_general error: cg_scheme != [weiss,delta]"
          end select
       endif
       !
       !
    case (1)
       if(cg_grad==0)then
          write(*,*) "                                                                                "
          write(*,*) "WARNING: analytic gradient not available with cg-method=1 (minimize f77 routine)"
          write(*,*) "         > we will force cg_grad=1 (so let the routine estimate the gradient)   "
          write(*,*) "                                                                                "
       endif
       select case (cg_scheme)
       case ("weiss")
          call fmin_cgminimize(array_bath,chi2_weiss_general,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,chi2_delta_general,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_general error: cg_scheme != [weiss,delta]"
       end select
       !
    end select
    !
    write(LOGfile,"(A,ES18.9,A,I5,A)")"chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,"  <--  All Orbs, All Spins"
    !
    suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
    unit=free_unit()
    open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
    write(unit,"(ES18.9,1x,I5)") chi,iter
    close(unit)
    !
    bath_(2:size(bath_))=array_bath
    call set_dmft_bath(bath_,dmft_bath) ! *** array_bath --> dmft_bath *** (per write fit result)
    call write_dmft_bath(dmft_bath,LOGfile)
    call save_dmft_bath(dmft_bath)
    !
    allocate(fgand(Nspin,Nspin,Norb,Norb,Ldelta))
    if(cg_scheme=='weiss')then
       fgand = g0and_bath_function(xi*Xdelta(:),dmft_bath)
    else
       fgand = delta_bath_function(xi*Xdelta(:),dmft_bath)
    endif
    call write_fit_result()
    deallocate(fgand)
    !
    call get_dmft_bath(dmft_bath,bath_)                ! ***  dmft_bath --> bath_ ***    (bath in output)
    call deallocate_dmft_bath(dmft_bath)
    deallocate(FGmatrix,Gdelta,Xdelta,Wdelta)
    deallocate(getIspin,getJspin)
    deallocate(getIorb,getJorb)
    deallocate(array_bath)
    !
  contains
    !
    subroutine write_fit_result()
      integer   :: i,j,s,l,iorb,jorb,ispin,jspin
      !
      do l=1,totNso
         iorb = getIorb(l)
         jorb = getJorb(l)
         ispin = getIspin(l)
         jspin = getJspin(l)
         suffix="_l"//reg(txtfy(iorb))//&
              "_m"//reg(txtfy(jorb))//&
              "_s"//reg(txtfy(ispin))//&
              "_r"//reg(txtfy(jspin))//reg(ed_file_suffix)
         unit=free_unit()
         if(cg_scheme=='weiss')then
            open(unit,file="fit_weiss"//reg(suffix)//".ed")
         else
            open(unit,file="fit_delta"//reg(suffix)//".ed")
         endif
         do i=1,Ldelta
            write(unit,"(5F24.15)")Xdelta(i),&
                 dimag(fg(ispin,jspin,iorb,jorb,i)),dimag(fgand(ispin,jspin,iorb,jorb,i)),&
                 dreal(fg(ispin,jspin,iorb,jorb,i)),dreal(fgand(ispin,jspin,iorb,jorb,i))
         enddo
         close(unit)
      enddo
    end subroutine write_fit_result
    !
  end subroutine chi2_fitgf_general










  !##################################################################
  ! THESE PROCEDURES EVALUATE THE \chi^2 FUNCTIONS TO MINIMIZE. 
  !##################################################################
  !
  !
  !+-----------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
  !+-----------------------------------------------------------------+
   function chi2_delta_general(a) result(chi2)
    real(8),dimension(:)                                         :: a
    real(8)                                                      :: chi2
    !
    select case(cg_norm)
     case ("elemental")
       chi2 = chi2_delta_general_elemental(a)
     case ("frobenius")
       chi2 = chi2_delta_general_frobenius(a)
     case default
       stop "chi2_fitgf_general error: cg_norm != [elemental,frobenius]"
    end select
    !
 end function chi2_delta_general


  function chi2_delta_general_elemental(a) result(chi2)
    real(8),dimension(:)                               :: a
    real(8)                                            :: chi2
    real(8),dimension(totNso)                          :: chi2_so
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
    real(8),dimension(Ldelta)                          :: Ctmp
    integer                                            :: i,l,iorb,jorb,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_general. a:",a
#endif
    !
    Delta = delta_general(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       Ctmp =  abs(Gdelta(l,:)-Delta(ispin,jspin,iorb,jorb,:))
       chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_so)
    chi2=chi2/Ldelta/totNso
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_general. Chi**2:",chi2
#endif
    !
  end function chi2_delta_general_elemental
  !
  
 !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
 function chi2_delta_general_frobenius(a) result(chi2)
    real(8),dimension(:)                                         :: a
    real(8)                                                      :: chi2
    real(8),dimension(Ldelta)                                    :: chi2_freq
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)           :: Delta
    complex(8),dimension(Nspin*Norb,Nspin*Norb)                  :: Delta_so
    integer                                                      :: l
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_delta_general_frobenius. a:",a
#endif
    !
    Delta = delta_general(a)
    !
    do l=1,Ldelta
       Delta_so    =  nn2so_reshape(delta(:,:,:,:,l) - FGmatrix(:,:,:,:,l),Nspin,Norb)
       chi2_freq(l) =  sqrt(trace(matmul(Delta_so,conjg(transpose(Delta_so)))))
    enddo
    !
    chi2 = sum(chi2_freq**cg_pow/Wdelta) !Weighted sum over matsubara frqs
    chi2 = chi2/Ldelta/(Nspin*Norb) !Normalization over {iw} and Nlso
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_general_frobenius. chi2:",chi2
#endif
    !
 end function chi2_delta_general_frobenius
  !
  !+--------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function.
  !+--------------------------------------------------------------+
   function chi2_weiss_general(a) result(chi2)
      real(8),dimension(:)                                         :: a
      real(8)                                                      :: chi2
      !
      select case(cg_norm)
       case ("elemental")
         chi2 = chi2_weiss_general_elemental(a)
       case ("frobenius")
         chi2 = chi2_weiss_general_frobenius(a)
       case default
         stop "chi2_fitgf_general error: cg_norm != [elemental,frobenius]"
      end select
      !
   end function chi2_weiss_general  

  function chi2_weiss_general_elemental(a) result(chi2)
   real(8),dimension(:)                               :: a
   real(8)                                            :: chi2
   real(8),dimension(totNso)                          :: chi2_so
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
   real(8),dimension(Ldelta)                          :: Ctmp
   integer                                            :: i,l,iorb,jorb,ispin,jspin
   !
#ifdef _DEBUG
   if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_general. a:",a
#endif
   !
   g0and = g0and_general(a)
   !
   do l=1,totNso
      iorb = getIorb(l)
      jorb = getJorb(l)
      ispin = getIspin(l)
      jspin = getJspin(l)
      !
      Ctmp = abs(Gdelta(l,:)-g0and(ispin,jspin,iorb,jorb,:))
      chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
   enddo
   !
   chi2=sum(chi2_so)
   chi2=chi2/Ldelta/totNso
#ifdef _DEBUG
   if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_general. Chi**2:",chi2
#endif
  !
  end function chi2_weiss_general_elemental

   !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
   function chi2_weiss_general_frobenius(a) result(chi2)
      real(8),dimension(:)                                         :: a
      real(8)                                                      :: chi2
      real(8),dimension(Ldelta)                                    :: chi2_freq
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)           :: g0and
      complex(8),dimension(Nspin*Norb,Nspin*Norb)                  :: Delta_so
      integer                                                      :: l
      !
#ifdef _DEBUG
      if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_general_frobenius. a:",a
#endif
      !
      g0and = g0and_general(a)
      !
      do l=1,Ldelta
         Delta_so    =  nn2so_reshape(g0and(:,:,:,:,l) - FGmatrix(:,:,:,:,l),Nspin,Norb)
         chi2_freq(l) =  sqrt(trace(matmul(Delta_so,conjg(transpose(Delta_so)))))
      enddo
      !
      chi2 = sum(chi2_freq**cg_pow/Wdelta) !Weighted sum over matsubara frqs
      chi2 = chi2/Ldelta/(Nspin*Norb) !Normalization over {iw} and Nlso
      !
#ifdef _DEBUG
      if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_general_frobenius. chi2:",chi2
#endif
      !
   end function chi2_weiss_general_frobenius



  !######################################################################
  ! THESE PROCEDURES EVALUATE THE >GRADIENTS< OF THE \chi^2 TO MINIMIZE. 
  !######################################################################
   function grad_chi2_delta_general(a) result(dchi2)
      real(8),dimension(:)                                         :: a
      real(8),dimension(size(a))                                   :: dchi2
      !
      select case(cg_norm)
       case ("elemental")
         dchi2 = grad_chi2_delta_general_elemental(a)
       case ("frobenius")
         dchi2 = grad_chi2_delta_general_frobenius(a)
       case default
         stop "chi2_fitgf_general error: cg_norm != [elemental,frobenius]"
      end select
      !
   end function grad_chi2_delta_general
  !
  !+---------------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of \Delta_Anderson function.
  !+---------------------------------------------------------------------+
  function grad_chi2_delta_general_elemental(a) result(dchi2)
    real(8),dimension(:)                                       :: a
    real(8),dimension(size(a))                                 :: dchi2
    real(8),dimension(totNso,size(a))                          :: df
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)         :: Delta
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
    complex(8),dimension(Ldelta)                               :: Ftmp
    real(8),dimension(Ldelta)                                  :: Ctmp
    integer                                                    :: i,j,l,iorb,jorb,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_general. a:",a
#endif
    !
    Delta  = delta_general(a)
    dDelta = grad_delta_general(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       Ftmp = Gdelta(l,:)-Delta(ispin,jspin,iorb,jorb,:)
       Ctmp = abs(Ftmp)**(cg_pow-2)
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Ftmp)*dreal(dDelta(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
               sum( dimag(Ftmp)*dimag(dDelta(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta )
       enddo
    enddo
    !
    dchi2 = -cg_pow*sum(df,1)/Ldelta/totNso     !sum over all orbital indices
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_general. dChi**2:",dchi2
#endif
    !
  end function grad_chi2_delta_general_elemental
  !
   !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
   function grad_chi2_delta_general_frobenius(a) result(dchi2)
      real(8),dimension(:)                                                 :: a
      real(8),dimension(size(a))                                           :: dchi2
      real(8),dimension(Ldelta,size(a))                                    :: df
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)                   :: Delta
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a))           :: dDelta
      complex(8),dimension(Ldelta)                                         :: Ftmp
      real(8),dimension(Ldelta,size(a))                                    :: dChi_freq
      integer                                                              :: i,j,idelta,iorb,jorb,ispin,jspin
      !
      Delta  = delta_general(a)
      dDelta = grad_delta_general(a)
      Ftmp=zero
      df=zero
      !
      do idelta=1,Ldelta
           do ispin=1,Nspin
              do jspin=1,Nspin
                 do iorb=1,Norb
                    do jorb=1,Norb
                       !
                       Ftmp(idelta) = Ftmp(idelta) + abs(Delta(ispin,jspin,iorb,jorb,idelta)-FGmatrix(ispin,jspin,iorb,jorb,idelta))**2
                       do j=1,size(a)
                          df(idelta,j) = df(idelta,j) + &
                             real(Delta(ispin,jspin,iorb,jorb,idelta) - FGmatrix(ispin,jspin,iorb,jorb,idelta)) * &
                             real(dDelta(ispin,jspin,iorb,jorb,idelta,j)) + &
                             imag(Delta(ispin,jspin,iorb,jorb,idelta) - FGmatrix(ispin,jspin,iorb,jorb,idelta)) * &
                             imag(dDelta(ispin,jspin,iorb,jorb,idelta,j))
                       enddo
                    enddo
                 enddo
              enddo
         enddo
         Ftmp(idelta) = cg_pow * (sqrt(Ftmp(idelta))**(cg_pow-2)) / Wdelta(idelta)
         dchi_freq(idelta,:) = Ftmp(idelta) * df(idelta,:)
      enddo
      !
      dchi2 = sum(dchi_freq,1)/Ldelta/(Nspin*Norb)
      !
#ifdef _DEBUG
      if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_general_frobenius. dChi2:",dchi2
#endif
      !
   end function grad_chi2_delta_general_frobenius
  !
  !+------------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of G0_Anderson function.
  !+------------------------------------------------------------------+
   function grad_chi2_weiss_general(a) result(dchi2)
      real(8),dimension(:)                                         :: a
      real(8),dimension(size(a))                                   :: dchi2
      !
      select case(cg_norm)
       case ("elemental")
         dchi2 = grad_chi2_weiss_general_elemental(a)
       case ("frobenius")
         dchi2 = grad_chi2_weiss_general_frobenius(a)
       case default
         stop "chi2_fitgf_general error: cg_norm != [elemental,frobenius]"
      end select
      !
   end function grad_chi2_weiss_general  

  function grad_chi2_weiss_general_elemental(a) result(dchi2)
    real(8),dimension(:)                                       :: a
    real(8),dimension(size(a))                                 :: dchi2
    real(8),dimension(totNso,size(a))                          :: df
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)         :: g0and
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dg0and
    complex(8),dimension(Ldelta)                               :: Ftmp
    real(8),dimension(Ldelta)                                  :: Ctmp
    integer                                                    :: i,j,l,iorb,jorb,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_general. a:",a
#endif
    !
    g0and  = g0and_general(a)
    dg0and = grad_g0and_general(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       Ftmp = Gdelta(l,:)-g0and(ispin,jspin,iorb,jorb,:)
       Ctmp = abs(Ftmp)**(cg_pow-2)
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Ftmp)*dreal(dg0and(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
               sum( dimag(Ftmp)*dimag(dg0and(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta )
       enddo
    enddo
    !
    dchi2 = -cg_pow*sum(df,1)/Ldelta/totNso     !sum over all orbital indices
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_general. dChi**2:",dchi2
#endif
    !
  end function grad_chi2_weiss_general_elemental

   !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
   function grad_chi2_weiss_general_frobenius(a) result(dchi2)
      real(8),dimension(:)                                                 :: a
      real(8),dimension(size(a))                                           :: dchi2
      real(8),dimension(Ldelta,size(a))                                    :: df
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)                   :: g0and
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a))           :: dg0and
      complex(8),dimension(Ldelta)                                         :: Ftmp
      real(8),dimension(Ldelta,size(a))                                    :: dChi_freq
      integer                                                              :: i,j,idelta,iorb,jorb,ispin,jspin
      !
      !
      print*, "                                                                     "
      print*, "WARNING: the analytic gradient of the Weiss field Frobenius distance "
      print*, "         has been found giving /QUALITATIVELY WRONG/ fitted spectra! "
      print*, "         > the issue has emerged in some Nlat=Nspin=Norb=1 test runs "
      print*, "         > while the bug is investigated please switch cg_scheme to  "
      print*, "           'delta' or cg_norm to 'elemental' (or give up analytic cg)"
      print*, "                                                                     "
      !
      !
      g0and  = g0and_general(a)
      dg0and = grad_g0and_general(a)
      Ftmp=zero
      df=zero
      !
      do idelta=1,Ldelta
         do ispin=1,Nspin
            do jspin=1,Nspin
               do iorb=1,Norb
                  do jorb=1,Norb
                     !
                     Ftmp(idelta) = Ftmp(idelta) + abs(g0and(ispin,jspin,iorb,jorb,idelta)-FGmatrix(ispin,jspin,iorb,jorb,idelta))**2
                     do j=1,size(a)
                        df(idelta,j) = df(idelta,j) + &
                           real(g0and(ispin,jspin,iorb,jorb,idelta) - FGmatrix(ispin,jspin,iorb,jorb,idelta)) * &
                           real(dg0and(ispin,jspin,iorb,jorb,idelta,j)) + &
                           imag(g0and(ispin,jspin,iorb,jorb,idelta) - FGmatrix(ispin,jspin,iorb,jorb,idelta)) * &
                           imag(dg0and(ispin,jspin,iorb,jorb,idelta,j))
                     enddo
                  enddo
               enddo
            enddo
         enddo
         Ftmp(idelta) = cg_pow * (sqrt(Ftmp(idelta))**(cg_pow-2)) / Wdelta(idelta)
         dchi_freq(idelta,:) = Ftmp(idelta) * df(idelta,:)
      enddo
      !
      dchi2 = sum(dchi_freq,1)/Ldelta/(Nspin*Norb)
      !
#ifdef _DEBUG
      if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_general_frobenius. dChi2:",dchi2
#endif
      !
   end function grad_chi2_weiss_general_frobenius








  !##################################################################
  ! THESE PROCEDURES EVALUATE THE ANDERSON FUNCTIONS:
  ! - \Delta (hybridization)
  ! - g0     (weiss field)
  !##################################################################
  !ACHTUNG! We use a direct dump of the array into the necessary element of the bath.
  ! rather than using aux functions in ED_BATH. This improves execution speed. 
  function delta_general(a) result(Delta)
    real(8),dimension(:)                               :: a
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
    integer                                            :: ispin,jspin,iorb,jorb,ibath
    integer                                            :: i,stride
    complex(8),dimension(Nspin*Norb,Nspin*Norb)        :: Haux,Htmp,Vk
    complex(8),dimension(Nspin,Nspin,Norb,Norb)        :: invH_knn
    real(8),dimension(Nbath,Nspin*Norb)                :: dummy_Vbath
    type(nsymm_vector),dimension(Nbath)                :: dummy_lambda
    !
    !Get Hs
    stride = 0
    do ibath=1,Nbath
       allocate(dummy_lambda(ibath)%element(Nlambdas))
       !
       dummy_vbath(ibath,:) = a(stride+1:stride+Norb*Nspin)
       stride = stride + Norb*Nspin
       dummy_lambda(ibath)%element=a(stride+1:stride+Nlambdas)
       stride=stride+Nlambdas
    enddo
    !
    Delta=zero
    do ibath=1,Nbath
       Vk = dzdiag(dummy_Vbath(ibath,:))
       invH_knn = Hgeneral_build(dummy_lambda(ibath)%element)
       Htmp     = nn2so_reshape( invH_knn,Nspin,Norb)
       do i=1,Ldelta
          Haux     = zeye(Nspin*Norb)*xi*Xdelta(i) - Htmp
          call inv(Haux)
          Haux = matmul(Vk,matmul(Haux,Vk))
          invH_knn = so2nn_reshape(Haux,Nspin,Norb)
          !
          Delta(:,:,:,:,i)=Delta(:,:,:,:,i) + invH_knn
       enddo
       deallocate(dummy_lambda(ibath)%element)
    enddo
    !
  end function delta_general
  !
  function g0and_general(a) result(G0and)
    real(8),dimension(:)                               :: a
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: G0and,Delta
    complex(8),dimension(Norb*Nspin,Norb*Nspin)        :: FGorb
    integer                                            :: i,Nso
    !
    Nso   = Norb*Nspin
    Delta = delta_general(a)
    do i=1,Ldelta
       FGorb = (xi*Xdelta(i)+xmu)*zeye(Nso) - nn2so_reshape(impHloc + Delta(:,:,:,:,i), Nspin,Norb)
       call inv(FGorb)
       G0and(:,:,:,:,i) = so2nn_reshape(FGorb,Nspin,Norb)
    enddo
  end function g0and_general
  !
  !
  !
  function grad_delta_general(a) result(dDelta)
    real(8),dimension(:)                                      :: a
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
    integer                                                   :: ispin,iorb,jorb,ibath
    integer                                                   :: i,k,ik,l,io,ig,counter
    complex(8),dimension(Nspin*Norb,Nspin*Norb)               :: H_reconstructed,Hbasis_so,Htmp,Vk
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Ldelta)        :: Haux
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)        :: invH_knn
    real(8),dimension(Nbath,Norb*Nspin)                       :: dummy_Vbath
    type(nsymm_vector),dimension(Nbath)                       :: dummy_lambda
    !
    !
    !Get Hs
    counter = 0
    do ibath=1,Nbath
       allocate(dummy_lambda(ibath)%element(Nlambdas))
       !
       dummy_vbath(ibath,:) = a(counter+1:counter+Nspin*Norb)
       counter = counter + Nspin*Norb
       dummy_lambda(ibath)%element=a(counter+1:counter+Nlambdas)
       counter = counter + Nlambdas
    enddo
    !
    dDelta=zero
    counter=0
    do ibath=1,Nbath
       Vk = dzdiag(dummy_vbath(ibath,:))
       H_reconstructed= nn2so_reshape( Hgeneral_build(dummy_lambda(ibath)%element) ,Nspin,Norb)
       do i=1,Ldelta
          Haux(:,:,i) = zeye(Nspin*Norb)*xi*Xdelta(i) - H_reconstructed
          call inv(Haux(:,:,i))
          !Derivate_V^ibath_g Haux(a,b) = Haux*V(a,b)*delta(a,g)+V*Haux(a,b)*delta(b,g)
          do ig=1,Nspin*Norb
             Htmp(:,:)  = 0.d0
             Htmp(:,ig) = Htmp(:,ig) + VH_prod(Vk,Haux(:,:,i),ig)
             Htmp(ig,:) = Htmp(ig,:) + HV_prod(Vk,Haux(:,:,i),ig)
             dDelta(:,:,:,:,i,counter+ig)=so2nn_reshape(Htmp,Nspin,Norb)
          enddo
       enddo
       counter = counter + Nspin*Norb

       !
       !Derivate_lambda_p
       do k=1,Nlambdas
          counter = counter + 1
          Hbasis_so=nn2so_reshape(Hgeneral_basis(k)%O,Nspin,Norb)
          do l=1,Ldelta
             Htmp = ((Haux(:,:,l) .x. Hbasis_so)) .x. Haux(:,:,l)
             Htmp = ((Vk .x. Htmp) .x. Vk)
             dDelta(:,:,:,:,l,counter)=so2nn_reshape(Htmp,Nspin,Norb)
          enddo
       enddo
       !
    enddo
  end function grad_delta_general
  !
  !
  function grad_g0and_general(a) result(dG0and)
    real(8),dimension(:)                                       :: a
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dG0and
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)         :: G0and
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
    complex(8),dimension(Nspin*Norb,Nspin*Norb)                :: dDelta_so,dG0and_so,G0and_so
    integer                                                    :: ispin,iorb,jorb
    integer                                                    :: ik,l
    !
    G0and  = g0and_general(a)
    dDelta = grad_delta_general(a)
    !
    dG0and = zero
    do l=1,Ldelta
       G0and_so=nn2so_reshape(g0and(:,:,:,:,l),Nspin,Norb)
       do ik=1,size(a)
          dDelta_so=nn2so_reshape(dDelta(:,:,:,:,l,ik),Nspin,Norb)
          ! dG0and_lso=matmul(-G0and_lso,dDelta_lso)
          ! dG0and_lso=matmul(dG0and_lso,G0and_lso)
          dG0and_so = (G0and_so .x. dDelta_so) .x. G0and_so
          dG0and(:,:,:,:,l,ik)=so2nn_reshape(dG0and_so,Nspin,Norb) !Check the sign, should it be +
       enddo
    enddo
    !
  end function grad_g0and_general


  !+-------------------------------------------------------------+
  !
  ! GENERALSUPERC UNTESTED
  !  
  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for GENERAL BATH SUPERC
  !MANY THINGS HERE NEED TO BE REVISED TO BE WRITTEN PROPERLY
  !THE STRUCTURE OF THIS PART OF THE CODE IS MESSY
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_general_superc(fg,bath_)
    complex(8),dimension(:,:,:,:,:,:)             :: fg ![2][Nspin][Nspin][Norb][Norb][Lmats]
    logical,dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb)    :: Hmask
    real(8),dimension(:),intent(inout)            :: bath_
    real(8),dimension(:),allocatable              :: array_bath
    integer                                       :: i,j,iorb,jorb,ispin,jspin,io,jo,ibath,inambu
    integer                                       :: iter,stride,counter,Asize
    real(8)                                       :: chi
    logical                                       :: check
    type(effective_bath)                          :: dmft_bath
    character(len=256)                            :: suffix
    integer                                       :: unit
    complex(8),dimension(:,:,:,:,:,:),allocatable :: fgand ![2][Nspin][][Norb][][Ldelta]  
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG chi2_fitgf_general_superc: Fit"
#endif
    !
    call assert_shape(fg,[2,Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf_general_superc","fg")
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_general_superc error: wrong bath dimensions"
    !
    call allocate_dmft_bath(dmft_bath)
    call set_dmft_bath(bath_,dmft_bath)
    allocate(array_bath(size(bath_)-1))
    Nlambdas  =bath_(1)
    array_bath=bath_(2:)
    !
    if(ed_all_g)then
       Hmask=.true.
       write(*,*) "!! WARNING !!"
       write(*,*) "USING ED_ALL_G=.true."
       write(*,*) "superc-general is able to choose which components should be fitted"
       write(*,*) "the use of ED_ALL_G=.false. is strongly suggested to avoid fit problems"
       ! Aren't we sure about hermiticity?
       ! -> Hmask=Hgeneral_mask(wdiag=.false.,uplo=.true.)
       !HERE WE ARE ON GENERAL SUPERC
       !WE ONLY WANT TO FIT (1,1) and (1,2) SPIN
       !SECTORS CORRESPONDING TO G0 AND F0
       !THEREFORE THE SECOND SPIN INDEX IS EQUIVALENT
       !TO THE INDEX SELECTING G0/F0
    else
       Hmask=Hgeneral_mask(wdiag=.true.,uplo=.false.)
    endif
    !
    ! Taking only upper-diagonal
    do ispin=1,Nspin*Nnambu
       do jspin=1,Nspin*Nnambu
          do iorb=1,Norb
             do jorb=1,Norb
                if(ispin>jspin .or. iorb>jorb)Hmask(ispin,jspin,iorb,jorb)=.false.
             enddo
          enddo
       enddo
    enddo
    !
    totNso=count(Hmask)
    !
    allocate(getInambu(totNso))
    allocate(getIspin(totNso),getJspin(totNso))
    allocate(getIorb(totNso) ,getJorb(totNso))
    counter=0
    !HERE NSPIN=1
    do ispin=1,Nspin*Nnambu
       do jspin=1,Nspin*Nnambu
          do iorb=1,Norb
             do jorb=1,Norb
                if (Hmask(ispin,jspin,iorb,jorb))then
                   counter=counter+1
                   getInambu(counter) = 1
                   if(ispin/=jspin)getInambu(counter)=2
                   getIspin(counter)  = 1
                   getIorb(counter)   = iorb
                   getJspin(counter)  = 1
                   getJorb(counter)   = jorb
                endif
             enddo
          enddo
       enddo
    enddo
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,6))Ldelta=size(fg,6)
    !
    allocate(Gdelta(totNso,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
    !
    select case(cg_weight)
    case default
       Wdelta=1d0
    case(2)
       Wdelta=1d0*arange(1,Ldelta)
    case(3)
       Wdelta=Xdelta
    end select
    !
    write(LOGfile,*)"Fitted functions",totNso
    do i=1,totNso
       Gdelta(i,1:Ldelta) = fg(getInambu(i),getIspin(i),getJspin(i),getIorb(i),getJorb(i),1:Ldelta)
    enddo
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG chi2_fitgf_general_superc: cg_method:"//str(cg_method)//&
         ", cg_grad:"//str(cg_grad)//&
         ", cg_scheme:"//str(cg_scheme)
#endif
    !
    select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
    case default
       if(cg_grad==0)then
          stop "ERROR: general-superc does not support minimization with analytical gradient (cg_grad=0). You MUST use numerical gradient (cg_grad=1)"
       else
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,chi2_weiss_general_superc,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,chi2_delta_general_superc,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_general_superc error: cg_scheme != [weiss,delta]"
          end select
       endif
       !
       !
    case (1)
       if(cg_grad==0)then
          write(*,*) "                                                                                "
          write(*,*) "WARNING: analytic gradient not available with cg-method=1 (minimize f77 routine)"
          write(*,*) "         > we will force cg_grad=1 (so let the routine estimate the gradient)   "
          write(*,*) "                                                                                "
       endif
       select case (cg_scheme)
       case ("weiss")
          call fmin_cgminimize(array_bath,chi2_weiss_general_superc,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,chi2_delta_general_superc,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_general_superc error: cg_scheme != [weiss,delta]"
       end select
       !
    end select
    !
    write(LOGfile,"(A,ES18.9,A,I5,A)")"chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,"  <--  All Orbs, All Spins"
    !
    suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
    unit=free_unit()
    open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
    write(unit,"(ES18.9,1x,I5)") chi,iter
    close(unit)
    !
    bath_(2:size(bath_))=array_bath
    call set_dmft_bath(bath_,dmft_bath) ! *** array_bath --> dmft_bath *** (per write fit result)
    call write_dmft_bath(dmft_bath,LOGfile)
    call save_dmft_bath(dmft_bath)
    !
    allocate(fgand(2,Nspin,Nspin,Norb,Norb,Ldelta))
    if(cg_scheme=='weiss')then
       fgand(1,:,:,:,:,:) = g0and_bath_function(xi*Xdelta(:),dmft_bath)
       fgand(2,:,:,:,:,:) = f0and_bath_function(xi*Xdelta(:),dmft_bath)
    else
       fgand(1,:,:,:,:,:) = delta_bath_function(xi*Xdelta(:),dmft_bath)
       fgand(2,:,:,:,:,:) =fdelta_bath_function(xi*Xdelta(:),dmft_bath)
    endif
    call write_fit_result()
    deallocate(fgand)
    !
    call get_dmft_bath(dmft_bath,bath_)                ! ***  dmft_bath --> bath_ ***    (bath in output)
    call deallocate_dmft_bath(dmft_bath)
    deallocate(Gdelta,Xdelta,Wdelta)
    deallocate(getIspin,getJspin)
    deallocate(getIorb,getJorb)
    deallocate(getInambu)
    deallocate(array_bath)
    !
  contains
    !
    subroutine write_fit_result()
      integer            :: i,j,s,l,iorb,jorb,ispin,jspin
      character(len=7)   :: ANtype
      !
      do l=1,totNso
         iorb = getIorb(l)
         jorb = getJorb(l)
         ispin = getIspin(l)
         jspin = getJspin(l)
         if(jspin==1)then
            ANtype="Weiss"
         else
            ANtype="AWeiss"
         endif
         suffix="_"//reg(ANtype)//&
              "_l"//reg(txtfy(iorb))//&
              "_m"//reg(txtfy(jorb))//&
              "_s"//reg(txtfy(ispin))//&
              reg(ed_file_suffix)
         unit=free_unit()
         if(cg_scheme=='weiss')then
            open(unit,file="fit_weiss"//reg(suffix)//".ed")
         else
            open(unit,file="fit_delta"//reg(suffix)//".ed")
         endif
         do i=1,Ldelta
            write(unit,"(5F24.15)")Xdelta(i),&
                 dimag(fg(jspin,ispin,ispin,iorb,jorb,i)),dimag(fgand(jspin,ispin,ispin,iorb,jorb,i)),&
                 dreal(fg(jspin,ispin,ispin,iorb,jorb,i)),dreal(fgand(jspin,ispin,ispin,iorb,jorb,i))
         enddo
         close(unit)
      enddo
    end subroutine write_fit_result
    !
  end subroutine chi2_fitgf_general_superc
  
  !##################################################################
  ! THESE PROCEDURES EVALUATE THE \chi^2 FUNCTIONS TO MINIMIZE. 
  !##################################################################
  !
  !
  !+-----------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
  !+-----------------------------------------------------------------+
  function chi2_delta_general_superc(a) result(chi2)
    real(8),dimension(:)                                 :: a
    real(8)                                              :: chi2
    real(8),dimension(totNso)                            :: chi2_so
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
    real(8),dimension(Ldelta)                            :: Ctmp
    integer                                              :: i,l,iorb,jorb,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_general_superc. a:",a
#endif
    !
    Delta = delta_general_superc(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       !HERE JSPIN GIVES G0/F0 SEE ABOVE
       Ctmp =  abs(Gdelta(l,:)-Delta(jspin,ispin,ispin,iorb,jorb,:))
       chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_so)
    chi2=chi2/Ldelta/totNso
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_general_superc. Chi**2:",chi2
#endif
    !
  end function chi2_delta_general_superc
  !
  !
  !+--------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function.
  !+--------------------------------------------------------------+
  function chi2_weiss_general_superc(a) result(chi2)
    real(8),dimension(:)                                 :: a
    real(8)                                              :: chi2
    real(8),dimension(totNso)                            :: chi2_so
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
    real(8),dimension(Ldelta)                            :: Ctmp
    integer                                              :: i,l,iorb,jorb,ispin,jspin,inambu
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_general_superc. a:",a
#endif
    !
    g0and = g0and_general_superc(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       inambu= getInambu(l)
       !
       Ctmp = abs(Gdelta(l,:)-g0and(inambu,ispin,jspin,iorb,jorb,:))
       chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_so)
    chi2=chi2/Ldelta/totNso
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_general_superc. Chi**2:",chi2
#endif
    !
  end function chi2_weiss_general_superc
  
  
  !##################################################################
  ! THESE PROCEDURES EVALUATE THE ANDERSON FUNCTIONS:
  ! - \Delta (hybridization)
  ! - g0     (weiss field)
  !##################################################################
  !ACHTUNG! We use a direct dump of the array into the necessary element of the bath.
  ! rather than using aux functions in ED_BATH. This improves execution speed. 
  function delta_general_superc(a) result(Delta)
    real(8),dimension(:)                                 :: a
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
    integer                                              :: ispin,jspin,iorb,jorb,ibath
    integer                                              :: i,stride,Nso
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)          :: Haux,Htmp,Vk
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb)          :: invH_knn
    real(8),dimension(Nbath,Nnambu*Nspin*Norb)           :: dummy_Vbath
    type(nsymm_vector),dimension(Nbath)                  :: dummy_lambda
    complex(8),dimension(Nnambu*Norb,Nnambu*Norb)        :: JJ
    !
    JJ=kron(pauli_sigma_z,zeye(Norb))
    !Get Hs
    Nso=Nspin*Norb
    stride = 0
    do ibath=1,Nbath
       allocate(dummy_lambda(ibath)%element(Nlambdas))
       !
       dummy_vbath(ibath,:Nso)    = a(stride+1:stride+Nso)
       dummy_vbath(ibath,Nso+1:)  =-a(stride+1:stride+Nso)
       stride = stride + Nso
       dummy_lambda(ibath)%element= a(stride+1:stride+Nlambdas)
       stride = stride+Nlambdas
    enddo
    !
    Delta=zero
    do ibath=1,Nbath
       Vk        = dzdiag(dummy_vbath(ibath,:))
       invH_knn  = Hgeneral_build(dummy_lambda(ibath)%element)
       Htmp      = nn2so_reshape( invH_knn,Nnambu*Nspin,Norb)
       do i=1,Ldelta
          Haux     = zeye(Nnambu*Nspin*Norb)*xi*Xdelta(i) - Htmp
          call inv(Haux)
          
          Haux     = matmul(matmul(Vk,Haux),Vk)
          invH_knn = so2nn_reshape(Haux,Nnambu*Nspin,Norb)
          !
          Delta(1,1,1,:,:,i)=Delta(1,1,1,:,:,i) + invH_knn(1,1,:,:) !normal
          Delta(2,1,1,:,:,i)=Delta(2,1,1,:,:,i) + invH_knn(1,2,:,:) !anomalous
       enddo
       deallocate(dummy_lambda(ibath)%element)
    enddo
    !
  end function delta_general_superc
  !
  function g0and_general_superc(a) result(G0and)
    real(8),dimension(:)                                 :: a
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Ldelta) :: G0and,Delta
    complex(8),dimension(Nnambu*Norb*Nspin,Nnambu*Norb*Nspin)          :: FGorb
    integer                                              :: i,Nson,Nso
    !
    Nso   = Norb*Nspin
    Nson  = Norb*Nspin*Nnambu
    Delta = delta_general_superc(a)
    !
    do i=1,Ldelta
       FGorb(      :Norb,      :Norb) = (xi*Xdelta(i)+xmu)*zeye(Nso) - nn2so_reshape(impHloc + Delta(1,:,:,:,:,i), Nspin,Norb)
       FGorb(      :Norb,Norb+1:    ) = - nn2so_reshape( Delta(2,:,:,:,:,i), Nspin,Norb)
       FGorb(Norb+1:    ,      :Norb) = - conjg(nn2so_reshape(Delta(2,:,:,:,:,i), Nspin,Norb))
       FGorb(Norb+1:    ,Norb+1:    ) = (xi*Xdelta(i)-xmu)*zeye(Nso) + conjg(nn2so_reshape(impHloc + Delta(1,:,:,:,:,i), Nspin,Norb))
       call inv(FGorb)
       G0and(1,:,:,:,:,i) = so2nn_reshape(FGorb(1:Norb,     1:Norb  ),Nspin,Norb)
       G0and(2,:,:,:,:,i) = so2nn_reshape(FGorb(1:Norb,Norb+1:2*Norb),Nspin,Norb)
    enddo
  end function g0and_general_superc

  !Diagonal complex matrix from real vector
  function dzdiag(xx) result(AA)
    real(8),dimension(:)                   :: xx
    complex(8),dimension(:,:),allocatable  :: AA
    integer                                :: NN,ii

    NN=size(xx,1)
    allocate(AA(NN,NN))
    AA=0.d0
    do ii=1,NN
       AA(ii,ii)=xx(ii)
    enddo
  end function dzdiag

  !Auxiliary function for general derivative in V_k
  function HV_prod(vv,HH,ig) result(AA)
    complex(8),dimension(:,:)             :: HH,vv
    complex(8),dimension(:,:),allocatable :: Atmp
    complex(8),dimension(:),allocatable   :: AA
    integer                               :: N,i,ig
    N=size(vv,1)
    allocate(Atmp(N,N),AA(N))
    Atmp=matmul(HH,vv)
    AA=Atmp(ig,:)   
  end function HV_prod
  function VH_prod(vv,HH,ig) result(AA)
    complex(8),dimension(:,:)             :: HH,vv
    complex(8),dimension(:,:),allocatable :: Atmp
    complex(8),dimension(:),allocatable   :: AA
    integer                               :: N,i,ig
    N=size(vv,1)
    allocate(Atmp(N,N),AA(N))
    Atmp=matmul(vv,HH)
    AA=Atmp(:,ig)   
  end function VH_prod

END MODULE ED_FIT_GENERAL
