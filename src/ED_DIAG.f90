MODULE ED_DIAG
  USE ED_INPUT_VARS
  USE ED_DIAG_NORMAL
  USE ED_DIAG_SUPERC
  USE ED_DIAG_NONSU2
  !
  implicit none
  private

  !>Diag hamiltonian
  public  :: diagonalize_impurity

contains

  subroutine  diagonalize_impurity()
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG diagonalize_impurity: Start digonalization"
#endif
    !
    write(LOGfile,"(A)")"Diagonalize impurity problem:"
    select case(ed_mode)
    case default  ;call diagonalize_impurity_normal()
    case("superc");call diagonalize_impurity_superc()
    case("nonsu2");call diagonalize_impurity_nonsu2()
    end select
#ifdef _DEBUG
    write(Logfile,"(A)")""
#endif
  end subroutine diagonalize_impurity

end MODULE ED_DIAG
