MODULE ED_HAMILTONIAN
  USE ED_HAMILTONIAN_NORMAL
  USE ED_HAMILTONIAN_SUPERC
  USE ED_HAMILTONIAN_NONSU2
  !
  implicit none
  private

  !>Build sparse hamiltonian of the sector
  public  :: build_Hv_sector_normal
  public  :: build_Hv_sector_superc
  public  :: build_Hv_sector_nonsu2
  public  :: delete_Hv_sector_normal
  public  :: delete_Hv_sector_superc
  public  :: delete_Hv_sector_nonsu2
  !
  public  :: vecDim_Hv_sector_normal
  public  :: vecDim_Hv_sector_superc
  public  :: vecDim_Hv_sector_nonsu2

  !> Tridiag sparse Hamiltonian of the sector
  public  :: tridiag_Hv_sector_normal
  public  :: tridiag_Hv_sector_superc
  public  :: tridiag_Hv_sector_nonsu2




end MODULE ED_HAMILTONIAN
