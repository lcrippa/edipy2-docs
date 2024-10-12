Environment initialization
===========================


.. function:: set_hloc(hloc=None,Nlat=None)

   This function sets the local Hamiltonian of the impurity problem. 
    
   :type hloc: np.array(dtype=complex)
   :param hloc: Local Hamiltonian matrix. This can have the following shapes:
   
    * :code:`[ed.Nspin*ed.Norb,ed.Nspin*ed.Norb]`: single-impurity case, 2-dimensional array
    * :code:`[ed.Nspin,ed.Nspin,ed.Norb,ed.Norb]`: single-impurity case, 4-dimensional array
    * :code:`[Nlat*ed.Nspin*ed.Norb,Nlat*ed.Nspin*ed.Norb]`: real-space DMFT case, 2-dimensional array.
    * :code:`[Nlat,Nspin*ed.Norb,ed.Nspin*ed.Norb]`: single-impurity case, 3-dimensional array.
    * :code:`[Nlat,ed.Nspin,ed.Nspin,ed.Norb,ed.Norb]`: single-impurity case, 5-dimensional array.
   
    The array is ordered in F convention inside the function.
    
   :type Nlat: np.array(dtype=int)
   :param Nlat: Number of inequivalent sites for real-space DMFT. The function will raise a ValueError if the dimensions of :code:`hloc` are inconsistent with the presence or absence of Nlat. 
    The EDIpack2 library will check the correctness of the dimensions of :code:`hloc` and terminate execution if inconsistent.
   
   :raise ValueError: If hloc is not provided (as in edipy1) or has the wrong shape
   
   :return: Nothing
   :rtype: None







