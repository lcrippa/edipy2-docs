.. _envinit:

Impurity Problem Initialization, Solution and Finalization
================================================================


.. function:: init_solver()

   This function initializes the ED environment for the impurity problem solution, and sets the bath reading it from the ``hamiltonian.restart`` file or initializing it in a symmetric way.
   The function can take different argument combinations.

   :type empty: None
   :param empty: If no argument is provided, a bath array is initialized automatically by the module for a single-site impurity problem with the specified bath geometry and number of levels.

    
   :type bath: np.array(dtype=float) **or** [float]
   :param bath: If a bath array is provided, it has to be a numpy array or a tuple of floats. It has to have a specific length, as specified by :func:`get_bath_dimension`
   
   
    The array is ordered in F convention inside the function.
    
    **Note**: if the dimension of the bath is at odds with the global system parameters, EDIpack2 will terminate execution with an error.
    
   :type dimension: int **or** [int]
   :param dimension: If an integer or a tuple of length 1 number is provided, a bath is initialized with that dimension
    
   :type dimensions: [int,int]
   :param dimensions: If a tuple of length 2 is provided, a bath is initialized with that shape. The first dimension is the number of inequivalent sites for real-space DMFT, the second is the bath size for each inequivalent site.
   
     
   :return: An array of floats that contains the bath parameters for the impurity problem. This is a required input of :command:`solve` and :command:`chi2_fitgf`. Its components are ordered differently depending on the bath geometry. They are (de)compactified for user interaction via :command:`bath_packaging`. Specific symmetrization operations are implemented and listed in the :ref:`bath` section.
   :rtype: np.array(dtype=float) 
    




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
    
    **Note**: the way the EDIpack2 library passes from 1 comulative to 2 or 3 running indices is, from slower to faster: ``lat``, ``spin``, ``orb``
    
   :type Nlat: int
   :param Nlat: Number of inequivalent sites for real-space DMFT. The function will raise a ValueError if the dimensions of ``hloc`` are inconsistent with the presence or absence of Nlat. 
    The EDIpack2 library will check the correctness of the dimensions of ``hloc`` and terminate execution if inconsistent.
   
   :raise ValueError: If hloc is not provided (as in edipy1) or has the wrong shape
   
   :return: Nothing
   :rtype: None


.. function:: solve(bath, sflag=True, iflag=True, fmpi=True, mpi_lanc=False)

   This function solves the impurity problem and calculates the observables, Green's function and self-energy.

   :type bath: np.array(dtype=float) 
   :param bath: The bath array returned by  :func:`init_solver`. If the bath dimensions are inconsistent with the global properties of the problem, EDIpack2 will exit with an error.
   
   :type sflag: bool
   :param sflag: for single-site DMFT, if :code:`False`, it disables the calculation of the Green's function and susceptibilities
   
   :type iflag: bool
   :param iflag: for real-space DMFT, if :code:`False`, it disables the calculation of the Green's function and susceptibilities
   
   :type fmpi: bool
   :param fmpi: if :code:`False`, for single-site DMFT, it disables MPI for the ED routine, if the communicator is used elsewhere
   
   :type mpi_lanc: bool
   :param mpi_lanc: if :code:`True`, for real-space DMFT sets the MPI parallelization for the ED routine. By default it is :code:`False`, and each inequivalent site is solved serially by a different core.
        
   :return: Nothing
   :rtype: None

.. function:: finalize_solver()

   This function cleans up the ED environment, deallocates the relevant arrays and makes a second call to :command:`init_solver` possible
           
   :return: Nothing
   :rtype: None

