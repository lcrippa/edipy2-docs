.. _bath:

Bath Manipulation
===========================

These functions manipulate the user-accessible bath array

.. function:: chi2_fitgf(*args,ispin=0,iorb=None,fmpi=True)

   This function fits the Weiss field or Hybridization function (delta) with a discretized version. The fit parameters are the bath parameters contained in the user-accessible array. Depending on the type of system we are considering (normal, superconductive, non-SU(2)) a different set of inputs has to be passed.
   
    
   :type args: [np.array(dtype=complex,np.array(dtype=complex),np.array(dtype=float)] **or** [np.array(dtype=complex,np.array(dtype=float)]
   :param args: The positional arguments are the function(s) to fit and the bath array. 
   
    If the system is not superconductive (:code:`ED_MODE=NORMAL` or :code:`ED_MODE=NONSU2`) the argumens are
   
    * :code:`g`: the function to fit
    * :code:`bath`: the bath
   
    If the system is superconductive (:code:`ED_MODE=SUPERC`) the arguments are

    * :code:`g`: the normal function to fit
    * :code:`f`: the anomalous function to fit
    * :code:`bath`: the bath 
   
    The dimensions of the previous arrays can vary:
   
    The dimension of :code:`bath` can be
 
    * :code:`Nb`: if single-site, the output of :func:`get_bath_gimension`
    * :code:`[Nlat,Nb]`: if real-space DMFT
   
    Accordingly, the dimension of g (and f) can be:
   
    * :code:`3`: in the single-site case,  an array of the shape :code:`[ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, Lfit]`. 
    * :code:`3`: in the real-space DMFT case, an array of the shape :code:`[Nlat*ed.Nspin*ed.Norb, Nlat*ed.Nspin*ed.Norb, Lfit]`
    * :code:`4`: in the real-space DMFT case, an array of the shape :code:`[Nlat, ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, Lfit]`
    * :code:`5`: in the single-site case, an array of the shape :code:`[ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, Lfit]`
    * :code:`6`: in the real-space DMFT case, an array of the shape :code:`[Nlat, ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, Lfit]`
   
    where :code:`Lfit` is a given number of frequencies.

   
 
   :type ispin: int 
   :param ispin: spin species to be fitted. If :code:`ed.Nspin=2`, the fitting function needs to be called twice. Only the corresponding elements of :code:`bath` will be updated each time
    
   :type iorb: int 
   :param iorb: the orbital to be fitted. If omitted, all orbitals will be fitted
   
   :type fmpi: bool 
   :param fmpi: flag to automatically do and broadcast the fit over MPI, if defined

   :raise ValueError: if the shapes of the positional arguments are incompatible
   :raise ValueError: if a number of positional arguments different from 2 or 3 are passed   
     
   :return: An array of floats that contains the bath parameters for the impurity problem. This is a required input of :func:`solve` and :func:`chi2_fitgf`. Its elements are ordered differently depending on the bath geometry. They are (de)compactified for user interaction via :func:`bath_packaging`. Specific symmetrization operations are implemented and listed in the :ref:`bath` section.
   :rtype: np.array(dtype=float) 





.. function:: get_bath_dimension()

   This function returns the correct dimension for the bath to be allocated (for each impurity) given the parameters of the system.
   
   :return: a number which is the dimension of the bath array for each impurity.
   :rtype: int  
