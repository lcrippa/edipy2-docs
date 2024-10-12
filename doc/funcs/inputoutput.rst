Input / Output Functions
===========================

These functions pass to the userspace the observables and response functions calculated and stored within the EDIpack2 library.

.. function:: get_sigma(ilat=None,ishape=None,axis="m",typ="n")

   This function gets from the EDIpack2 library the value of the self-energy calculated on the Matsubara or real-frequency axis, with parameters specified in the input file.
    
   :type ilat: int
   :param ilat: if the case of real-space DMFT, if only the self-energy of a specific inequivalent site is needed, this can be specified.
        
   :type ishape: int 
   :param ishape: this variable determines the shape of the returned array. Possible values:
   
    * :code:`None`: the same shape as :code:`Hloc` plus one axis for frequency 
    * :code:`3`: in the single-site case, it will return an array of the shape :code:`[ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, Nfrequencies=ed.Lmats/ed.Lreal]`. In the real-space DMFT case, it will return an array of the shape :code:`[Nlat*ed.Nspin*ed.Norb, Nlat*ed.Nspin*ed.Norb, Nfrequencies=Lmats/Lreal]`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, ValueError is returned.
    * :code:`4`: in the real-space DMFT case, it will return an array of the shape :code:`[Nlat, ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, Nfrequencies=ed.Lmats/ed.Lreal]`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, the output will have one dimension less.
    * :code:`5`: in the single-site case, it will return an array of the shape :code:`[ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, Nfrequencies=Lmats/Lreal]`.
    * :code:`6`: in the real-space DMFT case, it will return an array of the shape :code:`[Nlat, ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, Nfrequencies=ed.Lmats/ed.Lreal]`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, the output will have one dimension less.
   
    
   :type axis: str 
   :param axis: the axis along which to return the self-energy. Can be :code:`m` for Matsubara or :code:`r` for real. In the first case, the number of frequencies will be given by :code:`ed.Lmats`, in the second by :code:`ed.Lreal`.
   
   :type typ: str 
   :param typ: whether to return the normal or anomalous self-energy (for the superconducting case). Can be :code:`n` for normal or :code:`a` for anomalous.
   
   :raise ValueError: If :code:`ishape` is incompatible woth :code:`ilat` or not in the previous list.
   :raise ValueError: If :code:`axis` is not in the previous list.
     
   :return: An array of floats that contains the self-energy along the specific axis, with dimension set by :code:`ishape` and :code:`axis`.  
   :rtype: np.array(dtype=float) 



.. function:: build_sigma(zeta,ilat=None,ishape=None,typ="n")

   This function generates the self-energy for a user-chosen set of frequencies in the complex plane

   :type zeta: complex **or** [complex] **or** np.array(dtype=complex)
   :param zeta: user-defined array of frequencies in the whole complex plane.

    
   :type ilat: int
   :param ilat: if the case of real-space DMFT, if only the self-energy of a specific inequivalent site is needed, this can be specified.
        
   :type ishape: int 
   :param ishape: this variable determines the shape of the returned array. Possible values:
   
    * :code:`None`: the same shape as :code:`Hloc` plus one axis for frequency 
    * :code:`3`: in the single-site case, it will return an array of the shape :code:`[ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, len(zeta)]`. In the real-space DMFT case, it will return an array of the shape :code:`[Nlat*ed.Nspin*ed.Norb, Nlat*ed.Nspin*ed.Norb, len(zeta)]`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, ValueError is returned.
    * :code:`4`: in the real-space DMFT case, it will return an array of the shape :code:`[Nlat, ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, len(zeta)`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, the output will have one dimension less.
    * :code:`5`: in the single-site case, it will return an array of the shape :code:`[ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, len(zeta)]`.
    * :code:`6`: in the real-space DMFT case, it will return an array of the shape :code:`[Nlat, ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, len(zeta)]`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, the output will have one dimension less.
        
   :type typ: str 
   :param typ: whether to return the normal or anomalous self-energy (for the superconducting case). Can be :code:`n` for normal or :code:`a` for anomalous.
   
   :raise ValueError: If :code:`ishape` is incompatible woth :code:`ilat` or not in the previous list.
   :raise ValueError: If :code:`axis` is not in the previous list.
     
   :return: An array of floats that contains the self-energy along the specific axis, with dimension set by :code:`ishape` and :code:`zeta`.  
   :rtype: np.array(dtype=float) 




.. function:: get_gimp(self,ilat=None,ishape=None,axis="m",typ="n")

   This function gets from the EDIpack2 library the value of the Green's function calculated on the Matsubara or real-frequency axis, with parameters specified in the input file.
    
   :type ilat: int
   :param ilat: if the case of real-space DMFT, if only the Green's function of a specific inequivalent site is needed, this can be specified.
        
   :type ishape: int 
   :param ishape: this variable determines the shape of the returned array. Possible values:
   
    * :code:`None`: the same shape as :code:`Hloc` plus one axis for frequency 
    * :code:`3`: in the single-site case, it will return an array of the shape :code:`[ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, Nfrequencies=ed.Lmats/ed.Lreal]`. In the real-space DMFT case, it will return an array of the shape :code:`[Nlat*ed.Nspin*ed.Norb, Nlat*ed.Nspin*ed.Norb, Nfrequencies=Lmats/Lreal]`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, ValueError is returned.
    * :code:`4`: in the real-space DMFT case, it will return an array of the shape :code:`[Nlat, ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, Nfrequencies=ed.Lmats/ed.Lreal]`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, the output will have one dimension less.
    * :code:`5`: in the single-site case, it will return an array of the shape :code:`[ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, Nfrequencies=Lmats/Lreal]`.
    * :code:`6`: in the real-space DMFT case, it will return an array of the shape :code:`[Nlat, ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, Nfrequencies=ed.Lmats/ed.Lreal]`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, the output will have one dimension less.
   
    
   :type axis: str 
   :param axis: the axis along which to return the Green's functio. Can be :code:`m` for Matsubara or :code:`r` for real. In the first case, the number of frequencies will be given by :code:`ed.Lmats`, in the second by :code:`ed.Lreal`.
   
   :type typ: str 
   :param typ: whether to return the normal or anomalous Green's function (for the superconducting case). Can be :code:`n` for normal or :code:`a` for anomalous.
   
   :raise ValueError: If :code:`ishape` is incompatible woth :code:`ilat` or not in the previous list.
   :raise ValueError: If :code:`axis` is not in the previous list.
     
   :return: An array of floats that contains the Green's function along the specific axis, with dimension set by :code:`ishape` and :code:`axis`.  
   :rtype: np.array(dtype=float) 






.. function:: build_gimp(zeta,ilat=None,ishape=None,typ="n")

   This function generates the Green's function for a user-chosen set of frequencies in the complex plane

   :type zeta: complex **or** [complex] **or** np.array(dtype=complex)
   :param zeta: user-defined array of frequencies in the whole complex plane.

    
   :type ilat: int
   :param ilat: if the case of real-space DMFT, if only the Green's function of a specific inequivalent site is needed, this can be specified.
        
   :type ishape: int 
   :param ishape: this variable determines the shape of the returned array. Possible values:
   
    * :code:`None`: the same shape as :code:`Hloc` plus one axis for frequency 
    * :code:`3`: in the single-site case, it will return an array of the shape :code:`[ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, len(zeta)]`. In the real-space DMFT case, it will return an array of the shape :code:`[Nlat*ed.Nspin*ed.Norb, Nlat*ed.Nspin*ed.Norb, len(zeta)]`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, ValueError is returned.
    * :code:`4`: in the real-space DMFT case, it will return an array of the shape :code:`[Nlat, ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, len(zeta)`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, the output will have one dimension less.
    * :code:`5`: in the single-site case, it will return an array of the shape :code:`[ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, len(zeta)]`.
    * :code:`6`: in the real-space DMFT case, it will return an array of the shape :code:`[Nlat, ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, len(zeta)]`. :code:`Nlat` will be determined from the module by assessing the shape of Hloc. If :code:`ilat` is set, the output will have one dimension less.
        
   :type typ: str 
   :param typ: whether to return the normal or anomalous Green's function (for the superconducting case). Can be :code:`n` for normal or :code:`a` for anomalous.
   
   :raise ValueError: If :code:`ishape` is incompatible woth :code:`ilat` or not in the previous list.
   :raise ValueError: If :code:`axis` is not in the previous list.
     
   :return: An array of floats that contains the Green's function along the specific axis, with dimension set by :code:`ishape` and :code:`zeta`.  
   :rtype: np.array(dtype=float) 

