!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_subscf
!! NAME
!!  m_subscf
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2018- Xing Zhang
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_subscf

 use defs_basis
 use m_abicore
 use m_errors
 use defs_datatypes
 use defs_abitypes

 use m_cgwf
 use m_efield
 use m_gemm_nonlop

 use m_plowannier,       only : compute_oper_ks2sub
 use defs_wvltypes,      only : wvl_data
 use m_geometry,         only : metric,fixsym
 use m_crystal,          only : crystal_t
 use m_kg,               only : getph,getcut,mkkin,mkkpg
 use m_mkffnl,           only : mkffnl
 use m_pawtab,           only : pawtab_type,pawtab_get_lsize
 use m_paw_an,           only : paw_an_type,paw_an_init,paw_an_free,paw_an_nullify,paw_an_reset_flags
 use m_paw_ij,           only : paw_ij_type,paw_ij_init,paw_ij_free,paw_ij_nullify,paw_ij_reset_flags
 use m_pawfgr,           only : pawfgr_type
 use m_pawrad,           only : pawrad_type
 use m_pawang,           only : pawang_type
 use m_pawfgrtab,        only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free
 use m_paw_nhat,         only : nhatgrid,pawmknhat
 use m_paw_tools,        only : chkpawovlp
 use m_pawrhoij,         only : pawrhoij_type
 use m_paw_occupancies,  only : initrhoij
 use m_paw_denpot,       only : pawdenpot
 use m_pawdij,           only : pawdij, symdij,pawdijhat
 use m_energies,         only : energies_type, energies_init

 use m_symtk,            only : symmetrize_xred
 use m_spacepar,         only : setsym
 use m_drivexc,          only : check_kxc
 use m_setvtr,           only : setvtr
 use m_mkrho,            only : initro,mkrho
 use m_hamiltonian,      only : init_hamiltonian,destroy_hamiltonian, &
&                               load_spin_hamiltonian, load_k_hamiltonian,gs_hamiltonian_type
 use m_fock,             only : fock_type
 use m_electronpositron, only : electronpositron_type

 use m_fourier_interpol, only : transgrid
 use m_fft,              only : fftpac,fourdp

 use m_paw_dmft,         only : init_sc_dmft,destroy_sc_dmft,paw_dmft_type
 use m_cgprj,            only : ctocprj
 use m_paw_occupancies,  only : pawmkrhoij
 use m_paw_mkrho,        only : pawmkrho
 use m_pawcprj,          only : pawcprj_type, pawcprj_alloc,pawcprj_getdim

 implicit none

 private 

 public :: subscf_init
 public :: subscf_run
 public :: subscf_destroy

!!****t* m_subscf/subscf_type
!! NAME
!!  subscf_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! SOURCE
 type, public :: subscf_type

   type(datafiles_type),pointer :: dtfil => null()
   type(dataset_type),pointer :: dtset => null()
   type(crystal_t),pointer :: crystal => null()
   type(MPI_type),pointer :: mpi_enreg => null()
   type(pseudopotential_type),pointer :: psps => null()
   type(wvl_data),pointer :: wvl => null()

   type(pawtab_type), pointer :: pawtab(:) => null()
   type(pawfgr_type),pointer :: pawfgr => null()
   type(pawang_type),pointer :: pawang => null()
   type(pawrad_type), pointer :: pawrad(:) => null()
   type(pawrhoij_type),pointer :: pawrhoij(:) => null()

   integer, pointer :: kg(:,:) => null()
   integer, pointer :: npwarr(:) => null()
   real(dp), pointer :: ylm(:,:) => null()
   real(dp), pointer :: ylmgr(:,:,:) => null()

   integer,pointer :: mcprj => null()
   integer,pointer :: nfftf => null()
   integer,pointer :: mcg => null()
   real(dp), pointer :: ecore => null()
   real(dp), pointer :: cg(:,:) => null()
   real(dp), pointer :: occ(:) => null()

   real(dp),allocatable :: rhog(:,:),rhor(:,:)
   real(dp),allocatable :: eig_sub(:)
   complex(dpc),allocatable :: subham_sub(:,:)

   !useless
   type(electronpositron_type),pointer :: electronpositron=>null()
   integer :: pwind_alloc
   integer,allocatable :: pwind(:,:,:)
   real(dp),allocatable :: pwnsfac(:,:) 
   real(dp),allocatable :: taug(:,:),taur(:,:)

 end type subscf_type

contains 

!!****f* m_subscf/subscf_init
!! NAME
!! subscf_init
!!
!! FUNCTION
!! Initialize a subscf object
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! TODO
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine subscf_init(this,dtfil,dtset,psps,crystal,nfftf,pawtab,pawrad,pawang,pawfgr,mpi_enreg,&
&                      ylm,ylmgr,kg,cg,mcg,npwarr,mcprj,ecore,wvl,&
&                      occ)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subscf_init'
!End of the abilint section

 implicit none

 type(subscf_type),intent(inout):: this
 type(datafiles_type),intent(in),target :: dtfil
 type(dataset_type),intent(in),target :: dtset
 type(crystal_t),intent(in),target :: crystal
 type(MPI_type),intent(in),target :: mpi_enreg
 type(pseudopotential_type),intent(in),target :: psps
 type(pawtab_type), intent(in),target :: pawtab(psps%ntypat*psps%usepaw)
 type(pawfgr_type),intent(in),target :: pawfgr
 type(pawrad_type), intent(in),target :: pawrad(psps%ntypat*psps%usepaw)
 type(pawang_type),intent(in),target :: pawang
 type(wvl_data),intent(in),target :: wvl !useless

 real(dp), intent(in),target :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in),target :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in),target :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 integer, intent(in),target :: kg(3,dtset%mpw*dtset%mkmem)
 integer,intent(in),target :: nfftf,mcg,mcprj
 integer, intent(in),target :: npwarr(dtset%nkpt)
 real(dp), intent(in),target :: cg(2,mcg),ecore

 write(std_out,*) "Entering subscf_init"
 !initialize pointers
 this%dtfil=>dtfil
 this%dtset=>dtset
 this%psps=>psps
 this%mpi_enreg=>mpi_enreg
 this%crystal=>crystal

 this%wvl=>wvl

 this%pawtab=>pawtab
 this%pawrad=>pawrad
 this%pawang=>pawang
 this%pawfgr=>pawfgr
 this%npwarr=>npwarr

 this%occ=>occ

 this%nfftf=>nfftf
 this%mcprj=>mcprj
 this%mcg=>mcg
 this%cg=>cg
 this%kg=>kg
 this%ylm=>ylm
 this%ylmgr=>ylmgr

 this%ecore=>ecore

 !useless
 nullify (this%electronpositron)
 this%pwind_alloc = 1
 ABI_ALLOCATE(this%pwind,(this%pwind_alloc,2,3))
 ABI_ALLOCATE(this%pwnsfac,(2,this%pwind_alloc))
 this%pwind(:,:,:)=0
 this%pwnsfac(:,:)=zero

 ABI_DATATYPE_ALLOCATE(this%pawrhoij,(mpi_enreg%my_natom*psps%usepaw))

 ABI_ALLOCATE(this%rhor,(nfftf,dtset%nspden))
 ABI_ALLOCATE(this%rhog,(2,nfftf))

 !useless
 if(dtset%usekden.gt.0) then
   ABI_ALLOCATE(this%taur,(nfftf,dtset%nspden*dtset%usekden))
   ABI_ALLOCATE(this%taug,(2,nfftf*dtset%usekden))
 endif

end subroutine subscf_init



!!****f* m_subscf/subscf_destroy
!! NAME
!! subscf_destroy
!!
!! FUNCTION
!! destroy subscf object
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
subroutine subscf_destroy(this)

 implicit none

 type(subscf_type),intent(inout):: this

 this%dtfil=>null()
 this%dtset=>null()
 this%psps=>null()
 this%mpi_enreg=>null()
 this%crystal=>null()

 this%wvl=>null()

 this%pawtab=>null()
 this%pawrad=>null()
 this%pawang=>null()
 this%pawfgr=>null()
 this%npwarr=>null()

 this%occ=>null()

 this%nfftf=>null()
 this%mcprj=>null()
 this%mcg=>null()
 this%cg=>null()
 this%kg=>null()
 this%ylm=>null()
 this%ylmgr=>null()

 this%ecore=>null()

 this%electronpositron=>null()

 ABI_DEALLOCATE(this%pwind)
 ABI_DEALLOCATE(this%pwnsfac)
 if(this%dtset%usekden.gt.0) then
   ABI_DEALLOCATE(this%taur)
   ABI_DEALLOCATE(this%taug)
 endif

end subroutine subscf_destroy



!!****f* m_subscf/subscf_run
!! NAME
!! subscf_run
!!
!! FUNCTION
!! subscf driver
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! TODO
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine subscf_run(this,can2sub,dimsub)

 implicit none

 type(subscf_type),intent(inout):: this
 integer, intent(in) :: dimsub
 complex(dpc), intent(in) :: can2sub(dimsub,dimsub)

 call subscf_core(this,this%dtset,this%crystal,this%psps,this%pawtab,this%pawrad,this%pawang,this%pawfgr,this%mpi_enreg,&
&  this%nfftf,this%ecore,this%rhog,this%rhor,can2sub,dimsub)


end subroutine subscf_run


!!****f* m_subscf/subscf_core
!! NAME
!! subscf_core
!!
!! FUNCTION
!! main routine for scf calculation in subspace
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! TODO
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine subscf_core(this,dtset,crystal,psps,pawtab,pawrad,pawang,pawfgr,mpi_enreg,&
&                      nfftf,ecore,rhog,rhor,can2sub,dimsub)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subscf_core'
!End of the abilint section

 implicit none

 type(subscf_type),intent(inout) :: this
 type(dataset_type),intent(inout) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg
 type(crystal_t),intent(inout) :: crystal
 type(pseudopotential_type),intent(in) :: psps

 type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(inout) :: pawfgr
 type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)

 real(dp),intent(inout) :: rhog(2,nfftf)
 real(dp),intent(inout) :: rhor(nfftf,dtset%nspden)

 integer, intent(in) :: dimsub
 complex(dpc), intent(in) :: can2sub(dimsub,dimsub)

 integer,intent(inout) :: nfftf
 real(dp),intent(in) :: ecore

 integer :: ii,jj
 integer :: cplex,mgfftf,ipositron,optene,optres,ngrvdw,nkxc
 integer :: has_dijhat,has_dijnd,has_dijfock,has_dijU,has_vhartree
 integer :: iatom,istep,nstep,n1xccc,n3xccc,nzlmopt,istep_mix,option
 integer :: nhatgrdim,ider,idir,izero,ipert,moved_atm_inside,moved_rhor
 integer :: my_natom,usexcnhat,usefock,usecprj,forces_needed,stress_needed,use_hybcomp
 integer :: optcut,optgr0,optgr1,optgr2,optrad
 integer :: mband_cprj,my_nspinor,mcprj,ctocprj_choice,iorder_cprj
 real(dp) :: boxcut,compch_fft,compch_sph,gsqcut,ecut,ecutf,ucvol,ucvol_local
 real(dp) :: zion,vxcavg,hyb_mixing,hyb_mixing_sr

 logical :: dummy_nhatgr
 logical,save :: tfw_activated=.false.
 character(len=1500) :: message

 type(energies_type) :: energies

 !arrays
 integer :: ngfft(18),ngfftf(18),nfftot,nfftotf
 integer,allocatable :: dimcprj(:),dimcprj_srt(:),indsym(:,:,:),symrec(:,:,:),irrzon(:,:,:)
 real(dp),allocatable :: kxc(:,:),nhat(:,:),nhatgr(:,:,:),ph1d(:,:),ph1df(:,:)
 integer,allocatable :: l_size_atm(:)
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),strsxc(6),vpotzero(2)
 real(dp),allocatable :: vhartr(:),vpsp(:),vtrial(:,:),phnons(:,:,:),rhowfg(:,:),rhowfr(:,:)
 real(dp),allocatable :: vxc(:,:),vxc_hybcomp(:,:),vxctau(:,:,:),workr(:,:),xccc3d(:),ylmdiel(:,:)
 real(dp),allocatable :: grchempottn(:,:),grewtn(:,:),grhf(:,:),grnl(:),grvdw(:,:),grxc(:,:)

 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)

 type(fock_type),pointer :: fock

 type(pawcprj_type),pointer :: cprj(:,:)
 type(pawcprj_type),allocatable, target :: cprj_local(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_unsym(:)


 !useless
 integer :: use_sc_dmft
 type(paw_dmft_type) :: paw_dmft
 real(dp) :: qphon(3),rhopsg_dummy(0,0),rhopsr_dummy(0,0),rhor_dummy(0,0)

 write(std_out,*) "Entering subscf_core"

 ABI_ALLOCATE(this%eig_sub,(dimsub))
 ABI_ALLOCATE(this%subham_sub,(dimsub,dimsub))

 ipert=0;idir=0;cplex=1
 istep_mix=1
 hyb_mixing=zero;hyb_mixing_sr=zero
 ipositron = 0
 forces_needed = 0; stress_needed = 0
 moved_atm_inside = 0
 vpotzero(:)=zero

 nstep=dtset%nstep
 my_natom = mpi_enreg%my_natom
 ecut=dtset%ecut
 ecutf=ecut; if (psps%usepaw==1) ecutf=dtset%pawecutdg
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) ecutf=dtset%pawecutdg

 zion=zero
 do iatom=1,dtset%natom
   zion=zion+psps%ziontypat(dtset%typat(iatom))
 end do


 ngfft(:)=dtset%ngfft(:)
 if (psps%usepaw==1) then
   mgfftf=pawfgr%mgfft;ngfftf(:)=pawfgr%ngfft(:)
 else
   mgfftf=dtset%mgfft;ngfftf(:)=ngfft(:)
 end if

 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 nfftotf=product(ngfftf(1:3))

 call metric(gmet,gprimd,-1,rmet,crystal%rprimd,ucvol)
 ucvol_local = ucvol

 usecprj=0
 if (this%mcprj>0) then
  usecprj=1
 end if

 optres=merge(0,1,dtset%iscf<10)

 ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
 ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*dtset%natom))
 ABI_ALLOCATE(vhartr,(nfftf))
 ABI_ALLOCATE(vtrial,(nfftf,dtset%nspden))
 ABI_ALLOCATE(vpsp,(nfftf))
 ABI_ALLOCATE(vxc,(nfftf,dtset%nspden))
 ABI_ALLOCATE(vxctau,(nfftf,dtset%nspden,4*dtset%usekden))

 use_hybcomp=0
 if(mod(dtset%fockoptmix,100)==11)use_hybcomp=1
 ABI_ALLOCATE(vxc_hybcomp,(nfftf,dtset%nspden*use_hybcomp))

 call energies_init(energies)
 select case(dtset%usepotzero)
 case(0,1)
   energies%e_corepsp   = ecore / ucvol
   energies%e_corepspdc = zero
 case(2)
   energies%e_corepsp   = zero
   energies%e_corepspdc = zero
 end select


 usefock=dtset%usefock
 nullify(fock)
 if(usefock==1) MSG_ERROR('usefock == 1, NYI')

 strsxc=zero
 ABI_ALLOCATE(grchempottn,(3,dtset%natom))
 grchempottn(:,:)=zero
 ABI_ALLOCATE(grewtn,(3,dtset%natom))
 ngrvdw=0
 ABI_ALLOCATE(grvdw,(3,ngrvdw))

 nkxc=0
 if (dtset%iscf>0.and.modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79)) nkxc=2*min(dtset%nspden,2)-1
 if (nkxc>0) then
   MSG_ERROR('nkxc>0, NYI')
   call check_kxc(dtset%ixc,dtset%optdriver)
 end if
 ABI_ALLOCATE(kxc,(nfftf,nkxc))

 n1xccc=0;if (psps%n1xccc/=0) n1xccc=psps%n1xccc
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=nfftf
 ABI_ALLOCATE(xccc3d,(n3xccc))


 !symmetry related stuffs
 ABI_ALLOCATE(irrzon,(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(phnons,(2,nfftot**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(indsym,(4,dtset%nsym,dtset%natom))
 ABI_ALLOCATE(symrec,(3,3,dtset%nsym))
 irrzon(:,:,:)=0
 phnons(:,:,:)=zero
 indsym(:,:,:)=0
 symrec(:,:,:)=0
 if (dtset%nsym>1) then
   call setsym(indsym,irrzon,dtset%iscf,dtset%natom,&
&   nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
&   phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,crystal%xred)

!  Make sure dtset%iatfix does not break symmetry
   call fixsym(dtset%iatfix,indsym,dtset%natom,dtset%nsym)
 else
!  The symrec array is used by initberry even in case nsym = 1
   symrec(:,:,1) = 0
   symrec(1,1,1) = 1 ; symrec(2,2,1) = 1 ; symrec(3,3,1) = 1
 end if

!is this the right place?
 call symmetrize_xred(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,crystal%xred)
!Get cut-off for g-vectors
 if (psps%usepaw==1) then
   call wrtout(std_out,' FFT (fine) grid used in SCF cycle:','COLL')
 end if
!Compute large sphere cut-off (gsqcut):
 call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,k0,ngfftf)
!Compute structure factor phases:
 call getph(crystal%atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,crystal%xred)
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
   call getph(crystal%atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,crystal%xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if


 if(psps%usepaw==1) then
   ABI_ALLOCATE(nhat,(this%nfftf,dtset%nspden*psps%usepaw))
   if (nstep==0) nhat=zero
   ABI_DATATYPE_ALLOCATE(pawfgrtab,(my_natom))
   if (my_natom>0) then
     call pawtab_get_lsize(pawtab,l_size_atm,my_natom,dtset%typat,&
&      mpi_atmtab=mpi_enreg%my_atmtab)
     call pawfgrtab_init(pawfgrtab,cplex,l_size_atm,dtset%nspden,dtset%typat,&
&      mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
     ABI_DEALLOCATE(l_size_atm)
   end if
   compch_fft=-1.d5
   usexcnhat=maxval(pawtab(:)%usexcnhat)
   write(std_out,*) "usexcnhat = ",usexcnhat
   if (usexcnhat==0.and.dtset%ionmov==4.and.dtset%iscf<10) then
     message = 'You cannot simultaneously use ionmov=4 and such a PAW psp file!'
     MSG_ERROR(message)
   end if

   ABI_DATATYPE_ALLOCATE(paw_ij,(my_natom))
   ABI_DATATYPE_ALLOCATE(paw_an,(my_natom))
   call paw_an_nullify(paw_an)
   call paw_ij_nullify(paw_ij)
   has_dijhat=0;if (dtset%iscf==22) has_dijhat=1
   has_vhartree=0; if (dtset%prtvha > 0 .or. dtset%prtvclmb > 0) has_vhartree=1
   has_dijfock=0; if (usefock==1) has_dijfock=1
   has_dijnd=0;if(any(abs(dtset%nucdipmom)>tol8)) has_dijnd=1
   has_dijU=0; if (dtset%usepawu==5.or.dtset%usepawu==6) has_dijU=1
   call paw_an_init(paw_an,dtset%natom,dtset%ntypat,0,0,dtset%nspden,&
&   cplex,dtset%pawxcdev,dtset%typat,pawang,pawtab,has_vxc=1,has_vxc_ex=1,has_vhartree=has_vhartree,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call paw_ij_init(paw_ij,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,&
&   dtset%pawspnorb,dtset%natom,dtset%ntypat,dtset%typat,pawtab,&
&   has_dij=1,has_dijfock=has_dijfock,has_dijhartree=1,has_dijnd=has_dijnd,has_dijso=1,has_dijhat=has_dijhat,&
&   has_dijU=has_dijU,has_pawu_occ=1,has_exexch_pot=1,nucdipmom=dtset%nucdipmom,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

   compch_sph=-1.d5

!  Allocation of projected WF (optional)
   if (usecprj==1) then
     MSG_ERROR('NYI')
!     iorder_cprj=0
!     if (usefock==1) then
!       ctocprj_choice = 1
!       if (dtset%optforces == 1) then
!        ctocprj_choice = 2; ! ncpgr = 3 
!       end if
!     end if
   endif

 endif !end usepaw


 !prepare initial density
 write(std_out,*) "debug: dtfil%ireadwf,dtfil%ireadden:"
 write(std_out,*) this%dtfil%ireadwf, this%dtfil%ireadden

 use_sc_dmft=dtset%usedmft
 if(dtset%paral_kgb>0) use_sc_dmft=0
 call init_sc_dmft(dtset%nbandkss,dtset%dmftbandi,dtset%dmftbandf,dtset%dmft_read_occnd,dtset%mband,&
& dtset%nband,dtset%nkpt,dtset%nspden, &
& dtset%nspinor,dtset%nsppol,this%occ,dtset%usedmft,paw_dmft,use_sc_dmft,dtset%dmft_solv,mpi_enreg)

 if(this%dtfil%ireadwf/=0)then
!  Obtain the charge density from wfs that were read previously
!  Be careful: in PAW, rho does not include the compensation
!  density (to be added in scfcv.F90) !
!  tim_mkrho=1 ; mpi_enreg%paralbd=0
   if (psps%usepaw==1) then
     ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
     ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
     call mkrho(this%cg,dtset,gprimd,irrzon,this%kg,this%mcg,&
&      mpi_enreg,this%npwarr,this%occ,paw_dmft,phnons,rhowfg,rhowfr,crystal%rprimd,1,ucvol,this%wvl%den,this%wvl%wfs)
     call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,rhog,rhowfr,rhor)
   else
     call mkrho(this%cg,dtset,gprimd,irrzon,this%kg,this%mcg,&
&      mpi_enreg,this%npwarr,this%occ,paw_dmft,phnons,rhog,rhor,crystal%rprimd,1,ucvol,this%wvl%den,this%wvl%wfs)
   end if
 else if(this%dtfil%ireadwf==0)then
   call initro(crystal%atindx,dtset%densty,gmet,gsqcut,psps%usepaw,&
&     mgfftf,mpi_enreg,psps%mqgrid_vl,dtset%natom,crystal%nattyp,nfftf,&
&     ngfftf,dtset%nspden,psps%ntypat,dtset%paral_kgb,psps,pawtab,ph1df,&
&     psps%qgrid_vl,rhog,rhor,dtset%spinat,ucvol,psps%usepaw,&
&     dtset%ziontypat,dtset%znucl)
 endif

 if (psps%usepaw==1) then
   call initrhoij(dtset%pawcpxocc,dtset%lexexch,&
&   dtset%lpawu,mpi_enreg%my_natom,dtset%natom,dtset%nspden,dtset%nspinor,dtset%nsppol,&
&   dtset%ntypat,this%pawrhoij,dtset%pawspnorb,pawtab,dtset%spinat,dtset%typat,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

   ABI_ALLOCATE(dimcprj,(dtset%natom))
   ABI_ALLOCATE(dimcprj_srt,(dtset%natom))
   call pawcprj_getdim(dimcprj_srt,dtset%natom,crystal%nattyp,dtset%ntypat,dtset%typat,pawtab,'O')
   call pawcprj_getdim(dimcprj,dtset%natom,crystal%nattyp,dtset%ntypat,dtset%typat,pawtab,'R')

   my_nspinor=1
   mband_cprj=dtset%mband/mpi_enreg%nproc_band
   mcprj=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
   ABI_DATATYPE_ALLOCATE(cprj_local,(dtset%natom,mcprj))
   ctocprj_choice = 1
   call pawcprj_alloc(cprj_local,0,dimcprj_srt)
   cprj=> cprj_local
   iatom=0 ; iorder_cprj=0 !0: sorted; 1:un-sorted
   idir = 0
   call ctocprj(crystal%atindx,this%cg,ctocprj_choice,cprj_local,gmet,gprimd,&
&   iatom,idir,iorder_cprj,dtset%istwfk,this%kg,dtset%kptns,&
&   this%mcg,mcprj,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,&
&   dtset%mpw,dtset%natom,crystal%nattyp,dtset%nband,dtset%natom,ngfft,&
&   dtset%nkpt,dtset%nloalg,this%npwarr,dtset%nspinor,dtset%nsppol,&
&   dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&   ucvol,this%dtfil%unpaw,crystal%xred,this%ylm,this%ylmgr)

   pawrhoij_unsym=>this%pawrhoij
 
   call pawmkrhoij(crystal%atindx,crystal%atindx1,cprj,dimcprj_srt,dtset%istwfk,dtset%kptopt,&
&   dtset%mband,mband_cprj,mcprj,dtset%mkmem,mpi_enreg,dtset%natom,dtset%nband,dtset%nkpt,&
&   dtset%nspinor,dtset%nsppol,this%occ,dtset%paral_kgb,paw_dmft,dtset%pawprtvol,pawrhoij_unsym,&
&   this%dtfil%unpaw,dtset%usewvl,dtset%wtk)
 end if

 !scf iteration
 do istep=1,nstep

   if (istep==1) then
       !call symmetrize_xred(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,crystal%xred)

!      Get cut-off for g-vectors
!       if (psps%usepaw==1) then
!         call wrtout(std_out,' FFT (fine) grid used in SCF cycle:','COLL')
!       end if
!       call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,k0,ngfftf)

!      Compute structure factor phases and large sphere cut-off (gsqcut):
!       call getph(crystal%atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,crystal%xred)

!       if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
!         call getph(crystal%atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,crystal%xred)
!       else
!         ph1df(:,:)=ph1d(:,:)
!       end if

     if (psps%usepaw==1) then
!      Check for non-overlapping spheres
       call chkpawovlp(dtset%natom,psps%ntypat,dtset%pawovlp,pawtab,rmet,dtset%typat,crystal%xred)

!      Identify parts of the rectangular grid where the density has to be calculated
       optcut=0;optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
       if (forces_needed==1.or.(dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0)) then
         optgr1=dtset%pawstgylm
         if (stress_needed==1)  optrad=1
         if (dtset%pawprtwf==1) optrad=1
       end if

       call nhatgrid(crystal%atindx1,gmet,my_natom,dtset%natom,&
&        crystal%nattyp,ngfftf,psps%ntypat,optcut,optgr0,optgr1,optgr2,optrad,&
&        pawfgrtab,pawtab,crystal%rprimd,dtset%typat,ucvol,crystal%xred,&
&        comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&        comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)

!      Compute n_tild + n_hat
       qphon(:)=zero
       call pawmkrho(1,compch_fft,cplex,gprimd,idir,indsym,ipert,mpi_enreg,&
&        my_natom,dtset%natom,dtset%nspden,dtset%nsym,dtset%ntypat,dtset%paral_kgb,pawang,pawfgr,pawfgrtab,&
&        dtset%pawprtvol,this%pawrhoij,pawrhoij_unsym,pawtab,qphon,rhowfg,rhowfr,rhor,crystal%rprimd,dtset%symafm,&
&        symrec,dtset%typat,ucvol,dtset%usewvl,crystal%xred,rhog=rhog,pawnhat=nhat)
     end if
     write(std_out,*) "debug: compch_fft = ", compch_fft
     write(std_out,*) "debug: maxval(nhat) = ", maxval(nhat)
     
     moved_rhor = 0

   end if !end istep==1

   if (istep==1) then
!    PAW only: we sometimes have to compute compensation density
!    and eventually add it to density from WFs
     nhatgrdim=0
     dummy_nhatgr = .False.
!    This is executed only in the positron case.
!     if (psps%usepaw==1.and.(dtset%positron>=0.or.ipositron/=1) &
!&     .and.(usexcnhat==0 &
!&     .or.(dtset%xclevel==2.and.(this%dtfil%ireadwf/=0.or.this%dtfil%ireadden/=0)) &
!&     .or.(this%dtfil%ireadwf/=0.and.this%dtfil%ireadden==0))) then
!       nhatgrdim=0;if (dtset%xclevel==2) nhatgrdim=usexcnhat*dtset%pawnhatxc
!       ider=2*nhatgrdim;izero=0
!       if (nhatgrdim>0)   then
!         ABI_ALLOCATE(nhatgr,(cplex*nfftf,dtset%nspden,3*nhatgrdim))
!       else
!         ABI_ALLOCATE(nhatgr,(0,0,0))
!         dummy_nhatgr = .True.
!       end if
!       call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,dtset%natom,&
!&       nfftf,ngfftf,nhatgrdim,dtset%nspden,psps%ntypat,pawang,pawfgrtab,&
!&       nhatgr,nhat,this%pawrhoij,this%pawrhoij,pawtab,k0,crystal%rprimd,ucvol_local,dtset%usewvl,crystal%xred,&
!&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
!&       comm_fft=mpi_enreg%comm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0,&
!&       distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl)
!       if (this%dtfil%ireadwf/=0.and.this%dtfil%ireadden==0) then
!         rhor(:,:)=rhor(:,:)+nhat(:,:)
!         call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
!       end if
!     end if

     optene = 4 * optres
     if(dtset%iscf==-3) optene=4

     if (.not.allocated(nhatgr))  then
       ABI_ALLOCATE(nhatgr,(nfftf,dtset%nspden,3*nhatgrdim))
       dummy_nhatgr = .True.
     end if

     call setvtr(crystal%atindx1,dtset,energies,gmet,gprimd,grchempottn,grewtn,grvdw,gsqcut,&
&     istep,kxc,mgfftf,moved_atm_inside,moved_rhor,mpi_enreg,&
&     crystal%nattyp,nfftf,ngfftf,ngrvdw,nhat,nhatgr,nhatgrdim,nkxc,psps%ntypat,&
&     n1xccc,n3xccc,optene,pawrad,pawtab,ph1df,psps,rhog,rhor,rmet,crystal%rprimd,&
&     strsxc,ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,this%wvl,&
&     xccc3d,crystal%xred,electronpositron=this%electronpositron,&
&     taug=this%taug,taur=this%taur,vxc_hybcomp=vxc_hybcomp,vxctau=vxctau,add_tfw=tfw_activated)

     ! set the zero of the potentials here
     if(dtset%usepotzero==2) then
       vpsp(:) = vpsp(:) + ecore / ( zion * ucvol )
     end if

     if ((nhatgrdim>0.and.nstep>0).or.dummy_nhatgr) then
       ABI_DEALLOCATE(nhatgr)
     end if

   endif


!  ######################################################################
!  The following steps are done at every iteration
!  ----------------------------------------------------------------------
!  PAW: Compute energies and potentials in the augmentation regions (spheres)
!  Compute pseudopotential strengths (Dij quantities)
   if (psps%usepaw==1)then

!    Local exact exch.: impose occ. matrix if required
!     if (dtset%useexexch>0) then
!       call setrhoijpbe0(dtset,initialized0,istep,istep_mix,&
!&       spaceComm,my_natom,dtset%natom,dtset%ntypat,pawrhoij,pawtab,dtset%typat,&
!&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
!     end if

!    Computation of on-site densities/potentials/energies
     nzlmopt=0;if (istep_mix==2.and.dtset%pawnzlm>0) nzlmopt=-1
     if (istep_mix>2) nzlmopt=dtset%pawnzlm
     call paw_an_reset_flags(paw_an) ! Force the recomputation of on-site potentials
     call paw_ij_reset_flags(paw_ij,self_consistent=.true.) ! Force the recomputation of Dij
     option=0;if (dtset%iscf>0.and.dtset%iscf<10.and.nstep>0) option=1
     call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,ipert,dtset%ixc,my_natom,dtset%natom,&
&     dtset%nspden,psps%ntypat,dtset%nucdipmom,nzlmopt,option,paw_an,paw_an,paw_ij,pawang,dtset%pawprtvol,pawrad,&
&     this%pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,ucvol,psps%znuclpsp,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     hyb_mixing=hyb_mixing,hyb_mixing_sr=hyb_mixing_sr,&
&     electronpositron=this%electronpositron,vpotzero=vpotzero)

     write(std_out,*) "debug: compch_sph = ",compch_sph
     !do iatom=1,my_natom
     !   write(std_out,*) this%pawrhoij(iatom)%rhoijp(:,:)
     !enddo

!    Correct the average potential with the calculated constant vpotzero
!    Correct the total energies accordingly
!    vpotzero(1) = -beta/ucvol
!    vpotzero(2) = -1/ucvol sum_ij rho_ij gamma_ij
     write(message,'(a,f14.6,2x,f14.6)') &
&     ' average electrostatic smooth potential [Ha] , [eV]',SUM(vpotzero(:)),SUM(vpotzero(:))*Ha_eV
     call wrtout(std_out,message,'COLL')
     vtrial(:,:)=vtrial(:,:)+SUM(vpotzero(:))
     if(option/=1)then
!      Fix the direct total energy (non-zero only for charged systems)
       energies%e_paw=energies%e_paw-SUM(vpotzero(:))*dtset%charge
!      Fix the double counting total energy accordingly (for both charged AND
!      neutral systems)
       energies%e_pawdc=energies%e_pawdc-SUM(vpotzero(:))*zion+vpotzero(2)*dtset%charge
     end if


     call pawdij(cplex,dtset%enunit,gprimd,ipert,my_natom,dtset%natom,nfftf,nfftotf,&
&     dtset%nspden,psps%ntypat,paw_an,paw_ij,pawang,pawfgrtab,dtset%pawprtvol,&
&     pawrad,this%pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,k0,dtset%spnorbscl,&
&     ucvol_local,dtset%charge,vtrial,vxc,crystal%xred,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     mpi_comm_grid=mpi_enreg%comm_fft,&
&     hyb_mixing=hyb_mixing,hyb_mixing_sr=hyb_mixing_sr,&
&     nucdipmom=dtset%nucdipmom)

!    Symetrize Dij
     call symdij(gprimd,indsym,ipert,my_natom,dtset%natom,dtset%nsym,&
&     psps%ntypat,0,paw_ij,pawang,dtset%pawprtvol,pawtab,crystal%rprimd,dtset%symafm,symrec,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     if (has_dijhat==1) then
       call symdij(gprimd,indsym,ipert,my_natom,dtset%natom,dtset%nsym,&
&       psps%ntypat,1,paw_ij,pawang,dtset%pawprtvol,pawtab,crystal%rprimd,dtset%symafm,symrec,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     end if


   endif

   if(dtset%iscf>=0)then
     write(message, '(a,a,i4)' )ch10,' ITER STEP NUMBER  ',istep
     call wrtout(std_out,message,'COLL')
   end if

   call subscf_vtorho(this,dtset,psps,crystal,mpi_enreg,istep,pawfgr,paw_ij,&
&    this%kg,this%ylm,this%ylmgr,vtrial,energies,ph1d,fock,&
&    my_natom,dtset%natom,psps%ntypat,0,nfftf,gmet,gprimd,rmet,&
&    can2sub,dimsub)


   if (dtset%iscf>=10) then
     !compute total nenergy
     !check convergence
     !mix density
   endif

   !rhotov
   exit
 enddo

 call destroy_sc_dmft(paw_dmft)


end subroutine subscf_core




!!****f* m_subscf/subscf_vtorho
!! NAME
!! subscf_vtorho
!!
!! FUNCTION
!! potential to density
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! TODO
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine subscf_vtorho(this,dtset,psps,crystal,mpi_enreg,istep,pawfgr,paw_ij,&
& kg,ylm,ylmgr,vtrial,energies,ph1d,fock,my_natom,natom,ntypat,optforces,nfftf,&
& gmet,gprimd,rmet,can2sub,dimsub)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subscf_vtorho'
!End of the abilint section

 implicit none

 type(subscf_type),intent(inout) :: this
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(crystal_t),intent(in) :: crystal
 type(MPI_type),intent(in) :: mpi_enreg
 type(energies_type), intent(inout) :: energies
 type(fock_type),pointer, intent(inout) :: fock

 type(pawfgr_type), intent(in) :: pawfgr
 type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)

 integer, intent(in) :: istep,my_natom,natom,ntypat,nfftf,optforces
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 real(dp), intent(inout) :: vtrial(nfftf,dtset%nspden)
 real(dp), intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: gmet(3,3),gprimd(3,3),rmet(3,3)

 integer, intent(in) :: dimsub
 complex(dpc), intent(in) :: can2sub(dimsub,dimsub)

 integer :: usecprj_local,istwf_k
 integer :: isppol,ikg,ilm,nkpg,dimffnl,ider,idir
 integer :: ikpt_loc,ikpt,nkpt1,nband_k,my_ikpt
 integer :: n1,n2,n3,n4,n5,n6,npw_k

 !arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: rhodum(1),kpoint(3)
 real(dp), allocatable :: ylm_k(:,:),kinpw(:),kpg_k(:,:),subham(:)
 type(gs_hamiltonian_type) :: gs_hamk
 real(dp),allocatable :: cgrvtrial(:,:),vlocal(:,:,:,:),ffnl(:,:,:,:),ph3d(:,:,:),zshift(:)

 !useless
 integer :: mcgq,mkgq
 real(dp),allocatable :: cgq(:,:), pwnsfacq(:,:)


 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)

 usecprj_local=0;if (psps%usepaw==1) usecprj_local=1

 call init_hamiltonian(gs_hamk,psps,this%pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,natom,&
& dtset%typat,crystal%xred,dtset%nfft,dtset%mgfft,dtset%ngfft,crystal%rprimd,dtset%nloalg,&
& paw_ij=paw_ij,ph1d=ph1d,usecprj=usecprj_local,electronpositron=this%electronpositron,fock=fock,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
& nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)


 ABI_ALLOCATE(vlocal,(n4,n5,n6,gs_hamk%nvloc))

 nkpt1 = dtset%nkpt
 do isppol=1,dtset%nsppol
   ikpt_loc = 0
   ikg = 0

   ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
   call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
   call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal,2)
   ABI_DEALLOCATE(cgrvtrial)

   call load_spin_hamiltonian(gs_hamk,isppol,vlocal=vlocal,with_nonlocal=.true.)

   ikpt = 0
   do while (ikpt_loc < nkpt1)

     ikpt_loc = ikpt_loc + 1
     ikpt = ikpt_loc
     my_ikpt = mpi_enreg%my_kpttab(ikpt)
     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     npw_k=this%npwarr(ikpt)
     istwf_k=dtset%istwfk(ikpt)

     kpoint(:)=dtset%kptns(:,ikpt)

     ABI_ALLOCATE(zshift,(nband_k))
     zshift(:)=dtset%eshift 

     ABI_ALLOCATE(kg_k,(3,npw_k))
     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Set up remaining of the Hamiltonian

!    Compute (1/2) (2 Pi)**2 (k+G)**2:
     ABI_ALLOCATE(kinpw,(npw_k))
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,kpoint,npw_k,0,0)

!    Compute (k+G) vectors (only if useylm=1)
     nkpg=3*optforces*dtset%nloalg(3)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if ((mpi_enreg%paral_kgb/=1.or.istep<=1).and.nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if

!    Compute nonlocal form factors ffnl at all (k+G):
     ider=0;idir=0;dimffnl=1
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
     if (mpi_enreg%paral_kgb/=1.or.istep<=1) then
       call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&         gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&         npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&         psps%usepaw,psps%useylm,ylm_k,ylmgr)
     end if

!    compute and load nuclear dipole Hamiltonian at current k point
     if(any(abs(gs_hamk%nucdipmom)>0.0)) then
        MSG_ERROR('NYI')
!         if(allocated(nucdipmom_k)) then
!           ABI_DEALLOCATE(nucdipmom_k)
!         end if
!         ABI_ALLOCATE(nucdipmom_k,(npw_k*(npw_k+1)/2))
!         call mknucdipmom_k(gmet,kg_k,kpoint,natom,gs_hamk%nucdipmom,nucdipmom_k,npw_k,rprimd,ucvol,xred)
!         if(allocated(gs_hamk%nucdipmom_k)) then
!            ABI_DEALLOCATE(gs_hamk%nucdipmom_k)
!         end if
!         ABI_ALLOCATE(gs_hamk%nucdipmom_k,(npw_k*(npw_k+1)/2))
!         call load_k_hamiltonian(gs_hamk,nucdipmom_k=nucdipmom_k)
     end if

!    Load k-dependent part in the Hamiltonian datastructure
!       - Compute 3D phase factors
!       - Prepare various tabs in case of band-FFT parallelism
!       - Load k-dependent quantities in the Hamiltonian
     ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))

     call load_k_hamiltonian(gs_hamk,kpt_k=dtset%kptns(:,ikpt),istwf_k=istwf_k,npw_k=npw_k,&
&         kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,&
&         compute_ph3d=(mpi_enreg%paral_kgb/=1.or.istep<=1),&
&         compute_gbound=(mpi_enreg%paral_kgb/=1))

     ! Setup gemm_nonlop
     if (gemm_nonlop_use_gemm) then
         !set the global variable indicating to gemm_nonlop where to get its data from
         gemm_nonlop_ikpt_this_proc_being_treated = my_ikpt
         if (istep <= 1) then
           !Init the arrays
           call make_gemm_nonlop(my_ikpt,gs_hamk%npw_fft_k,gs_hamk%lmnmax, &
&           gs_hamk%ntypat, gs_hamk%indlmn, gs_hamk%nattyp, gs_hamk%istwf_k, gs_hamk%ucvol, gs_hamk%ffnl_k,&
&           gs_hamk%ph3d_k)
         end if
     end if

     !useless
     mcgq=1 ; mkgq=1
     ABI_ALLOCATE(cgq,(2,mcgq))
     ABI_ALLOCATE(pwnsfacq,(2,mkgq))

     call subscf_mkham(this,dtset,mpi_enreg,gs_hamk,isppol,ikpt,nband_k,&
&      cgq,pwnsfacq,mcgq,mkgq,zshift,can2sub,dimsub)

   enddo
 enddo


end subroutine subscf_vtorho


!!****f* m_subscf/subscf_mkham
!! NAME
!! subscf_mkham
!!
!! FUNCTION
!! build subspace hamiltonian
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! TODO
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine subscf_mkham(this,dtset,mpi_enreg,gs_hamk,isppol,ikpt,nband_k,&
&                       cgq,pwnsfacq,mcgq,mkgq,zshift,can2sub,dimsub)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subscf_mkham'
!End of the abilint section

 implicit none

 type(subscf_type),intent(inout):: this
 type(gs_hamiltonian_type), intent(inout) :: gs_hamk
 integer,intent(in) :: isppol,ikpt,nband_k
 real(dp),intent(in) :: zshift(nband_k)

 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in) :: mpi_enreg

 integer, intent(in) :: dimsub
 complex(dpc), intent(in) :: can2sub(dimsub,dimsub)

 integer :: npw_k,mgsc

 integer :: use_subovl = 0
 integer :: berryopt = 0, nline = 0
 integer :: chkexit = 0, quit = 0

 integer :: icg = 0, igsc = 0, my_nspinor = 1, inonsc = 1
 integer :: wfoptalg, ierr, iband,ii,isubh

 real(dp),allocatable :: gsc(:,:),subham(:)
 real(dp),allocatable :: subovl(:),subvnl(:),resid_k(:)
 complex(dpc),allocatable::subham_full(:,:)

 real(dp), allocatable :: rwork(:)
 complex(dpc), allocatable :: zwork(:)
 integer :: lwork,info


 !useless stuffs
 integer, intent(in) :: mcgq,mkgq
 real(dp), intent(in) :: cgq(2,mcgq),pwnsfacq(2,mkgq)
 real(dp) :: dphase_k(3)
 type(efield_type) :: dtefield

 ABI_ALLOCATE(subovl,(nband_k*(nband_k+1)*use_subovl))
 ABI_ALLOCATE(subvnl,(nband_k*(nband_k+1)*(1-gs_hamk%usepaw)))
 ABI_ALLOCATE(resid_k,(nband_k))

 wfoptalg=mod(dtset%wfoptalg,100)

 npw_k=this%npwarr(ikpt)
 mgsc=nband_k*npw_k*my_nspinor*gs_hamk%usepaw
 ABI_STAT_ALLOCATE(gsc,(2,mgsc), ierr)
 ABI_CHECK(ierr==0, "out of memory in gsc")
 gsc=zero

 ABI_ALLOCATE(subham,(nband_k*(nband_k+1)))
 subham(:) = zero

 call cgwf(berryopt,this%cg,cgq,chkexit,dtset%cpus,dphase_k,dtefield,this%dtfil%filnam_ds(1),&
&          gsc,gs_hamk,icg,igsc,ikpt,inonsc,isppol,dtset%mband,this%mcg,mcgq,mgsc,mkgq,&
&          mpi_enreg,dtset%mpw,nband_k,dtset%nbdblock,dtset%nkpt,nline,npw_k,this%npwarr,my_nspinor,&
&          dtset%nsppol,dtset%ortalg,dtset%prtvol,this%pwind,this%pwind_alloc,this%pwnsfac,pwnsfacq,quit,resid_k,&
&          subham,subovl,subvnl,dtset%tolrde,dtset%tolwfr,use_subovl,wfoptalg,zshift)


 ABI_ALLOCATE(subham_full,(nband_k,nband_k))
 isubh = 1
 do iband=1,nband_k
   do ii=1,iband
     if(ii/=iband)then
       subham_full(ii,iband)=cmplx(subham(isubh),subham(isubh+1),kind=dp)
       subham_full(iband,ii)=cmplx(subham(isubh),-subham(isubh+1),kind=dp)
     else
       if(abs(subham(isubh+1)).gt.1.d-8) MSG_ERROR('hamiltonian is not hermitian!') 
       subham_full(ii,iband)=cmplx(subham(isubh),0.0,kind=dp) 
     endif
     isubh=isubh+2
   enddo
 enddo

 call compute_oper_ks2sub(subham_full,this%subham_sub,can2sub,nband_k,dimsub)

 ABI_ALLOCATE(rwork,(3*dimsub-2))
 lwork = 65*dimsub ! Value to optimize speed of the diagonalization
 ABI_ALLOCATE(zwork,(lwork))
 call zheev('v','u',dimsub,this%subham_sub,dimsub,this%eig_sub,zwork,lwork,rwork,info)
 write(std_out,*) "debug: eig_sub"
 write(std_out,*) this%eig_sub(:)


 ABI_DEALLOCATE(subham)
 ABI_DEALLOCATE(subham_full) 
 ABI_DEALLOCATE(rwork)
 ABI_DEALLOCATE(zwork)

end subroutine subscf_mkham


end module m_subscf
