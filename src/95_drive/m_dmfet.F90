!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dmfet
!! NAME
!!  m_dmfet
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


module m_dmfet

 use defs_basis
 use m_abicore
 use m_errors
 use defs_datatypes
 use defs_abitypes

 use m_crystal,          only : crystal_t
 use defs_wvltypes,      only : wvl_data
 use m_plowannier

 use m_cgprj,            only : ctocprj
 use m_kg,               only : getph
 use m_pawang,           only : pawang_type
 use m_pawtab,           only : pawtab_type
 use m_pawfgr,           only : pawfgr_type
 use m_pawcprj,          only : pawcprj_type, pawcprj_getdim, pawcprj_alloc
 use m_pawrad,           only : pawrad_type

 use m_subscf

 implicit none

 private

 public :: dmfet_init
 public :: dmfet_run

 type, public :: dmfet_type

   integer,pointer :: nfftf => null()
   integer,pointer :: mcg => null()
   real(dp), pointer :: cg(:,:) => null()

   type(datafiles_type),pointer :: dtfil => null()
   type(dataset_type),pointer :: dtset => null()
   type(pseudopotential_type),pointer :: psps => null()
   type(crystal_t),pointer :: crystal => null()
   type(MPI_type),pointer :: mpi_enreg => null()

   integer, pointer :: kg(:,:) => null()
   integer, pointer :: npwarr(:) => null()
   type(pawtab_type), pointer :: pawtab(:) => null()
   type(pawfgr_type),pointer :: pawfgr => null()
   type(pawrad_type), pointer :: pawrad(:) => null()
   type(pawang_type),pointer :: pawang => null()

   real(dp), pointer :: ylm(:,:) => null()
   real(dp), pointer :: ylmgr(:,:,:) => null()

   real(dp), pointer :: eigen(:) => null()
   real(dp), pointer :: occ(:) => null()
   real(dp) :: e_fermie,ecore

   complex(dpc),allocatable :: can2sub(:,:)
   real(dp), allocatable ::  occ_wan(:)
   integer :: n_canonical, dim_all, dim_sub, n_frozen

   type(wvl_data),pointer :: wvl => null()

 end type dmfet_type


contains
!!****f* m_dmfet/dmfet_init
!! NAME
!! dmfet_init
!!
!! FUNCTION
!! Initialize a dmfet object
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
subroutine dmfet_init(this,crystal,dtfil,dtset,psps,mpi_enreg,&
& kg,nfftf,pawtab,pawrad,pawang,pawfgr,npwarr,ylm,ylmgr,mcg,cg,eigen,occ,e_fermie,ecore,wvl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_init'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this
 type(datafiles_type),intent(in),target :: dtfil
 type(dataset_type),intent(in),target :: dtset
 type(pseudopotential_type),intent(in),target :: psps
 type(crystal_t),intent(in),target :: crystal
 type(MPI_type),intent(in),target :: mpi_enreg
 type(wvl_data),intent(in),target :: wvl !useless
 integer, intent(in),target :: npwarr(dtset%nkpt)

 type(pawfgr_type),intent(in),target :: pawfgr

 integer,intent(in),target :: nfftf,mcg
 real(dp), intent(in),target :: cg(2,mcg)
 integer, intent(in),target :: kg(3,dtset%mpw*dtset%mkmem)
 type(pawtab_type), intent(in),target :: pawtab(psps%ntypat*psps%usepaw)
 type(pawrad_type), intent(in),target :: pawrad(psps%ntypat*psps%usepaw)
 type(pawang_type), intent(in),target :: pawang
 real(dp), intent(in),target :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in),target :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)

 real(dp), intent(in),target :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in),target :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

 real(dp), intent(in) :: e_fermie,ecore

 logical :: DEBUG=.FALSE.

 DBG_ENTER("COLL")

 this%dtfil=>dtfil
 this%dtset=>dtset
 this%psps=>psps
 this%crystal=>crystal
 this%mpi_enreg=>mpi_enreg
 this%nfftf=>nfftf

 this%kg=>kg
 this%npwarr=>npwarr
 this%mcg=>mcg
 this%cg=>cg

 this%pawfgr=>pawfgr
 this%pawtab=>pawtab
 this%pawrad=>pawrad
 this%pawang=>pawang

 this%ylm=>ylm
 this%ylmgr=>ylmgr

 this%wvl=>wvl

 this%eigen=>eigen
 this%occ=>occ
 this%e_fermie = e_fermie
 this%ecore = ecore

 this%n_canonical = 0
 this%dim_all = 0
 this%dim_sub = 0
 this%n_frozen = 0

 DBG_EXIT("COLL")

end subroutine dmfet_init




!!****f* m_dmfet/dmfet_subspac
!! NAME
!! dmfet_subspac
!!
!! FUNCTION
!! Construct the subspace for a dmfet calculation
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
subroutine dmfet_subspac(this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_subspac'
!End of the abilint section

 type(dmfet_type), intent(inout) :: this

 type(plowannier_type) :: wan
 type(pawcprj_type),allocatable, target :: cprj_local(:,:)
 type(pawcprj_type),pointer :: cprj(:,:)
 integer,allocatable :: dimcprj(:)
 real(dp),allocatable :: ph1d(:,:)


 integer :: ngfft(18)
 integer :: usecprj,mband_cprj,mcprj,my_nspinor
 integer :: iatom,idir,iorder_cprj,ctocprj_choice,ncpgr


 ABI_ALLOCATE(dimcprj,(this%dtset%natom))
 call pawcprj_getdim(dimcprj,this%dtset%natom,this%crystal%nattyp,this%dtset%ntypat,this%dtset%typat,this%pawtab,'R')

 ngfft(:)=this%dtset%ngfft(:)
 ABI_ALLOCATE(ph1d,(2,3*(2*this%dtset%mgfft+1)*this%dtset%natom))
 call getph(this%crystal%atindx,this%dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,this%crystal%xred)

 usecprj=1
 my_nspinor=max(1,this%dtset%nspinor/this%mpi_enreg%nproc_spinor)
 mband_cprj=this%dtset%mband/this%mpi_enreg%nproc_band
 mcprj=my_nspinor*mband_cprj*this%dtset%mkmem*this%dtset%nsppol
 ABI_DATATYPE_ALLOCATE(cprj_local,(this%dtset%natom,mcprj))
 ncpgr = 0 ; ctocprj_choice = 1
 call pawcprj_alloc(cprj_local,ncpgr,dimcprj)
 cprj=> cprj_local
 iatom=0 ; iorder_cprj=1 ! cprj are not ordered
 idir = 0
 call ctocprj(this%crystal%atindx,this%cg,ctocprj_choice,cprj_local,this%crystal%gmet,this%crystal%gprimd,&
&   iatom,idir,iorder_cprj,this%dtset%istwfk,this%kg,this%dtset%kptns,&
&   this%mcg,mcprj,this%dtset%mgfft,this%dtset%mkmem,this%mpi_enreg,this%psps%mpsang,&
&   this%dtset%mpw,this%dtset%natom,this%crystal%nattyp,this%dtset%nband,this%dtset%natom,ngfft,&
&   this%dtset%nkpt,this%dtset%nloalg,this%npwarr,this%dtset%nspinor,this%dtset%nsppol,&
&   this%dtset%ntypat,this%dtset%paral_kgb,ph1d,this%psps,this%crystal%rmet,this%dtset%typat,&
&   this%crystal%ucvol,this%dtfil%unpaw,this%crystal%xred,this%ylm,this%ylmgr)


 call init_plowannier(this%dtset,wan)
 this%dim_all = wan%size_wan
 this%n_canonical = wan%bandf_wan-wan%bandi_wan+1

 call compute_coeff_plowannier(this%crystal,cprj,dimcprj,this%dtset,this%eigen,this%e_fermie,&
&     this%mpi_enreg,this%occ,wan,this%pawtab,this%psps,usecprj,this%dtfil%unpaw,this%pawrad,this%dtfil)

 call dmfet_wan2sub(this,wan)

 call destroy_plowannier(wan)


end subroutine dmfet_subspac


!!****f* m_dmfet/dmfet_wan2sub
!! NAME
!! dmfet_wan2sub
!!
!! FUNCTION
!! Construct subspace orbitals using Wannier orbitals
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
subroutine dmfet_wan2sub(this,wan)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_wan2sub'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this
 type(plowannier_type),intent(inout) :: wan

 type(operwan_type), allocatable :: operwan(:,:,:)
 complex(dpc), allocatable :: operks(:,:,:,:)

 real(dp), allocatable :: onedm_wan(:,:,:,:), loc2sub(:,:)
 !complex(dpc), allocatable :: uuT(:,:,:), uTu(:,:,:)
 complex(dpc), allocatable :: opersub(:,:)
 integer:: dim_sub, n_full

 integer :: nbands,isppol,iband1,ibandc,ikpt
 integer :: ii,jj

 real(dp), allocatable :: eig(:), rwork(:)
 complex(dpc), allocatable :: zwork(:)
 integer :: lwork,info



 nbands = this%n_canonical

 !Inialize an empty Wannier operator
 ABI_DATATYPE_ALLOCATE(operwan,(wan%nkpt,wan%natom_wan,wan%natom_wan))
 call initialize_operwan(wan,operwan)

 !Creation of the KS occupation operator
 ABI_ALLOCATE(operks,(wan%nkpt,nbands,nbands,wan%nsppol))
 operks = czero
 do isppol = 1,wan%nsppol
   do iband1 = 1,nbands
     ibandc = iband1 + wan%bandi_wan - 1
     do ikpt = 1,wan%nkpt
       operks(ikpt,iband1,iband1,isppol) = this%occ(((ikpt-1)*this%dtset%mband+ibandc+(isppol-1)*wan%nkpt*this%dtset%mband))
     end do
   end do
 end do

 !compute the occupation in wannier basis
 do ikpt = 1,wan%nkpt
   call compute_oper_ks2wan(wan,operks,operwan,ikpt)
 end do
 ABI_DEALLOCATE(operks)

 ABI_ALLOCATE(onedm_wan,(wan%size_wan,wan%size_wan,wan%nkpt,wan%nsppol))
 call oper_to_matrix(wan,operwan,onedm_wan)
 ABI_DEALLOCATE(operwan)
!   write(std_out,*)"debug:1pdm_wan"
!   do ii=1,wan%size_wan
!      write(std_out,*) onedm_wan(ii,:,1,1) 
!   enddo

 ABI_ALLOCATE(this%occ_wan,(wan%size_wan))
 ABI_ALLOCATE(loc2sub,(wan%size_wan,wan%size_wan))
 call build_subspace(onedm_wan,loc2sub,this%occ_wan,6,wan%size_wan,dim_sub,n_full)
 ABI_ALLOCATE(this%can2sub,(wan%size_wan,wan%size_wan))
 call canonical_to_sub(wan,loc2sub,this%can2sub,wan%size_wan,1,1)

!   ABI_ALLOCATE(uuT,(wan%size_wan,wan%size_wan,wan%nkpt))
!   ABI_ALLOCATE(uTu,(wan%bandf_wan-wan%bandi_wan+1,wan%bandf_wan-wan%bandi_wan+1,wan%nkpt))
!   call test_unitary(wan,uuT,uTu,1)
!   write(std_out,*)"debug:uuT"
!   do ii=1,wan%size_wan
!      write(std_out,*) uuT(ii,:,1)
!   enddo

!   write(std_out,*)"debug:uTu"
!   do ii=1,wan%bandf_wan-wan%bandi_wan+1
!      write(std_out,*) uTu(ii,:,1)
!   enddo

 ABI_DEALLOCATE(loc2sub)
 ABI_DEALLOCATE(onedm_wan)

 this%dim_sub = dim_sub
 this%n_frozen = n_full


 ! Creation of the KS occupation operator
 ABI_ALLOCATE(operks,(wan%nkpt,nbands,nbands,wan%nsppol))
 operks = czero
 do isppol = 1,wan%nsppol
   do iband1 = 1,nbands
     ibandc = iband1 + wan%bandi_wan - 1
     do ikpt = 1,wan%nkpt
       operks(ikpt,iband1,iband1,isppol) = this%eigen(((ikpt-1)*this%dtset%mband+ibandc+(isppol-1)*wan%nkpt*this%dtset%mband))
     end do
   end do
 end do

 ABI_ALLOCATE(opersub,(wan%size_wan,wan%size_wan))
 call compute_oper_ks2sub(operks(1,:,:,1),opersub,this%can2sub,nbands,wan%size_wan)

 ABI_ALLOCATE(eig,(wan%size_wan))
 ABI_ALLOCATE(rwork,(3*wan%size_wan-2))
 lwork = 65*wan%size_wan ! Value to optimize speed of the diagonalization
 ABI_ALLOCATE(zwork,(lwork))
 call zheev('v','u',wan%size_wan,opersub,wan%size_wan,eig,zwork,lwork,rwork,info)
 write(ab_out,*) "debug: eigen_sub"
 write(ab_out,*) eig(:)

 ABI_DEALLOCATE(operks)
 ABI_DEALLOCATE(eig)
 ABI_DEALLOCATE(rwork)
 ABI_DEALLOCATE(zwork)
 ABI_DEALLOCATE(opersub)

end subroutine dmfet_wan2sub



!!****f* m_dmfet/dmfet_core
!! NAME
!! dmfet_core
!!
!! FUNCTION
!! Main function for performing dmfet calculations using subspace orbitals
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
subroutine dmfet_core(this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_core'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this
 type(subscf_type) :: totscf_args 

 integer :: mcprj

 mcprj = 0  !do not save cproj

 !total scf calc in subspace
 call subscf_init(totscf_args,this%dtfil,this%dtset,this%psps,this%crystal,&
& this%nfftf,this%pawtab,this%pawrad,this%pawang,this%pawfgr,this%mpi_enreg,&
& this%ylm,this%ylmgr,this%kg,this%cg,this%mcg,&
& this%npwarr,mcprj,this%ecore,this%wvl,&
& this%occ,0)

 call subscf_run(totscf_args,this%can2sub,this%dim_all)


end subroutine dmfet_core



!!****f* m_dmfet/dmfet_run
!! NAME
!! dmfet_run
!!
!! FUNCTION
!! Run a dmfet calculation
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
subroutine dmfet_run(this)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmfet_run'
!End of the abilint section

 implicit none

 type(dmfet_type), intent(inout) :: this

 character(len=500) :: message


 write(message,'(2a,i3)') ch10,&
&   '================================================================================='
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 write(message,'(2a,i3)') ch10,&
&   '==                           Welcome to sDMFET module                         =='
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')


 call dmfet_subspac(this)
 call dmfet_core(this)

end subroutine dmfet_run


end module m_dmfet
