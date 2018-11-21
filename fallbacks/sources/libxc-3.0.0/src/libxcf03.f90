!! Copyright (C) 2016 Micael Oliveira
!! All rights reserved.
!!
!! This file is dual-licensed under a GPL and a BSD license
!!
!! GPL License:
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! BSD License:
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!! 1. Redistributions of source code must retain the above copyright
!! notice, this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above
!! copyright notice, this list of conditions and the following
!! disclaimer in the documentation and/or other materials provided
!! with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its
!! contributors may be used to endorse or promote products derived
!! from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
!! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
!! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
!! STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
!! OF THE POSSIBILITY OF SUCH DAMAGE.
!! $Id: libxc_master.F03 12341 2016-04-20 18:50:29Z dstrubbe $
module xc_f03_lib_m
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: &
    ! version
    xc_f03_version, &
    xc_f03_version_string, &
    ! info
    xc_f03_func_info_t, &
    xc_f03_func_info_get_number, &
    xc_f03_func_info_get_kind, &
    xc_f03_func_info_get_name, &
    xc_f03_func_info_get_family, &
    xc_f03_func_info_get_flags, &
    xc_f03_func_info_get_refs, &
    ! func
    xc_f03_func_t, &
    xc_f03_func_init, &
    xc_f03_func_end, &
    xc_f03_func_get_info, &
    xc_f03_functional_get_name, &
    xc_f03_functional_get_number, &
    xc_f03_family_from_id, &
    ! lda
    xc_f03_lda, &
    xc_f03_lda_exc, &
    xc_f03_lda_exc_vxc, &
    xc_f03_lda_vxc, &
    xc_f03_lda_fxc, &
    xc_f03_lda_kxc, &
    xc_f03_lda_x_set_params, &
    xc_f03_lda_x_1d_set_params, &
    xc_f03_lda_c_1d_csc_set_params, &
    xc_f03_lda_c_xalpha_set_params, &
    xc_f03_lda_c_2d_prm_set_params, &
    xc_f03_lda_c_vwn_set_params, &
    xc_f03_lda_xc_ksdt_set_params, &
    ! gga
    xc_f03_gga, &
    xc_f03_gga_exc, &
    xc_f03_gga_exc_vxc, &
    xc_f03_gga_vxc, &
    xc_f03_gga_fxc, &
    xc_f03_gga_kxc, &
    xc_f03_gga_lb_modified, &
    xc_f03_gga_x_b86_set_params, &
    xc_f03_gga_x_b88_set_params, &
    xc_f03_gga_x_pbe_set_params, &
    xc_f03_gga_c_pbe_set_params, &
    xc_f03_gga_x_pw91_set_params, &
    xc_f03_gga_x_pw91_set_params2, &
    xc_f03_gga_x_rpbe_set_params, &
    xc_f03_gga_x_optx_set_params, &
    xc_f03_gga_c_lyp_set_params, &
    xc_f03_gga_lb_set_params, &
    xc_f03_gga_k_tflw_set_params, &
    xc_f03_gga_x_2d_b88_set_params, &
    xc_f03_gga_x_wpbeh_set_params, &
    xc_f03_gga_x_hjs_set_params, &
    xc_f03_gga_x_ityh_set_params, &
    xc_f03_gga_x_sfat_set_params, &
    xc_f03_gga_x_ssb_sw_set_params, &
    xc_f03_gga_x_kt_set_params, &
    xc_f03_gga_x_lambda_set_params, &
    xc_f03_gga_ak13_get_asymptotic, &
    xc_f03_hyb_exx_coef, &
    xc_f03_hyb_cam_coef, &
    xc_f03_nlc_coef, &
    xc_f03_hyb_gga_xc_hse_set_params, &
    xc_f03_hyb_gga_xc_pbeh_set_params, &
    ! mgga
    xc_f03_mgga, &
    xc_f03_mgga_exc, &
    xc_f03_mgga_exc_vxc, &
    xc_f03_mgga_vxc, &
    xc_f03_mgga_fxc, &
    xc_f03_mgga_x_tb09_set_params, &
    xc_f03_mgga_x_tpss_set_params, &
    xc_f03_mgga_c_bc95_set_params, &
    xc_f03_mgga_c_pkzb_set_params
  integer(c_int), parameter, public :: &
    XC_UNPOLARIZED = 1, & ! Spin unpolarized
    XC_POLARIZED = 2 ! Spin polarized
    integer(c_int), parameter, public :: &
    XC_NON_RELATIVISTIC = 0, & ! Functional includes or not relativistic
    XC_RELATIVISTIC = 1 ! corrections. Only available in some functionals.
  ! Kinds
  integer(c_int), parameter, public :: &
    XC_EXCHANGE = 0, &
    XC_CORRELATION = 1, &
    XC_EXCHANGE_CORRELATION = 2, &
    XC_KINETIC = 3
  ! Families of xc functionals
  integer(c_int), parameter, public :: &
    XC_FAMILY_UNKNOWN = -1, &
    XC_FAMILY_NONE = 0, &
    XC_FAMILY_LDA = 1, &
    XC_FAMILY_GGA = 2, &
    XC_FAMILY_MGGA = 4, &
    XC_FAMILY_LCA = 8, &
    XC_FAMILY_OEP = 16, &
    XC_FAMILY_HYB_GGA = 32, &
    XC_FAMILY_HYB_MGGA = 64
  integer(c_int), parameter, public :: &
    XC_FLAGS_HAVE_EXC = 1, &
    XC_FLAGS_HAVE_VXC = 2, &
    XC_FLAGS_HAVE_FXC = 4, &
    XC_FLAGS_HAVE_KXC = 8, &
    XC_FLAGS_HAVE_LXC = 16, &
    XC_FLAGS_1D = 32, &
    XC_FLAGS_2D = 64, &
    XC_FLAGS_3D = 128, &
    XC_FLAGS_HYB_CAM = 256, &
    XC_FLAGS_HYB_CAMY = 512, &
    XC_FLAGS_VV10 = 1024, &
    XC_FLAGS_HYB_LC = 2048, &
    XC_FLAGS_HYB_LCY = 4096, &
    XC_FLAGS_STABLE = 8192, &
    XC_FLAGS_DEVELOPMENT = 16384
  integer(c_int), parameter, public :: &
    XC_TAU_EXPLICIT = 0, &
    XC_TAU_EXPANSION = 1
  ! List of functionals
  integer(c_int), parameter, public :: XC_LDA_X = 1 ! Exchange
  integer(c_int), parameter, public :: XC_LDA_C_WIGNER = 2 ! Wigner parametrization
  integer(c_int), parameter, public :: XC_LDA_C_RPA = 3 ! Random Phase Approximation
  integer(c_int), parameter, public :: XC_LDA_C_HL = 4 ! Hedin & Lundqvist
  integer(c_int), parameter, public :: XC_LDA_C_GL = 5 ! Gunnarson & Lundqvist
  integer(c_int), parameter, public :: XC_LDA_C_XALPHA = 6 ! Slater Xalpha
  integer(c_int), parameter, public :: XC_LDA_C_VWN = 7 ! Vosko, Wilk, & Nusair (5)
  integer(c_int), parameter, public :: XC_LDA_C_VWN_RPA = 8 ! Vosko, Wilk, & Nusair (RPA)
  integer(c_int), parameter, public :: XC_LDA_C_PZ = 9 ! Perdew & Zunger
  integer(c_int), parameter, public :: XC_LDA_C_PZ_MOD = 10 ! Perdew & Zunger (Modified)
  integer(c_int), parameter, public :: XC_LDA_C_OB_PZ = 11 ! Ortiz & Ballone (PZ)
  integer(c_int), parameter, public :: XC_LDA_C_PW = 12 ! Perdew & Wang
  integer(c_int), parameter, public :: XC_LDA_C_PW_MOD = 13 ! Perdew & Wang (Modified)
  integer(c_int), parameter, public :: XC_LDA_C_OB_PW = 14 ! Ortiz & Ballone (PW)
  integer(c_int), parameter, public :: XC_LDA_C_2D_AMGB = 15 ! Attaccalite et al
  integer(c_int), parameter, public :: XC_LDA_C_2D_PRM = 16 ! Pittalis, Rasanen & Marques correlation in 2D
  integer(c_int), parameter, public :: XC_LDA_C_vBH = 17 ! von Barth & Hedin
  integer(c_int), parameter, public :: XC_LDA_C_1D_CSC = 18 ! Casula, Sorella, and Senatore 1D correlation
  integer(c_int), parameter, public :: XC_LDA_X_2D = 19 ! Exchange in 2D
  integer(c_int), parameter, public :: XC_LDA_XC_TETER93 = 20 ! Teter 93 parametrization
  integer(c_int), parameter, public :: XC_LDA_X_1D = 21 ! Exchange in 1D
  integer(c_int), parameter, public :: XC_LDA_C_ML1 = 22 ! Modified LSD (version 1) of Proynov and Salahub
  integer(c_int), parameter, public :: XC_LDA_C_ML2 = 23 ! Modified LSD (version 2) of Proynov and Salahub
  integer(c_int), parameter, public :: XC_LDA_C_GOMBAS = 24 ! Gombas parametrization
  integer(c_int), parameter, public :: XC_LDA_C_PW_RPA = 25 ! Perdew & Wang fit of the RPA
  integer(c_int), parameter, public :: XC_LDA_C_1D_LOOS = 26 ! P-F Loos correlation LDA
  integer(c_int), parameter, public :: XC_LDA_C_RC04 = 27 ! Ragot-Cortona
  integer(c_int), parameter, public :: XC_LDA_C_VWN_1 = 28 ! Vosko, Wilk, & Nusair (1)
  integer(c_int), parameter, public :: XC_LDA_C_VWN_2 = 29 ! Vosko, Wilk, & Nusair (2)
  integer(c_int), parameter, public :: XC_LDA_C_VWN_3 = 30 ! Vosko, Wilk, & Nusair (3)
  integer(c_int), parameter, public :: XC_LDA_C_VWN_4 = 31 ! Vosko, Wilk, & Nusair (4)
  integer(c_int), parameter, public :: XC_LDA_XC_ZLP = 43 ! Zhao, Levy & Parr, Eq. (20)
  integer(c_int), parameter, public :: XC_LDA_K_TF = 50 ! Thomas-Fermi kinetic energy functional
  integer(c_int), parameter, public :: XC_LDA_K_LP = 51 ! Lee and Parr Gaussian ansatz
  integer(c_int), parameter, public :: XC_LDA_XC_KSDT = 259 ! Karasiev et al. parametrization
  integer(c_int), parameter, public :: XC_GGA_X_GAM = 32 ! GAM functional from Minnesota
  integer(c_int), parameter, public :: XC_GGA_C_GAM = 33 ! GAM functional from Minnesota
  integer(c_int), parameter, public :: XC_GGA_X_HCTH_A = 34 ! HCTH-A
  integer(c_int), parameter, public :: XC_GGA_X_EV93 = 35 ! Engel and Vosko
  integer(c_int), parameter, public :: XC_GGA_X_BGCP = 38 ! Burke, Cancio, Gould, and Pittalis
  integer(c_int), parameter, public :: XC_GGA_C_BGCP = 39 ! Burke, Cancio, Gould, and Pittalis
  integer(c_int), parameter, public :: XC_GGA_X_LAMBDA_OC2_N = 40 ! lambda_OC2(N) version of PBE
  integer(c_int), parameter, public :: XC_GGA_X_B86_R = 41 ! Revised Becke 86 Xalpha,beta,gamma (with mod. grad. correction)
  integer(c_int), parameter, public :: XC_GGA_X_LAMBDA_CH_N = 44 ! lambda_CH(N) version of PBE
  integer(c_int), parameter, public :: XC_GGA_X_LAMBDA_LO_N = 45 ! lambda_LO(N) version of PBE
  integer(c_int), parameter, public :: XC_GGA_X_HJS_B88_V2 = 46 ! HJS screened exchange corrected B88 version
  integer(c_int), parameter, public :: XC_GGA_C_Q2D = 47 ! Chiodo et al
  integer(c_int), parameter, public :: XC_GGA_X_Q2D = 48 ! Chiodo et al
  integer(c_int), parameter, public :: XC_GGA_X_PBE_MOL = 49 ! Del Campo, Gazquez, Trickey and Vela (PBE-like)
  integer(c_int), parameter, public :: XC_GGA_K_TFVW = 52 ! Thomas-Fermi plus von Weiszaecker correction
  integer(c_int), parameter, public :: XC_GGA_K_REVAPBEINT = 53 ! interpolated version of REVAPBE
  integer(c_int), parameter, public :: XC_GGA_K_APBEINT = 54 ! interpolated version of APBE
  integer(c_int), parameter, public :: XC_GGA_K_REVAPBE = 55 ! revised APBE
  integer(c_int), parameter, public :: XC_GGA_X_AK13 = 56 ! Armiento & Kuemmel 2013
  integer(c_int), parameter, public :: XC_GGA_K_MEYER = 57 ! Meyer, Wang, and Young
  integer(c_int), parameter, public :: XC_GGA_X_LV_RPW86 = 58 ! Berland and Hyldgaard
  integer(c_int), parameter, public :: XC_GGA_X_PBE_TCA = 59 ! PBE revised by Tognetti et al
  integer(c_int), parameter, public :: XC_GGA_X_PBEINT = 60 ! PBE for hybrid interfaces
  integer(c_int), parameter, public :: XC_GGA_C_ZPBEINT = 61 ! spin-dependent gradient correction to PBEint
  integer(c_int), parameter, public :: XC_GGA_C_PBEINT = 62 ! PBE for hybrid interfaces
  integer(c_int), parameter, public :: XC_GGA_C_ZPBESOL = 63 ! spin-dependent gradient correction to PBEsol
  integer(c_int), parameter, public :: XC_GGA_XC_OPBE_D = 65 ! oPBE_D functional of Goerigk and Grimme
  integer(c_int), parameter, public :: XC_GGA_XC_OPWLYP_D = 66 ! oPWLYP-D functional of Goerigk and Grimme
  integer(c_int), parameter, public :: XC_GGA_XC_OBLYP_D = 67 ! oBLYP-D functional of Goerigk and Grimme
  integer(c_int), parameter, public :: XC_GGA_X_VMT84_GE = 68 ! VMT{8,4} with constraint satisfaction with mu = mu_GE
  integer(c_int), parameter, public :: XC_GGA_X_VMT84_PBE = 69 ! VMT{8,4} with constraint satisfaction with mu = mu_PBE
  integer(c_int), parameter, public :: XC_GGA_X_VMT_GE = 70 ! Vela, Medel, and Trickey with mu = mu_GE
  integer(c_int), parameter, public :: XC_GGA_X_VMT_PBE = 71 ! Vela, Medel, and Trickey with mu = mu_PBE
  integer(c_int), parameter, public :: XC_GGA_C_N12_SX = 79 ! N12-SX functional from Minnesota
  integer(c_int), parameter, public :: XC_GGA_C_N12 = 80 ! N12 functional from Minnesota
  integer(c_int), parameter, public :: XC_GGA_X_N12 = 82 ! N12 functional from Minnesota
  integer(c_int), parameter, public :: XC_GGA_C_REGTPSS = 83 ! Regularized TPSS correlation (ex-VPBE)
  integer(c_int), parameter, public :: XC_GGA_C_OP_XALPHA = 84 ! one-parameter progressive functional (XALPHA version)
  integer(c_int), parameter, public :: XC_GGA_C_OP_G96 = 85 ! one-parameter progressive functional (G96 version)
  integer(c_int), parameter, public :: XC_GGA_C_OP_PBE = 86 ! one-parameter progressive functional (PBE version)
  integer(c_int), parameter, public :: XC_GGA_C_OP_B88 = 87 ! one-parameter progressive functional (B88 version)
  integer(c_int), parameter, public :: XC_GGA_C_FT97 = 88 ! Filatov & Thiel correlation
  integer(c_int), parameter, public :: XC_GGA_C_SPBE = 89 ! PBE correlation to be used with the SSB exchange
  integer(c_int), parameter, public :: XC_GGA_X_SSB_SW = 90 ! Swarta, Sola and Bickelhaupt correction to PBE
  integer(c_int), parameter, public :: XC_GGA_X_SSB = 91 ! Swarta, Sola and Bickelhaupt
  integer(c_int), parameter, public :: XC_GGA_X_SSB_D = 92 ! Swarta, Sola and Bickelhaupt dispersion
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_407P = 93 ! HCTH/407+
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_P76 = 94 ! HCTH p=7/6
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_P14 = 95 ! HCTH p=1/4
  integer(c_int), parameter, public :: XC_GGA_XC_B97_GGA1 = 96 ! Becke 97 GGA-1
  integer(c_int), parameter, public :: XC_GGA_C_HCTH_A = 97 ! HCTH-A
  integer(c_int), parameter, public :: XC_GGA_X_BPCCAC = 98 ! BPCCAC (GRAC for the energy)
  integer(c_int), parameter, public :: XC_GGA_C_REVTCA = 99 ! Tognetti, Cortona, Adamo (revised)
  integer(c_int), parameter, public :: XC_GGA_C_TCA = 100 ! Tognetti, Cortona, Adamo
  integer(c_int), parameter, public :: XC_GGA_X_PBE = 101 ! Perdew, Burke & Ernzerhof exchange
  integer(c_int), parameter, public :: XC_GGA_X_PBE_R = 102 ! Perdew, Burke & Ernzerhof exchange (revised)
  integer(c_int), parameter, public :: XC_GGA_X_B86 = 103 ! Becke 86 Xalpha,beta,gamma
  integer(c_int), parameter, public :: XC_GGA_X_HERMAN = 104 ! Herman et al original GGA
  integer(c_int), parameter, public :: XC_GGA_X_B86_MGC = 105 ! Becke 86 Xalpha,beta,gamma (with mod. grad. correction)
  integer(c_int), parameter, public :: XC_GGA_X_B88 = 106 ! Becke 88
  integer(c_int), parameter, public :: XC_GGA_X_G96 = 107 ! Gill 96
  integer(c_int), parameter, public :: XC_GGA_X_PW86 = 108 ! Perdew & Wang 86
  integer(c_int), parameter, public :: XC_GGA_X_PW91 = 109 ! Perdew & Wang 91
  integer(c_int), parameter, public :: XC_GGA_X_OPTX = 110 ! Handy & Cohen OPTX 01
  integer(c_int), parameter, public :: XC_GGA_X_DK87_R1 = 111 ! dePristo & Kress 87 (version R1)
  integer(c_int), parameter, public :: XC_GGA_X_DK87_R2 = 112 ! dePristo & Kress 87 (version R2)
  integer(c_int), parameter, public :: XC_GGA_X_LG93 = 113 ! Lacks & Gordon 93
  integer(c_int), parameter, public :: XC_GGA_X_FT97_A = 114 ! Filatov & Thiel 97 (version A)
  integer(c_int), parameter, public :: XC_GGA_X_FT97_B = 115 ! Filatov & Thiel 97 (version B)
  integer(c_int), parameter, public :: XC_GGA_X_PBE_SOL = 116 ! Perdew, Burke & Ernzerhof exchange (solids)
  integer(c_int), parameter, public :: XC_GGA_X_RPBE = 117 ! Hammer, Hansen & Norskov (PBE-like)
  integer(c_int), parameter, public :: XC_GGA_X_WC = 118 ! Wu & Cohen
  integer(c_int), parameter, public :: XC_GGA_X_MPW91 = 119 ! Modified form of PW91 by Adamo & Barone
  integer(c_int), parameter, public :: XC_GGA_X_AM05 = 120 ! Armiento & Mattsson 05 exchange
  integer(c_int), parameter, public :: XC_GGA_X_PBEA = 121 ! Madsen (PBE-like)
  integer(c_int), parameter, public :: XC_GGA_X_MPBE = 122 ! Adamo & Barone modification to PBE
  integer(c_int), parameter, public :: XC_GGA_X_XPBE = 123 ! xPBE reparametrization by Xu & Goddard
  integer(c_int), parameter, public :: XC_GGA_X_2D_B86_MGC = 124 ! Becke 86 MGC for 2D systems
  integer(c_int), parameter, public :: XC_GGA_X_BAYESIAN = 125 ! Bayesian best fit for the enhancement factor
  integer(c_int), parameter, public :: XC_GGA_X_PBE_JSJR = 126 ! JSJR reparametrization by Pedroza, Silva & Capelle
  integer(c_int), parameter, public :: XC_GGA_X_2D_B88 = 127 ! Becke 88 in 2D
  integer(c_int), parameter, public :: XC_GGA_X_2D_B86 = 128 ! Becke 86 Xalpha,beta,gamma
  integer(c_int), parameter, public :: XC_GGA_X_2D_PBE = 129 ! Perdew, Burke & Ernzerhof exchange in 2D
  integer(c_int), parameter, public :: XC_GGA_C_PBE = 130 ! Perdew, Burke & Ernzerhof correlation
  integer(c_int), parameter, public :: XC_GGA_C_LYP = 131 ! Lee, Yang & Parr
  integer(c_int), parameter, public :: XC_GGA_C_P86 = 132 ! Perdew 86
  integer(c_int), parameter, public :: XC_GGA_C_PBE_SOL = 133 ! Perdew, Burke & Ernzerhof correlation SOL
  integer(c_int), parameter, public :: XC_GGA_C_PW91 = 134 ! Perdew & Wang 91
  integer(c_int), parameter, public :: XC_GGA_C_AM05 = 135 ! Armiento & Mattsson 05 correlation
  integer(c_int), parameter, public :: XC_GGA_C_XPBE = 136 ! xPBE reparametrization by Xu & Goddard
  integer(c_int), parameter, public :: XC_GGA_C_LM = 137 ! Langreth and Mehl correlation
  integer(c_int), parameter, public :: XC_GGA_C_PBE_JRGX = 138 ! JRGX reparametrization by Pedroza, Silva & Capelle
  integer(c_int), parameter, public :: XC_GGA_X_OPTB88_VDW = 139 ! Becke 88 reoptimized to be used with vdW functional of Dion et al
  integer(c_int), parameter, public :: XC_GGA_X_PBEK1_VDW = 140 ! PBE reparametrization for vdW
  integer(c_int), parameter, public :: XC_GGA_X_OPTPBE_VDW = 141 ! PBE reparametrization for vdW
  integer(c_int), parameter, public :: XC_GGA_X_RGE2 = 142 ! Regularized PBE
  integer(c_int), parameter, public :: XC_GGA_C_RGE2 = 143 ! Regularized PBE
  integer(c_int), parameter, public :: XC_GGA_X_RPW86 = 144 ! refitted Perdew & Wang 86
  integer(c_int), parameter, public :: XC_GGA_X_KT1 = 145 ! Keal and Tozer version 1
  integer(c_int), parameter, public :: XC_GGA_XC_KT2 = 146 ! Keal and Tozer version 2
  integer(c_int), parameter, public :: XC_GGA_C_WL = 147 ! Wilson & Levy
  integer(c_int), parameter, public :: XC_GGA_C_WI = 148 ! Wilson & Ivanov
  integer(c_int), parameter, public :: XC_GGA_X_MB88 = 149 ! Modified Becke 88 for proton transfer
  integer(c_int), parameter, public :: XC_GGA_X_SOGGA = 150 ! Second-order generalized gradient approximation
  integer(c_int), parameter, public :: XC_GGA_X_SOGGA11 = 151 ! Second-order generalized gradient approximation 2011
  integer(c_int), parameter, public :: XC_GGA_C_SOGGA11 = 152 ! Second-order generalized gradient approximation 2011
  integer(c_int), parameter, public :: XC_GGA_C_WI0 = 153 ! Wilson & Ivanov initial version
  integer(c_int), parameter, public :: XC_GGA_XC_TH1 = 154 ! Tozer and Handy v. 1
  integer(c_int), parameter, public :: XC_GGA_XC_TH2 = 155 ! Tozer and Handy v. 2
  integer(c_int), parameter, public :: XC_GGA_XC_TH3 = 156 ! Tozer and Handy v. 3
  integer(c_int), parameter, public :: XC_GGA_XC_TH4 = 157 ! Tozer and Handy v. 4
  integer(c_int), parameter, public :: XC_GGA_X_C09X = 158 ! C09x to be used with the VdW of Rutgers-Chalmers
  integer(c_int), parameter, public :: XC_GGA_C_SOGGA11_X = 159 ! To be used with HYB_GGA_X_SOGGA11_X
  integer(c_int), parameter, public :: XC_GGA_X_LB = 160 ! van Leeuwen & Baerends
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_93 = 161 ! HCTH functional fitted to 93 molecules
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_120 = 162 ! HCTH functional fitted to 120 molecules
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_147 = 163 ! HCTH functional fitted to 147 molecules
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_407 = 164 ! HCTH functional fitted to 407 molecules
  integer(c_int), parameter, public :: XC_GGA_XC_EDF1 = 165 ! Empirical functionals from Adamson, Gill, and Pople
  integer(c_int), parameter, public :: XC_GGA_XC_XLYP = 166 ! XLYP functional
  integer(c_int), parameter, public :: XC_GGA_XC_B97_D = 170 ! Grimme functional to be used with C6 vdW term
  integer(c_int), parameter, public :: XC_GGA_XC_PBE1W = 173 ! Functionals fitted for water
  integer(c_int), parameter, public :: XC_GGA_XC_MPWLYP1W = 174 ! Functionals fitted for water
  integer(c_int), parameter, public :: XC_GGA_XC_PBELYP1W = 175 ! Functionals fitted for water
  integer(c_int), parameter, public :: XC_GGA_X_LBM = 182 ! van Leeuwen & Baerends modified
  integer(c_int), parameter, public :: XC_GGA_X_OL2 = 183 ! Exchange form based on Ou-Yang and Levy v.2
  integer(c_int), parameter, public :: XC_GGA_X_APBE = 184 ! mu fixed from the semiclassical neutral atom
  integer(c_int), parameter, public :: XC_GGA_K_APBE = 185 ! mu fixed from the semiclassical neutral atom
  integer(c_int), parameter, public :: XC_GGA_C_APBE = 186 ! mu fixed from the semiclassical neutral atom
  integer(c_int), parameter, public :: XC_GGA_K_TW1 = 187 ! Tran and Wesolowski set 1 (Table II)
  integer(c_int), parameter, public :: XC_GGA_K_TW2 = 188 ! Tran and Wesolowski set 2 (Table II)
  integer(c_int), parameter, public :: XC_GGA_K_TW3 = 189 ! Tran and Wesolowski set 3 (Table II)
  integer(c_int), parameter, public :: XC_GGA_K_TW4 = 190 ! Tran and Wesolowski set 4 (Table II)
  integer(c_int), parameter, public :: XC_GGA_X_HTBS = 191 ! Haas, Tran, Blaha, and Schwarz
  integer(c_int), parameter, public :: XC_GGA_X_AIRY = 192 ! Constantin et al based on the Airy gas
  integer(c_int), parameter, public :: XC_GGA_X_LAG = 193 ! Local Airy Gas
  integer(c_int), parameter, public :: XC_GGA_XC_MOHLYP = 194 ! Functional for organometallic chemistry
  integer(c_int), parameter, public :: XC_GGA_XC_MOHLYP2 = 195 ! Functional for barrier heights
  integer(c_int), parameter, public :: XC_GGA_XC_TH_FL = 196 ! Tozer and Handy v. FL
  integer(c_int), parameter, public :: XC_GGA_XC_TH_FC = 197 ! Tozer and Handy v. FC
  integer(c_int), parameter, public :: XC_GGA_XC_TH_FCFO = 198 ! Tozer and Handy v. FCFO
  integer(c_int), parameter, public :: XC_GGA_XC_TH_FCO = 199 ! Tozer and Handy v. FCO
  integer(c_int), parameter, public :: XC_GGA_C_OPTC = 200 ! Optimized correlation functional of Cohen and Handy
  integer(c_int), parameter, public :: XC_GGA_C_PBELOC = 246 ! Semilocal dynamical correlation
  integer(c_int), parameter, public :: XC_GGA_XC_VV10 = 255 ! Vydrov and Van Voorhis
  integer(c_int), parameter, public :: XC_GGA_C_PBEFE = 258 ! PBE for formation energies
  integer(c_int), parameter, public :: XC_GGA_C_OP_PW91 = 262 ! one-parameter progressive functional (PW91 version)
  integer(c_int), parameter, public :: XC_GGA_X_PBEFE = 265 ! PBE for formation energies
  integer(c_int), parameter, public :: XC_GGA_X_CAP = 270 ! Correct Asymptotic Potential
  integer(c_int), parameter, public :: XC_GGA_K_VW = 500 ! von Weiszaecker functional
  integer(c_int), parameter, public :: XC_GGA_K_GE2 = 501 ! Second-order gradient expansion (l = 1/9)
  integer(c_int), parameter, public :: XC_GGA_K_GOLDEN = 502 ! TF-lambda-vW form by Golden (l = 13/45)
  integer(c_int), parameter, public :: XC_GGA_K_YT65 = 503 ! TF-lambda-vW form by Yonei and Tomishima (l = 1/5)
  integer(c_int), parameter, public :: XC_GGA_K_BALTIN = 504 ! TF-lambda-vW form by Baltin (l = 5/9)
  integer(c_int), parameter, public :: XC_GGA_K_LIEB = 505 ! TF-lambda-vW form by Lieb (l = 0.185909191)
  integer(c_int), parameter, public :: XC_GGA_K_ABSP1 = 506 ! gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)]
  integer(c_int), parameter, public :: XC_GGA_K_ABSP2 = 507 ! gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)]
  integer(c_int), parameter, public :: XC_GGA_K_GR = 508 ! gamma-TFvW form by Gazquez and Robles
  integer(c_int), parameter, public :: XC_GGA_K_LUDENA = 509 ! gamma-TFvW form by Ludena
  integer(c_int), parameter, public :: XC_GGA_K_GP85 = 510 ! gamma-TFvW form by Ghosh and Parr
  integer(c_int), parameter, public :: XC_GGA_K_PEARSON = 511 ! Pearson
  integer(c_int), parameter, public :: XC_GGA_K_OL1 = 512 ! Ou-Yang and Levy v.1
  integer(c_int), parameter, public :: XC_GGA_K_OL2 = 513 ! Ou-Yang and Levy v.2
  integer(c_int), parameter, public :: XC_GGA_K_FR_B88 = 514 ! Fuentealba & Reyes (B88 version)
  integer(c_int), parameter, public :: XC_GGA_K_FR_PW86 = 515 ! Fuentealba & Reyes (PW86 version)
  integer(c_int), parameter, public :: XC_GGA_K_DK = 516 ! DePristo and Kress
  integer(c_int), parameter, public :: XC_GGA_K_PERDEW = 517 ! Perdew
  integer(c_int), parameter, public :: XC_GGA_K_VSK = 518 ! Vitos, Skriver, and Kollar
  integer(c_int), parameter, public :: XC_GGA_K_VJKS = 519 ! Vitos, Johansson, Kollar, and Skriver
  integer(c_int), parameter, public :: XC_GGA_K_ERNZERHOF = 520 ! Ernzerhof
  integer(c_int), parameter, public :: XC_GGA_K_LC94 = 521 ! Lembarki & Chermette
  integer(c_int), parameter, public :: XC_GGA_K_LLP = 522 ! Lee, Lee & Parr
  integer(c_int), parameter, public :: XC_GGA_K_THAKKAR = 523 ! Thakkar 1992
  integer(c_int), parameter, public :: XC_GGA_X_WPBEH = 524 ! short-range version of the PBE
  integer(c_int), parameter, public :: XC_GGA_X_HJS_PBE = 525 ! HJS screened exchange PBE version
  integer(c_int), parameter, public :: XC_GGA_X_HJS_PBE_SOL = 526 ! HJS screened exchange PBE_SOL version
  integer(c_int), parameter, public :: XC_GGA_X_HJS_B88 = 527 ! HJS screened exchange B88 version
  integer(c_int), parameter, public :: XC_GGA_X_HJS_B97X = 528 ! HJS screened exchange B97x version
  integer(c_int), parameter, public :: XC_GGA_X_ITYH = 529 ! short-range recipe for exchange GGA functionals
  integer(c_int), parameter, public :: XC_GGA_X_SFAT = 530 ! short-range recipe for exchange GGA functionals
  integer(c_int), parameter, public :: XC_HYB_GGA_X_N12_SX = 81 ! N12-SX functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97_1p = 266 ! version of B97 by Cohen and Handy
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3PW91 = 401 ! The original (ACM) hybrid of Becke
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3LYP = 402 ! The (in)famous B3LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3P86 = 403 ! Perdew 86 hybrid similar to B3PW91
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_O3LYP = 404 ! hybrid using the optx functional
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_mPW1K = 405 ! mixture of mPW91 and PW91 optimized for kinetics
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_PBEH = 406 ! aka PBE0 or PBE1PBE
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97 = 407 ! Becke 97
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97_1 = 408 ! Becke 97-1
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97_2 = 410 ! Becke 97-2
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_X3LYP = 411 ! hybrid by Xu and Goddard
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B1WC = 412 ! Becke 1-parameter mixture of WC and PBE
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97_K = 413 ! Boese-Martin for Kinetics
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97_3 = 414 ! Becke 97-3
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MPW3PW = 415 ! mixture with the mPW functional
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B1LYP = 416 ! Becke 1-parameter mixture of B88 and LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B1PW91 = 417 ! Becke 1-parameter mixture of B88 and PW91
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_mPW1PW = 418 ! Becke 1-parameter mixture of mPW91 and PW91
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MPW3LYP = 419 ! mixture of mPW and LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_1a = 420 ! Schmider-Becke 98 parameterization 1a
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_1b = 421 ! Schmider-Becke 98 parameterization 1b
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_1c = 422 ! Schmider-Becke 98 parameterization 1c
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_2a = 423 ! Schmider-Becke 98 parameterization 2a
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_2b = 424 ! Schmider-Becke 98 parameterization 2b
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_2c = 425 ! Schmider-Becke 98 parameterization 2c
  integer(c_int), parameter, public :: XC_HYB_GGA_X_SOGGA11_X = 426 ! Hybrid based on SOGGA11 form
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HSE03 = 427 ! the 2003 version of the screened hybrid HSE
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HSE06 = 428 ! the 2006 version of the screened hybrid HSE
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HJS_PBE = 429 ! HJS hybrid screened exchange PBE version
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HJS_PBE_SOL = 430 ! HJS hybrid screened exchange PBE_SOL version
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HJS_B88 = 431 ! HJS hybrid screened exchange B88 version
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HJS_B97X = 432 ! HJS hybrid screened exchange B97x version
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAM_B3LYP = 433 ! CAM version of B3LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_TUNED_CAM_B3LYP = 434 ! CAM version of B3LYP tuned for excitations
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_BHANDH = 435 ! Becke half-and-half
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_BHANDHLYP = 436 ! Becke half-and-half with B88 exchange
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MB3LYP_RC04 = 437 ! B3LYP with RC04 LDA
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MPWLYP1M = 453 ! MPW with 1 par. for metals/LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_REVB3LYP = 454 ! Revised B3LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAMY_BLYP = 455 ! BLYP with yukawa screening
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_PBE0_13 = 456 ! PBE0-1/3
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3LYPs = 459 ! B3LYP* functional
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_WB97 = 463 ! Chai and Head-Gordon
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_WB97X = 464 ! Chai and Head-Gordon
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LRC_WPBEH = 465 ! Long-range corrected functional by Rorhdanz et al
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_WB97X_V = 466 ! Mardirossian and Head-Gordon
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LCY_PBE = 467 ! PBE with yukawa screening
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LCY_BLYP = 468 ! BLYP with yukawa screening
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LC_VV10 = 469 ! Vydrov and Van Voorhis
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAMY_B3LYP = 470 ! B3LYP with Yukawa screening
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_WB97X_D = 471 ! Chai and Head-Gordon
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HPBEINT = 472 ! hPBEint
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LRC_WPBE = 473 ! Long-range corrected functional by Rorhdanz et al
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3LYP5 = 475 ! B3LYP with VWN functional 5 instead of RPA
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_EDF2 = 476 ! Empirical functional from Lin, George and Gill
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAP0 = 477 ! Correct Asymptotic Potential hybrid
  integer(c_int), parameter, public :: XC_MGGA_C_DLDF = 37 ! Dispersionless Density Functional
  integer(c_int), parameter, public :: XC_MGGA_XC_ZLP = 42 ! Zhao, Levy & Parr, Eq. (21)
  integer(c_int), parameter, public :: XC_MGGA_XC_OTPSS_D = 64 ! oTPSS_D functional of Goerigk and Grimme
  integer(c_int), parameter, public :: XC_MGGA_C_CS = 72 ! Colle and Salvetti
  integer(c_int), parameter, public :: XC_MGGA_C_MN12_SX = 73 ! Worker for MN12-SX functional
  integer(c_int), parameter, public :: XC_MGGA_C_MN12_L = 74 ! MN12-L functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M11_L = 75 ! M11-L functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M11 = 76 ! Worker for M11 functional
  integer(c_int), parameter, public :: XC_MGGA_C_M08_SO = 77 ! Worker for M08-SO functional
  integer(c_int), parameter, public :: XC_MGGA_C_M08_HX = 78 ! Worker for M08-HX functional
  integer(c_int), parameter, public :: XC_MGGA_X_LTA = 201 ! Local tau approximation of Ernzerhof & Scuseria
  integer(c_int), parameter, public :: XC_MGGA_X_TPSS = 202 ! Perdew, Tao, Staroverov & Scuseria exchange
  integer(c_int), parameter, public :: XC_MGGA_X_M06_L = 203 ! M06-Local functional of Minnesota
  integer(c_int), parameter, public :: XC_MGGA_X_GVT4 = 204 ! GVT4 from Van Voorhis and Scuseria
  integer(c_int), parameter, public :: XC_MGGA_X_TAU_HCTH = 205 ! tau-HCTH from Boese and Handy
  integer(c_int), parameter, public :: XC_MGGA_X_BR89 = 206 ! Becke-Roussel 89
  integer(c_int), parameter, public :: XC_MGGA_X_BJ06 = 207 ! Becke & Johnson correction to Becke-Roussel 89
  integer(c_int), parameter, public :: XC_MGGA_X_TB09 = 208 ! Tran & Blaha correction to Becke & Johnson
  integer(c_int), parameter, public :: XC_MGGA_X_RPP09 = 209 ! Rasanen, Pittalis, and Proetto correction to Becke & Johnson
  integer(c_int), parameter, public :: XC_MGGA_X_2D_PRHG07 = 210 ! Pittalis, Rasanen, Helbig, Gross Exchange Functional
  integer(c_int), parameter, public :: XC_MGGA_X_2D_PRHG07_PRP10 = 211 ! PRGH07 with PRP10 correction
  integer(c_int), parameter, public :: XC_MGGA_X_REVTPSS = 212 ! revised Perdew, Tao, Staroverov & Scuseria exchange
  integer(c_int), parameter, public :: XC_MGGA_X_PKZB = 213 ! Perdew, Kurth, Zupan, and Blaha
  integer(c_int), parameter, public :: XC_MGGA_X_M05 = 214 ! Worker for M05 functional
  integer(c_int), parameter, public :: XC_MGGA_X_M05_2X = 215 ! Worker for M05-2X functional
  integer(c_int), parameter, public :: XC_MGGA_X_M06_HF = 216 ! Worker for M06-HF functional
  integer(c_int), parameter, public :: XC_MGGA_X_M06 = 217 ! Worker for M06 functional
  integer(c_int), parameter, public :: XC_MGGA_X_M06_2X = 218 ! Worker for M06-2X functional
  integer(c_int), parameter, public :: XC_MGGA_X_M08_HX = 219 ! Worker for M08-HX functional
  integer(c_int), parameter, public :: XC_MGGA_X_M08_SO = 220 ! Worker for M08-SO functional
  integer(c_int), parameter, public :: XC_MGGA_X_MS0 = 221 ! MS exchange of Sun, Xiao, and Ruzsinszky
  integer(c_int), parameter, public :: XC_MGGA_X_MS1 = 222 ! MS1 exchange of Sun, et al
  integer(c_int), parameter, public :: XC_MGGA_X_MS2 = 223 ! MS2 exchange of Sun, et al
  integer(c_int), parameter, public :: XC_MGGA_X_M11 = 225 ! Worker for M11 functional
  integer(c_int), parameter, public :: XC_MGGA_X_M11_L = 226 ! M11-L functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_X_MN12_L = 227 ! MN12-L functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_CC06 = 229 ! Cancio and Chou 2006
  integer(c_int), parameter, public :: XC_MGGA_X_MK00 = 230 ! Exchange for accurate virtual orbital energies
  integer(c_int), parameter, public :: XC_MGGA_C_TPSS = 231 ! Perdew, Tao, Staroverov & Scuseria correlation
  integer(c_int), parameter, public :: XC_MGGA_C_VSXC = 232 ! VSxc from Van Voorhis and Scuseria (correlation part)
  integer(c_int), parameter, public :: XC_MGGA_C_M06_L = 233 ! M06-Local functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M06_HF = 234 ! Worker for M06-HF functional
  integer(c_int), parameter, public :: XC_MGGA_C_M06 = 235 ! Worker for M06 functional
  integer(c_int), parameter, public :: XC_MGGA_C_M06_2X = 236 ! Worker for M06-2X functional
  integer(c_int), parameter, public :: XC_MGGA_C_M05 = 237 ! Worker for M05 functional
  integer(c_int), parameter, public :: XC_MGGA_C_M05_2X = 238 ! Worker for M05-2X functional
  integer(c_int), parameter, public :: XC_MGGA_C_PKZB = 239 ! Perdew, Kurth, Zupan, and Blaha
  integer(c_int), parameter, public :: XC_MGGA_C_BC95 = 240 ! Becke correlation 95
  integer(c_int), parameter, public :: XC_MGGA_C_REVTPSS = 241 ! revised TPSS correlation
  integer(c_int), parameter, public :: XC_MGGA_XC_TPSSLYP1W = 242 ! Functionals fitted for water
  integer(c_int), parameter, public :: XC_MGGA_X_MK00B = 243 ! Exchange for accurate virtual orbital energies (v. B)
  integer(c_int), parameter, public :: XC_MGGA_X_BLOC = 244 ! functional with balanced localization
  integer(c_int), parameter, public :: XC_MGGA_X_MODTPSS = 245 ! Modified Perdew, Tao, Staroverov & Scuseria exchange
  integer(c_int), parameter, public :: XC_MGGA_C_TPSSLOC = 247 ! Semilocal dynamical correlation
  integer(c_int), parameter, public :: XC_MGGA_X_MBEEF = 249 ! mBEEF exchange
  integer(c_int), parameter, public :: XC_MGGA_X_MBEEFVDW = 250 ! mBEEF-vdW exchange
  integer(c_int), parameter, public :: XC_MGGA_XC_B97M_V = 254 ! Mardirossian and Head-Gordon
  integer(c_int), parameter, public :: XC_MGGA_X_MVS = 257 ! MVS exchange of Sun, Perdew, and Ruzsinszky
  integer(c_int), parameter, public :: XC_MGGA_X_MN15_L = 260 ! MN15-L functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_MN15_L = 261 ! MN15-L functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_X_SCAN = 263 ! SCAN exchange of Sun, Ruzsinszky, and Perdew
  integer(c_int), parameter, public :: XC_MGGA_C_SCAN = 267 ! SCAN correlation
  integer(c_int), parameter, public :: XC_MGGA_C_MN15 = 269 ! MN15 functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_DLDF = 36 ! Dispersionless Density Functional
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_MS2H = 224 ! MS2 hybrid exchange of Sun, et al
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_MN12_SX = 248 ! MN12-SX hybrid functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_SCAN0 = 264 ! SCAN hybrid
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_MN15 = 268 ! MN15 functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_M05 = 438 ! M05 functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_M05_2X = 439 ! M05-2X functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_B88B95 = 440 ! Mixture of B88 with BC95 (B1B95)
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_B86B95 = 441 ! Mixture of B86 with BC95
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_PW86B95 = 442 ! Mixture of PW86 with BC95
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_BB1K = 443 ! Mixture of B88 with BC95 from Zhao and Truhlar
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_M06_HF = 444 ! M06-HF functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_MPW1B95 = 445 ! Mixture of mPW91 with BC95 from Zhao and Truhlar
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_MPWB1K = 446 ! Mixture of mPW91 with BC95 for kinetics
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_X1B95 = 447 ! Mixture of X with BC95
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_XB1K = 448 ! Mixture of X with BC95 for kinetics
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_M06 = 449 ! M06 functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_M06_2X = 450 ! M06-2X functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_PW6B95 = 451 ! Mixture of PW91 with BC95 from Zhao and Truhlar
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_PWB6K = 452 ! Mixture of PW91 with BC95 from Zhao and Truhlar for kinetics
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_TPSSH = 457 ! TPSS hybrid
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_REVTPSSH = 458 ! revTPSS hybrid
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_M08_HX = 460 ! M08-HX functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_M08_SO = 461 ! M08-SO functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_M11 = 462 ! M11 functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_MVSH = 474 ! MVS hybrid
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_WB97M_V = 531 ! Mardirossian and Head-Gordon
  ! These are old names kept for compatibility, and that should disappear soon
  integer(c_int), parameter, public :: &
    XC_GGA_C_VPBE = 83, &
    XC_GGA_XC_LB = 160, &
    XC_GGA_K_ABSR1 = 506, &
    XC_GGA_K_ABSR2 = 507
  !----------------------------------------------------------------
  interface
    subroutine xc_version(major, minor, micro) bind(c)
      import
      integer(c_int), intent(out) :: major, minor, micro
    end subroutine xc_version
    type(c_ptr) function xc_version_string() bind(c)
      import
    end function xc_version_string
  end interface
  !----------------------------------------------------------------
  type :: xc_f03_func_info_t
    private
    type(c_ptr) :: ptr = C_NULL_PTR
  end type xc_f03_func_info_t
  interface
    integer(c_int) function xc_func_info_get_number(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_number
    integer(c_int) function xc_func_info_get_kind(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_kind
    type(c_ptr) function xc_func_info_get_name(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_name
    integer(c_int) function xc_func_info_get_family(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_family
    integer(c_int) function xc_func_info_get_flags(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_flags
    type(c_ptr) function xc_func_info_get_ref(info, number) bind(c)
      import
      type(c_ptr), value :: info
      integer(c_int), value :: number
    end function xc_func_info_get_ref
  end interface
  !----------------------------------------------------------------
  type :: xc_f03_func_t
    private
    type(c_ptr) :: ptr = C_NULL_PTR
  end type xc_f03_func_t
  interface
    type(c_ptr) function xc_func_alloc() bind(c)
      import
    end function xc_func_alloc
    integer(c_int) function xc_func_init(p, functional, nspin) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: functional, nspin
    end function xc_func_init
    subroutine xc_func_end(p) bind(c)
      import
      type(c_ptr), value :: p
    end subroutine xc_func_end
    subroutine xc_func_free(p) bind(c)
      import
      type(c_ptr), value :: p
    end subroutine xc_func_free
    type(c_ptr) function xc_func_get_info(p) bind(c)
      import
      type(c_ptr), value :: p
    end function xc_func_get_info
    type(c_ptr) function xc_functional_get_name(number) bind(c)
      import
      integer(c_int), value :: number
    end function xc_functional_get_name
    integer(c_int) function xc_functional_get_number(func_string) bind(c)
      import
      character(kind=c_char) :: func_string(*)
    end function xc_functional_get_number
    integer(c_int) function xc_family_from_id(id, family, number) bind(c)
      import
      integer(c_int), value :: id
      type(c_ptr), value :: family, number
    end function xc_family_from_id
  end interface
  ! LDAs
  !----------------------------------------------------------------
  interface
    subroutine xc_lda(p, np, rho, zk, vrho, v2rho2, v3rho3) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*)
      real(c_double), intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*)
    end subroutine xc_lda
    subroutine xc_lda_exc(p, np, rho, zk) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*)
      real(c_double), intent(out) :: zk(*)
    end subroutine xc_lda_exc
    subroutine xc_lda_exc_vxc(p, np, rho, zk, vrho) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*)
      real(c_double), intent(out) :: zk(*), vrho(*)
    end subroutine xc_lda_exc_vxc
    subroutine xc_lda_vxc(p, np, rho, vrho) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*)
      real(c_double), intent(out) :: vrho(*)
    end subroutine xc_lda_vxc
    subroutine xc_lda_fxc(p, np, rho, v2rho2) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*)
      real(c_double), intent(out) :: v2rho2(*)
    end subroutine xc_lda_fxc
    subroutine xc_lda_kxc(p, np, rho, v3rho3) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*)
      real(c_double), intent(out) :: v3rho3(*)
    end subroutine xc_lda_kxc
  end interface
  interface
    subroutine xc_lda_x_set_params(p, alpha, relativistic, omega) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: alpha, omega
      integer(c_int), value :: relativistic
    end subroutine xc_lda_x_set_params
    subroutine xc_lda_x_1d_set_params(p, interaction, bb) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: interaction
      real(c_double), value :: bb
    end subroutine xc_lda_x_1d_set_params
    subroutine xc_lda_c_1d_csc_set_params(p, interaction, bb) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: interaction
      real(c_double), value :: bb
    end subroutine xc_lda_c_1d_csc_set_params
    subroutine xc_lda_c_xalpha_set_params(p, alpha) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: alpha
    end subroutine xc_lda_c_xalpha_set_params
    subroutine xc_lda_c_2d_prm_set_params(p, N) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: N
    end subroutine xc_lda_c_2d_prm_set_params
    subroutine xc_lda_c_vwn_set_params(p, spin_interpolation) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: spin_interpolation
    end subroutine xc_lda_c_vwn_set_params
    subroutine xc_lda_xc_ksdt_set_params(p, t) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: t
    end subroutine xc_lda_xc_ksdt_set_params
  end interface
  ! GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_gga(p, np, rho, sigma, zk, vrho, vsigma, &
      v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*)
      real(c_double), intent(out) :: zk(*), vrho(*), vsigma(*)
      real(c_double), intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
      real(c_double), intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    end subroutine xc_gga
    subroutine xc_gga_exc(p, np, rho, sigma, zk) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*)
      real(c_double), intent(out) :: zk(*)
    end subroutine xc_gga_exc
    subroutine xc_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*)
      real(c_double), intent(out) :: zk(*), vrho(*), vsigma(*)
    end subroutine xc_gga_exc_vxc
    subroutine xc_gga_vxc(p, np, rho, sigma, vrho, vsigma) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*)
      real(c_double), intent(out) :: vrho(*), vsigma(*)
    end subroutine xc_gga_vxc
    subroutine xc_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*)
      real(c_double), intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    end subroutine xc_gga_fxc
    subroutine xc_gga_kxc(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*)
      real(c_double), intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    end subroutine xc_gga_kxc
  end interface
  interface
    subroutine xc_gga_lb_modified(p, np, rho, grho, r, dedd) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), grho(*)
      real(c_double), value :: r
      real(c_double), intent(out) :: dedd(*)
    end subroutine xc_gga_lb_modified
    subroutine xc_gga_x_b86_set_params(p, beta, gamma, omega) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: beta, gamma, omega
    end subroutine xc_gga_x_b86_set_params
    subroutine xc_gga_x_b88_set_params(p, beta, gamma) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: beta, gamma
    end subroutine xc_gga_x_b88_set_params
    subroutine xc_gga_x_pbe_set_params(p, kappa, mu) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: kappa, mu
    end subroutine xc_gga_x_pbe_set_params
    subroutine xc_gga_c_pbe_set_params(p, beta) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: beta
    end subroutine xc_gga_c_pbe_set_params
    subroutine xc_gga_x_pw91_set_params(p, a, b, c, d, f, alpha, expo) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: a, b, c, d, f, alpha, expo
    end subroutine xc_gga_x_pw91_set_params
    subroutine xc_gga_x_pw91_set_params2(p, bt, alpha, expo) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: bt, alpha, expo
    end subroutine xc_gga_x_pw91_set_params2
    subroutine xc_gga_x_rpbe_set_params(p, kappa, mu) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: kappa, mu
    end subroutine xc_gga_x_rpbe_set_params
    subroutine xc_gga_x_optx_set_params(p, a, b, gamma) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: a, b, gamma
    end subroutine xc_gga_x_optx_set_params
    subroutine xc_gga_c_lyp_set_params(p, a, b, c, d) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: a, b, c, d
    end subroutine xc_gga_c_lyp_set_params
    subroutine xc_gga_lb_set_params(p, modified, threshold, ip, qtot) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: modified
      real(c_double), value :: threshold, ip, qtot
    end subroutine xc_gga_lb_set_params
    subroutine xc_gga_k_tflw_set_params(p, gamma, lambda, n) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: gamma, lambda
      integer(c_int), value :: n
    end subroutine xc_gga_k_tflw_set_params
    subroutine xc_gga_x_2d_b88_set_params(p, beta) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: beta
    end subroutine xc_gga_x_2d_b88_set_params
    subroutine xc_gga_x_wpbeh_set_params(p, omega) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: omega
    end subroutine xc_gga_x_wpbeh_set_params
    subroutine xc_gga_x_hjs_set_params(p, omega) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: omega
    end subroutine xc_gga_x_hjs_set_params
    subroutine xc_gga_x_ityh_set_params(p, func_id, omega) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: func_id
      real(c_double), value :: omega
    end subroutine xc_gga_x_ityh_set_params
    subroutine xc_gga_x_sfat_set_params(p, func_id, omega) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: func_id
      real(c_double), value :: omega
    end subroutine xc_gga_x_sfat_set_params
    subroutine xc_gga_x_ssb_sw_set_params(p, a, b, c, d, e) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: a, b, c, d, e
    end subroutine xc_gga_x_ssb_sw_set_params
    subroutine xc_gga_x_kt_set_params(p, gamma, delta) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: gamma, delta
    end subroutine xc_gga_x_kt_set_params
    subroutine xc_gga_x_lambda_set_params(p, n) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: n
    end subroutine xc_gga_x_lambda_set_params
    real(c_double) function xc_gga_ak13_get_asymptotic(homo) bind(c)
      import
      real(c_double), value :: homo
    end function xc_gga_ak13_get_asymptotic
  end interface
  interface
    real(c_double) function xc_hyb_exx_coef(p) bind(c)
      import
      type(c_ptr), value :: p
    end function xc_hyb_exx_coef
    subroutine xc_hyb_cam_coef(p, omega, alpha, beta) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), intent(out) :: omega, alpha, beta
    end subroutine xc_hyb_cam_coef
    subroutine xc_nlc_coef(p, nlc_b, nlc_c) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), intent(out) :: nlc_b, nlc_c
    end subroutine xc_nlc_coef
    subroutine xc_hyb_gga_xc_hse_set_params(p, alpha, omega) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: alpha, omega
    end subroutine xc_hyb_gga_xc_hse_set_params
    subroutine xc_hyb_gga_xc_pbeh_set_params(p, alpha) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: alpha
    end subroutine xc_hyb_gga_xc_pbeh_set_params
  end interface
  ! the meta-GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double), intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
      real(c_double), intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), v2rholapl(*), &
                                     v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)
    end subroutine xc_mgga
    subroutine xc_mgga_exc(p, np, rho, sigma, lapl, tau, zk) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double), intent(out) :: zk(*)
    end subroutine xc_mgga_exc
    subroutine xc_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double), intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    end subroutine xc_mgga_exc_vxc
    subroutine xc_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double), intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
    end subroutine xc_mgga_vxc
    subroutine xc_mgga_fxc(p, np, rho, sigma, lapl, tau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau) bind(c)
      import
      type(c_ptr), value :: p
      integer(c_int), value :: np
      real(c_double), intent(in) :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double), intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), v2rholapl(*), &
                                     v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)
    end subroutine xc_mgga_fxc
  end interface
  interface
    subroutine xc_mgga_x_tb09_set_params(p, c) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: c
    end subroutine xc_mgga_x_tb09_set_params
    subroutine xc_mgga_x_tpss_set_params(p, b, c, e, kappa, mu) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: b, c, e, kappa, mu
    end subroutine xc_mgga_x_tpss_set_params
    subroutine xc_mgga_c_bc95_set_params(p, css, copp) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: css, copp
    end subroutine xc_mgga_c_bc95_set_params
    subroutine xc_mgga_c_pkzb_set_params(p, beta, d, c0_0, c0_1, c0_2, c0_3) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: beta, d, c0_0, c0_1, c0_2, c0_3
    end subroutine xc_mgga_c_pkzb_set_params
  end interface
contains
  !----------------------------------------------------------------
  subroutine xc_f03_version(major, minor, micro)
    integer, intent(out) :: major, minor, micro
    call xc_version(major, minor, micro)
  end subroutine xc_f03_version
  subroutine xc_f03_version_string(version)
    character(len=*), intent(out) :: version
    type(c_ptr) :: c_version
    c_version = xc_version_string()
    call c_to_f_string_ptr(c_version, version)
  end subroutine xc_f03_version_string
  !----------------------------------------------------------------
  integer function xc_f03_func_info_get_number(info) result(number)
    type(xc_f03_func_info_t), intent(in) :: info
    number = xc_func_info_get_number(info%ptr)
  end function xc_f03_func_info_get_number
  integer function xc_f03_func_info_get_kind(info) result(kind)
    type(xc_f03_func_info_t), intent(in) :: info
    kind = xc_func_info_get_kind(info%ptr)
  end function xc_f03_func_info_get_kind
  character(len=128) function xc_f03_func_info_get_name(info) result(name)
    type(xc_f03_func_info_t), intent(in) :: info
    call c_to_f_string_ptr(xc_func_info_get_name(info%ptr), name)
  end function xc_f03_func_info_get_name
  integer function xc_f03_func_info_get_family(info) result(family)
    type(xc_f03_func_info_t), intent(in) :: info
    family = xc_func_info_get_family(info%ptr)
  end function xc_f03_func_info_get_family
  integer function xc_f03_func_info_get_flags(info) result(flags)
    type(xc_f03_func_info_t), intent(in) :: info
    flags = xc_func_info_get_flags(info%ptr)
  end function xc_f03_func_info_get_flags
  character(len=120) function xc_f03_func_info_get_refs(info, number) result(ref)
    type(xc_f03_func_info_t), intent(in) :: info
    integer, intent(inout) :: number ! number of the reference. Must be 0 in the first call
    type(c_ptr) :: c_ref
    c_ref = xc_func_info_get_ref(info%ptr, number)
    if (c_associated(c_ref)) then
      call c_to_f_string_ptr(c_ref, ref)
      number = number + 1
      if (.not. c_associated(xc_func_info_get_ref(info%ptr, number))) number = -1
    end if
  end function xc_f03_func_info_get_refs
  !----------------------------------------------------------------
  subroutine xc_f03_func_init(p, functional, nspin)
    type(xc_f03_func_t), intent(inout) :: p
    integer, intent(in) :: functional
    integer, intent(in) :: nspin
    integer :: ierr
    p%ptr = xc_func_alloc()
    ierr = xc_func_init(p%ptr, functional, nspin)
  end subroutine xc_f03_func_init
  subroutine xc_f03_func_end(p)
    type(xc_f03_func_t), intent(inout) :: p
    call xc_func_end(p%ptr)
    call xc_func_free(p%ptr)
  end subroutine xc_f03_func_end
  type(xc_f03_func_info_t) function xc_f03_func_get_info(p) result(info)
    type(xc_f03_func_t), intent(in) :: p
    info%ptr = xc_func_get_info(p%ptr)
  end function xc_f03_func_get_info
  character(len=128) function xc_f03_functional_get_name(number) result(name)
    integer, intent(in) :: number
    call c_to_f_string_ptr(xc_functional_get_name(number), name)
  end function xc_f03_functional_get_name
  integer function xc_f03_functional_get_number(func_string) result(number)
    character(len=*), intent(in) :: func_string
    number = xc_functional_get_number(f_to_c_string(func_string))
  end function xc_f03_functional_get_number
  integer function xc_f03_family_from_id(id, family, number)
    integer, intent(in) :: id
    integer, intent(out), optional, target :: family, number
    type(c_ptr) c_family, c_number
    integer, pointer :: f_family, f_number
    if (present(family)) then
      f_family => family
      call c_f_pointer(c_family, f_family)
    else
      c_family = C_NULL_PTR
    end if
    if (present(number)) then
      f_number => number
      call c_f_pointer(c_number, f_number)
    else
      c_number = C_NULL_PTR
    end if
    xc_f03_family_from_id = xc_family_from_id(id, c_family, c_number)
  end function xc_f03_family_from_id
  ! LDAs
  !----------------------------------------------------------------
  subroutine xc_f03_lda(p, np, rho, zk, vrho, v2rho2, v3rho3)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*)
    real(c_double), intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*)
    call xc_lda(p%ptr, np, rho, zk, vrho, v2rho2, v3rho3)
  end subroutine xc_f03_lda
  subroutine xc_f03_lda_exc(p, np, rho, zk)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*)
    real(c_double), intent(out) :: zk(*)
    call xc_lda_exc(p%ptr, np, rho, zk)
  end subroutine xc_f03_lda_exc
  subroutine xc_f03_lda_exc_vxc(p, np, rho, zk, vrho)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*)
    real(c_double), intent(out) :: zk(*), vrho(*)
    call xc_lda_exc_vxc(p%ptr, np, rho, zk, vrho)
  end subroutine xc_f03_lda_exc_vxc
  subroutine xc_f03_lda_vxc(p, np, rho, vrho)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*)
    real(c_double), intent(out) :: vrho(*)
    call xc_lda_vxc(p%ptr, np, rho, vrho)
  end subroutine xc_f03_lda_vxc
  subroutine xc_f03_lda_fxc(p, np, rho, v2rho2)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*)
    real(c_double), intent(out) :: v2rho2(*)
    call xc_lda_fxc(p%ptr, np, rho, v2rho2)
  end subroutine xc_f03_lda_fxc
  subroutine xc_f03_lda_kxc(p, np, rho, v3rho3)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*)
    real(c_double), intent(out) :: v3rho3(*)
    call xc_lda_kxc(p%ptr, np, rho, v3rho3)
  end subroutine xc_f03_lda_kxc
  subroutine xc_f03_lda_x_set_params(p, alpha, relativistic, omega)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: alpha, omega
    integer, intent(in) :: relativistic
    call xc_lda_x_set_params(p%ptr, alpha, relativistic, omega)
  end subroutine xc_f03_lda_x_set_params
  subroutine xc_f03_lda_x_1d_set_params(p, interaction, bb)
    type(xc_f03_func_t), intent(inout) :: p
    integer, intent(in) :: interaction
    real(c_double), intent(in) :: bb
    call xc_lda_x_1d_set_params(p%ptr, interaction, bb)
  end subroutine xc_f03_lda_x_1d_set_params
  subroutine xc_f03_lda_c_1d_csc_set_params(p, interaction, bb)
    type(xc_f03_func_t), intent(inout) :: p
    integer, intent(in) :: interaction
    real(c_double), intent(in) :: bb
    call xc_lda_c_1d_csc_set_params(p%ptr, interaction, bb)
  end subroutine xc_f03_lda_c_1d_csc_set_params
  subroutine xc_f03_lda_c_xalpha_set_params(p, alpha)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: alpha
    call xc_lda_c_xalpha_set_params(p%ptr, alpha)
  end subroutine xc_f03_lda_c_xalpha_set_params
  subroutine xc_f03_lda_c_2d_prm_set_params(p, n)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: n
    call xc_lda_c_2d_prm_set_params(p%ptr, n)
  end subroutine xc_f03_lda_c_2d_prm_set_params
  subroutine xc_f03_lda_c_vwn_set_params(p, spin_interpolation)
    type(xc_f03_func_t), intent(inout) :: p
    integer, intent(in) :: spin_interpolation
    call xc_lda_c_vwn_set_params(p%ptr, spin_interpolation)
  end subroutine xc_f03_lda_c_vwn_set_params
  subroutine xc_f03_lda_xc_ksdt_set_params(p, t)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: t
    call xc_lda_xc_ksdt_set_params(p%ptr, t)
  end subroutine xc_f03_lda_xc_ksdt_set_params
  ! GGAs
  !----------------------------------------------------------------
  subroutine xc_f03_gga(p, np, rho, sigma, zk, vrho, vsigma, &
    v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), sigma(*)
    real(c_double), intent(out) :: zk(*), vrho(*), vsigma(*)
    real(c_double), intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    real(c_double), intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    call xc_gga(p%ptr, np, rho, sigma, zk, vrho, vsigma, &
      v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
  end subroutine xc_f03_gga
  subroutine xc_f03_gga_exc(p, np, rho, sigma, zk)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), sigma(*)
    real(c_double), intent(out) :: zk(*)
    call xc_gga_exc(p%ptr, np, rho, sigma, zk)
  end subroutine xc_f03_gga_exc
  subroutine xc_f03_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), sigma(*)
    real(c_double), intent(out) :: zk(*), vrho(*), vsigma(*)
    call xc_gga_exc_vxc(p%ptr, np, rho, sigma, zk, vrho, vsigma)
  end subroutine xc_f03_gga_exc_vxc
  subroutine xc_f03_gga_vxc(p, np, rho, sigma, vrho, vsigma)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), sigma(*)
    real(c_double), intent(out) :: vrho(*), vsigma(*)
    call xc_gga_vxc(p%ptr, np, rho, sigma, vrho, vsigma)
  end subroutine xc_f03_gga_vxc
  subroutine xc_f03_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), sigma(*)
    real(c_double), intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    call xc_gga_fxc(p%ptr, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
  end subroutine xc_f03_gga_fxc
  subroutine xc_f03_gga_kxc(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), sigma(*)
    real(c_double), intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    call xc_gga_kxc(p%ptr, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
  end subroutine xc_f03_gga_kxc
  subroutine xc_f03_gga_lb_modified(p, np, rho, grho, r, dedd)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), grho(*)
    real(c_double), intent(in) :: r
    real(c_double), intent(out) :: dedd(*)
    call xc_gga_lb_modified(p%ptr, np, rho, grho, r, dedd)
  end subroutine xc_f03_gga_lb_modified
  subroutine xc_f03_gga_x_b86_set_params(p, beta, gamma, omega)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: beta, gamma, omega
    call xc_gga_x_b86_set_params(p%ptr, beta, gamma, omega)
  end subroutine xc_f03_gga_x_b86_set_params
  subroutine xc_f03_gga_x_b88_set_params(p, beta, gamma)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: beta, gamma
    call xc_gga_x_b88_set_params(p%ptr, beta, gamma)
  end subroutine xc_f03_gga_x_b88_set_params
  subroutine xc_f03_gga_x_pbe_set_params(p, kappa, mu)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: kappa, mu
    call xc_gga_x_pbe_set_params(p%ptr, kappa, mu)
  end subroutine xc_f03_gga_x_pbe_set_params
  subroutine xc_f03_gga_c_pbe_set_params(p, beta)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: beta
    call xc_gga_c_pbe_set_params(p%ptr, beta)
  end subroutine xc_f03_gga_c_pbe_set_params
  subroutine xc_f03_gga_x_pw91_set_params(p, a, b, c, d, f, alpha, expo)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: a, b, c, d, f, alpha, expo
    call xc_gga_x_pw91_set_params(p%ptr, a, b, c, d, f, alpha, expo)
  end subroutine xc_f03_gga_x_pw91_set_params
  subroutine xc_f03_gga_x_pw91_set_params2(p, bt, alpha, expo)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: bt, alpha, expo
    call xc_gga_x_pw91_set_params2(p%ptr, bt, alpha, expo)
  end subroutine xc_f03_gga_x_pw91_set_params2
  subroutine xc_f03_gga_x_rpbe_set_params(p, kappa, mu)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: kappa, mu
    call xc_gga_x_rpbe_set_params(p%ptr, kappa, mu)
  end subroutine xc_f03_gga_x_rpbe_set_params
  subroutine xc_f03_gga_x_optx_set_params(p, a, b, gamma)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: a, b, gamma
    call xc_gga_x_optx_set_params(p%ptr, a, b, gamma)
  end subroutine xc_f03_gga_x_optx_set_params
  subroutine xc_f03_gga_c_lyp_set_params(p, a, b, c, d)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: a, b, c, d
    call xc_gga_c_lyp_set_params(p%ptr, a, b, c, d)
  end subroutine xc_f03_gga_c_lyp_set_params
  subroutine xc_f03_gga_lb_set_params(p, modified, threshold, ip, qtot)
    type(xc_f03_func_t), intent(inout) :: p
    integer, intent(in) :: modified
    real(c_double), intent(in) :: threshold, ip, qtot
    call xc_gga_lb_set_params(p%ptr, modified, threshold, ip, qtot)
  end subroutine xc_f03_gga_lb_set_params
  subroutine xc_f03_gga_k_tflw_set_params(p, gamma, lambda, n)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: gamma, lambda
    integer, intent(in) :: n
    call xc_gga_k_tflw_set_params(p%ptr, gamma, lambda, n)
  end subroutine xc_f03_gga_k_tflw_set_params
  subroutine xc_f03_gga_x_2d_b88_set_params(p, beta)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: beta
    call xc_gga_x_2d_b88_set_params(p%ptr, beta)
  end subroutine xc_f03_gga_x_2d_b88_set_params
  subroutine xc_f03_gga_x_wpbeh_set_params(p, omega)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: omega
    call xc_gga_x_wpbeh_set_params(p%ptr, omega)
  end subroutine xc_f03_gga_x_wpbeh_set_params
  subroutine xc_f03_gga_x_hjs_set_params(p, omega)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: omega
    call xc_gga_x_hjs_set_params(p%ptr, omega)
  end subroutine xc_f03_gga_x_hjs_set_params
  subroutine xc_f03_gga_x_ityh_set_params(p, func_id, omega)
    type(xc_f03_func_t), intent(inout) :: p
    integer, intent(in) :: func_id
    real(c_double), intent(in) :: omega
    call xc_gga_x_ityh_set_params(p%ptr, func_id, omega)
  end subroutine xc_f03_gga_x_ityh_set_params
  subroutine xc_f03_gga_x_sfat_set_params(p, func_id, omega)
    type(xc_f03_func_t), intent(inout) :: p
    integer, intent(in) :: func_id
    real(c_double), intent(in) :: omega
    call xc_gga_x_sfat_set_params(p%ptr, func_id, omega)
  end subroutine xc_f03_gga_x_sfat_set_params
  subroutine xc_f03_gga_x_ssb_sw_set_params(p, a, b, c, d, e)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: a, b, c, d, e
    call xc_gga_x_ssb_sw_set_params(p%ptr, a, b, c, d, e)
  end subroutine xc_f03_gga_x_ssb_sw_set_params
  subroutine xc_f03_gga_x_kt_set_params(p, gamma, delta)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: gamma, delta
    call xc_gga_x_kt_set_params(p%ptr, gamma, delta)
  end subroutine xc_f03_gga_x_kt_set_params
  subroutine xc_f03_gga_x_lambda_set_params(p, n)
    type(xc_f03_func_t), intent(inout) :: p
    integer, intent(in) :: n
    call xc_gga_x_lambda_set_params(p%ptr, n)
  end subroutine xc_f03_gga_x_lambda_set_params
  real(c_double) function xc_f03_gga_ak13_get_asymptotic(homo) result(asymptotic)
    real(c_double), intent(in) :: homo
    asymptotic = xc_gga_ak13_get_asymptotic(homo)
  end function xc_f03_gga_ak13_get_asymptotic
  real(c_double) function xc_f03_hyb_exx_coef(p) result(coef)
    type(xc_f03_func_t), intent(in) :: p
    coef = xc_hyb_exx_coef(p%ptr)
  end function xc_f03_hyb_exx_coef
  subroutine xc_f03_hyb_cam_coef(p, omega, alpha, beta)
    type(xc_f03_func_t), intent(in) :: p
    real(c_double), intent(out) :: omega, alpha, beta
    call xc_hyb_cam_coef(p%ptr, omega, alpha, beta)
  end subroutine xc_f03_hyb_cam_coef
  subroutine xc_f03_nlc_coef(p, nlc_b, nlc_c)
    type(xc_f03_func_t), intent(in) :: p
    real(c_double), intent(out) :: nlc_b, nlc_c
    call xc_nlc_coef(p%ptr, nlc_b, nlc_c)
  end subroutine xc_f03_nlc_coef
  subroutine xc_f03_hyb_gga_xc_hse_set_params(p, alpha, omega)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: alpha, omega
    call xc_hyb_gga_xc_hse_set_params(p%ptr, alpha, omega)
  end subroutine xc_f03_hyb_gga_xc_hse_set_params
  subroutine xc_f03_hyb_gga_xc_pbeh_set_params(p, alpha)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: alpha
    call xc_hyb_gga_xc_pbeh_set_params(p%ptr, alpha)
  end subroutine xc_f03_hyb_gga_xc_pbeh_set_params
  ! the meta-GGAs
  !----------------------------------------------------------------
  subroutine xc_f03_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
    v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
    v2sigmalapl, v2sigmatau, v2lapltau)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double), intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(c_double), intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), v2rholapl(*), &
                                      v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)
    call xc_mgga(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)
  end subroutine xc_f03_mgga
  subroutine xc_f03_mgga_exc(p, np, rho, sigma, lapl, tau, zk)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double), intent(out) :: zk(*)
    call xc_mgga_exc(p%ptr, np, rho, sigma, lapl, tau, zk)
  end subroutine xc_f03_mgga_exc
  subroutine xc_f03_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double), intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    call xc_mgga_exc_vxc(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
  end subroutine xc_f03_mgga_exc_vxc
  subroutine xc_f03_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double), intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
    call xc_mgga_vxc(p%ptr, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
  end subroutine xc_f03_mgga_vxc
  subroutine xc_f03_mgga_fxc(p, np, rho, sigma, lapl, tau, &
    v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
    v2sigmalapl, v2sigmatau, v2lapltau)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real(c_double), intent(in) :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double), intent(out) :: v2rho2(*), v2sigma2(*), v2lapl2(*), v2tau2(*), v2rhosigma(*), &
                                      v2rholapl(*), v2rhotau(*), v2sigmalapl(*), v2sigmatau(*), v2lapltau(*)
    call xc_mgga_fxc(p%ptr, np, rho, sigma, lapl, tau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)
  end subroutine xc_f03_mgga_fxc
  subroutine xc_f03_mgga_x_tb09_set_params(p, c)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: c
    call xc_mgga_x_tb09_set_params(p%ptr, c)
  end subroutine xc_f03_mgga_x_tb09_set_params
  subroutine xc_f03_mgga_x_tpss_set_params(p, b, c, e, kappa, mu)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: b, c, e, kappa, mu
    call xc_mgga_x_tpss_set_params(p%ptr, b, c, e, kappa, mu)
  end subroutine xc_f03_mgga_x_tpss_set_params
  subroutine xc_f03_mgga_c_bc95_set_params(p, css, copp)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: css, copp
    call xc_mgga_c_bc95_set_params(p%ptr, css, copp)
  end subroutine xc_f03_mgga_c_bc95_set_params
  subroutine xc_f03_mgga_c_pkzb_set_params(p, beta, d, c0_0, c0_1, c0_2, c0_3)
    type(xc_f03_func_t), intent(inout) :: p
    real(c_double), intent(in) :: beta, d, c0_0, c0_1, c0_2, c0_3
    call xc_mgga_c_pkzb_set_params(p%ptr, beta, d, c0_0, c0_1, c0_2, c0_3)
  end subroutine xc_f03_mgga_c_pkzb_set_params
  ! Helper functions to convert between C and Fortran strings
  ! Based on the routines by Joseph M. Krahn
  function f_to_c_string(f_string) result(c_string)
    character(len=*), intent(in) :: f_string
    character(kind=c_char,len=1) :: c_string(len_trim(f_string)+1)
    integer :: i, strlen
    strlen = len_trim(f_string)
    forall (i=1:strlen)
      c_string(i) = f_string(i:i)
    end forall
    c_string(strlen+1) = C_NULL_CHAR
  end function f_to_c_string
  subroutine c_to_f_string(c_string, f_string)
    character(kind=c_char,len=1), intent(in) :: c_string(*)
    character(len=*), intent(out) :: f_string
    integer :: i
    i = 1
    do while(c_string(i) /= C_NULL_CHAR .and. i <= len(f_string))
      f_string(i:i) = c_string(i)
      i = i + 1
    end do
    if (i < len(f_string)) f_string(i:) = ' '
  end subroutine c_to_f_string
  subroutine c_to_f_string_ptr(c_string, f_string)
    type(c_ptr), intent(in) :: c_string
    character(len=*), intent(out) :: f_string
    character(len=1, kind=c_char), pointer :: p_chars(:)
    integer :: i
    if (.not. c_associated(c_string)) then
      f_string = ' '
    else
      call c_f_pointer(c_string, p_chars, [huge(0)])
      i = 1
      do while(p_chars(i) /= C_NULL_CHAR .and. i <= len(f_string))
        f_string(i:i) = p_chars(i)
        i = i + 1
      end do
      if (i < len(f_string)) f_string(i:) = ' '
    end if
  end subroutine c_to_f_string_ptr
end module xc_f03_lib_m
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
