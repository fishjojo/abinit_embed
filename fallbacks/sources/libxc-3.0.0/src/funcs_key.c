#include "util.h"

XC(functional_key_t) XC(functional_keys)[] = {
{"lda_x", 1},
{"lda_c_wigner", 2},
{"lda_c_rpa", 3},
{"lda_c_hl", 4},
{"lda_c_gl", 5},
{"lda_c_xalpha", 6},
{"lda_c_vwn", 7},
{"lda_c_vwn_rpa", 8},
{"lda_c_pz", 9},
{"lda_c_pz_mod", 10},
{"lda_c_ob_pz", 11},
{"lda_c_pw", 12},
{"lda_c_pw_mod", 13},
{"lda_c_ob_pw", 14},
{"lda_c_2d_amgb", 15},
{"lda_c_2d_prm", 16},
{"lda_c_vbh", 17},
{"lda_c_1d_csc", 18},
{"lda_x_2d", 19},
{"lda_xc_teter93", 20},
{"lda_x_1d", 21},
{"lda_c_ml1", 22},
{"lda_c_ml2", 23},
{"lda_c_gombas", 24},
{"lda_c_pw_rpa", 25},
{"lda_c_1d_loos", 26},
{"lda_c_rc04", 27},
{"lda_c_vwn_1", 28},
{"lda_c_vwn_2", 29},
{"lda_c_vwn_3", 30},
{"lda_c_vwn_4", 31},
{"lda_xc_zlp", 43},
{"lda_k_tf", 50},
{"lda_k_lp", 51},
{"lda_xc_ksdt", 259},
{"gga_x_gam", 32},
{"gga_c_gam", 33},
{"gga_x_hcth_a", 34},
{"gga_x_ev93", 35},
{"gga_x_bgcp", 38},
{"gga_c_bgcp", 39},
{"gga_x_lambda_oc2_n", 40},
{"gga_x_b86_r", 41},
{"gga_x_lambda_ch_n", 44},
{"gga_x_lambda_lo_n", 45},
{"gga_x_hjs_b88_v2", 46},
{"gga_c_q2d", 47},
{"gga_x_q2d", 48},
{"gga_x_pbe_mol", 49},
{"gga_k_tfvw", 52},
{"gga_k_revapbeint", 53},
{"gga_k_apbeint", 54},
{"gga_k_revapbe", 55},
{"gga_x_ak13", 56},
{"gga_k_meyer", 57},
{"gga_x_lv_rpw86", 58},
{"gga_x_pbe_tca", 59},
{"gga_x_pbeint", 60},
{"gga_c_zpbeint", 61},
{"gga_c_pbeint", 62},
{"gga_c_zpbesol", 63},
{"gga_xc_opbe_d", 65},
{"gga_xc_opwlyp_d", 66},
{"gga_xc_oblyp_d", 67},
{"gga_x_vmt84_ge", 68},
{"gga_x_vmt84_pbe", 69},
{"gga_x_vmt_ge", 70},
{"gga_x_vmt_pbe", 71},
{"gga_c_n12_sx", 79},
{"gga_c_n12", 80},
{"gga_x_n12", 82},
{"gga_c_regtpss", 83},
{"gga_c_op_xalpha", 84},
{"gga_c_op_g96", 85},
{"gga_c_op_pbe", 86},
{"gga_c_op_b88", 87},
{"gga_c_ft97", 88},
{"gga_c_spbe", 89},
{"gga_x_ssb_sw", 90},
{"gga_x_ssb", 91},
{"gga_x_ssb_d", 92},
{"gga_xc_hcth_407p", 93},
{"gga_xc_hcth_p76", 94},
{"gga_xc_hcth_p14", 95},
{"gga_xc_b97_gga1", 96},
{"gga_c_hcth_a", 97},
{"gga_x_bpccac", 98},
{"gga_c_revtca", 99},
{"gga_c_tca", 100},
{"gga_x_pbe", 101},
{"gga_x_pbe_r", 102},
{"gga_x_b86", 103},
{"gga_x_herman", 104},
{"gga_x_b86_mgc", 105},
{"gga_x_b88", 106},
{"gga_x_g96", 107},
{"gga_x_pw86", 108},
{"gga_x_pw91", 109},
{"gga_x_optx", 110},
{"gga_x_dk87_r1", 111},
{"gga_x_dk87_r2", 112},
{"gga_x_lg93", 113},
{"gga_x_ft97_a", 114},
{"gga_x_ft97_b", 115},
{"gga_x_pbe_sol", 116},
{"gga_x_rpbe", 117},
{"gga_x_wc", 118},
{"gga_x_mpw91", 119},
{"gga_x_am05", 120},
{"gga_x_pbea", 121},
{"gga_x_mpbe", 122},
{"gga_x_xpbe", 123},
{"gga_x_2d_b86_mgc", 124},
{"gga_x_bayesian", 125},
{"gga_x_pbe_jsjr", 126},
{"gga_x_2d_b88", 127},
{"gga_x_2d_b86", 128},
{"gga_x_2d_pbe", 129},
{"gga_c_pbe", 130},
{"gga_c_lyp", 131},
{"gga_c_p86", 132},
{"gga_c_pbe_sol", 133},
{"gga_c_pw91", 134},
{"gga_c_am05", 135},
{"gga_c_xpbe", 136},
{"gga_c_lm", 137},
{"gga_c_pbe_jrgx", 138},
{"gga_x_optb88_vdw", 139},
{"gga_x_pbek1_vdw", 140},
{"gga_x_optpbe_vdw", 141},
{"gga_x_rge2", 142},
{"gga_c_rge2", 143},
{"gga_x_rpw86", 144},
{"gga_x_kt1", 145},
{"gga_xc_kt2", 146},
{"gga_c_wl", 147},
{"gga_c_wi", 148},
{"gga_x_mb88", 149},
{"gga_x_sogga", 150},
{"gga_x_sogga11", 151},
{"gga_c_sogga11", 152},
{"gga_c_wi0", 153},
{"gga_xc_th1", 154},
{"gga_xc_th2", 155},
{"gga_xc_th3", 156},
{"gga_xc_th4", 157},
{"gga_x_c09x", 158},
{"gga_c_sogga11_x", 159},
{"gga_x_lb", 160},
{"gga_xc_hcth_93", 161},
{"gga_xc_hcth_120", 162},
{"gga_xc_hcth_147", 163},
{"gga_xc_hcth_407", 164},
{"gga_xc_edf1", 165},
{"gga_xc_xlyp", 166},
{"gga_xc_b97_d", 170},
{"gga_xc_pbe1w", 173},
{"gga_xc_mpwlyp1w", 174},
{"gga_xc_pbelyp1w", 175},
{"gga_x_lbm", 182},
{"gga_x_ol2", 183},
{"gga_x_apbe", 184},
{"gga_k_apbe", 185},
{"gga_c_apbe", 186},
{"gga_k_tw1", 187},
{"gga_k_tw2", 188},
{"gga_k_tw3", 189},
{"gga_k_tw4", 190},
{"gga_x_htbs", 191},
{"gga_x_airy", 192},
{"gga_x_lag", 193},
{"gga_xc_mohlyp", 194},
{"gga_xc_mohlyp2", 195},
{"gga_xc_th_fl", 196},
{"gga_xc_th_fc", 197},
{"gga_xc_th_fcfo", 198},
{"gga_xc_th_fco", 199},
{"gga_c_optc", 200},
{"gga_c_pbeloc", 246},
{"gga_xc_vv10", 255},
{"gga_c_pbefe", 258},
{"gga_c_op_pw91", 262},
{"gga_x_pbefe", 265},
{"gga_x_cap", 270},
{"gga_k_vw", 500},
{"gga_k_ge2", 501},
{"gga_k_golden", 502},
{"gga_k_yt65", 503},
{"gga_k_baltin", 504},
{"gga_k_lieb", 505},
{"gga_k_absp1", 506},
{"gga_k_absp2", 507},
{"gga_k_gr", 508},
{"gga_k_ludena", 509},
{"gga_k_gp85", 510},
{"gga_k_pearson", 511},
{"gga_k_ol1", 512},
{"gga_k_ol2", 513},
{"gga_k_fr_b88", 514},
{"gga_k_fr_pw86", 515},
{"gga_k_dk", 516},
{"gga_k_perdew", 517},
{"gga_k_vsk", 518},
{"gga_k_vjks", 519},
{"gga_k_ernzerhof", 520},
{"gga_k_lc94", 521},
{"gga_k_llp", 522},
{"gga_k_thakkar", 523},
{"gga_x_wpbeh", 524},
{"gga_x_hjs_pbe", 525},
{"gga_x_hjs_pbe_sol", 526},
{"gga_x_hjs_b88", 527},
{"gga_x_hjs_b97x", 528},
{"gga_x_ityh", 529},
{"gga_x_sfat", 530},
{"hyb_gga_x_n12_sx", 81},
{"hyb_gga_xc_b97_1p", 266},
{"hyb_gga_xc_b3pw91", 401},
{"hyb_gga_xc_b3lyp", 402},
{"hyb_gga_xc_b3p86", 403},
{"hyb_gga_xc_o3lyp", 404},
{"hyb_gga_xc_mpw1k", 405},
{"hyb_gga_xc_pbeh", 406},
{"hyb_gga_xc_b97", 407},
{"hyb_gga_xc_b97_1", 408},
{"hyb_gga_xc_b97_2", 410},
{"hyb_gga_xc_x3lyp", 411},
{"hyb_gga_xc_b1wc", 412},
{"hyb_gga_xc_b97_k", 413},
{"hyb_gga_xc_b97_3", 414},
{"hyb_gga_xc_mpw3pw", 415},
{"hyb_gga_xc_b1lyp", 416},
{"hyb_gga_xc_b1pw91", 417},
{"hyb_gga_xc_mpw1pw", 418},
{"hyb_gga_xc_mpw3lyp", 419},
{"hyb_gga_xc_sb98_1a", 420},
{"hyb_gga_xc_sb98_1b", 421},
{"hyb_gga_xc_sb98_1c", 422},
{"hyb_gga_xc_sb98_2a", 423},
{"hyb_gga_xc_sb98_2b", 424},
{"hyb_gga_xc_sb98_2c", 425},
{"hyb_gga_x_sogga11_x", 426},
{"hyb_gga_xc_hse03", 427},
{"hyb_gga_xc_hse06", 428},
{"hyb_gga_xc_hjs_pbe", 429},
{"hyb_gga_xc_hjs_pbe_sol", 430},
{"hyb_gga_xc_hjs_b88", 431},
{"hyb_gga_xc_hjs_b97x", 432},
{"hyb_gga_xc_cam_b3lyp", 433},
{"hyb_gga_xc_tuned_cam_b3lyp", 434},
{"hyb_gga_xc_bhandh", 435},
{"hyb_gga_xc_bhandhlyp", 436},
{"hyb_gga_xc_mb3lyp_rc04", 437},
{"hyb_gga_xc_mpwlyp1m", 453},
{"hyb_gga_xc_revb3lyp", 454},
{"hyb_gga_xc_camy_blyp", 455},
{"hyb_gga_xc_pbe0_13", 456},
{"hyb_gga_xc_b3lyps", 459},
{"hyb_gga_xc_wb97", 463},
{"hyb_gga_xc_wb97x", 464},
{"hyb_gga_xc_lrc_wpbeh", 465},
{"hyb_gga_xc_wb97x_v", 466},
{"hyb_gga_xc_lcy_pbe", 467},
{"hyb_gga_xc_lcy_blyp", 468},
{"hyb_gga_xc_lc_vv10", 469},
{"hyb_gga_xc_camy_b3lyp", 470},
{"hyb_gga_xc_wb97x_d", 471},
{"hyb_gga_xc_hpbeint", 472},
{"hyb_gga_xc_lrc_wpbe", 473},
{"hyb_gga_xc_b3lyp5", 475},
{"hyb_gga_xc_edf2", 476},
{"hyb_gga_xc_cap0", 477},
{"mgga_c_dldf", 37},
{"mgga_xc_zlp", 42},
{"mgga_xc_otpss_d", 64},
{"mgga_c_cs", 72},
{"mgga_c_mn12_sx", 73},
{"mgga_c_mn12_l", 74},
{"mgga_c_m11_l", 75},
{"mgga_c_m11", 76},
{"mgga_c_m08_so", 77},
{"mgga_c_m08_hx", 78},
{"mgga_x_lta", 201},
{"mgga_x_tpss", 202},
{"mgga_x_m06_l", 203},
{"mgga_x_gvt4", 204},
{"mgga_x_tau_hcth", 205},
{"mgga_x_br89", 206},
{"mgga_x_bj06", 207},
{"mgga_x_tb09", 208},
{"mgga_x_rpp09", 209},
{"mgga_x_2d_prhg07", 210},
{"mgga_x_2d_prhg07_prp10", 211},
{"mgga_x_revtpss", 212},
{"mgga_x_pkzb", 213},
{"mgga_x_m05", 214},
{"mgga_x_m05_2x", 215},
{"mgga_x_m06_hf", 216},
{"mgga_x_m06", 217},
{"mgga_x_m06_2x", 218},
{"mgga_x_m08_hx", 219},
{"mgga_x_m08_so", 220},
{"mgga_x_ms0", 221},
{"mgga_x_ms1", 222},
{"mgga_x_ms2", 223},
{"mgga_x_m11", 225},
{"mgga_x_m11_l", 226},
{"mgga_x_mn12_l", 227},
{"mgga_c_cc06", 229},
{"mgga_x_mk00", 230},
{"mgga_c_tpss", 231},
{"mgga_c_vsxc", 232},
{"mgga_c_m06_l", 233},
{"mgga_c_m06_hf", 234},
{"mgga_c_m06", 235},
{"mgga_c_m06_2x", 236},
{"mgga_c_m05", 237},
{"mgga_c_m05_2x", 238},
{"mgga_c_pkzb", 239},
{"mgga_c_bc95", 240},
{"mgga_c_revtpss", 241},
{"mgga_xc_tpsslyp1w", 242},
{"mgga_x_mk00b", 243},
{"mgga_x_bloc", 244},
{"mgga_x_modtpss", 245},
{"mgga_c_tpssloc", 247},
{"mgga_x_mbeef", 249},
{"mgga_x_mbeefvdw", 250},
{"mgga_xc_b97m_v", 254},
{"mgga_x_mvs", 257},
{"mgga_x_mn15_l", 260},
{"mgga_c_mn15_l", 261},
{"mgga_x_scan", 263},
{"mgga_c_scan", 267},
{"mgga_c_mn15", 269},
{"hyb_mgga_x_dldf", 36},
{"hyb_mgga_x_ms2h", 224},
{"hyb_mgga_x_mn12_sx", 248},
{"hyb_mgga_x_scan0", 264},
{"hyb_mgga_x_mn15", 268},
{"hyb_mgga_xc_m05", 438},
{"hyb_mgga_xc_m05_2x", 439},
{"hyb_mgga_xc_b88b95", 440},
{"hyb_mgga_xc_b86b95", 441},
{"hyb_mgga_xc_pw86b95", 442},
{"hyb_mgga_xc_bb1k", 443},
{"hyb_mgga_xc_m06_hf", 444},
{"hyb_mgga_xc_mpw1b95", 445},
{"hyb_mgga_xc_mpwb1k", 446},
{"hyb_mgga_xc_x1b95", 447},
{"hyb_mgga_xc_xb1k", 448},
{"hyb_mgga_xc_m06", 449},
{"hyb_mgga_xc_m06_2x", 450},
{"hyb_mgga_xc_pw6b95", 451},
{"hyb_mgga_xc_pwb6k", 452},
{"hyb_mgga_xc_tpssh", 457},
{"hyb_mgga_xc_revtpssh", 458},
{"hyb_mgga_xc_m08_hx", 460},
{"hyb_mgga_xc_m08_so", 461},
{"hyb_mgga_xc_m11", 462},
{"hyb_mgga_x_mvsh", 474},
{"hyb_mgga_xc_wb97m_v", 531},
{"", -1}
};
