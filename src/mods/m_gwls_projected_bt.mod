  W  ±   k820309    l          18.0        ³Ø[                                                                                                          
       m_gwls_Projected_BT.F90 M_GWLS_PROJECTED_BT              COMPUTE_PROJECTED_BT_SHIFT_LANCZOS COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED BT_LSTERNHEIMER PROJECTED_BT_A_MATRIX PROJECTED_BT_L_MATRIX PROJECTED_BT_BETA_MATRIX BT_LSTERNHEIMER_LANCZOS PROJECTED_BT_A_MATRIX_LANCZOS PROJECTED_BT_L_MATRIX_LANCZOS PROJECTED_BT_BETA_MATRIX_LANCZOS BT_LSTERNHEIMER_MODEL_LANCZOS PROJECTED_BT_A_MATRIX_MODEL_LANCZOS PROJECTED_BT_L_MATRIX_MODEL_LANCZOS PROJECTED_BT_BETA_MATRIX_MODEL_LANCZOS                      @                              
                            @                              
                            @                              
                @          @                              
                            @                              
                            @                              
                            @                              
                            @                              
                            @                         	     
                                                      
     
                                                           
                            @                              
                        @                               '                                      @  @               E                '              B      #COMM_WORLD    #ME    #NPROC    #ME_G0    #COMM_ATOM    #NPROC_ATOM    #MY_NATOM    #MY_ATMTAB    #PARAL_PERT    #COMM_PERT    #COMM_CELL_PERT    #ME_PERT    #NPROC_PERT    #DISTRB_PERT    #PARAL_IMG    #MY_NIMAGE    #COMM_IMG    #ME_IMG     #NPROC_IMG !   #DISTRB_IMG "   #MY_IMGTAB #   #COMM_CELL $   #ME_CELL %   #NPROC_CELL &   #COMM_FFT '   #ME_FFT (   #NPROC_FFT )   #DISTRIBFFT *   #PARALBD :   #COMM_BAND ;   #ME_BAND <   #NPROC_BAND =   #PARAL_SPINOR >   #COMM_SPINOR ?   #ME_SPINOR @   #NPROC_SPINOR A   #COMM_KPT B   #ME_KPT C   #NPROC_KPT D   #PROC_DISTRB E   #MY_ISPPOLTAB F   #PARAL_KGB G   #BANDPP H   #COMM_BANDSPINORFFT I   #COMM_BANDFFT J   #COMM_KPTBAND K   #COMM_SPINORFFT L   #COMM_BANDSPINOR M   #MY_KGTAB N   #MY_KPTTAB O   #PW_UNBAL_THRESH P   #KPTDSTRB Q   #KPT_LOC2FBZ_SP R   #KPT_LOC2IBZ_SP S   #MKMEM T   #PARAL_HF U   #COMM_HF V   #ME_HF W   #NPROC_HF X   #DISTRB_HF Y   #COMM_WVL Z   #ME_WVL [   #NPROC_WVL \   #NSCATTERARR ]   #NGATHERARR ^   #NGFFT3_IONIC _                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        &                                                                   ?              y                                                                                                h       	                                                       l       
                                                       p                                                              t                                                              x                                                                                                &                                                                                            È                                                              Ì                                                              Ð                                                               Ô                                                         !     Ø                                                       "            à                             &                                                                                     #            (                            &                                                                                       $     p                                                        %     t                                                        &     x                                                        '     |                                                        (                                                             )                                                           *     h                  #DISTRIBFFT_TYPE +                           ?              y#DISTRIBFFT_TYPE +                                                                 @  @               E           +     'h                   #N2_COARSE ,   #N2_FINE -   #TAB_FFTWF2_DISTRIB .   #TAB_FFTDP2_DISTRIB /   #TAB_FFTDP3_DISTRIB 0   #TAB_FFTWF2DG_DISTRIB 1   #TAB_FFTDP2DG_DISTRIB 2   #TAB_FFTDP3DG_DISTRIB 3   #TAB_FFTWF2_LOCAL 4   #TAB_FFTDP2_LOCAL 5   #TAB_FFTDP3_LOCAL 6   #TAB_FFTWF2DG_LOCAL 7   #TAB_FFTDP2DG_LOCAL 8   #TAB_FFTDP3DG_LOCAL 9                $                              ,                                          >                                                      0                 $                              -                                         >                                                      0               $                              .                                         &                                                      $                              /            P                             &                                                      $                              0                                         &                                                      $                              1            à                             &                                                      $                              2            (                            &                                                      $                              3            p                            &                                                      $                              4            ¸             	               &                                                      $                              5                          
               &                                                      $                              6            H                            &                                                      $                              7                                        &                                                      $                              8            Ø                            &                                                      $                              9                                         &                                                                                       :                                                             ;                                                             <                                                             =                                                              >            !                                                  ?     ¤      "                                                  @     ¨      #                                                  A     ¬      $                                                  B     °      %                                                  C     ´      &                                                  D     ¸      '                                                E            À             (               &                   &                   &                                                                                       F            8             )     p          p            p                                                                      G     @      *                                                  H     D      +                                                  I     H      ,                                                  J     L      -                                                  K     P      .                                                  L     T      /                                                  M     X      0                                                N            `             1               &                   &                                                                                     O            À             2               &                                                                                      P           3   
                                             Q                         4               &                   &                   &                                                                                     R                         5               &                   &                   &                                                                                     S                          6               &                   &                   &                                                                                     T            x             7               &                                                                                       U     À      8                                                  V     Ä      9                                                  W     È      :                                                  X     Ì      ;                                                Y            Ð             <               &                   &                   &                                                                                       Z     H      =                                                  [     L      >                                                  \     P      ?                                               ]            X             @              &                   &                                                                   ?              y                                                                                        ^            ¸             A              &                   &                                                                   ?              y                                                                                           _           B                     @ @                          `     '                    #MPI_VAL a                                               a                                  @ @                          b     '                    #MPI_VAL c                                               c                                  @ @                          d     '                    #MPI_VAL e                                               e                                  @ @                          f     '                    #MPI_VAL g                                               g                                  @ @                          h     '                    #MPI_VAL i                                               i                                  @ @                          j     '                    #MPI_VAL k                                               k                                  @ @                          l     '                    #MPI_VAL m                                               m                                  @ @                          n     '                    #MPI_VAL o                                               o                                  @ @                          p     '                    #MPI_VAL q                                               q                                  @ @                          r     '                    #MPI_VAL s                                               s                             @@                               t                                                        u                                                        v                        @                               w                                                     x                   
                &                                                                                       y                                 	             ð?        (1.0,0.0)                                            z                                 	                     ð?(0.0,1.0)                                            {                                 	                       (0.0,0.0)#         @                                  |                    #KMAX }   #PREC ~             
                                  }                     
                                  ~                       @                                           #MPI_TYPE    #         @                                                     #NPW_K t   #PSIK    #PSIK_ALLTOALL    #DIRECTION             
                                                    
 4    p          p          5 r u       p          5 r u                              
                                                    
 5    p          p          5 r v       p          5 r v                               
                                             #         @                                                     #NPW_G    #NPW_K t   #SEED_VECTOR    #RIGHT_VEC                @                                                 
                                                         p          5 r        5 r                                                                                        p          5 r        5 r                                @                                                                 &                   &                                                                                                                       &                   &                                           #         @                                                      #N    #MATRIX              
                                                      
                                                           p        5 O p        p          5 O p          5 O p            5 O p          5 O p                          #         @                                                                                                                                                                                           &                   &                                                                                                                       &                   &                                                                                                                       &                   &                                                                                                                                                                               &                   &                                                                                                                       &                   &                                                                                                                       &                   &                                                                                                                                                                               &                   &                                                                                                                       &                   &                                                                                                                       &                   &                                           #         @                                                       #NFREQ    #LIST_EXTERNAL_OMEGA    #LMAX    #MODIFIED_LBASIS    #KMAX_NUMERIC     #NPT_GAUSS ¡   #DIELECTRIC_ARRAY ¢   #ARRAY_INTEGRAND £             
                                                      
                                                     
 
   p          5  p        r        5  p        r                                
@ @                                                   
@ @                                                        p        5 r t   p          5 r t     5  p        r        5 r t     5  p        r                                
@ @                                                     
                                  ¡                    
                                 ¢                            p        5  p        r    p        5  p        r    p          5  p        r      5  p        r       5  p        r ¡   n                                       1    5  p        r      5  p        r       5  p        r ¡   n                                      1                                    D                                £                           p         5  p        r ¡   n                                           1p           5  p        r ¡   n                                      1  5  p        r         5  p        r ¡   n                                      1  5  p        r                                #         @                                   ¤                    #NFREQ ¥   #LIST_EXTERNAL_OMEGA ¦   #LMAX §   #BLOCKSIZE_EPS ¨   #MODEL_LANCZOS_VECTOR_BELONGS_TO_THIS_NODE ©   #MODEL_LANCZOS_VECTOR_INDEX ª   #MODIFIED_LBASIS «   #KMAX_NUMERIC ¬   #NPT_GAUSS ­   #DIELECTRIC_ARRAY ®   #ARRAY_INTEGRAND ¯             
                                  ¥                    
                                 ¦                    
    p          5  p        r ¥       5  p        r ¥                               
@ @                               §                     
                                  ¨                    
                                  ©                     "   p          5  p        r §       5  p        r §                              
                                  ª                     #   p          5  p        r §       5  p        r §                              
@ @                              «                           p        5 r t   p          5 r t     5  p        r §       5 r t     5  p        r §                               
@ @                               ¬                     
                                  ­                    
                                 ®                     !       p        5  p        r ¨   p        5  p        r §   p          5  p        r §     5  p        r ¨      5  p 	       r ­   n                                       1    5  p        r §     5  p        r ¨      5  p 	       r ­   n                                      1                                    D @                              ¯                     $      p         5  p 	       r ­   n                                           1p           5  p 	       r ­   n                                      1  5  p        r ¥        5  p 	       r ­   n                                      1  5  p        r ¥                                            @                                                4      fn#fn )   Ô   ¹  b   uapp(M_GWLS_PROJECTED_BT      @   J  M_GWLS_UTILITY    Í  @   J  M_GWLS_WF !     @   J  M_GWLS_TIMINGLOG #   M  @   j  M_GWLS_HAMILTONIAN #     @   J  M_GWLS_LINEQSOLVER !   Í  @   J  M_GWLS_GWLANCZOS $     @   J  M_GWLS_LANCZOSBASIS '   M  @   J  M_GWLS_DIELECTRICARRAY )     @   J  M_GWLS_LANCZOSRESOLVENTS    Í  @   J  DEFS_BASIS      @   J  M_ABICORE    M  @   J  M_XMPI '     P       #UNLPOLY+ISO_C_BINDING '   Ý  A     MPI_TYPE+DEFS_ABITYPES 2   
  H   a   MPI_TYPE%COMM_WORLD+DEFS_ABITYPES *   f
  H   a   MPI_TYPE%ME+DEFS_ABITYPES -   ®
  H   a   MPI_TYPE%NPROC+DEFS_ABITYPES -   ö
  H   a   MPI_TYPE%ME_G0+DEFS_ABITYPES 1   >  H   a   MPI_TYPE%COMM_ATOM+DEFS_ABITYPES 2     H   a   MPI_TYPE%NPROC_ATOM+DEFS_ABITYPES 0   Î  H   a   MPI_TYPE%MY_NATOM+DEFS_ABITYPES 1     ô   a   MPI_TYPE%MY_ATMTAB+DEFS_ABITYPES 2   
  H   a   MPI_TYPE%PARAL_PERT+DEFS_ABITYPES 1   R  H   a   MPI_TYPE%COMM_PERT+DEFS_ABITYPES 6     H   a   MPI_TYPE%COMM_CELL_PERT+DEFS_ABITYPES /   â  H   a   MPI_TYPE%ME_PERT+DEFS_ABITYPES 2   *  H   a   MPI_TYPE%NPROC_PERT+DEFS_ABITYPES 3   r     a   MPI_TYPE%DISTRB_PERT+DEFS_ABITYPES 1     H   a   MPI_TYPE%PARAL_IMG+DEFS_ABITYPES 1   N  H   a   MPI_TYPE%MY_NIMAGE+DEFS_ABITYPES 0     H   a   MPI_TYPE%COMM_IMG+DEFS_ABITYPES .   Þ  H   a   MPI_TYPE%ME_IMG+DEFS_ABITYPES 1   &  H   a   MPI_TYPE%NPROC_IMG+DEFS_ABITYPES 2   n     a   MPI_TYPE%DISTRB_IMG+DEFS_ABITYPES 1        a   MPI_TYPE%MY_IMGTAB+DEFS_ABITYPES 1     H   a   MPI_TYPE%COMM_CELL+DEFS_ABITYPES /   Þ  H   a   MPI_TYPE%ME_CELL+DEFS_ABITYPES 2   &  H   a   MPI_TYPE%NPROC_CELL+DEFS_ABITYPES 0   n  H   a   MPI_TYPE%COMM_FFT+DEFS_ABITYPES .   ¶  H   a   MPI_TYPE%ME_FFT+DEFS_ABITYPES 1   þ  H   a   MPI_TYPE%NPROC_FFT+DEFS_ABITYPES 2   F  Ú   a   MPI_TYPE%DISTRIBFFT+DEFS_ABITYPES -           DISTRIBFFT_TYPE+M_DISTRIBFFT 7   ¬  ¥   a   DISTRIBFFT_TYPE%N2_COARSE+M_DISTRIBFFT 5   Q  ¥   a   DISTRIBFFT_TYPE%N2_FINE+M_DISTRIBFFT @   ö     a   DISTRIBFFT_TYPE%TAB_FFTWF2_DISTRIB+M_DISTRIBFFT @        a   DISTRIBFFT_TYPE%TAB_FFTDP2_DISTRIB+M_DISTRIBFFT @        a   DISTRIBFFT_TYPE%TAB_FFTDP3_DISTRIB+M_DISTRIBFFT B   ²     a   DISTRIBFFT_TYPE%TAB_FFTWF2DG_DISTRIB+M_DISTRIBFFT B   F     a   DISTRIBFFT_TYPE%TAB_FFTDP2DG_DISTRIB+M_DISTRIBFFT B   Ú     a   DISTRIBFFT_TYPE%TAB_FFTDP3DG_DISTRIB+M_DISTRIBFFT >   n     a   DISTRIBFFT_TYPE%TAB_FFTWF2_LOCAL+M_DISTRIBFFT >        a   DISTRIBFFT_TYPE%TAB_FFTDP2_LOCAL+M_DISTRIBFFT >        a   DISTRIBFFT_TYPE%TAB_FFTDP3_LOCAL+M_DISTRIBFFT @   *     a   DISTRIBFFT_TYPE%TAB_FFTWF2DG_LOCAL+M_DISTRIBFFT @   ¾     a   DISTRIBFFT_TYPE%TAB_FFTDP2DG_LOCAL+M_DISTRIBFFT @   R     a   DISTRIBFFT_TYPE%TAB_FFTDP3DG_LOCAL+M_DISTRIBFFT /   æ  H   a   MPI_TYPE%PARALBD+DEFS_ABITYPES 1   .  H   a   MPI_TYPE%COMM_BAND+DEFS_ABITYPES /   v  H   a   MPI_TYPE%ME_BAND+DEFS_ABITYPES 2   ¾  H   a   MPI_TYPE%NPROC_BAND+DEFS_ABITYPES 4     H   a   MPI_TYPE%PARAL_SPINOR+DEFS_ABITYPES 3   N  H   a   MPI_TYPE%COMM_SPINOR+DEFS_ABITYPES 1     H   a   MPI_TYPE%ME_SPINOR+DEFS_ABITYPES 4   Þ  H   a   MPI_TYPE%NPROC_SPINOR+DEFS_ABITYPES 0   &   H   a   MPI_TYPE%COMM_KPT+DEFS_ABITYPES .   n   H   a   MPI_TYPE%ME_KPT+DEFS_ABITYPES 1   ¶   H   a   MPI_TYPE%NPROC_KPT+DEFS_ABITYPES 3   þ   Ä   a   MPI_TYPE%PROC_DISTRB+DEFS_ABITYPES 4   Â!     a   MPI_TYPE%MY_ISPPOLTAB+DEFS_ABITYPES 1   ^"  H   a   MPI_TYPE%PARAL_KGB+DEFS_ABITYPES .   ¦"  H   a   MPI_TYPE%BANDPP+DEFS_ABITYPES :   î"  H   a   MPI_TYPE%COMM_BANDSPINORFFT+DEFS_ABITYPES 4   6#  H   a   MPI_TYPE%COMM_BANDFFT+DEFS_ABITYPES 4   ~#  H   a   MPI_TYPE%COMM_KPTBAND+DEFS_ABITYPES 6   Æ#  H   a   MPI_TYPE%COMM_SPINORFFT+DEFS_ABITYPES 7   $  H   a   MPI_TYPE%COMM_BANDSPINOR+DEFS_ABITYPES 0   V$  ¬   a   MPI_TYPE%MY_KGTAB+DEFS_ABITYPES 1   %     a   MPI_TYPE%MY_KPTTAB+DEFS_ABITYPES 7   %  H   a   MPI_TYPE%PW_UNBAL_THRESH+DEFS_ABITYPES 0   Þ%  Ä   a   MPI_TYPE%KPTDSTRB+DEFS_ABITYPES 6   ¢&  Ä   a   MPI_TYPE%KPT_LOC2FBZ_SP+DEFS_ABITYPES 6   f'  Ä   a   MPI_TYPE%KPT_LOC2IBZ_SP+DEFS_ABITYPES -   *(     a   MPI_TYPE%MKMEM+DEFS_ABITYPES 0   ¾(  H   a   MPI_TYPE%PARAL_HF+DEFS_ABITYPES /   )  H   a   MPI_TYPE%COMM_HF+DEFS_ABITYPES -   N)  H   a   MPI_TYPE%ME_HF+DEFS_ABITYPES 0   )  H   a   MPI_TYPE%NPROC_HF+DEFS_ABITYPES 1   Þ)  Ä   a   MPI_TYPE%DISTRB_HF+DEFS_ABITYPES 0   ¢*  H   a   MPI_TYPE%COMM_WVL+DEFS_ABITYPES .   ê*  H   a   MPI_TYPE%ME_WVL+DEFS_ABITYPES 1   2+  H   a   MPI_TYPE%NPROC_WVL+DEFS_ABITYPES 3   z+    a   MPI_TYPE%NSCATTERARR+DEFS_ABITYPES 2   ,    a   MPI_TYPE%NGATHERARR+DEFS_ABITYPES 4   -  H   a   MPI_TYPE%NGFFT3_IONIC+DEFS_ABITYPES (   Ú-  ]      MPI_GROUP+MPI_CONSTANTS 0   7.  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS '   .  ]      MPI_INFO+MPI_CONSTANTS /   Ü.  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS *   $/  ]      MPI_MESSAGE+MPI_CONSTANTS 2   /  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS &   É/  ]      MPI_WIN+MPI_CONSTANTS .   &0  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS +   n0  ]      MPI_DATATYPE+MPI_CONSTANTS 3   Ë0  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS *   1  ]      MPI_REQUEST+MPI_CONSTANTS 2   p1  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS '   ¸1  ]      MPI_FILE+MPI_CONSTANTS /   2  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS -   ]2  ]      MPI_ERRHANDLER+MPI_CONSTANTS 5   º2  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '   3  ]      MPI_COMM+MPI_CONSTANTS /   _3  H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS %   §3  ]      MPI_OP+MPI_CONSTANTS -   4  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS )   L4  @       NPW_K+M_GWLS_HAMILTONIAN *   4  @       NPW_KB+M_GWLS_HAMILTONIAN )   Ì4  @       NPW_G+M_GWLS_HAMILTONIAN -   5  @       BLOCKSIZE+M_GWLS_HAMILTONIAN 2   L5         LIST_OMEGA+M_GWLS_DIELECTRICARRAY '   Ø5  y       CMPLX_1+M_GWLS_UTILITY '   Q6  y       CMPLX_I+M_GWLS_UTILITY '   Ê6  y       CMPLX_0+M_GWLS_UTILITY A   C7  \       SETUP_LANCZOSRESOLVENTS+M_GWLS_LANCZOSRESOLVENTS F   7  @   a   SETUP_LANCZOSRESOLVENTS%KMAX+M_GWLS_LANCZOSRESOLVENTS F   ß7  @   a   SETUP_LANCZOSRESOLVENTS%PREC+M_GWLS_LANCZOSRESOLVENTS -   8  N       MPI_ENREG+M_GWLS_HAMILTONIAN 7   m8         WF_BLOCK_DISTRIBUTE+M_GWLS_HAMILTONIAN <   ì8  ´   a   WF_BLOCK_DISTRIBUTE%PSIK+M_GWLS_HAMILTONIAN E    9  ´   a   WF_BLOCK_DISTRIBUTE%PSIK_ALLTOALL+M_GWLS_HAMILTONIAN A   T:  @   a   WF_BLOCK_DISTRIBUTE%DIRECTION+M_GWLS_HAMILTONIAN ^   :  ~       COMPUTE_RESOLVENT_COLUMN_SHIFT_LANCZOS_RIGHT_VECTORS+M_GWLS_LANCZOSRESOLVENTS )   ;  @     NPW_G+M_GWLS_HAMILTONIAN j   R;     a   COMPUTE_RESOLVENT_COLUMN_SHIFT_LANCZOS_RIGHT_VECTORS%SEED_VECTOR+M_GWLS_LANCZOSRESOLVENTS h   æ;     a   COMPUTE_RESOLVENT_COLUMN_SHIFT_LANCZOS_RIGHT_VECTORS%RIGHT_VEC+M_GWLS_LANCZOSRESOLVENTS 5   z<  ¤       LR_M_MATRIX+M_GWLS_LANCZOSRESOLVENTS 8   =  ¤       HAMILTONIAN_QK+M_GWLS_LANCZOSRESOLVENTS ?   Â=  [       INVERT_GENERAL_MATRIX+M_GWLS_LANCZOSRESOLVENTS A   >  @   a   INVERT_GENERAL_MATRIX%N+M_GWLS_LANCZOSRESOLVENTS F   ]>  ü   a   INVERT_GENERAL_MATRIX%MATRIX+M_GWLS_LANCZOSRESOLVENTS C   Y?  H       CLEANUP_LANCZOSRESOLVENTS+M_GWLS_LANCZOSRESOLVENTS     ¡?  @       BT_LSTERNHEIMER &   á?  ¤       PROJECTED_BT_A_MATRIX &   @  ¤       PROJECTED_BT_L_MATRIX )   )A  ¤       PROJECTED_BT_BETA_MATRIX (   ÍA  @       BT_LSTERNHEIMER_LANCZOS .   B  ¤       PROJECTED_BT_A_MATRIX_LANCZOS .   ±B  ¤       PROJECTED_BT_L_MATRIX_LANCZOS 1   UC  ¤       PROJECTED_BT_BETA_MATRIX_LANCZOS .   ùC  @       BT_LSTERNHEIMER_MODEL_LANCZOS 4   9D  ¤       PROJECTED_BT_A_MATRIX_MODEL_LANCZOS 4   ÝD  ¤       PROJECTED_BT_L_MATRIX_MODEL_LANCZOS 7   E  ¤       PROJECTED_BT_BETA_MATRIX_MODEL_LANCZOS 3   %F  ×       COMPUTE_PROJECTED_BT_SHIFT_LANCZOS 9   üF  @   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS%NFREQ G   <G  ´   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS%LIST_EXTERNAL_OMEGA 8   ðG  @   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS%LMAX C   0H  ô   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS%MODIFIED_LBASIS @   $I  @   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS%KMAX_NUMERIC =   dI  @   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS%NPT_GAUSS D   ¤I    a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS%DIELECTRIC_ARRAY C   ªK  Ï  a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS%ARRAY_INTEGRAND ?   yM  9      COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED E   ²N  @   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED%NFREQ S   òN  ´   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED%LIST_EXTERNAL_OMEGA D   ¦O  @   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED%LMAX M   æO  @   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED%BLOCKSIZE_EPS i   &P  ´   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED%MODEL_LANCZOS_VECTOR_BELONGS_TO_THIS_NODE Z   ÚP  ´   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED%MODEL_LANCZOS_VECTOR_INDEX O   Q  ô   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED%MODIFIED_LBASIS L   R  @   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED%KMAX_NUMERIC I   ÂR  @   a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED%NPT_GAUSS P   S    a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED%DIELECTRIC_ARRAY O   U  Ï  a   COMPUTE_PROJECTED_BT_SHIFT_LANCZOS_DISTRIBUTED%ARRAY_INTEGRAND 1   ×V  @      LR_KMAX+M_GWLS_LANCZOSRESOLVENTS 