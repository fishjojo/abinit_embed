  5d  é   k820309    l          18.0        ±²Ø[                                                                                                          
       defs_rectypes.F90 DEFS_RECTYPES                                                     
                         @                               
                                                           
       MPI_TYPE                      @                              
       PAWFGR_TYPE                                                        
    #GROUPEQ    #INFOEQ    #MESSAGEEQ    #WINEQ    #DATATYPEEQ 	   #REQUESTEQ 
   #FILEEQ    #ERRHANDLEREQ    #COMMEQ    #OPEQ                                                           
    #GROUPNEQ    #INFONEQ    #MESSAGENEQ    #WINNEQ    #DATATYPENEQ    #REQUESTNEQ    #FILENEQ    #ERRHANDLERNEQ    #COMMNEQ    #OPNEQ                      @              A                '              B      #COMM_WORLD    #ME    #NPROC    #ME_G0    #COMM_ATOM    #NPROC_ATOM    #MY_NATOM     #MY_ATMTAB !   #PARAL_PERT "   #COMM_PERT #   #COMM_CELL_PERT $   #ME_PERT %   #NPROC_PERT &   #DISTRB_PERT '   #PARAL_IMG (   #MY_NIMAGE )   #COMM_IMG *   #ME_IMG +   #NPROC_IMG ,   #DISTRB_IMG -   #MY_IMGTAB .   #COMM_CELL /   #ME_CELL 0   #NPROC_CELL 1   #COMM_FFT 2   #ME_FFT 3   #NPROC_FFT 4   #DISTRIBFFT 5   #PARALBD E   #COMM_BAND F   #ME_BAND G   #NPROC_BAND H   #PARAL_SPINOR I   #COMM_SPINOR J   #ME_SPINOR K   #NPROC_SPINOR L   #COMM_KPT M   #ME_KPT N   #NPROC_KPT O   #PROC_DISTRB P   #MY_ISPPOLTAB Q   #PARAL_KGB R   #BANDPP S   #COMM_BANDSPINORFFT T   #COMM_BANDFFT U   #COMM_KPTBAND V   #COMM_SPINORFFT W   #COMM_BANDSPINOR X   #MY_KGTAB Y   #MY_KPTTAB Z   #PW_UNBAL_THRESH [   #KPTDSTRB \   #KPT_LOC2FBZ_SP ]   #KPT_LOC2IBZ_SP ^   #MKMEM _   #PARAL_HF `   #COMM_HF a   #ME_HF b   #NPROC_HF c   #DISTRB_HF d   #COMM_WVL e   #ME_WVL f   #NPROC_WVL g   #NSCATTERARR h   #NGATHERARR i   #NGFFT3_IONIC j                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                !                                         &                                                                   °              y                                                                                           "     h       	                                                  #     l       
                                                  $     p                                                         %     t                                                         &     x                                                       '                                         &                                                                                       (     È                                                         )     Ì                                                         *     Ð                                                         +     Ô                                                         ,     Ø                                                       -            à                             &                                                                                     .            (                            &                                                                                       /     p                                                        0     t                                                        1     x                                                        2     |                                                        3                                                             4                                                           5     h                  #DISTRIBFFT_TYPE 6                           °              y#DISTRIBFFT_TYPE 6                                                                    @              E           6     'h                   #N2_COARSE 7   #N2_FINE 8   #TAB_FFTWF2_DISTRIB 9   #TAB_FFTDP2_DISTRIB :   #TAB_FFTDP3_DISTRIB ;   #TAB_FFTWF2DG_DISTRIB <   #TAB_FFTDP2DG_DISTRIB =   #TAB_FFTDP3DG_DISTRIB >   #TAB_FFTWF2_LOCAL ?   #TAB_FFTDP2_LOCAL @   #TAB_FFTDP3_LOCAL A   #TAB_FFTWF2DG_LOCAL B   #TAB_FFTDP2DG_LOCAL C   #TAB_FFTDP3DG_LOCAL D                $                              7                                          °                                                      0                 $                              8                                         °                                                      0               $                              9                                         &                                                      $                              :            P                             &                                                      $                              ;                                         &                                                      $                              <            à                             &                                                      $                              =            (                            &                                                      $                              >            p                            &                                                      $                              ?            ¸             	               &                                                      $                              @                          
               &                                                      $                              A            H                            &                                                      $                              B                                        &                                                      $                              C            Ø                            &                                                      $                              D                                         &                                                                                       E                                                             F                                                             G                                                             H                                                              I            !                                                  J     ¤      "                                                  K     ¨      #                                                  L     ¬      $                                                  M     °      %                                                  N     ´      &                                                  O     ¸      '                                                P            À             (               &                   &                   &                                                                                       Q            8             )     p          p            p                                                                      R     @      *                                                  S     D      +                                                  T     H      ,                                                  U     L      -                                                  V     P      .                                                  W     T      /                                                  X     X      0                                                Y            `             1               &                   &                                                                                     Z            À             2               &                                                                                      [           3   
                                             \                         4               &                   &                   &                                                                                     ]                         5               &                   &                   &                                                                                     ^                          6               &                   &                   &                                                                                     _            x             7               &                                                                                       `     À      8                                                  a     Ä      9                                                  b     È      :                                                  c     Ì      ;                                                d            Ð             <               &                   &                   &                                                                                       e     H      =                                                  f     L      >                                                  g     P      ?                                               h            X             @              &                   &                                                                   °              y                                                                                        i            ¸             A              &                   &                                                                   °              y                                                                                           j           B                        @                          k     '8             	      #MGFFT l   #NFFT m   #MGFFTC n   #NFFTC o   #USEFINEGRID p   #COATOFIN q   #FINTOCOA r   #NGFFT s   #NGFFTC t                 $                              l                                 $                              m                                $                              n                                $                              o                                $                              p                              $                              q                                         &                                                       $                              r            `                             &                                                         $                              s            ¨                   p          p            p                                        $                              t            ð              	     p          p            p                                           @                          u     '                    #MPI_VAL v                                               v                                    @                          w     '                    #MPI_VAL x                                               x                                    @                          y     '                    #MPI_VAL z                                               z                                    @                          {     '                    #MPI_VAL |                                               |                                    @                          }     '                    #MPI_VAL ~                                               ~                                    @                               '                    #MPI_VAL                                                                                    @                               '                    #MPI_VAL                                                                                    @                               '                    #MPI_VAL                                                                                    @                               '                    #MPI_VAL                                                                                    @                               '                    #MPI_VAL                                                                                     @                               '                                                                                                                              %         @                                                          #LHS    #RHS              
                                                     #MPI_GROUP u             
                                                     #MPI_GROUP u   %         @                                                          #LHS    #RHS              
                                                     #MPI_INFO w             
                                                     #MPI_INFO w   %         @                                                          #LHS    #RHS              
                                                     #MPI_MESSAGE y             
                                                     #MPI_MESSAGE y   %         @                                                          #LHS    #RHS              
                                                     #MPI_WIN {             
                                                     #MPI_WIN {   %         @                               	                           #LHS    #RHS              
                                                     #MPI_DATATYPE }             
                                                     #MPI_DATATYPE }   %         @                               
                           #LHS    #RHS              
                                                     #MPI_REQUEST              
                                                     #MPI_REQUEST    %         @                                                          #LHS    #RHS              
                                                     #MPI_FILE              
                                                     #MPI_FILE    %         @                                                          #LHS    #RHS              
                                                     #MPI_ERRHANDLER              
                                                     #MPI_ERRHANDLER    %         @                                                          #LHS    #RHS              
                                                     #MPI_COMM              
                                                     #MPI_COMM    %         @                                                          #LHS    #RHS              
                                                     #MPI_OP              
                                                     #MPI_OP    %         @                                                          #LHS    #RHS               
                                                     #MPI_GROUP u             
                                                      #MPI_GROUP u   %         @                                                          #LHS ¡   #RHS ¢             
                                  ¡                   #MPI_INFO w             
                                  ¢                   #MPI_INFO w   %         @                                                          #LHS £   #RHS ¤             
                                  £                   #MPI_MESSAGE y             
                                  ¤                   #MPI_MESSAGE y   %         @                                                          #LHS ¥   #RHS ¦             
                                  ¥                   #MPI_WIN {             
                                  ¦                   #MPI_WIN {   %         @                                                          #LHS §   #RHS ¨             
                                  §                   #MPI_DATATYPE }             
                                  ¨                   #MPI_DATATYPE }   %         @                                                          #LHS ©   #RHS ª             
                                  ©                   #MPI_REQUEST              
                                  ª                   #MPI_REQUEST    %         @                                                          #LHS «   #RHS ¬             
                                  «                   #MPI_FILE              
                                  ¬                   #MPI_FILE    %         @                                                          #LHS ­   #RHS ®             
                                  ­                   #MPI_ERRHANDLER              
                                  ®                   #MPI_ERRHANDLER    %         @                                                          #LHS ¯   #RHS °             
                                  ¯                   #MPI_COMM              
                                  °                   #MPI_COMM    %         @                                                          #LHS ±   #RHS ²             
                                  ±                   #MPI_OP              
                                  ²                   #MPI_OP                      @                          ³     '                    #X ´   #Y µ   #Z ¶                                               ´                                                               µ                                                              ¶                                    @              @           ·     'È                    #UCVOL ¸   #GCART ¹   #TR º   #RMET »                                              ¸                
                                             ¹                                         &                   &                                                                                      º            h                 
  p          p            p                                                                     »     	                        
  p          p          p            p          p                                            @              @           ¼     'X                   #LMNMAX ½   #INTLEN ¾   #NPSP ¿   #NLPSP À   #PSPINFO Á   #INDLMN Â   #RADII Ã   #PROJEC Ä   #MAT_EXP_PSP_NL Å   #EIVAL Æ   #EIVEC Ç                                               ½                                                               ¾                                                              ¿                                                              À                                                            Á                                         &                   &                                                                                     Â            p                             &                   &                   &                                                                                    Ã            è                 
            &                   &                                                                                    Ä            H                
            &                   &                   &                                                                                    Å            À             	   
            &                   &                   &                   &                                                                                    Æ            P             
   
            &                   &                   &                                                                                    Ç            È                
            &                   &                   &                   &                                                             @              @           È     '¸                    #MIN_PT É   #MAX_PT Ê   #NTRANCHE Ë   #NPT Ì   #DISPLS Í   #VCOUNT Î   #PT0 Ï   #PT1 Ð                                               É                                                               Ê                                                              Ë                                                              Ì                                                            Í                                         &                                                                                     Î            X                             &                                                                                       Ï                           #VEC_INT ³                                               Ð            ¬              #VEC_INT ³                     @                          Ñ     '                   #NPTREC Ò   #MAP Ó   #PAR Ô                                               Ò                                                             Ó                                         &                                                                                       Ô     ¸       P              #RECPARALL_TYPE È                     @               Á           Õ     'ð                   #QUITREC Ö   #MIN_NREC ×   #NFFTREC Ø   #NGPU Ù   #GPUDEVICE Ú   #LOAD Û   #TP Ü   #DEBUG Ý   #TRONC Þ   #NGFFTREC ß   #EFERMI à   #ZT_P á   #PAR â   #PAWFGR ã   #MPI ä   #NL å   #INF æ                                               Ö                                                               ×                                                              Ø                                                              Ù                                                              Ú                                                              Û                                                              Ü                                                              Ý                                                              Þ             	                                                  ß            $              
     p          p            p                                                                     à     p          
                                            á            x                 
            &                   &                                                                                       â     ¸       Ø              #RECPARALL_TYPE È                                               ã     8                   #PAWFGR_TYPE k                                              ä            È             #MPI_TYPE                                                å     X      Ð             #NLPSPREC_TYPE ¼                                               æ     È       (             #METRICREC_TYPE ·          (      fn#fn    È   @   J   DEFS_BASIS      @   J   M_ABICORE    H  I   J  DEFS_ABITYPES      L   J  M_PAWFGR #   Ý  Æ      i@+MPI_CONSTANTS #   £  Ð      i@+MPI_CONSTANTS '   s  A      MPI_TYPE+DEFS_ABITYPES 2   ´  H   a   MPI_TYPE%COMM_WORLD+DEFS_ABITYPES *   ü  H   a   MPI_TYPE%ME+DEFS_ABITYPES -   D  H   a   MPI_TYPE%NPROC+DEFS_ABITYPES -     H   a   MPI_TYPE%ME_G0+DEFS_ABITYPES 1   Ô  H   a   MPI_TYPE%COMM_ATOM+DEFS_ABITYPES 2   	  H   a   MPI_TYPE%NPROC_ATOM+DEFS_ABITYPES 0   d	  H   a   MPI_TYPE%MY_NATOM+DEFS_ABITYPES 1   ¬	  ô   a   MPI_TYPE%MY_ATMTAB+DEFS_ABITYPES 2    
  H   a   MPI_TYPE%PARAL_PERT+DEFS_ABITYPES 1   è
  H   a   MPI_TYPE%COMM_PERT+DEFS_ABITYPES 6   0  H   a   MPI_TYPE%COMM_CELL_PERT+DEFS_ABITYPES /   x  H   a   MPI_TYPE%ME_PERT+DEFS_ABITYPES 2   À  H   a   MPI_TYPE%NPROC_PERT+DEFS_ABITYPES 3        a   MPI_TYPE%DISTRB_PERT+DEFS_ABITYPES 1     H   a   MPI_TYPE%PARAL_IMG+DEFS_ABITYPES 1   ä  H   a   MPI_TYPE%MY_NIMAGE+DEFS_ABITYPES 0   ,  H   a   MPI_TYPE%COMM_IMG+DEFS_ABITYPES .   t  H   a   MPI_TYPE%ME_IMG+DEFS_ABITYPES 1   ¼  H   a   MPI_TYPE%NPROC_IMG+DEFS_ABITYPES 2        a   MPI_TYPE%DISTRB_IMG+DEFS_ABITYPES 1        a   MPI_TYPE%MY_IMGTAB+DEFS_ABITYPES 1   ,  H   a   MPI_TYPE%COMM_CELL+DEFS_ABITYPES /   t  H   a   MPI_TYPE%ME_CELL+DEFS_ABITYPES 2   ¼  H   a   MPI_TYPE%NPROC_CELL+DEFS_ABITYPES 0     H   a   MPI_TYPE%COMM_FFT+DEFS_ABITYPES .   L  H   a   MPI_TYPE%ME_FFT+DEFS_ABITYPES 1     H   a   MPI_TYPE%NPROC_FFT+DEFS_ABITYPES 2   Ü  Ú   a   MPI_TYPE%DISTRIBFFT+DEFS_ABITYPES -   ¶        DISTRIBFFT_TYPE+M_DISTRIBFFT 7   B  ¥   a   DISTRIBFFT_TYPE%N2_COARSE+M_DISTRIBFFT 5   ç  ¥   a   DISTRIBFFT_TYPE%N2_FINE+M_DISTRIBFFT @        a   DISTRIBFFT_TYPE%TAB_FFTWF2_DISTRIB+M_DISTRIBFFT @         a   DISTRIBFFT_TYPE%TAB_FFTDP2_DISTRIB+M_DISTRIBFFT @   ´     a   DISTRIBFFT_TYPE%TAB_FFTDP3_DISTRIB+M_DISTRIBFFT B   H     a   DISTRIBFFT_TYPE%TAB_FFTWF2DG_DISTRIB+M_DISTRIBFFT B   Ü     a   DISTRIBFFT_TYPE%TAB_FFTDP2DG_DISTRIB+M_DISTRIBFFT B   p     a   DISTRIBFFT_TYPE%TAB_FFTDP3DG_DISTRIB+M_DISTRIBFFT >        a   DISTRIBFFT_TYPE%TAB_FFTWF2_LOCAL+M_DISTRIBFFT >        a   DISTRIBFFT_TYPE%TAB_FFTDP2_LOCAL+M_DISTRIBFFT >   ,     a   DISTRIBFFT_TYPE%TAB_FFTDP3_LOCAL+M_DISTRIBFFT @   À     a   DISTRIBFFT_TYPE%TAB_FFTWF2DG_LOCAL+M_DISTRIBFFT @   T     a   DISTRIBFFT_TYPE%TAB_FFTDP2DG_LOCAL+M_DISTRIBFFT @   è     a   DISTRIBFFT_TYPE%TAB_FFTDP3DG_LOCAL+M_DISTRIBFFT /   |  H   a   MPI_TYPE%PARALBD+DEFS_ABITYPES 1   Ä  H   a   MPI_TYPE%COMM_BAND+DEFS_ABITYPES /     H   a   MPI_TYPE%ME_BAND+DEFS_ABITYPES 2   T  H   a   MPI_TYPE%NPROC_BAND+DEFS_ABITYPES 4     H   a   MPI_TYPE%PARAL_SPINOR+DEFS_ABITYPES 3   ä  H   a   MPI_TYPE%COMM_SPINOR+DEFS_ABITYPES 1   ,  H   a   MPI_TYPE%ME_SPINOR+DEFS_ABITYPES 4   t  H   a   MPI_TYPE%NPROC_SPINOR+DEFS_ABITYPES 0   ¼  H   a   MPI_TYPE%COMM_KPT+DEFS_ABITYPES .     H   a   MPI_TYPE%ME_KPT+DEFS_ABITYPES 1   L  H   a   MPI_TYPE%NPROC_KPT+DEFS_ABITYPES 3     Ä   a   MPI_TYPE%PROC_DISTRB+DEFS_ABITYPES 4   X     a   MPI_TYPE%MY_ISPPOLTAB+DEFS_ABITYPES 1   ô  H   a   MPI_TYPE%PARAL_KGB+DEFS_ABITYPES .   <   H   a   MPI_TYPE%BANDPP+DEFS_ABITYPES :      H   a   MPI_TYPE%COMM_BANDSPINORFFT+DEFS_ABITYPES 4   Ì   H   a   MPI_TYPE%COMM_BANDFFT+DEFS_ABITYPES 4   !  H   a   MPI_TYPE%COMM_KPTBAND+DEFS_ABITYPES 6   \!  H   a   MPI_TYPE%COMM_SPINORFFT+DEFS_ABITYPES 7   ¤!  H   a   MPI_TYPE%COMM_BANDSPINOR+DEFS_ABITYPES 0   ì!  ¬   a   MPI_TYPE%MY_KGTAB+DEFS_ABITYPES 1   "     a   MPI_TYPE%MY_KPTTAB+DEFS_ABITYPES 7   ,#  H   a   MPI_TYPE%PW_UNBAL_THRESH+DEFS_ABITYPES 0   t#  Ä   a   MPI_TYPE%KPTDSTRB+DEFS_ABITYPES 6   8$  Ä   a   MPI_TYPE%KPT_LOC2FBZ_SP+DEFS_ABITYPES 6   ü$  Ä   a   MPI_TYPE%KPT_LOC2IBZ_SP+DEFS_ABITYPES -   À%     a   MPI_TYPE%MKMEM+DEFS_ABITYPES 0   T&  H   a   MPI_TYPE%PARAL_HF+DEFS_ABITYPES /   &  H   a   MPI_TYPE%COMM_HF+DEFS_ABITYPES -   ä&  H   a   MPI_TYPE%ME_HF+DEFS_ABITYPES 0   ,'  H   a   MPI_TYPE%NPROC_HF+DEFS_ABITYPES 1   t'  Ä   a   MPI_TYPE%DISTRB_HF+DEFS_ABITYPES 0   8(  H   a   MPI_TYPE%COMM_WVL+DEFS_ABITYPES .   (  H   a   MPI_TYPE%ME_WVL+DEFS_ABITYPES 1   È(  H   a   MPI_TYPE%NPROC_WVL+DEFS_ABITYPES 3   )    a   MPI_TYPE%NSCATTERARR+DEFS_ABITYPES 2   *    a   MPI_TYPE%NGATHERARR+DEFS_ABITYPES 4   (+  H   a   MPI_TYPE%NGFFT3_IONIC+DEFS_ABITYPES %   p+  À       PAWFGR_TYPE+M_PAWFGR +   0,  H   a   PAWFGR_TYPE%MGFFT+M_PAWFGR *   x,  H   a   PAWFGR_TYPE%NFFT+M_PAWFGR ,   À,  H   a   PAWFGR_TYPE%MGFFTC+M_PAWFGR +   -  H   a   PAWFGR_TYPE%NFFTC+M_PAWFGR 1   P-  H   a   PAWFGR_TYPE%USEFINEGRID+M_PAWFGR .   -     a   PAWFGR_TYPE%COATOFIN+M_PAWFGR .   ,.     a   PAWFGR_TYPE%FINTOCOA+M_PAWFGR +   À.     a   PAWFGR_TYPE%NGFFT+M_PAWFGR ,   \/     a   PAWFGR_TYPE%NGFFTC+M_PAWFGR (   ø/  ]       MPI_GROUP+MPI_CONSTANTS 0   U0  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS '   0  ]       MPI_INFO+MPI_CONSTANTS /   ú0  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS *   B1  ]       MPI_MESSAGE+MPI_CONSTANTS 2   1  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS &   ç1  ]       MPI_WIN+MPI_CONSTANTS .   D2  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS +   2  ]       MPI_DATATYPE+MPI_CONSTANTS 3   é2  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS *   13  ]       MPI_REQUEST+MPI_CONSTANTS 2   3  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS '   Ö3  ]       MPI_FILE+MPI_CONSTANTS /   34  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS -   {4  ]       MPI_ERRHANDLER+MPI_CONSTANTS 5   Ø4  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '    5  ]       MPI_COMM+MPI_CONSTANTS /   }5  H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS %   Å5  ]       MPI_OP+MPI_CONSTANTS -   "6  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS '   j6  P       #UNLPOLY+ISO_C_BINDING    º6  p       DP+DEFS_BASIS &   *7  b       GROUPEQ+MPI_CONSTANTS *   7  W   a   GROUPEQ%LHS+MPI_CONSTANTS *   ã7  W   a   GROUPEQ%RHS+MPI_CONSTANTS %   :8  b       INFOEQ+MPI_CONSTANTS )   8  V   a   INFOEQ%LHS+MPI_CONSTANTS )   ò8  V   a   INFOEQ%RHS+MPI_CONSTANTS (   H9  b       MESSAGEEQ+MPI_CONSTANTS ,   ª9  Y   a   MESSAGEEQ%LHS+MPI_CONSTANTS ,   :  Y   a   MESSAGEEQ%RHS+MPI_CONSTANTS $   \:  b       WINEQ+MPI_CONSTANTS (   ¾:  U   a   WINEQ%LHS+MPI_CONSTANTS (   ;  U   a   WINEQ%RHS+MPI_CONSTANTS )   h;  b       DATATYPEEQ+MPI_CONSTANTS -   Ê;  Z   a   DATATYPEEQ%LHS+MPI_CONSTANTS -   $<  Z   a   DATATYPEEQ%RHS+MPI_CONSTANTS (   ~<  b       REQUESTEQ+MPI_CONSTANTS ,   à<  Y   a   REQUESTEQ%LHS+MPI_CONSTANTS ,   9=  Y   a   REQUESTEQ%RHS+MPI_CONSTANTS %   =  b       FILEEQ+MPI_CONSTANTS )   ô=  V   a   FILEEQ%LHS+MPI_CONSTANTS )   J>  V   a   FILEEQ%RHS+MPI_CONSTANTS +    >  b       ERRHANDLEREQ+MPI_CONSTANTS /   ?  \   a   ERRHANDLEREQ%LHS+MPI_CONSTANTS /   ^?  \   a   ERRHANDLEREQ%RHS+MPI_CONSTANTS %   º?  b       COMMEQ+MPI_CONSTANTS )   @  V   a   COMMEQ%LHS+MPI_CONSTANTS )   r@  V   a   COMMEQ%RHS+MPI_CONSTANTS #   È@  b       OPEQ+MPI_CONSTANTS '   *A  T   a   OPEQ%LHS+MPI_CONSTANTS '   ~A  T   a   OPEQ%RHS+MPI_CONSTANTS '   ÒA  b       GROUPNEQ+MPI_CONSTANTS +   4B  W   a   GROUPNEQ%LHS+MPI_CONSTANTS +   B  W   a   GROUPNEQ%RHS+MPI_CONSTANTS &   âB  b       INFONEQ+MPI_CONSTANTS *   DC  V   a   INFONEQ%LHS+MPI_CONSTANTS *   C  V   a   INFONEQ%RHS+MPI_CONSTANTS )   ðC  b       MESSAGENEQ+MPI_CONSTANTS -   RD  Y   a   MESSAGENEQ%LHS+MPI_CONSTANTS -   «D  Y   a   MESSAGENEQ%RHS+MPI_CONSTANTS %   E  b       WINNEQ+MPI_CONSTANTS )   fE  U   a   WINNEQ%LHS+MPI_CONSTANTS )   »E  U   a   WINNEQ%RHS+MPI_CONSTANTS *   F  b       DATATYPENEQ+MPI_CONSTANTS .   rF  Z   a   DATATYPENEQ%LHS+MPI_CONSTANTS .   ÌF  Z   a   DATATYPENEQ%RHS+MPI_CONSTANTS )   &G  b       REQUESTNEQ+MPI_CONSTANTS -   G  Y   a   REQUESTNEQ%LHS+MPI_CONSTANTS -   áG  Y   a   REQUESTNEQ%RHS+MPI_CONSTANTS &   :H  b       FILENEQ+MPI_CONSTANTS *   H  V   a   FILENEQ%LHS+MPI_CONSTANTS *   òH  V   a   FILENEQ%RHS+MPI_CONSTANTS ,   HI  b       ERRHANDLERNEQ+MPI_CONSTANTS 0   ªI  \   a   ERRHANDLERNEQ%LHS+MPI_CONSTANTS 0   J  \   a   ERRHANDLERNEQ%RHS+MPI_CONSTANTS &   bJ  b       COMMNEQ+MPI_CONSTANTS *   ÄJ  V   a   COMMNEQ%LHS+MPI_CONSTANTS *   K  V   a   COMMNEQ%RHS+MPI_CONSTANTS $   pK  b       OPNEQ+MPI_CONSTANTS (   ÒK  T   a   OPNEQ%LHS+MPI_CONSTANTS (   &L  T   a   OPNEQ%RHS+MPI_CONSTANTS    zL  e       VEC_INT    ßL  H   a   VEC_INT%X    'M  H   a   VEC_INT%Y    oM  H   a   VEC_INT%Z    ·M  x       METRICREC_TYPE %   /N  H   a   METRICREC_TYPE%UCVOL %   wN  ¬   a   METRICREC_TYPE%GCART "   #O     a   METRICREC_TYPE%TR $   ¿O  ¼   a   METRICREC_TYPE%RMET    {P  ×       NLPSPREC_TYPE %   RQ  H   a   NLPSPREC_TYPE%LMNMAX %   Q  H   a   NLPSPREC_TYPE%INTLEN #   âQ  H   a   NLPSPREC_TYPE%NPSP $   *R  H   a   NLPSPREC_TYPE%NLPSP &   rR  ¬   a   NLPSPREC_TYPE%PSPINFO %   S  Ä   a   NLPSPREC_TYPE%INDLMN $   âS  ¬   a   NLPSPREC_TYPE%RADII %   T  Ä   a   NLPSPREC_TYPE%PROJEC -   RU  Ü   a   NLPSPREC_TYPE%MAT_EXP_PSP_NL $   .V  Ä   a   NLPSPREC_TYPE%EIVAL $   òV  Ü   a   NLPSPREC_TYPE%EIVEC    ÎW  ©       RECPARALL_TYPE &   wX  H   a   RECPARALL_TYPE%MIN_PT &   ¿X  H   a   RECPARALL_TYPE%MAX_PT (   Y  H   a   RECPARALL_TYPE%NTRANCHE #   OY  H   a   RECPARALL_TYPE%NPT &   Y     a   RECPARALL_TYPE%DISPLS &   +Z     a   RECPARALL_TYPE%VCOUNT #   ¿Z  ]   a   RECPARALL_TYPE%PT0 #   [  ]   a   RECPARALL_TYPE%PT1    y[  n       RECGPU_TYPE #   ç[  H   a   RECGPU_TYPE%NPTREC     /\     a   RECGPU_TYPE%MAP     Ã\  d   a   RECGPU_TYPE%PAR    ']        RECURSION_TYPE '   3^  H   a   RECURSION_TYPE%QUITREC (   {^  H   a   RECURSION_TYPE%MIN_NREC '   Ã^  H   a   RECURSION_TYPE%NFFTREC $   _  H   a   RECURSION_TYPE%NGPU )   S_  H   a   RECURSION_TYPE%GPUDEVICE $   _  H   a   RECURSION_TYPE%LOAD "   ã_  H   a   RECURSION_TYPE%TP %   +`  H   a   RECURSION_TYPE%DEBUG %   s`  H   a   RECURSION_TYPE%TRONC (   »`     a   RECURSION_TYPE%NGFFTREC &   Wa  H   a   RECURSION_TYPE%EFERMI $   a  ¬   a   RECURSION_TYPE%ZT_P #   Kb  d   a   RECURSION_TYPE%PAR &   ¯b  a   a   RECURSION_TYPE%PAWFGR #   c  ^   a   RECURSION_TYPE%MPI "   nc  c   a   RECURSION_TYPE%NL #   Ñc  d   a   RECURSION_TYPE%INF 