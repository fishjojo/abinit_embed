  A  �   k820309    l          18.0        ��[                                                                                                          
       m_pred_diisrelax.F90 M_PRED_DIISRELAX              PRED_DIISRELAX                                                     
                                                           
                            @                              
                            @                              
                                                           
                            @                              
       XCART2XRED XRED2XCART                      @                              
       HESSINIT HESSUPDT                   @                                '�             *      #DELAYPERM 	   #DIISMEMORY 
   #GOPRECON    #JELLSLAB    #NATOM    #NCONEQ    #PH_FREEZ_DISP_ADDSTRAIN    #PH_FREEZ_DISP_OPTION    #PH_FREEZ_DISP_NAMPL    #PH_NQSHIFT    #NNOS    #NSYM    #NTYPAT    #OPTCELL    #RESTARTXF    #SIGNPERM    #IONMOV    #BMASS    #DTION    #FRICTION    #MDWALL    #NOSEINERT    #STRPRECON    #VIS     #IATFIX !   #SYMAFM "   #SYMREL #   #TNONS $   #TYPAT %   #PRTATLIST &   #PH_NGQPT '   #PH_FREEZ_DISP_AMPL (   #PH_QSHIFT )   #AMU_CURR *   #AMASS +   #GOPRECPRM ,   #MDTEMP -   #STRTARGET .   #QMASS /   #ZNUCL 0   #FNAMEABI_HES 1   #FILNAM_DS 2                � $                              	                                � $                              
                               � $                                                             � $                                                            � $                                                             � $                                                             � $                                                             � $                                                             � $                                           	                   � $                                   $       
                   � $                                   (                          � $                                   ,                          � $                                   0                          � $                                   4                          � $                                   8                          � $                                   <                          � $                                   @                          � $                                  H          
                � $                                  P          
                � $                                  X          
                � $                                  `          
                � $                                  h          
                � $                                  p          
                � $                                   x          
               �$                              !            �                             &                   &                                                       �$                              "            �                             &                                                       �$                              #            (                            &                   &                   &                                                       �$                             $            �                
            &                   &                                                       �$                              %                                         &                                                       �$                              &            H                            &                                                       �$                              '            �                            &                                                       �$                             (            �                 
            &                   &                                                       �$                             )            8             !   
            &                   &                                                       �$                             *            �             "   
            &                                                       �$                             +            �             #   
            &                                                       �$                             ,            (             $   
            &                                                       �$                             -            p             %   
            &                                                       �$                             .            �             &   
            &                                                       �$                             /                          '   
            &                                                       �$                             0            H             (   
            &                                                        �$                             1           �      )       .            �$                             2            �            *               &                                                                     @               A           3     'H                   #IHIST 4   #MXHIST 5   #ISVUSED 6   #ISARUSED 7   #ACELL 8   #RPRIMD 9   #XRED :   #FCART ;   #STRTEN <   #VEL =   #VEL_CELL >   #ETOT ?   #EKIN @   #ENTROPY A   #TIME B               � $                              4                                          �                                                      0                � $                              5                                         �                                                      0                 � $                              6                               � $                              7                            � $                             8                             
            &                   &                                                     � $                             9            p                 
            &                   &                   &                                                     � $                             :            �                 
            &                   &                   &                                                     � $                             ;            `                
            &                   &                   &                                                     � $                             <            �             	   
            &                   &                                                     � $                             =            8             
   
            &                   &                   &                                                     � $                             >            �                
            &                   &                   &                                                     � $                             ?            (                
            &                                                     � $                             @            p                
            &                                                     � $                             A            �                
            &                                                     � $                             B                             
            &                                                          @ @                          C     '                    #MPI_VAL D                �                               D                                  @ @                          E     '                    #MPI_VAL F                �                               F                                  @ @                          G     '                    #MPI_VAL H                �                               H                                  @ @                          I     '                    #MPI_VAL J                �                               J                                  @ @                          K     '                    #MPI_VAL L                �                               L                                  @ @                          M     '                    #MPI_VAL N                �                               N                                  @ @                          O     '                    #MPI_VAL P                �                               P                                  @ @                          Q     '                    #MPI_VAL R                �                               R                                  @ @                          S     '                    #MPI_VAL T                �                               T                                  @ @                          U     '                    #MPI_VAL V                �                               V                                    @                          W     '                        #         @                                  X                    #NATOM Y   #RPRIMD Z   #XCART [   #XRED \             
                                  Y                     
                                 Z     	              
 N   p          p          p            p          p                                   
                                 [                    
 O   p          p          5 O p            p          5 O p                                                                   \                    
 P    p          p          5 O p            p          5 O p                          #         @                                  ]                    #NATOM ^   #RPRIMD _   #XCART `   #XRED a             
                                  ^                     
                                 _     	              
 R   p          p          p            p          p                                                                   `                    
 T    p          p          5 O p            p          5 O p                                   
                                 a                    
 S   p          p          5 O p            p          5 O p                          #         @                                  b                    #AB_MOVER c   #HESSIN d   #INIT_MATRIX e   #NDIM f   #UCVOL g             
                                  c     �             #ABIMOVER                                             d                    
       p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                    
                                 e     	              
    p          p          p            p          p                                    
                                  f                     
                                 g     
      #         @                                  h                 	   #HESSIN i   #IATFIX j   #NATOM k   #NDIM l   #VIN m   #VIN_PREV n   #VOUT o   #VOUT_PREV p   #NIMAGE q            
                                i                    
       p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                   
                                  j                        p          p          5 O p            p          5 O p                                    
                                  k                     
                                  l                    
                                 m                    
    p          5 O p            5 O p                                   
                                 n                    
    p          5 O p            5 O p                                   
                                 o                    
    p          5 O p            5 O p                                   
                                 p                    
    p          5 O p            5 O p                                    
                                 q                                                      r            #         @                                  s                    #ACELL t   #HIST u   #NATOM v   #RPRIMD w   #XRED x   #ZDEBUG y                                             t                   
     p          p            p                                    
                                  u     H             #ABIHIST 3             
                                  v                                                     w     	              
     p          p          p            p          p                                                                   x                    
     p          p          5 O p            p          5 O p                                    
                                  y           %         @                                z                           #HIST {   #STEP |             
                                  {     H             #ABIHIST 3             
                                  |           #         @                                  }                    #ACELL ~   #HIST    #NATOM �   #RPRIMD �   #XRED �   #ZDEBUG �             
                                 ~                   
    p          p            p                                    
                                      H              #ABIHIST 3             
                                  �                     
                                 �     	              
    p          p          p            p          p                                   
                                 �                    
    p          p          5 O p            p          5 O p                                    
                                  �           #         @                                   �                    #AB_MOVER �   #HIST �   #ITIME �   #NTIME �   #ZDEBUG �   #IEXIT �            
  @                               �     �             #ABIMOVER              
D `                               �     H              #ABIHIST 3             
                                  �                     
                                  �                     
  @                               �                     
                                  �              �   .      fn#fn &   �      b   uapp(M_PRED_DIISRELAX    �   @   J  DEFS_BASIS    -  @   J  M_ABICORE    m  @   J  M_ABIMOVER    �  @   J  M_ABIHIST $   �  @   J  M_LINALG_INTERFACES    -  V   J  M_GEOMETRY    �  R   J  M_BFGS $   �  �      ABIMOVER+M_ABIMOVER .   |  H   a   ABIMOVER%DELAYPERM+M_ABIMOVER /   �  H   a   ABIMOVER%DIISMEMORY+M_ABIMOVER -     H   a   ABIMOVER%GOPRECON+M_ABIMOVER -   T  H   a   ABIMOVER%JELLSLAB+M_ABIMOVER *   �  H   a   ABIMOVER%NATOM+M_ABIMOVER +   �  H   a   ABIMOVER%NCONEQ+M_ABIMOVER <   ,  H   a   ABIMOVER%PH_FREEZ_DISP_ADDSTRAIN+M_ABIMOVER 9   t  H   a   ABIMOVER%PH_FREEZ_DISP_OPTION+M_ABIMOVER 8   �  H   a   ABIMOVER%PH_FREEZ_DISP_NAMPL+M_ABIMOVER /     H   a   ABIMOVER%PH_NQSHIFT+M_ABIMOVER )   L  H   a   ABIMOVER%NNOS+M_ABIMOVER )   �  H   a   ABIMOVER%NSYM+M_ABIMOVER +   �  H   a   ABIMOVER%NTYPAT+M_ABIMOVER ,   $	  H   a   ABIMOVER%OPTCELL+M_ABIMOVER .   l	  H   a   ABIMOVER%RESTARTXF+M_ABIMOVER -   �	  H   a   ABIMOVER%SIGNPERM+M_ABIMOVER +   �	  H   a   ABIMOVER%IONMOV+M_ABIMOVER *   D
  H   a   ABIMOVER%BMASS+M_ABIMOVER *   �
  H   a   ABIMOVER%DTION+M_ABIMOVER -   �
  H   a   ABIMOVER%FRICTION+M_ABIMOVER +     H   a   ABIMOVER%MDWALL+M_ABIMOVER .   d  H   a   ABIMOVER%NOSEINERT+M_ABIMOVER .   �  H   a   ABIMOVER%STRPRECON+M_ABIMOVER (   �  H   a   ABIMOVER%VIS+M_ABIMOVER +   <  �   a   ABIMOVER%IATFIX+M_ABIMOVER +   �  �   a   ABIMOVER%SYMAFM+M_ABIMOVER +   |  �   a   ABIMOVER%SYMREL+M_ABIMOVER *   @  �   a   ABIMOVER%TNONS+M_ABIMOVER *   �  �   a   ABIMOVER%TYPAT+M_ABIMOVER .   �  �   a   ABIMOVER%PRTATLIST+M_ABIMOVER -     �   a   ABIMOVER%PH_NGQPT+M_ABIMOVER 7   �  �   a   ABIMOVER%PH_FREEZ_DISP_AMPL+M_ABIMOVER .   T  �   a   ABIMOVER%PH_QSHIFT+M_ABIMOVER -      �   a   ABIMOVER%AMU_CURR+M_ABIMOVER *   �  �   a   ABIMOVER%AMASS+M_ABIMOVER .   (  �   a   ABIMOVER%GOPRECPRM+M_ABIMOVER +   �  �   a   ABIMOVER%MDTEMP+M_ABIMOVER .   P  �   a   ABIMOVER%STRTARGET+M_ABIMOVER *   �  �   a   ABIMOVER%QMASS+M_ABIMOVER *   x  �   a   ABIMOVER%ZNUCL+M_ABIMOVER 1     P   a   ABIMOVER%FNAMEABI_HES+M_ABIMOVER .   \  �   a   ABIMOVER%FILNAM_DS+M_ABIMOVER "   �  �       ABIHIST+M_ABIHIST (   �  �   a   ABIHIST%IHIST+M_ABIHIST )   �  �   a   ABIHIST%MXHIST+M_ABIHIST *   >  H   a   ABIHIST%ISVUSED+M_ABIHIST +   �  H   a   ABIHIST%ISARUSED+M_ABIHIST (   �  �   a   ABIHIST%ACELL+M_ABIHIST )   z  �   a   ABIHIST%RPRIMD+M_ABIHIST '   >  �   a   ABIHIST%XRED+M_ABIHIST (     �   a   ABIHIST%FCART+M_ABIHIST )   �  �   a   ABIHIST%STRTEN+M_ABIHIST &   r  �   a   ABIHIST%VEL+M_ABIHIST +   6  �   a   ABIHIST%VEL_CELL+M_ABIHIST '   �  �   a   ABIHIST%ETOT+M_ABIHIST '   �  �   a   ABIHIST%EKIN+M_ABIHIST *   "   �   a   ABIHIST%ENTROPY+M_ABIHIST '   �   �   a   ABIHIST%TIME+M_ABIHIST (   J!  ]      MPI_GROUP+MPI_CONSTANTS 0   �!  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS '   �!  ]      MPI_INFO+MPI_CONSTANTS /   L"  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS *   �"  ]      MPI_MESSAGE+MPI_CONSTANTS 2   �"  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS &   9#  ]      MPI_WIN+MPI_CONSTANTS .   �#  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS +   �#  ]      MPI_DATATYPE+MPI_CONSTANTS 3   ;$  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS *   �$  ]      MPI_REQUEST+MPI_CONSTANTS 2   �$  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS '   (%  ]      MPI_FILE+MPI_CONSTANTS /   �%  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS -   �%  ]      MPI_ERRHANDLER+MPI_CONSTANTS 5   *&  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '   r&  ]      MPI_COMM+MPI_CONSTANTS /   �&  H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS %   '  ]      MPI_OP+MPI_CONSTANTS -   t'  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS '   �'  P       #UNLPOLY+ISO_C_BINDING &   (  t       XCART2XRED+M_GEOMETRY ,   �(  @   a   XCART2XRED%NATOM+M_GEOMETRY -   �(  �   a   XCART2XRED%RPRIMD+M_GEOMETRY ,   t)  �   a   XCART2XRED%XCART+M_GEOMETRY +   8*  �   a   XCART2XRED%XRED+M_GEOMETRY &   �*  t       XRED2XCART+M_GEOMETRY ,   p+  @   a   XRED2XCART%NATOM+M_GEOMETRY -   �+  �   a   XRED2XCART%RPRIMD+M_GEOMETRY ,   d,  �   a   XRED2XCART%XCART+M_GEOMETRY +   (-  �   a   XRED2XCART%XRED+M_GEOMETRY     �-  �       HESSINIT+M_BFGS )   t.  V   a   HESSINIT%AB_MOVER+M_BFGS '   �.  �   a   HESSINIT%HESSIN+M_BFGS ,   �/  �   a   HESSINIT%INIT_MATRIX+M_BFGS %   z0  @   a   HESSINIT%NDIM+M_BFGS &   �0  @   a   HESSINIT%UCVOL+M_BFGS     �0  �       HESSUPDT+M_BFGS '   �1  �   a   HESSUPDT%HESSIN+M_BFGS '   �2  �   a   HESSUPDT%IATFIX+M_BFGS &   k3  @   a   HESSUPDT%NATOM+M_BFGS %   �3  @   a   HESSUPDT%NDIM+M_BFGS $   �3  �   a   HESSUPDT%VIN+M_BFGS )   �4  �   a   HESSUPDT%VIN_PREV+M_BFGS %   35  �   a   HESSUPDT%VOUT+M_BFGS *   �5  �   a   HESSUPDT%VOUT_PREV+M_BFGS '   {6  @   a   HESSUPDT%NIMAGE+M_BFGS #   �6  @       STD_OUT+DEFS_BASIS #   �6  �       HIST2VAR+M_ABIHIST )   �7  �   a   HIST2VAR%ACELL+M_ABIHIST (   8  U   a   HIST2VAR%HIST+M_ABIHIST )   n8  @   a   HIST2VAR%NATOM+M_ABIHIST *   �8  �   a   HIST2VAR%RPRIMD+M_ABIHIST (   b9  �   a   HIST2VAR%XRED+M_ABIHIST *   &:  @   a   HIST2VAR%ZDEBUG+M_ABIHIST ,   f:  d       ABIHIST_FINDINDEX+M_ABIHIST 1   �:  U   a   ABIHIST_FINDINDEX%HIST+M_ABIHIST 1   ;  @   a   ABIHIST_FINDINDEX%STEP+M_ABIHIST #   _;  �       VAR2HIST+M_ABIHIST )   �;  �   a   VAR2HIST%ACELL+M_ABIHIST (   }<  U   a   VAR2HIST%HIST+M_ABIHIST )   �<  @   a   VAR2HIST%NATOM+M_ABIHIST *   =  �   a   VAR2HIST%RPRIMD+M_ABIHIST (   �=  �   a   VAR2HIST%XRED+M_ABIHIST *   �>  @   a   VAR2HIST%ZDEBUG+M_ABIHIST    �>  �       PRED_DIISRELAX (   W?  V   a   PRED_DIISRELAX%AB_MOVER $   �?  U   a   PRED_DIISRELAX%HIST %   @  @   a   PRED_DIISRELAX%ITIME %   B@  @   a   PRED_DIISRELAX%NTIME &   �@  @   a   PRED_DIISRELAX%ZDEBUG %   �@  @   a   PRED_DIISRELAX%IEXIT 