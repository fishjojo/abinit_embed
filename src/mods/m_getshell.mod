  N'  Z   k820309    l          18.0        3��[                                                                                                          
       m_getshell.F90 M_GETSHELL              GETSHELL                                                     
                                                           
                            @                              
                            @                              
                                                           
                            @                              
       GETKGRID                @ @                               '                    #MPI_VAL                 �                                                                 @ @                          	     '                    #MPI_VAL 
                �                               
                                  @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                   @                               '                        #         @                                                      #CHKSYMBREAK    #IOUT    #ISCF    #KPT     #KPTOPT !   #KPTRLATT "   #KPTRLEN #   #MSYM $   #NKPT %   #NKPT_COMPUTED &   #NSHIFTK '   #NSYM (   #RPRIMD )   #SHIFTK *   #SYMAFM +   #SYMREL ,   #VACUUM -   #WTK .   #FULLBZ /   #NKPTHF 0   #KPTHF 1   #DOWNSAMPLING 2             
                                                       
                                                       
                                                      
                                                     
 3    p          p          5 O p 	           p          5 O p 	                                   
                                  !                     
                                 "     	               0    p          p          p            p          p                                                                    #     
                 
                                  $                     
                                  %                     
                                 &                      
                                 '                      
                                  (                     
                                 )     	              
 1   p          p          p            p          p                                    
                                *     v             
 2    p          p          p �           p          p �                                  
                                  +                     ,   p          5 O p            5 O p                                   
                                  ,                     -   p          p          p          5 O p            p          p          5 O p                                    
                                  -                    .   p          p            p                                   
                                .                    
 4    p          5 O p 	           5 O p 	                                                                 /                   
 5              &                   &                                                     
                                 0                                                      1                   
 6              &                   &                                                     
                                 2                    /   p          p            p                                                                       3                                                      1#         @                                  4                    #MESSAGE 5   #LEVEL 6   #MODE_PARAL 7   #FILE 8   #LINE 9   #NODUMP :   #NOSTOP ;             
                                5                    1           
                                6                    1           
                                7                    1           
                               8                    1           
                                 9                     
                                 :                     
                                 ;                      @                               <            #         @                                  =                    #UNIT >   #MSG ?   #MODE_PARAL @   #DO_FLUSH A             
                                  >                     
                                ?                    1           
                               @                    1           
                                 A                      @                               B            #         @                                  C                     n                         	              Cabi_cabort                    #         @                                   D                    #GMET E   #KNEIGH F   #KG_NEIGH H   #KPTINDEX I   #KPTOPT K   #KPTRLATT L   #KPT2 M   #KPT3 N   #MKMEM O   #MKMEM_MAX P   #MVWTK Q   #NKPT2 G   #NKPT3 J   #NNEIGH R   #NSHIFTK S   #RMET T   #RPRIMD U   #SHIFTK V   #WTK2 W   #COMM X             
                                 E     	              
    p          p          p            p          p                                   D                                 F                         p          p          5 � p        r G       p          5 � p        r G                              D                                 H                             p        5 � p        r G   p        p        p          p          5 � p        r G     p            p          5 � p        r G     p                                   D                                 I                         p          p          5 � p        r J       p          5 � p        r J                               
                                  K                     
D @                               L     	                   p          p          p            p          p                                   
                                 M                    
    p          p          5 � p        r G       p          5 � p        r G                              D @                              N                    
     p          p          5 � p        r J       p          5 � p        r J                               
                                  O                     D @                               P                     D                                Q                    
     p          p          5 � p        r G       p          5 � p        r G                               
                                  G                     
  @                               J                     D                                 R                      
D @                               S                      
                                 T     	              
    p          p          p            p          p                                    
  @                              U     	              
    p          p          p            p          p                                   
                                 V                    
 	   p          p          5 � p        r S       p          5 � p        r S                              
                                 W                    
 
   p          5 � p        r G       5 � p        r G                               
  @                               X              �   "      fn#fn     �      b   uapp(M_GETSHELL    �   @   J  DEFS_BASIS      @   J  M_ABICORE    [  @   J  M_XMPI    �  @   J  M_ERRORS $   �  @   J  M_LINALG_INTERFACES      I   J  M_KPTS (   d  ]      MPI_GROUP+MPI_CONSTANTS 0   �  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS '   	  ]      MPI_INFO+MPI_CONSTANTS /   f  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS *   �  ]      MPI_MESSAGE+MPI_CONSTANTS 2     H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS &   S  ]      MPI_WIN+MPI_CONSTANTS .   �  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS +   �  ]      MPI_DATATYPE+MPI_CONSTANTS 3   U  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS *   �  ]      MPI_REQUEST+MPI_CONSTANTS 2   �  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS '   B  ]      MPI_FILE+MPI_CONSTANTS /   �  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS -   �  ]      MPI_ERRHANDLER+MPI_CONSTANTS 5   D  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '   �  ]      MPI_COMM+MPI_CONSTANTS /   �  H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS %   1  ]      MPI_OP+MPI_CONSTANTS -   �  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS '   �  P       #UNLPOLY+ISO_C_BINDING     &	  U      GETKGRID+M_KPTS ,   {
  @   a   GETKGRID%CHKSYMBREAK+M_KPTS %   �
  @   a   GETKGRID%IOUT+M_KPTS %   �
  @   a   GETKGRID%ISCF+M_KPTS $   ;  �   a   GETKGRID%KPT+M_KPTS '   �  @   a   GETKGRID%KPTOPT+M_KPTS )   ?  �   a   GETKGRID%KPTRLATT+M_KPTS (   �  @   a   GETKGRID%KPTRLEN+M_KPTS %   3  @   a   GETKGRID%MSYM+M_KPTS %   s  @   a   GETKGRID%NKPT+M_KPTS .   �  @   a   GETKGRID%NKPT_COMPUTED+M_KPTS (   �  @   a   GETKGRID%NSHIFTK+M_KPTS %   3  @   a   GETKGRID%NSYM+M_KPTS '   s  �   a   GETKGRID%RPRIMD+M_KPTS '   '  �   a   GETKGRID%SHIFTK+M_KPTS '   �  �   a   GETKGRID%SYMAFM+M_KPTS '     �   a   GETKGRID%SYMREL+M_KPTS '   c  �   a   GETKGRID%VACUUM+M_KPTS $   �  �   a   GETKGRID%WTK+M_KPTS '   �  �   a   GETKGRID%FULLBZ+M_KPTS '   ?  @   a   GETKGRID%NKPTHF+M_KPTS &     �   a   GETKGRID%KPTHF+M_KPTS -   #  �   a   GETKGRID%DOWNSAMPLING+M_KPTS "   �  q       XMPI_PARAL+M_XMPI "   (  �       MSG_HNDL+M_ERRORS *   �  L   a   MSG_HNDL%MESSAGE+M_ERRORS (     L   a   MSG_HNDL%LEVEL+M_ERRORS -   \  L   a   MSG_HNDL%MODE_PARAL+M_ERRORS '   �  L   a   MSG_HNDL%FILE+M_ERRORS '   �  @   a   MSG_HNDL%LINE+M_ERRORS )   4  @   a   MSG_HNDL%NODUMP+M_ERRORS )   t  @   a   MSG_HNDL%NOSTOP+M_ERRORS #   �  @       STD_OUT+DEFS_BASIS $   �  y       WRTOUT+M_SPECIALMSG )   m  @   a   WRTOUT%UNIT+M_SPECIALMSG (   �  L   a   WRTOUT%MSG+M_SPECIALMSG /   �  L   a   WRTOUT%MODE_PARAL+M_SPECIALMSG -   E  @   a   WRTOUT%DO_FLUSH+M_SPECIALMSG "   �  @       AB_OUT+DEFS_BASIS $   �  �       ABI_CABORT+M_ERRORS    \  2      GETSHELL    �  �   a   GETSHELL%GMET     B  �   a   GETSHELL%KNEIGH "     D  a   GETSHELL%KG_NEIGH "   Z  �   a   GETSHELL%KPTINDEX     .  @   a   GETSHELL%KPTOPT "   n  �   a   GETSHELL%KPTRLATT    "   �   a   GETSHELL%KPT2    �   �   a   GETSHELL%KPT3    �!  @   a   GETSHELL%MKMEM #   
"  @   a   GETSHELL%MKMEM_MAX    J"  �   a   GETSHELL%MVWTK    #  @   a   GETSHELL%NKPT2    ^#  @   a   GETSHELL%NKPT3     �#  @   a   GETSHELL%NNEIGH !   �#  @   a   GETSHELL%NSHIFTK    $  �   a   GETSHELL%RMET     �$  �   a   GETSHELL%RPRIMD     �%  �   a   GETSHELL%SHIFTK    Z&  �   a   GETSHELL%WTK2    '  @   a   GETSHELL%COMM 