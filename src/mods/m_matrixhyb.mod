  4  =   k820309    l          18.0        ���[                                                                                                          
       m_MatrixHyb.F90 M_MATRIXHYB       	       MATRIXHYB_INIT MATRIXHYB_SETSIZE MATRIXHYB_CLEAR MATRIXHYB_ASSIGN MATRIXHYB_INVERSE MATRIXHYB_GETDET MATRIXHYB_PRINT MATRIXHYB_DESTROY MATRIXHYB                                                     
                      @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL 	                �                               	                                  @ @                          
     '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                                 @ @                               '                    #MPI_VAL                 �                                                  #         @                                                      #MESSAGE              
                                                    1 #         @                                                      #MESSAGE              
                                                    1 #         @                                                       n                         �              Cabi_cabort                                      @               A                'H                   #SIZE    #TAIL    #ITECH    #WMAX    #MAT     #MAT_TAU !   #MAT_OMEGA "                � D                                                              � $                                                            � D                                                                                                                ��������                         � D                                                           � $                                                            
            &                   &                                                      � $                              !            p                             &                   &                                                      � $                             "            �                             &                   &                   &                                           #         @                                  #                    #THIS $   #ITECH %   #SIZE &   #WMAX '             
D                                 $     H              #MATRIXHYB              
                                  %                     
 @                               &                     
 @                               '           #         @                                  (                    #THIS )   #NEW_TAIL *             
D @                               )     H              #MATRIXHYB              
                                  *           #         @                                   +                    #THIS ,             
D                                 ,     H              #MATRIXHYB    #         @                                   -                    #THIS .   #MATRIX /             
D @                               .     H              #MATRIXHYB              
                                  /     H             #MATRIXHYB    #         @                                   0                    #THIS 1   #DETERMINANT 2             
D @                               1     H              #MATRIXHYB              F @                               2     
       #         @                                   3                    #THIS 4   #DET 5             
D @                               4     H              #MATRIXHYB              D @                               5     
       #         @                                   6                    #THIS 7   #OSTREAM 8   #OPT_PRINT 9             
                                  7     H             #MATRIXHYB              
 @                               8                     
 @                               9           #         @                                   :                    #THIS ;             
D                                 ;     H              #MATRIXHYB       �   $      fn#fn !   �   �   b   uapp(M_MATRIXHYB    e  @   J  M_GLOBAL (   �  ]      MPI_GROUP+MPI_CONSTANTS 0     H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS '   J  ]      MPI_INFO+MPI_CONSTANTS /   �  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS *   �  ]      MPI_MESSAGE+MPI_CONSTANTS 2   L  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS &   �  ]      MPI_WIN+MPI_CONSTANTS .   �  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS +   9  ]      MPI_DATATYPE+MPI_CONSTANTS 3   �  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS *   �  ]      MPI_REQUEST+MPI_CONSTANTS 2   ;  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS '   �  ]      MPI_FILE+MPI_CONSTANTS /   �  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS -   (  ]      MPI_ERRHANDLER+MPI_CONSTANTS 5   �  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '   �  ]      MPI_COMM+MPI_CONSTANTS /   *  H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS %   r  ]      MPI_OP+MPI_CONSTANTS -   �  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS      U       ERROR+M_GLOBAL '   l  L   a   ERROR%MESSAGE+M_GLOBAL !   �  U       WARNALL+M_GLOBAL )   	  L   a   WARNALL%MESSAGE+M_GLOBAL $   Y	  �       ABI_CABORT+M_ERRORS    �	  �       MATRIXHYB    �
  H   !   MATRIXHYB%SIZE    �
  H   a   MATRIXHYB%TAIL       �   !   MATRIXHYB%ITECH    �  H   !   MATRIXHYB%WMAX    
  �   a   MATRIXHYB%MAT "   �  �   a   MATRIXHYB%MAT_TAU $   b  �   a   MATRIXHYB%MAT_OMEGA    &  q       MATRIXHYB_INIT $   �  W   a   MATRIXHYB_INIT%THIS %   �  @   a   MATRIXHYB_INIT%ITECH $   .  @   a   MATRIXHYB_INIT%SIZE $   n  @   a   MATRIXHYB_INIT%WMAX "   �  `       MATRIXHYB_SETSIZE '     W   a   MATRIXHYB_SETSIZE%THIS +   e  @   a   MATRIXHYB_SETSIZE%NEW_TAIL     �  R       MATRIXHYB_CLEAR %   �  W   a   MATRIXHYB_CLEAR%THIS !   N  ^       MATRIXHYB_ASSIGN &   �  W   a   MATRIXHYB_ASSIGN%THIS (     W   a   MATRIXHYB_ASSIGN%MATRIX "   Z  c       MATRIXHYB_INVERSE '   �  W   a   MATRIXHYB_INVERSE%THIS .     @   a   MATRIXHYB_INVERSE%DETERMINANT !   T  [       MATRIXHYB_GETDET &   �  W   a   MATRIXHYB_GETDET%THIS %     @   a   MATRIXHYB_GETDET%DET     F  n       MATRIXHYB_PRINT %   �  W   a   MATRIXHYB_PRINT%THIS (     @   a   MATRIXHYB_PRINT%OSTREAM *   K  @   a   MATRIXHYB_PRINT%OPT_PRINT "   �  R       MATRIXHYB_DESTROY '   �  W   a   MATRIXHYB_DESTROY%THIS 