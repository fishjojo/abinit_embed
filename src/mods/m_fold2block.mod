  ,  w   k820309    l          18.0        ��[                                                                                                          
       m_fold2block.F90 M_FOLD2BLOCK                                                     
                                                           
                            @                              
                        @                               '                    #MPI_VAL                 �                                                                   @                               '                    #MPI_VAL                 �                                                                   @                               '                    #MPI_VAL 	                �                               	                                    @                          
     '                    #MPI_VAL                 �                                                                   @                               '                    #MPI_VAL                 �                                                                   @                               '                    #MPI_VAL                 �                                                                   @                               '                    #MPI_VAL                 �                                                                   @                               '                    #MPI_VAL                 �                                                                   @                               '                    #MPI_VAL                 �                                                                   @                               '                    #MPI_VAL                 �                                                                                                                                                                                                                                                          264                                                       #         @                                                      #MESSAGE    #LEVEL    #MODE_PARAL    #FILE    #LINE     #NODUMP !   #NOSTOP "             
                                                    1           
                                                    1           
                                                    1           
                                                   1           
                                                       
                                 !                     
                                 "           %         @                               #                           #LHS $   #RHS %             
                                  $                   #MPI_GROUP              
                                  %                   #MPI_GROUP    %         @                               &                           #LHS '   #RHS (             
                                  '                   #MPI_INFO              
                                  (                   #MPI_INFO    %         @                               )                           #LHS *   #RHS +             
                                  *                   #MPI_MESSAGE              
                                  +                   #MPI_MESSAGE    %         @                               ,                           #LHS -   #RHS .             
                                  -                   #MPI_WIN 
             
                                  .                   #MPI_WIN 
   %         @                               /                           #LHS 0   #RHS 1             
                                  0                   #MPI_DATATYPE              
                                  1                   #MPI_DATATYPE    %         @                               2                           #LHS 3   #RHS 4             
                                  3                   #MPI_REQUEST              
                                  4                   #MPI_REQUEST    %         @                               5                           #LHS 6   #RHS 7             
                                  6                   #MPI_FILE              
                                  7                   #MPI_FILE    %         @                               8                           #LHS 9   #RHS :             
                                  9                   #MPI_ERRHANDLER              
                                  :                   #MPI_ERRHANDLER    %         @                               ;                           #LHS <   #RHS =             
                                  <                   #MPI_COMM              
                                  =                   #MPI_COMM    %         @                               >                           #LHS ?   #RHS @             
                                  ?                   #MPI_OP              
                                  @                   #MPI_OP    %         @                               A                           #LHS B   #RHS C             
                                  B                   #MPI_GROUP              
                                  C                   #MPI_GROUP    %         @                               D                           #LHS E   #RHS F             
                                  E                   #MPI_INFO              
                                  F                   #MPI_INFO    %         @                               G                           #LHS H   #RHS I             
                                  H                   #MPI_MESSAGE              
                                  I                   #MPI_MESSAGE    %         @                               J                           #LHS K   #RHS L             
                                  K                   #MPI_WIN 
             
                                  L                   #MPI_WIN 
   %         @                               M                           #LHS N   #RHS O             
                                  N                   #MPI_DATATYPE              
                                  O                   #MPI_DATATYPE    %         @                               P                           #LHS Q   #RHS R             
                                  Q                   #MPI_REQUEST              
                                  R                   #MPI_REQUEST    %         @                               S                           #LHS T   #RHS U             
                                  T                   #MPI_FILE              
                                  U                   #MPI_FILE    %         @                               V                           #LHS W   #RHS X             
                                  W                   #MPI_ERRHANDLER              
                                  X                   #MPI_ERRHANDLER    %         @                               Y                           #LHS Z   #RHS [             
                                  Z                   #MPI_COMM              
                                  [                   #MPI_COMM    %         @                               \                           #LHS ]   #RHS ^             
                                  ]                   #MPI_OP              
                                  ^                   #MPI_OP    #         @                                  _                     n                         �              Cabi_cabort                    #         @                                   `                    #FX a   #FY b   #FZ c   #VECTOR d   #COEFC f   #NV e   #WEIGHTS g             
  @                               a                     
  @                               b                     
  @                               c                    
                                  d                        p          p          5 � p        r e       p          5 � p        r e                              
                                 f                    
    p          p          5 � p        r e       p          5 � p        r e                               
                                  e                    
D                                g                    
     p            5 � p        r a   5 � p        r b   5 � p        r c         5 � p        r a   5 � p        r b   5 � p        r c                     #         @                                   h                    #XX i   #YY j   #ZZ k   #FX l   #FY m   #FZ n   #NKVAL o             
                                 i     
                
                                 j     
                
                                 k     
                
                                  l                     
                                  m                     
                                  n                    
D                                o                    
     p          p            5 � p        r l   5 � p        r m   5 � p        r n       p            5 � p        r l   5 � p        r m   5 � p        r n                     #         @                                   p                    #FOLDS q   #FNAME r             
D                                 q                    	    p          p            p                                    
D @                              r                           #         @                                   s                    #IKPT t   #NKPT u   #KPT v             
                                  t                     
                                  u                     
                                 v                   
    p          p            p                             �   &      fn#fn    �   @   J   DEFS_BASIS      @   J   M_ABICORE    F  @   J   M_ERRORS (   �  ]       MPI_GROUP+MPI_CONSTANTS 0   �  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS '   +  ]       MPI_INFO+MPI_CONSTANTS /   �  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS *   �  ]       MPI_MESSAGE+MPI_CONSTANTS 2   -  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS &   u  ]       MPI_WIN+MPI_CONSTANTS .   �  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS +     ]       MPI_DATATYPE+MPI_CONSTANTS 3   w  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS *   �  ]       MPI_REQUEST+MPI_CONSTANTS 2     H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS '   d  ]       MPI_FILE+MPI_CONSTANTS /   �  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS -   	  ]       MPI_ERRHANDLER+MPI_CONSTANTS 5   f  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '   �  ]       MPI_COMM+MPI_CONSTANTS /     H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS %   S  ]       MPI_OP+MPI_CONSTANTS -   �  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS    �  p       DP+DEFS_BASIS !   h  s       FNLEN+DEFS_BASIS #   �  @       STD_OUT+DEFS_BASIS "   	  �       MSG_HNDL+M_ERRORS *   �	  L   a   MSG_HNDL%MESSAGE+M_ERRORS (   
  L   a   MSG_HNDL%LEVEL+M_ERRORS -   O
  L   a   MSG_HNDL%MODE_PARAL+M_ERRORS '   �
  L   a   MSG_HNDL%FILE+M_ERRORS '   �
  @   a   MSG_HNDL%LINE+M_ERRORS )   '  @   a   MSG_HNDL%NODUMP+M_ERRORS )   g  @   a   MSG_HNDL%NOSTOP+M_ERRORS &   �  b       GROUPEQ+MPI_CONSTANTS *   	  W   a   GROUPEQ%LHS+MPI_CONSTANTS *   `  W   a   GROUPEQ%RHS+MPI_CONSTANTS %   �  b       INFOEQ+MPI_CONSTANTS )     V   a   INFOEQ%LHS+MPI_CONSTANTS )   o  V   a   INFOEQ%RHS+MPI_CONSTANTS (   �  b       MESSAGEEQ+MPI_CONSTANTS ,   '  Y   a   MESSAGEEQ%LHS+MPI_CONSTANTS ,   �  Y   a   MESSAGEEQ%RHS+MPI_CONSTANTS $   �  b       WINEQ+MPI_CONSTANTS (   ;  U   a   WINEQ%LHS+MPI_CONSTANTS (   �  U   a   WINEQ%RHS+MPI_CONSTANTS )   �  b       DATATYPEEQ+MPI_CONSTANTS -   G  Z   a   DATATYPEEQ%LHS+MPI_CONSTANTS -   �  Z   a   DATATYPEEQ%RHS+MPI_CONSTANTS (   �  b       REQUESTEQ+MPI_CONSTANTS ,   ]  Y   a   REQUESTEQ%LHS+MPI_CONSTANTS ,   �  Y   a   REQUESTEQ%RHS+MPI_CONSTANTS %     b       FILEEQ+MPI_CONSTANTS )   q  V   a   FILEEQ%LHS+MPI_CONSTANTS )   �  V   a   FILEEQ%RHS+MPI_CONSTANTS +     b       ERRHANDLEREQ+MPI_CONSTANTS /     \   a   ERRHANDLEREQ%LHS+MPI_CONSTANTS /   �  \   a   ERRHANDLEREQ%RHS+MPI_CONSTANTS %   7  b       COMMEQ+MPI_CONSTANTS )   �  V   a   COMMEQ%LHS+MPI_CONSTANTS )   �  V   a   COMMEQ%RHS+MPI_CONSTANTS #   E  b       OPEQ+MPI_CONSTANTS '   �  T   a   OPEQ%LHS+MPI_CONSTANTS '   �  T   a   OPEQ%RHS+MPI_CONSTANTS '   O  b       GROUPNEQ+MPI_CONSTANTS +   �  W   a   GROUPNEQ%LHS+MPI_CONSTANTS +     W   a   GROUPNEQ%RHS+MPI_CONSTANTS &   _  b       INFONEQ+MPI_CONSTANTS *   �  V   a   INFONEQ%LHS+MPI_CONSTANTS *     V   a   INFONEQ%RHS+MPI_CONSTANTS )   m  b       MESSAGENEQ+MPI_CONSTANTS -   �  Y   a   MESSAGENEQ%LHS+MPI_CONSTANTS -   (  Y   a   MESSAGENEQ%RHS+MPI_CONSTANTS %   �  b       WINNEQ+MPI_CONSTANTS )   �  U   a   WINNEQ%LHS+MPI_CONSTANTS )   8  U   a   WINNEQ%RHS+MPI_CONSTANTS *   �  b       DATATYPENEQ+MPI_CONSTANTS .   �  Z   a   DATATYPENEQ%LHS+MPI_CONSTANTS .   I  Z   a   DATATYPENEQ%RHS+MPI_CONSTANTS )   �  b       REQUESTNEQ+MPI_CONSTANTS -     Y   a   REQUESTNEQ%LHS+MPI_CONSTANTS -   ^  Y   a   REQUESTNEQ%RHS+MPI_CONSTANTS &   �  b       FILENEQ+MPI_CONSTANTS *     V   a   FILENEQ%LHS+MPI_CONSTANTS *   o  V   a   FILENEQ%RHS+MPI_CONSTANTS ,   �  b       ERRHANDLERNEQ+MPI_CONSTANTS 0   '  \   a   ERRHANDLERNEQ%LHS+MPI_CONSTANTS 0   �  \   a   ERRHANDLERNEQ%RHS+MPI_CONSTANTS &   �  b       COMMNEQ+MPI_CONSTANTS *   A  V   a   COMMNEQ%LHS+MPI_CONSTANTS *   �  V   a   COMMNEQ%RHS+MPI_CONSTANTS $   �  b       OPNEQ+MPI_CONSTANTS (   O   T   a   OPNEQ%LHS+MPI_CONSTANTS (   �   T   a   OPNEQ%RHS+MPI_CONSTANTS $   �   �       ABI_CABORT+M_ERRORS    �!  �       SORTC    "  @   a   SORTC%FX    Z"  @   a   SORTC%FY    �"  @   a   SORTC%FZ    �"  �   a   SORTC%VECTOR    �#  �   a   SORTC%COEFC    �$  @   a   SORTC%NV    �$  4  a   SORTC%WEIGHTS    �%  �       NEWK    y&  @   a   NEWK%XX    �&  @   a   NEWK%YY    �&  @   a   NEWK%ZZ    9'  @   a   NEWK%FX    y'  @   a   NEWK%FY    �'  @   a   NEWK%FZ    �'  T  a   NEWK%NKVAL    M)  ^       GETARGS    �)  �   a   GETARGS%FOLDS    ?*  P   a   GETARGS%FNAME    �*  e       PROGRESS    �*  @   a   PROGRESS%IKPT    4+  @   a   PROGRESS%NKPT    t+  �   a   PROGRESS%KPT 