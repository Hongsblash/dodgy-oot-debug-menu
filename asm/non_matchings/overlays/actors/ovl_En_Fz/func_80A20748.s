glabel func_80A20748
/* 00988 80A20748 8C8F0004 */  lw      $t7, 0x0004($a0)           ## 00000004
/* 0098C 80A2074C 2401FFFE */  addiu   $at, $zero, 0xFFFE         ## $at = FFFFFFFE
/* 00990 80A20750 3C1980A2 */  lui     $t9, %hi(func_80A20774)    ## $t9 = 80A20000
/* 00994 80A20754 240E0002 */  addiu   $t6, $zero, 0x0002         ## $t6 = 00000002
/* 00998 80A20758 27390774 */  addiu   $t9, $t9, %lo(func_80A20774) ## $t9 = 80A20774
/* 0099C 80A2075C 01E1C024 */  and     $t8, $t7, $at              
/* 009A0 80A20760 A08E0260 */  sb      $t6, 0x0260($a0)           ## 00000260
/* 009A4 80A20764 A0800248 */  sb      $zero, 0x0248($a0)         ## 00000248
/* 009A8 80A20768 AC980004 */  sw      $t8, 0x0004($a0)           ## 00000004
/* 009AC 80A2076C 03E00008 */  jr      $ra                        
/* 009B0 80A20770 AC99014C */  sw      $t9, 0x014C($a0)           ## 0000014C


