glabel func_80B26B18
/* 00468 80B26B18 8C8F0004 */  lw      $t7, 0x0004($a0)           ## 00000004
/* 0046C 80B26B1C 9099040D */  lbu     $t9, 0x040D($a0)           ## 0000040D
/* 00470 80B26B20 2401FFFE */  addiu   $at, $zero, 0xFFFE         ## $at = FFFFFFFE
/* 00474 80B26B24 3C0980B2 */  lui     $t1, %hi(func_80B27318)    ## $t1 = 80B20000
/* 00478 80B26B28 240E0014 */  addiu   $t6, $zero, 0x0014         ## $t6 = 00000014
/* 0047C 80B26B2C 25297318 */  addiu   $t1, $t1, %lo(func_80B27318) ## $t1 = 80B27318
/* 00480 80B26B30 01E1C024 */  and     $t8, $t7, $at              
/* 00484 80B26B34 03214024 */  and     $t0, $t9, $at              
/* 00488 80B26B38 A08E0194 */  sb      $t6, 0x0194($a0)           ## 00000194
/* 0048C 80B26B3C AC980004 */  sw      $t8, 0x0004($a0)           ## 00000004
/* 00490 80B26B40 A088040D */  sb      $t0, 0x040D($a0)           ## 0000040D
/* 00494 80B26B44 03E00008 */  jr      $ra                        
/* 00498 80B26B48 AC890190 */  sw      $t1, 0x0190($a0)           ## 00000190
