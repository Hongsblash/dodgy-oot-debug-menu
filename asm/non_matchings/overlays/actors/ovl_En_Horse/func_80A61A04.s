glabel func_80A61A04
/* 06714 80A61A04 8C8F01F0 */  lw      $t7, 0x01F0($a0)           ## 000001F0
/* 06718 80A61A08 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 0671C 80A61A0C 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 06720 80A61A10 240E0013 */  addiu   $t6, $zero, 0x0013         ## $t6 = 00000013
/* 06724 80A61A14 01E1C025 */  or      $t8, $t7, $at              ## $t8 = 00010000
/* 06728 80A61A18 AC8E014C */  sw      $t6, 0x014C($a0)           ## 0000014C
/* 0672C 80A61A1C AC9801F0 */  sw      $t8, 0x01F0($a0)           ## 000001F0
/* 06730 80A61A20 03E00008 */  jr      $ra                        
/* 06734 80A61A24 E4840068 */  swc1    $f4, 0x0068($a0)           ## 00000068


