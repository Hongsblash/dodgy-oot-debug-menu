glabel func_80B9E280
/* 01070 80B9E280 3C0E0501 */  lui     $t6, 0x0501                ## $t6 = 05010000
/* 01074 80B9E284 3C0F80BA */  lui     $t7, %hi(func_80B9E2A8)    ## $t7 = 80BA0000
/* 01078 80B9E288 25CE44B0 */  addiu   $t6, $t6, 0x44B0           ## $t6 = 050144B0
/* 0107C 80B9E28C 25EFE2A8 */  addiu   $t7, $t7, %lo(func_80B9E2A8) ## $t7 = 80B9E2A8
/* 01080 80B9E290 A080017C */  sb      $zero, 0x017C($a0)         ## 0000017C
/* 01084 80B9E294 A080017D */  sb      $zero, 0x017D($a0)         ## 0000017D
/* 01088 80B9E298 A080017E */  sb      $zero, 0x017E($a0)         ## 0000017E
/* 0108C 80B9E29C AC8E0174 */  sw      $t6, 0x0174($a0)           ## 00000174
/* 01090 80B9E2A0 03E00008 */  jr      $ra                        
/* 01094 80B9E2A4 AC8F0164 */  sw      $t7, 0x0164($a0)           ## 00000164
