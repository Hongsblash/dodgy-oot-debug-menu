glabel func_808B8EE0
/* 005D0 808B8EE0 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 005D4 808B8EE4 3C0E808C */  lui     $t6, %hi(func_808B8F08)    ## $t6 = 808C0000
/* 005D8 808B8EE8 25CE8F08 */  addiu   $t6, $t6, %lo(func_808B8F08) ## $t6 = 808B8F08
/* 005DC 808B8EEC AC8E0164 */  sw      $t6, 0x0164($a0)           ## 00000164
/* 005E0 808B8EF0 A4800032 */  sh      $zero, 0x0032($a0)         ## 00000032
/* 005E4 808B8EF4 E4800068 */  swc1    $f0, 0x0068($a0)           ## 00000068
/* 005E8 808B8EF8 E4800064 */  swc1    $f0, 0x0064($a0)           ## 00000064
/* 005EC 808B8EFC E4800060 */  swc1    $f0, 0x0060($a0)           ## 00000060
/* 005F0 808B8F00 03E00008 */  jr      $ra                        
/* 005F4 808B8F04 E480005C */  swc1    $f0, 0x005C($a0)           ## 0000005C


