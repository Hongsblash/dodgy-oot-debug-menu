glabel func_80B9C110
/* 00F90 80B9C110 AFA50004 */  sw      $a1, 0x0004($sp)           
/* 00F94 80B9C114 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 00F98 80B9C118 948E0168 */  lhu     $t6, 0x0168($a0)           ## 00000168
/* 00F9C 80B9C11C 3C01BF80 */  lui     $at, 0xBF80                ## $at = BF800000
/* 00FA0 80B9C120 44812000 */  mtc1    $at, $f4                   ## $f4 = -1.00
/* 00FA4 80B9C124 3C1880BA */  lui     $t8, %hi(func_80B9C14C)    ## $t8 = 80BA0000
/* 00FA8 80B9C128 2718C14C */  addiu   $t8, $t8, %lo(func_80B9C14C) ## $t8 = 80B9C14C
/* 00FAC 80B9C12C 35CF0008 */  ori     $t7, $t6, 0x0008           ## $t7 = 00000008
/* 00FB0 80B9C130 A48F0168 */  sh      $t7, 0x0168($a0)           ## 00000168
/* 00FB4 80B9C134 AC980164 */  sw      $t8, 0x0164($a0)           ## 00000164
/* 00FB8 80B9C138 E4800064 */  swc1    $f0, 0x0064($a0)           ## 00000064
/* 00FBC 80B9C13C E4800060 */  swc1    $f0, 0x0060($a0)           ## 00000060
/* 00FC0 80B9C140 E480005C */  swc1    $f0, 0x005C($a0)           ## 0000005C
/* 00FC4 80B9C144 03E00008 */  jr      $ra                        
/* 00FC8 80B9C148 E484006C */  swc1    $f4, 0x006C($a0)           ## 0000006C


