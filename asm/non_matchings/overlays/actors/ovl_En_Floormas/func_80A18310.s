glabel func_80A18310
/* 00E00 80A18310 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00E04 80A18314 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00E08 80A18318 00803025 */  or      $a2, $a0, $zero            ## $a2 = 00000000
/* 00E0C 80A1831C 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 00E10 80A18320 24A50EA4 */  addiu   $a1, $a1, 0x0EA4           ## $a1 = 06000EA4
/* 00E14 80A18324 AFA60018 */  sw      $a2, 0x0018($sp)           
/* 00E18 80A18328 0C02947A */  jal     SkelAnimeChangeAnimationDefaultStop              
/* 00E1C 80A1832C 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 00E20 80A18330 8FA60018 */  lw      $a2, 0x0018($sp)           
/* 00E24 80A18334 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 00E28 80A18338 3C0F80A2 */  lui     $t7, %hi(func_80A19C6C)    ## $t7 = 80A20000
/* 00E2C 80A1833C 84CE00B6 */  lh      $t6, 0x00B6($a2)           ## 000000B6
/* 00E30 80A18340 25EF9C6C */  addiu   $t7, $t7, %lo(func_80A19C6C) ## $t7 = 80A19C6C
/* 00E34 80A18344 ACCF0190 */  sw      $t7, 0x0190($a2)           ## 00000190
/* 00E38 80A18348 E4C00068 */  swc1    $f0, 0x0068($a2)           ## 00000068
/* 00E3C 80A1834C E4C00060 */  swc1    $f0, 0x0060($a2)           ## 00000060
/* 00E40 80A18350 A4CE0032 */  sh      $t6, 0x0032($a2)           ## 00000032
/* 00E44 80A18354 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 00E48 80A18358 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 00E4C 80A1835C 03E00008 */  jr      $ra                        
/* 00E50 80A18360 00000000 */  nop


