glabel func_80862154
/* 02B04 80862154 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 02B08 80862158 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 02B0C 8086215C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 02B10 80862160 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 02B14 80862164 3C050601 */  lui     $a1, 0x0601                ## $a1 = 06010000
/* 02B18 80862168 24A58604 */  addiu   $a1, $a1, 0x8604           ## $a1 = 06008604
/* 02B1C 8086216C 0C02947A */  jal     SkelAnimeChangeAnimationDefaultStop              
/* 02B20 80862170 24840188 */  addiu   $a0, $a0, 0x0188           ## $a0 = 00000188
/* 02B24 80862174 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 02B28 80862178 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 02B2C 8086217C 2405383A */  addiu   $a1, $zero, 0x383A         ## $a1 = 0000383A
/* 02B30 80862180 3C01C000 */  lui     $at, 0xC000                ## $at = C0000000
/* 02B34 80862184 44812000 */  mtc1    $at, $f4                   ## $f4 = -2.00
/* 02B38 80862188 240E0008 */  addiu   $t6, $zero, 0x0008         ## $t6 = 00000008
/* 02B3C 8086218C A20E07C8 */  sb      $t6, 0x07C8($s0)           ## 000007C8
/* 02B40 80862190 240F0008 */  addiu   $t7, $zero, 0x0008         ## $t7 = 00000008
/* 02B44 80862194 E6040068 */  swc1    $f4, 0x0068($s0)           ## 00000068
/* 02B48 80862198 AFAF0010 */  sw      $t7, 0x0010($sp)           
/* 02B4C 8086219C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 02B50 808621A0 24054000 */  addiu   $a1, $zero, 0x4000         ## $a1 = 00004000
/* 02B54 808621A4 240600FF */  addiu   $a2, $zero, 0x00FF         ## $a2 = 000000FF
/* 02B58 808621A8 0C00D09B */  jal     func_8003426C              
/* 02B5C 808621AC 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 02B60 808621B0 3C058086 */  lui     $a1, %hi(func_808621D4)    ## $a1 = 80860000
/* 02B64 808621B4 24A521D4 */  addiu   $a1, $a1, %lo(func_808621D4) ## $a1 = 808621D4
/* 02B68 808621B8 0C217D94 */  jal     func_8085F650              
/* 02B6C 808621BC 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 02B70 808621C0 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 02B74 808621C4 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 02B78 808621C8 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 02B7C 808621CC 03E00008 */  jr      $ra                        
/* 02B80 808621D0 00000000 */  nop


