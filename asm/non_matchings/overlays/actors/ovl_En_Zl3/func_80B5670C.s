glabel func_80B5670C
/* 0335C 80B5670C 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 03360 80B56710 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 03364 80B56714 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 03368 80B56718 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0336C 80B5671C 0C2D5378 */  jal     func_80B54DE0              
/* 03370 80B56720 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 03374 80B56724 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03378 80B56728 0C2D4D9B */  jal     func_80B5366C              
/* 0337C 80B5672C 8FA50024 */  lw      $a1, 0x0024($sp)           
/* 03380 80B56730 0C2D4D33 */  jal     func_80B534CC              
/* 03384 80B56734 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03388 80B56738 0C2D4E53 */  jal     func_80B5394C              
/* 0338C 80B5673C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03390 80B56740 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03394 80B56744 0C2D56DE */  jal     func_80B55B78              
/* 03398 80B56748 00402825 */  or      $a1, $v0, $zero            ## $a1 = 00000000
/* 0339C 80B5674C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 033A0 80B56750 0C2D592A */  jal     func_80B564A8              
/* 033A4 80B56754 8FA50024 */  lw      $a1, 0x0024($sp)           
/* 033A8 80B56758 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 033AC 80B5675C 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 033B0 80B56760 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 033B4 80B56764 03E00008 */  jr      $ra                        
/* 033B8 80B56768 00000000 */  nop


