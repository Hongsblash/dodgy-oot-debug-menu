glabel func_80862578
/* 02F28 80862578 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 02F2C 8086257C AFB00020 */  sw      $s0, 0x0020($sp)           
/* 02F30 80862580 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 02F34 80862584 3C01C080 */  lui     $at, 0xC080                ## $at = C0800000
/* 02F38 80862588 44813000 */  mtc1    $at, $f6                   ## $f6 = -4.00
/* 02F3C 8086258C 908F07E2 */  lbu     $t7, 0x07E2($a0)           ## 000007E2
/* 02F40 80862590 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 02F44 80862594 240E000B */  addiu   $t6, $zero, 0x000B         ## $t6 = 0000000B
/* 02F48 80862598 2401000E */  addiu   $at, $zero, 0x000E         ## $at = 0000000E
/* 02F4C 8086259C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 02F50 808625A0 A08E07C8 */  sb      $t6, 0x07C8($a0)           ## 000007C8
/* 02F54 808625A4 A08007DE */  sb      $zero, 0x07DE($a0)         ## 000007DE
/* 02F58 808625A8 A0800808 */  sb      $zero, 0x0808($a0)         ## 00000808
/* 02F5C 808625AC E4860068 */  swc1    $f6, 0x0068($a0)           ## 00000068
/* 02F60 808625B0 15E10009 */  bne     $t7, $at, .L808625D8       
/* 02F64 808625B4 E48401A4 */  swc1    $f4, 0x01A4($a0)           ## 000001A4
/* 02F68 808625B8 24180050 */  addiu   $t8, $zero, 0x0050         ## $t8 = 00000050
/* 02F6C 808625BC AFB80010 */  sw      $t8, 0x0010($sp)           
/* 02F70 808625C0 24058000 */  addiu   $a1, $zero, 0x8000         ## $a1 = FFFF8000
/* 02F74 808625C4 24060078 */  addiu   $a2, $zero, 0x0078         ## $a2 = 00000078
/* 02F78 808625C8 0C00D09B */  jal     func_8003426C              
/* 02F7C 808625CC 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 02F80 808625D0 10000014 */  beq     $zero, $zero, .L80862624   
/* 02F84 808625D4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
.L808625D8:
/* 02F88 808625D8 24190050 */  addiu   $t9, $zero, 0x0050         ## $t9 = 00000050
/* 02F8C 808625DC AFB90010 */  sw      $t9, 0x0010($sp)           
/* 02F90 808625E0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 02F94 808625E4 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 02F98 808625E8 24060078 */  addiu   $a2, $zero, 0x0078         ## $a2 = 00000078
/* 02F9C 808625EC 0C00D09B */  jal     func_8003426C              
/* 02FA0 808625F0 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 02FA4 808625F4 920807E2 */  lbu     $t0, 0x07E2($s0)           ## 000007E2
/* 02FA8 808625F8 2401000F */  addiu   $at, $zero, 0x000F         ## $at = 0000000F
/* 02FAC 808625FC 26040188 */  addiu   $a0, $s0, 0x0188           ## $a0 = 00000188
/* 02FB0 80862600 15010004 */  bne     $t0, $at, .L80862614       
/* 02FB4 80862604 3C050601 */  lui     $a1, 0x0601                ## $a1 = 06010000
/* 02FB8 80862608 24090024 */  addiu   $t1, $zero, 0x0024         ## $t1 = 00000024
/* 02FBC 8086260C 10000004 */  beq     $zero, $zero, .L80862620   
/* 02FC0 80862610 A60907E0 */  sh      $t1, 0x07E0($s0)           ## 000007E0
.L80862614:
/* 02FC4 80862614 24A58604 */  addiu   $a1, $a1, 0x8604           ## $a1 = 06008604
/* 02FC8 80862618 0C0294A7 */  jal     func_800A529C              
/* 02FCC 8086261C 24060000 */  addiu   $a2, $zero, 0x0000         ## $a2 = 00000000
.L80862620:
/* 02FD0 80862620 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
.L80862624:
/* 02FD4 80862624 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 02FD8 80862628 2405389E */  addiu   $a1, $zero, 0x389E         ## $a1 = 0000389E
/* 02FDC 8086262C 3C058086 */  lui     $a1, %hi(func_80862650)    ## $a1 = 80860000
/* 02FE0 80862630 24A52650 */  addiu   $a1, $a1, %lo(func_80862650) ## $a1 = 80862650
/* 02FE4 80862634 0C217D94 */  jal     func_8085F650              
/* 02FE8 80862638 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 02FEC 8086263C 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 02FF0 80862640 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 02FF4 80862644 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 02FF8 80862648 03E00008 */  jr      $ra                        
/* 02FFC 8086264C 00000000 */  nop


