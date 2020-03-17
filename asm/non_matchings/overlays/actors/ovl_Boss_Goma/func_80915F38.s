glabel func_80915F38
/* 00528 80915F38 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 0052C 80915F3C AFB00028 */  sw      $s0, 0x0028($sp)           
/* 00530 80915F40 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00534 80915F44 AFBF002C */  sw      $ra, 0x002C($sp)           
/* 00538 80915F48 3C040601 */  lui     $a0, 0x0601                ## $a0 = 06010000
/* 0053C 80915F4C AFA50034 */  sw      $a1, 0x0034($sp)           
/* 00540 80915F50 0C028800 */  jal     SkelAnime_GetFrameCount
              
/* 00544 80915F54 2484B2FC */  addiu   $a0, $a0, 0xB2FC           ## $a0 = 0600B2FC
/* 00548 80915F58 44822000 */  mtc1    $v0, $f4                   ## $f4 = 0.00
/* 0054C 80915F5C 3C01C000 */  lui     $at, 0xC000                ## $at = C0000000
/* 00550 80915F60 44814000 */  mtc1    $at, $f8                   ## $f8 = -2.00
/* 00554 80915F64 468021A0 */  cvt.s.w $f6, $f4                   
/* 00558 80915F68 3C050601 */  lui     $a1, 0x0601                ## $a1 = 06010000
/* 0055C 80915F6C 240E0002 */  addiu   $t6, $zero, 0x0002         ## $t6 = 00000002
/* 00560 80915F70 AFAE0014 */  sw      $t6, 0x0014($sp)           
/* 00564 80915F74 24A5B2FC */  addiu   $a1, $a1, 0xB2FC           ## $a1 = 0600B2FC
/* 00568 80915F78 2604014C */  addiu   $a0, $s0, 0x014C           ## $a0 = 0000014C
/* 0056C 80915F7C E7A60010 */  swc1    $f6, 0x0010($sp)           
/* 00570 80915F80 3C063F80 */  lui     $a2, 0x3F80                ## $a2 = 3F800000
/* 00574 80915F84 24070000 */  addiu   $a3, $zero, 0x0000         ## $a3 = 00000000
/* 00578 80915F88 0C029468 */  jal     SkelAnime_ChangeAnimation
              
/* 0057C 80915F8C E7A80018 */  swc1    $f8, 0x0018($sp)           
/* 00580 80915F90 8E080004 */  lw      $t0, 0x0004($s0)           ## 00000004
/* 00584 80915F94 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 00588 80915F98 3C0F8091 */  lui     $t7, %hi(func_80917D98)    ## $t7 = 80910000
/* 0058C 80915F9C 2401FFFA */  addiu   $at, $zero, 0xFFFA         ## $at = FFFFFFFA
/* 00590 80915FA0 25EF7D98 */  addiu   $t7, $t7, %lo(func_80917D98) ## $t7 = 80917D98
/* 00594 80915FA4 24180001 */  addiu   $t8, $zero, 0x0001         ## $t8 = 00000001
/* 00598 80915FA8 241904B0 */  addiu   $t9, $zero, 0x04B0         ## $t9 = 000004B0
/* 0059C 80915FAC 3C041001 */  lui     $a0, 0x1001                ## $a0 = 10010000
/* 005A0 80915FB0 01014824 */  and     $t1, $t0, $at              
/* 005A4 80915FB4 AE0F0190 */  sw      $t7, 0x0190($s0)           ## 00000190
/* 005A8 80915FB8 A61801BE */  sh      $t8, 0x01BE($s0)           ## 000001BE
/* 005AC 80915FBC A60001C0 */  sh      $zero, 0x01C0($s0)         ## 000001C0
/* 005B0 80915FC0 A60001C2 */  sh      $zero, 0x01C2($s0)         ## 000001C2
/* 005B4 80915FC4 A61901D2 */  sh      $t9, 0x01D2($s0)           ## 000001D2
/* 005B8 80915FC8 A60001D0 */  sh      $zero, 0x01D0($s0)         ## 000001D0
/* 005BC 80915FCC AE090004 */  sw      $t1, 0x0004($s0)           ## 00000004
/* 005C0 80915FD0 348400FF */  ori     $a0, $a0, 0x00FF           ## $a0 = 100100FF
/* 005C4 80915FD4 E6000068 */  swc1    $f0, 0x0068($s0)           ## 00000068
/* 005C8 80915FD8 0C03E803 */  jal     Audio_SetBGM
              
/* 005CC 80915FDC E60000C4 */  swc1    $f0, 0x00C4($s0)           ## 000000C4
/* 005D0 80915FE0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 005D4 80915FE4 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 005D8 80915FE8 24053812 */  addiu   $a1, $zero, 0x3812         ## $a1 = 00003812
/* 005DC 80915FEC 8FBF002C */  lw      $ra, 0x002C($sp)           
/* 005E0 80915FF0 8FB00028 */  lw      $s0, 0x0028($sp)           
/* 005E4 80915FF4 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000
/* 005E8 80915FF8 03E00008 */  jr      $ra                        
/* 005EC 80915FFC 00000000 */  nop


