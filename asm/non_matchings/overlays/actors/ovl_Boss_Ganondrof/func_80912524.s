glabel func_80912524
/* 01EE4 80912524 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 01EE8 80912528 00803825 */  or      $a3, $a0, $zero            ## $a3 = 00000000
/* 01EEC 8091252C AFBF0014 */  sw      $ra, 0x0014($sp)
/* 01EF0 80912530 3C040601 */  lui     $a0, 0x0601                ## $a0 = 06010000
/* 01EF4 80912534 AFA5001C */  sw      $a1, 0x001C($sp)
/* 01EF8 80912538 248429E0 */  addiu   $a0, $a0, 0x29E0           ## $a0 = 060129E0
/* 01EFC 8091253C 0C028800 */  jal     SkelAnime_GetFrameCount

/* 01F00 80912540 AFA70018 */  sw      $a3, 0x0018($sp)
/* 01F04 80912544 44822000 */  mtc1    $v0, $f4                   ## $f4 = 0.00
/* 01F08 80912548 8FA70018 */  lw      $a3, 0x0018($sp)
/* 01F0C 8091254C 3C050601 */  lui     $a1, 0x0601                ## $a1 = 06010000
/* 01F10 80912550 468021A0 */  cvt.s.w $f6, $f4
/* 01F14 80912554 24A529E0 */  addiu   $a1, $a1, 0x29E0           ## $a1 = 060129E0
/* 01F18 80912558 3C06C040 */  lui     $a2, 0xC040                ## $a2 = C0400000
/* 01F1C 8091255C 24E4014C */  addiu   $a0, $a3, 0x014C           ## $a0 = 0000014C
/* 01F20 80912560 0C0294D3 */  jal     SkelAnime_ChangeAnimationTransitionRepeat
/* 01F24 80912564 E4E601D0 */  swc1    $f6, 0x01D0($a3)           ## 000001D0
/* 01F28 80912568 8FA70018 */  lw      $a3, 0x0018($sp)
/* 01F2C 8091256C 3C0E8091 */  lui     $t6, %hi(func_80912594)    ## $t6 = 80910000
/* 01F30 80912570 25CE2594 */  addiu   $t6, $t6, %lo(func_80912594) ## $t6 = 80912594
/* 01F34 80912574 240F0014 */  addiu   $t7, $zero, 0x0014         ## $t7 = 00000014
/* 01F38 80912578 ACEE0190 */  sw      $t6, 0x0190($a3)           ## 00000190
/* 01F3C 8091257C A4EF01BC */  sh      $t7, 0x01BC($a3)           ## 000001BC
/* 01F40 80912580 A4E001A2 */  sh      $zero, 0x01A2($a3)         ## 000001A2
/* 01F44 80912584 8FBF0014 */  lw      $ra, 0x0014($sp)
/* 01F48 80912588 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 01F4C 8091258C 03E00008 */  jr      $ra
/* 01F50 80912590 00000000 */  nop


