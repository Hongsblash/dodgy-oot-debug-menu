glabel func_80A004BC
/* 0051C 80A004BC 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00520 80A004C0 AFBF0014 */  sw      $ra, 0x0014($sp)
/* 00524 80A004C4 848E00B6 */  lh      $t6, 0x00B6($a0)           ## 000000B6
/* 00528 80A004C8 3C0180A0 */  lui     $at, %hi(D_80A019F0)       ## $at = 80A00000
/* 0052C 80A004CC A48E0196 */  sh      $t6, 0x0196($a0)           ## 00000196
/* 00530 80A004D0 AFA40018 */  sw      $a0, 0x0018($sp)
/* 00534 80A004D4 0C041184 */  jal     cosf

/* 00538 80A004D8 C42C19F0 */  lwc1    $f12, %lo(D_80A019F0)($at)
/* 0053C 80A004DC 3C0140A0 */  lui     $at, 0x40A0                ## $at = 40A00000
/* 00540 80A004E0 44812000 */  mtc1    $at, $f4                   ## $f4 = 5.00
/* 00544 80A004E4 8FA70018 */  lw      $a3, 0x0018($sp)
/* 00548 80A004E8 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 0054C 80A004EC 46040182 */  mul.s   $f6, $f0, $f4
/* 00550 80A004F0 C4E80028 */  lwc1    $f8, 0x0028($a3)           ## 00000028
/* 00554 80A004F4 24A50FC0 */  addiu   $a1, $a1, 0x0FC0           ## $a1 = 06000FC0
/* 00558 80A004F8 3C06C0A0 */  lui     $a2, 0xC0A0                ## $a2 = C0A00000
/* 0055C 80A004FC 24E4014C */  addiu   $a0, $a3, 0x014C           ## $a0 = 0000014C
/* 00560 80A00500 46083280 */  add.s   $f10, $f6, $f8
/* 00564 80A00504 0C0294D3 */  jal     SkelAnime_ChangeAnimationTransitionRate
/* 00568 80A00508 E4EA0280 */  swc1    $f10, 0x0280($a3)          ## 00000280
/* 0056C 80A0050C 8FA70018 */  lw      $a3, 0x0018($sp)
/* 00570 80A00510 3C1880A0 */  lui     $t8, %hi(func_80A00C70)    ## $t8 = 80A00000
/* 00574 80A00514 240F003C */  addiu   $t7, $zero, 0x003C         ## $t7 = 0000003C
/* 00578 80A00518 27180C70 */  addiu   $t8, $t8, %lo(func_80A00C70) ## $t8 = 80A00C70
/* 0057C 80A0051C A4EF0194 */  sh      $t7, 0x0194($a3)           ## 00000194
/* 00580 80A00520 ACF80190 */  sw      $t8, 0x0190($a3)           ## 00000190
/* 00584 80A00524 8FBF0014 */  lw      $ra, 0x0014($sp)
/* 00588 80A00528 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 0058C 80A0052C 03E00008 */  jr      $ra
/* 00590 80A00530 00000000 */  nop


