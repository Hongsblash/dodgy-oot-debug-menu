glabel func_8083F570
/* 0D360 8083F570 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 0D364 8083F574 AFBF002C */  sw      $ra, 0x002C($sp)
/* 0D368 8083F578 AFB10028 */  sw      $s1, 0x0028($sp)
/* 0D36C 8083F57C AFB00024 */  sw      $s0, 0x0024($sp)
/* 0D370 8083F580 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 0D374 8083F584 C4800838 */  lwc1    $f0, 0x0838($a0)           ## 00000838
/* 0D378 8083F588 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0D37C 8083F58C 00A08825 */  or      $s1, $a1, $zero            ## $s1 = 00000000
/* 0D380 8083F590 46001032 */  c.eq.s  $f2, $f0
/* 0D384 8083F594 00000000 */  nop
/* 0D388 8083F598 4503005F */  bc1tl   .L8083F718
/* 0D38C 8083F59C 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 0D390 8083F5A0 948E0088 */  lhu     $t6, 0x0088($a0)           ## 00000088
/* 0D394 8083F5A4 3C188085 */  lui     $t8, %hi(D_808535F0)       ## $t8 = 80850000
/* 0D398 8083F5A8 31CF0008 */  andi    $t7, $t6, 0x0008           ## $t7 = 00000000
/* 0D39C 8083F5AC 51E0005A */  beql    $t7, $zero, .L8083F718
/* 0D3A0 8083F5B0 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 0D3A4 8083F5B4 8F1835F0 */  lw      $t8, %lo(D_808535F0)($t8)
/* 0D3A8 8083F5B8 33190030 */  andi    $t9, $t8, 0x0030           ## $t9 = 00000000
/* 0D3AC 8083F5BC 53200056 */  beql    $t9, $zero, .L8083F718
/* 0D3B0 8083F5C0 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 0D3B4 8083F5C4 848800B6 */  lh      $t0, 0x00B6($a0)           ## 000000B6
/* 0D3B8 8083F5C8 8489007E */  lh      $t1, 0x007E($a0)           ## 0000007E
/* 0D3BC 8083F5CC 4602003C */  c.lt.s  $f0, $f2
/* 0D3C0 8083F5D0 34018000 */  ori     $at, $zero, 0x8000         ## $at = 00008000
/* 0D3C4 8083F5D4 01091023 */  subu    $v0, $t0, $t1
/* 0D3C8 8083F5D8 00021400 */  sll     $v0, $v0, 16
/* 0D3CC 8083F5DC 45000004 */  bc1f    .L8083F5F0
/* 0D3D0 8083F5E0 00021403 */  sra     $v0, $v0, 16
/* 0D3D4 8083F5E4 00411021 */  addu    $v0, $v0, $at
/* 0D3D8 8083F5E8 00021400 */  sll     $v0, $v0, 16
/* 0D3DC 8083F5EC 00021403 */  sra     $v0, $v0, 16
.L8083F5F0:
/* 0D3E0 8083F5F0 04400003 */  bltz    $v0, .L8083F600
/* 0D3E4 8083F5F4 00021823 */  subu    $v1, $zero, $v0
/* 0D3E8 8083F5F8 10000001 */  beq     $zero, $zero, .L8083F600
/* 0D3EC 8083F5FC 00401825 */  or      $v1, $v0, $zero            ## $v1 = 00000000
.L8083F600:
/* 0D3F0 8083F600 28614001 */  slti    $at, $v1, 0x4001
/* 0D3F4 8083F604 14200043 */  bne     $at, $zero, .L8083F714
/* 0D3F8 8083F608 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 0D3FC 8083F60C 3C068085 */  lui     $a2, %hi(func_8084C81C)    ## $a2 = 80850000
/* 0D400 8083F610 24C6C81C */  addiu   $a2, $a2, %lo(func_8084C81C) ## $a2 = 8084C81C
/* 0D404 8083F614 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0D408 8083F618 0C20D716 */  jal     func_80835C58
/* 0D40C 8083F61C 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 0D410 8083F620 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 0D414 8083F624 C6040838 */  lwc1    $f4, 0x0838($s0)           ## 00000838
/* 0D418 8083F628 34018000 */  ori     $at, $zero, 0x8000         ## $at = 00008000
/* 0D41C 8083F62C 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0D420 8083F630 4604103C */  c.lt.s  $f2, $f4
/* 0D424 8083F634 3C060400 */  lui     $a2, 0x0400                ## $a2 = 04000000
/* 0D428 8083F638 3C040400 */  lui     $a0, 0x0400                ## $a0 = 04000000
/* 0D42C 8083F63C 45020014 */  bc1fl   .L8083F690
/* 0D430 8083F640 860C007E */  lh      $t4, 0x007E($s0)           ## 0000007E
/* 0D434 8083F644 860A007E */  lh      $t2, 0x007E($s0)           ## 0000007E
/* 0D438 8083F648 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 0D43C 8083F64C 24C62700 */  addiu   $a2, $a2, 0x2700           ## $a2 = 04002700
/* 0D440 8083F650 01415821 */  addu    $t3, $t2, $at
/* 0D444 8083F654 0C20C899 */  jal     func_80832264
/* 0D448 8083F658 A60B00B6 */  sh      $t3, 0x00B6($s0)           ## 000000B6
/* 0D44C 8083F65C 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 0D450 8083F660 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0D454 8083F664 0C20CBD5 */  jal     func_80832F54
/* 0D458 8083F668 2406009D */  addiu   $a2, $zero, 0x009D         ## $a2 = 0000009D
/* 0D45C 8083F66C 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 0D460 8083F670 24052581 */  addiu   $a1, $zero, 0x2581         ## $a1 = 00002581
/* 0D464 8083F674 240603E7 */  addiu   $a2, $zero, 0x03E7         ## $a2 = 000003E7
/* 0D468 8083F678 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 0D46C 8083F67C 0C02003E */  jal     func_800800F8
/* 0D470 8083F680 AFA00010 */  sw      $zero, 0x0010($sp)
/* 0D474 8083F684 1000001E */  beq     $zero, $zero, .L8083F700
/* 0D478 8083F688 860E00B6 */  lh      $t6, 0x00B6($s0)           ## 000000B6
/* 0D47C 8083F68C 860C007E */  lh      $t4, 0x007E($s0)           ## 0000007E
.L8083F690:
/* 0D480 8083F690 24842708 */  addiu   $a0, $a0, 0x2708           ## $a0 = 00002708
/* 0D484 8083F694 0C028800 */  jal     SkelAnime_GetFrameCount

/* 0D488 8083F698 A60C00B6 */  sh      $t4, 0x00B6($s0)           ## 000000B6
/* 0D48C 8083F69C 44823000 */  mtc1    $v0, $f6                   ## $f6 = 0.00
/* 0D490 8083F6A0 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 0D494 8083F6A4 3C060400 */  lui     $a2, 0x0400                ## $a2 = 04000000
/* 0D498 8083F6A8 46803220 */  cvt.s.w $f8, $f6
/* 0D49C 8083F6AC 240D0002 */  addiu   $t5, $zero, 0x0002         ## $t5 = 00000002
/* 0D4A0 8083F6B0 AFAD0018 */  sw      $t5, 0x0018($sp)
/* 0D4A4 8083F6B4 24C62708 */  addiu   $a2, $a2, 0x2708           ## $a2 = 04002708
/* 0D4A8 8083F6B8 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 0D4AC 8083F6BC 260501B4 */  addiu   $a1, $s0, 0x01B4           ## $a1 = 000001B4
/* 0D4B0 8083F6C0 E7A80010 */  swc1    $f8, 0x0010($sp)
/* 0D4B4 8083F6C4 3C07BF80 */  lui     $a3, 0xBF80                ## $a3 = BF800000
/* 0D4B8 8083F6C8 E7A20014 */  swc1    $f2, 0x0014($sp)
/* 0D4BC 8083F6CC 0C028FC2 */  jal     SkelAnime_LinkChangeAnimation
/* 0D4C0 8083F6D0 E7A2001C */  swc1    $f2, 0x001C($sp)
/* 0D4C4 8083F6D4 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 0D4C8 8083F6D8 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0D4CC 8083F6DC 0C20CBD5 */  jal     func_80832F54
/* 0D4D0 8083F6E0 2406009D */  addiu   $a2, $zero, 0x009D         ## $a2 = 0000009D
/* 0D4D4 8083F6E4 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 0D4D8 8083F6E8 24052582 */  addiu   $a1, $zero, 0x2582         ## $a1 = 00002582
/* 0D4DC 8083F6EC 240603E7 */  addiu   $a2, $zero, 0x03E7         ## $a2 = 000003E7
/* 0D4E0 8083F6F0 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 0D4E4 8083F6F4 0C02003E */  jal     func_800800F8
/* 0D4E8 8083F6F8 AFA00010 */  sw      $zero, 0x0010($sp)
/* 0D4EC 8083F6FC 860E00B6 */  lh      $t6, 0x00B6($s0)           ## 000000B6
.L8083F700:
/* 0D4F0 8083F700 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 0D4F4 8083F704 0C20C884 */  jal     func_80832210
/* 0D4F8 8083F708 A60E083C */  sh      $t6, 0x083C($s0)           ## 0000083C
/* 0D4FC 8083F70C 10000002 */  beq     $zero, $zero, .L8083F718
/* 0D500 8083F710 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L8083F714:
/* 0D504 8083F714 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L8083F718:
/* 0D508 8083F718 8FBF002C */  lw      $ra, 0x002C($sp)
/* 0D50C 8083F71C 8FB00024 */  lw      $s0, 0x0024($sp)
/* 0D510 8083F720 8FB10028 */  lw      $s1, 0x0028($sp)
/* 0D514 8083F724 03E00008 */  jr      $ra
/* 0D518 8083F728 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000


