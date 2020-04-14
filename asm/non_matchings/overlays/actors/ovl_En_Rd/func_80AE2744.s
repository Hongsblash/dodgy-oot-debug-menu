glabel func_80AE2744
/* 00344 80AE2744 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 00348 80AE2748 AFB00020 */  sw      $s0, 0x0020($sp)
/* 0034C 80AE274C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00350 80AE2750 AFBF0024 */  sw      $ra, 0x0024($sp)
/* 00354 80AE2754 24840188 */  addiu   $a0, $a0, 0x0188           ## $a0 = 00000188
/* 00358 80AE2758 AFA50034 */  sw      $a1, 0x0034($sp)
/* 0035C 80AE275C 0C02927F */  jal     SkelAnime_FrameUpdateMatrix

/* 00360 80AE2760 AFA4002C */  sw      $a0, 0x002C($sp)
/* 00364 80AE2764 2604030E */  addiu   $a0, $s0, 0x030E           ## $a0 = 0000030E
/* 00368 80AE2768 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 0036C 80AE276C 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 00370 80AE2770 24070064 */  addiu   $a3, $zero, 0x0064         ## $a3 = 00000064
/* 00374 80AE2774 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS

/* 00378 80AE2778 AFA00010 */  sw      $zero, 0x0010($sp)
/* 0037C 80AE277C 26040310 */  addiu   $a0, $s0, 0x0310           ## $a0 = 00000310
/* 00380 80AE2780 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 00384 80AE2784 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 00388 80AE2788 24070064 */  addiu   $a3, $zero, 0x0064         ## $a3 = 00000064
/* 0038C 80AE278C 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS

/* 00390 80AE2790 AFA00010 */  sw      $zero, 0x0010($sp)
/* 00394 80AE2794 860E001C */  lh      $t6, 0x001C($s0)           ## 0000001C
/* 00398 80AE2798 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
/* 0039C 80AE279C 55C1001C */  bnel    $t6, $at, .L80AE2810
/* 003A0 80AE27A0 860F030C */  lh      $t7, 0x030C($s0)           ## 0000030C
/* 003A4 80AE27A4 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 003A8 80AE27A8 C60601A0 */  lwc1    $f6, 0x01A0($s0)           ## 000001A0
/* 003AC 80AE27AC 46062032 */  c.eq.s  $f4, $f6
/* 003B0 80AE27B0 00000000 */  nop
/* 003B4 80AE27B4 45020016 */  bc1fl   .L80AE2810
/* 003B8 80AE27B8 860F030C */  lh      $t7, 0x030C($s0)           ## 0000030C
/* 003BC 80AE27BC 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 003C0 80AE27C0 00000000 */  nop
/* 003C4 80AE27C4 3C013F00 */  lui     $at, 0x3F00                ## $at = 3F000000
/* 003C8 80AE27C8 44814000 */  mtc1    $at, $f8                   ## $f8 = 0.50
/* 003CC 80AE27CC 8FA4002C */  lw      $a0, 0x002C($sp)
/* 003D0 80AE27D0 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 003D4 80AE27D4 4600403E */  c.le.s  $f8, $f0
/* 003D8 80AE27D8 00000000 */  nop
/* 003DC 80AE27DC 45000007 */  bc1f    .L80AE27FC
/* 003E0 80AE27E0 00000000 */  nop
/* 003E4 80AE27E4 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 003E8 80AE27E8 24A55D98 */  addiu   $a1, $a1, 0x5D98           ## $a1 = 06005D98
/* 003EC 80AE27EC 0C0294BE */  jal     SkelAnime_ChangeAnimDefaultRepeat
/* 003F0 80AE27F0 8FA4002C */  lw      $a0, 0x002C($sp)
/* 003F4 80AE27F4 10000018 */  beq     $zero, $zero, .L80AE2858
/* 003F8 80AE27F8 8E0A0118 */  lw      $t2, 0x0118($s0)           ## 00000118
.L80AE27FC:
/* 003FC 80AE27FC 0C0294BE */  jal     SkelAnime_ChangeAnimDefaultRepeat
/* 00400 80AE2800 24A557AC */  addiu   $a1, $a1, 0x57AC           ## $a1 = 000057AC
/* 00404 80AE2804 10000014 */  beq     $zero, $zero, .L80AE2858
/* 00408 80AE2808 8E0A0118 */  lw      $t2, 0x0118($s0)           ## 00000118
/* 0040C 80AE280C 860F030C */  lh      $t7, 0x030C($s0)           ## 0000030C
.L80AE2810:
/* 00410 80AE2810 25F8FFFF */  addiu   $t8, $t7, 0xFFFF           ## $t8 = FFFFFFFF
/* 00414 80AE2814 A618030C */  sh      $t8, 0x030C($s0)           ## 0000030C
/* 00418 80AE2818 8619030C */  lh      $t9, 0x030C($s0)           ## 0000030C
/* 0041C 80AE281C 5720000E */  bnel    $t9, $zero, .L80AE2858
/* 00420 80AE2820 8E0A0118 */  lw      $t2, 0x0118($s0)           ## 00000118
/* 00424 80AE2824 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 00428 80AE2828 00000000 */  nop
/* 0042C 80AE282C 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 00430 80AE2830 44811000 */  mtc1    $at, $f2                   ## $f2 = 10.00
/* 00434 80AE2834 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 00438 80AE2838 46020282 */  mul.s   $f10, $f0, $f2
/* 0043C 80AE283C E60401A0 */  swc1    $f4, 0x01A0($s0)           ## 000001A0
/* 00440 80AE2840 46025400 */  add.s   $f16, $f10, $f2
/* 00444 80AE2844 4600848D */  trunc.w.s $f18, $f16
/* 00448 80AE2848 44099000 */  mfc1    $t1, $f18
/* 0044C 80AE284C 00000000 */  nop
/* 00450 80AE2850 A609030C */  sh      $t1, 0x030C($s0)           ## 0000030C
/* 00454 80AE2854 8E0A0118 */  lw      $t2, 0x0118($s0)           ## 00000118
.L80AE2858:
/* 00458 80AE2858 51400011 */  beql    $t2, $zero, .L80AE28A0
/* 0045C 80AE285C 920D0305 */  lbu     $t5, 0x0305($s0)           ## 00000305
/* 00460 80AE2860 920B0305 */  lbu     $t3, 0x0305($s0)           ## 00000305
/* 00464 80AE2864 55600034 */  bnel    $t3, $zero, .L80AE2938
/* 00468 80AE2868 8FB90034 */  lw      $t9, 0x0034($sp)
/* 0046C 80AE286C 860C001C */  lh      $t4, 0x001C($s0)           ## 0000001C
/* 00470 80AE2870 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
/* 00474 80AE2874 11810005 */  beq     $t4, $at, .L80AE288C
/* 00478 80AE2878 00000000 */  nop
/* 0047C 80AE287C 0C2B8C77 */  jal     func_80AE31DC
/* 00480 80AE2880 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00484 80AE2884 1000002C */  beq     $zero, $zero, .L80AE2938
/* 00488 80AE2888 8FB90034 */  lw      $t9, 0x0034($sp)
.L80AE288C:
/* 0048C 80AE288C 0C2B8E4B */  jal     func_80AE392C
/* 00490 80AE2890 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00494 80AE2894 10000028 */  beq     $zero, $zero, .L80AE2938
/* 00498 80AE2898 8FB90034 */  lw      $t9, 0x0034($sp)
/* 0049C 80AE289C 920D0305 */  lbu     $t5, 0x0305($s0)           ## 00000305
.L80AE28A0:
/* 004A0 80AE28A0 51A0000C */  beql    $t5, $zero, .L80AE28D4
/* 004A4 80AE28A4 3C014316 */  lui     $at, 0x4316                ## $at = 43160000
/* 004A8 80AE28A8 860E001C */  lh      $t6, 0x001C($s0)           ## 0000001C
/* 004AC 80AE28AC 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
/* 004B0 80AE28B0 11C10005 */  beq     $t6, $at, .L80AE28C8
/* 004B4 80AE28B4 00000000 */  nop
/* 004B8 80AE28B8 0C2B8DEF */  jal     func_80AE37BC
/* 004BC 80AE28BC 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 004C0 80AE28C0 10000004 */  beq     $zero, $zero, .L80AE28D4
/* 004C4 80AE28C4 3C014316 */  lui     $at, 0x4316                ## $at = 43160000
.L80AE28C8:
/* 004C8 80AE28C8 0C2B8E4B */  jal     func_80AE392C
/* 004CC 80AE28CC 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 004D0 80AE28D0 3C014316 */  lui     $at, 0x4316                ## $at = 43160000
.L80AE28D4:
/* 004D4 80AE28D4 44813000 */  mtc1    $at, $f6                   ## $f6 = 150.00
/* 004D8 80AE28D8 C6080090 */  lwc1    $f8, 0x0090($s0)           ## 00000090
/* 004DC 80AE28DC A2000305 */  sb      $zero, 0x0305($s0)         ## 00000305
/* 004E0 80AE28E0 4606403E */  c.le.s  $f8, $f6
/* 004E4 80AE28E4 00000000 */  nop
/* 004E8 80AE28E8 45020013 */  bc1fl   .L80AE2938
/* 004EC 80AE28EC 8FB90034 */  lw      $t9, 0x0034($sp)
/* 004F0 80AE28F0 0C00B779 */  jal     func_8002DDE4
/* 004F4 80AE28F4 8FA40034 */  lw      $a0, 0x0034($sp)
/* 004F8 80AE28F8 5040000F */  beql    $v0, $zero, .L80AE2938
/* 004FC 80AE28FC 8FB90034 */  lw      $t9, 0x0034($sp)
/* 00500 80AE2900 860F001C */  lh      $t7, 0x001C($s0)           ## 0000001C
/* 00504 80AE2904 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
/* 00508 80AE2908 11E10008 */  beq     $t7, $at, .L80AE292C
/* 0050C 80AE290C 00000000 */  nop
/* 00510 80AE2910 92180305 */  lbu     $t8, 0x0305($s0)           ## 00000305
/* 00514 80AE2914 17000005 */  bne     $t8, $zero, .L80AE292C
/* 00518 80AE2918 00000000 */  nop
/* 0051C 80AE291C 0C2B8DEF */  jal     func_80AE37BC
/* 00520 80AE2920 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00524 80AE2924 10000004 */  beq     $zero, $zero, .L80AE2938
/* 00528 80AE2928 8FB90034 */  lw      $t9, 0x0034($sp)
.L80AE292C:
/* 0052C 80AE292C 0C2B8E4B */  jal     func_80AE392C
/* 00530 80AE2930 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00534 80AE2934 8FB90034 */  lw      $t9, 0x0034($sp)
.L80AE2938:
/* 00538 80AE2938 3C080001 */  lui     $t0, 0x0001                ## $t0 = 00010000
/* 0053C 80AE293C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00540 80AE2940 01194021 */  addu    $t0, $t0, $t9
/* 00544 80AE2944 8D081DE4 */  lw      $t0, 0x1DE4($t0)           ## 00011DE4
/* 00548 80AE2948 3109005F */  andi    $t1, $t0, 0x005F           ## $t1 = 00000000
/* 0054C 80AE294C 55200004 */  bnel    $t1, $zero, .L80AE2960
/* 00550 80AE2950 8FBF0024 */  lw      $ra, 0x0024($sp)
/* 00554 80AE2954 0C00BE0A */  jal     Audio_PlayActorSound2

/* 00558 80AE2958 240538E4 */  addiu   $a1, $zero, 0x38E4         ## $a1 = 000038E4
/* 0055C 80AE295C 8FBF0024 */  lw      $ra, 0x0024($sp)
.L80AE2960:
/* 00560 80AE2960 8FB00020 */  lw      $s0, 0x0020($sp)
/* 00564 80AE2964 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000
/* 00568 80AE2968 03E00008 */  jr      $ra
/* 0056C 80AE296C 00000000 */  nop
