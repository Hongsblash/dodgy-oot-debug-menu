glabel func_808627C4
/* 03174 808627C4 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 03178 808627C8 AFB00020 */  sw      $s0, 0x0020($sp)
/* 0317C 808627CC 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 03180 808627D0 AFBF0024 */  sw      $ra, 0x0024($sp)
/* 03184 808627D4 AFA5002C */  sw      $a1, 0x002C($sp)
/* 03188 808627D8 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 0318C 808627DC 0C00CEAE */  jal     func_80033AB8
/* 03190 808627E0 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 03194 808627E4 10400005 */  beq     $v0, $zero, .L808627FC
/* 03198 808627E8 26040188 */  addiu   $a0, $s0, 0x0188           ## $a0 = 00000188
/* 0319C 808627EC 0C2183B0 */  jal     func_80860EC0
/* 031A0 808627F0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 031A4 808627F4 10000030 */  beq     $zero, $zero, .L808628B8
/* 031A8 808627F8 8FBF0024 */  lw      $ra, 0x0024($sp)
.L808627FC:
/* 031AC 808627FC 3C050601 */  lui     $a1, 0x0601                ## $a1 = 06010000
/* 031B0 80862800 24A5E2B0 */  addiu   $a1, $a1, 0xE2B0           ## $a1 = 0600E2B0
/* 031B4 80862804 0C0294D3 */  jal     SkelAnime_ChangeAnimTransitionRepeat
/* 031B8 80862808 3C06C000 */  lui     $a2, 0xC000                ## $a2 = C0000000
/* 031BC 8086280C 8605008A */  lh      $a1, 0x008A($s0)           ## 0000008A
/* 031C0 80862810 240E0001 */  addiu   $t6, $zero, 0x0001         ## $t6 = 00000001
/* 031C4 80862814 AFAE0010 */  sw      $t6, 0x0010($sp)
/* 031C8 80862818 260400B6 */  addiu   $a0, $s0, 0x00B6           ## $a0 = 000000B6
/* 031CC 8086281C 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 031D0 80862820 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS

/* 031D4 80862824 24070FA0 */  addiu   $a3, $zero, 0x0FA0         ## $a3 = 00000FA0
/* 031D8 80862828 8FAF002C */  lw      $t7, 0x002C($sp)
/* 031DC 8086282C 3C180001 */  lui     $t8, 0x0001                ## $t8 = 00010000
/* 031E0 80862830 3C014080 */  lui     $at, 0x4080                ## $at = 40800000
/* 031E4 80862834 030FC021 */  addu    $t8, $t8, $t7
/* 031E8 80862838 8F181DE4 */  lw      $t8, 0x1DE4($t8)           ## 00011DE4
/* 031EC 8086283C 33190001 */  andi    $t9, $t8, 0x0001           ## $t9 = 00000000
/* 031F0 80862840 53200006 */  beql    $t9, $zero, .L8086285C
/* 031F4 80862844 44813000 */  mtc1    $at, $f6                   ## $f6 = 4.00
/* 031F8 80862848 3C01C080 */  lui     $at, 0xC080                ## $at = C0800000
/* 031FC 8086284C 44812000 */  mtc1    $at, $f4                   ## $f4 = -4.00
/* 03200 80862850 10000004 */  beq     $zero, $zero, .L80862864
/* 03204 80862854 E6040068 */  swc1    $f4, 0x0068($s0)           ## 00000068
/* 03208 80862858 44813000 */  mtc1    $at, $f6                   ## $f6 = -4.00
.L8086285C:
/* 0320C 8086285C 00000000 */  nop
/* 03210 80862860 E6060068 */  swc1    $f6, 0x0068($s0)           ## 00000068
.L80862864:
/* 03214 80862864 860800B6 */  lh      $t0, 0x00B6($s0)           ## 000000B6
/* 03218 80862868 25093FFF */  addiu   $t1, $t0, 0x3FFF           ## $t1 = 00003FFF
/* 0321C 8086286C 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 03220 80862870 A6090032 */  sh      $t1, 0x0032($s0)           ## 00000032
/* 03224 80862874 3C0141A0 */  lui     $at, 0x41A0                ## $at = 41A00000
/* 03228 80862878 44811000 */  mtc1    $at, $f2                   ## $f2 = 20.00
/* 0322C 8086287C 240C0018 */  addiu   $t4, $zero, 0x0018         ## $t4 = 00000018
/* 03230 80862880 3C058086 */  lui     $a1, %hi(func_808628C8)    ## $a1 = 80860000
/* 03234 80862884 46020202 */  mul.s   $f8, $f0, $f2
/* 03238 80862888 A20C07C8 */  sb      $t4, 0x07C8($s0)           ## 000007C8
/* 0323C 8086288C 24A528C8 */  addiu   $a1, $a1, %lo(func_808628C8) ## $a1 = 808628C8
/* 03240 80862890 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03244 80862894 46024280 */  add.s   $f10, $f8, $f2
/* 03248 80862898 4600540D */  trunc.w.s $f16, $f10
/* 0324C 8086289C 440B8000 */  mfc1    $t3, $f16
/* 03250 808628A0 0C217D94 */  jal     EnTest_SetupAction
/* 03254 808628A4 AE0B07E8 */  sw      $t3, 0x07E8($s0)           ## 000007E8
/* 03258 808628A8 44809000 */  mtc1    $zero, $f18                ## $f18 = 0.00
/* 0325C 808628AC 00000000 */  nop
/* 03260 808628B0 E61207EC */  swc1    $f18, 0x07EC($s0)          ## 000007EC
/* 03264 808628B4 8FBF0024 */  lw      $ra, 0x0024($sp)
.L808628B8:
/* 03268 808628B8 8FB00020 */  lw      $s0, 0x0020($sp)
/* 0326C 808628BC 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 03270 808628C0 03E00008 */  jr      $ra
/* 03274 808628C4 00000000 */  nop
