glabel func_80B07734
/* 001F4 80B07734 27BDFFC0 */  addiu   $sp, $sp, 0xFFC0           ## $sp = FFFFFFC0
/* 001F8 80B07738 AFBF002C */  sw      $ra, 0x002C($sp)
/* 001FC 80B0773C AFB00028 */  sw      $s0, 0x0028($sp)
/* 00200 80B07740 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00204 80B07744 0C2C1DAC */  jal     func_80B076B0
/* 00208 80B07748 AFA50044 */  sw      $a1, 0x0044($sp)
/* 0020C 80B0774C 920E02A4 */  lbu     $t6, 0x02A4($s0)           ## 000002A4
/* 00210 80B07750 8FA40044 */  lw      $a0, 0x0044($sp)
/* 00214 80B07754 3C190001 */  lui     $t9, 0x0001                ## $t9 = 00010000
/* 00218 80B07758 000E7900 */  sll     $t7, $t6,  4
/* 0021C 80B0775C 01EE7821 */  addu    $t7, $t7, $t6
/* 00220 80B07760 000F7880 */  sll     $t7, $t7,  2
/* 00224 80B07764 008FC021 */  addu    $t8, $a0, $t7
/* 00228 80B07768 0338C821 */  addu    $t9, $t9, $t8
/* 0022C 80B0776C 8F3917B4 */  lw      $t9, 0x17B4($t9)           ## 000117B4
/* 00230 80B07770 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 00234 80B07774 3C0680B1 */  lui     $a2, %hi(D_80B0B43C)       ## $a2 = 80B10000
/* 00238 80B07778 03214021 */  addu    $t0, $t9, $at
/* 0023C 80B0777C 3C018016 */  lui     $at, 0x8016                ## $at = 80160000
/* 00240 80B07780 AC286FC0 */  sw      $t0, 0x6FC0($at)           ## 80166FC0
/* 00244 80B07784 8609001C */  lh      $t1, 0x001C($s0)           ## 0000001C
/* 00248 80B07788 26050198 */  addiu   $a1, $s0, 0x0198           ## $a1 = 00000198
/* 0024C 80B0778C 260B01DC */  addiu   $t3, $s0, 0x01DC           ## $t3 = 000001DC
/* 00250 80B07790 00095080 */  sll     $t2, $t1,  2
/* 00254 80B07794 00CA3021 */  addu    $a2, $a2, $t2
/* 00258 80B07798 260C023C */  addiu   $t4, $s0, 0x023C           ## $t4 = 0000023C
/* 0025C 80B0779C 240D0010 */  addiu   $t5, $zero, 0x0010         ## $t5 = 00000010
/* 00260 80B077A0 AFAD0018 */  sw      $t5, 0x0018($sp)
/* 00264 80B077A4 AFAC0014 */  sw      $t4, 0x0014($sp)
/* 00268 80B077A8 8CC6B43C */  lw      $a2, %lo(D_80B0B43C)($a2)
/* 0026C 80B077AC AFAB0010 */  sw      $t3, 0x0010($sp)
/* 00270 80B077B0 AFA50034 */  sw      $a1, 0x0034($sp)
/* 00274 80B077B4 0C0291BE */  jal     SkelAnime_InitSV
/* 00278 80B077B8 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 0027C 80B077BC 860E001C */  lh      $t6, 0x001C($s0)           ## 0000001C
/* 00280 80B077C0 3C0580B1 */  lui     $a1, %hi(D_80B0B454)       ## $a1 = 80B10000
/* 00284 80B077C4 8FA40034 */  lw      $a0, 0x0034($sp)
/* 00288 80B077C8 000E7880 */  sll     $t7, $t6,  2
/* 0028C 80B077CC 00AF2821 */  addu    $a1, $a1, $t7
/* 00290 80B077D0 0C0294BE */  jal     SkelAnime_ChangeAnimDefaultRepeat
/* 00294 80B077D4 8CA5B454 */  lw      $a1, %lo(D_80B0B454)($a1)
/* 00298 80B077D8 8618001C */  lh      $t8, 0x001C($s0)           ## 0000001C
/* 0029C 80B077DC 3C0880B1 */  lui     $t0, %hi(D_80B0B484)       ## $t0 = 80B10000
/* 002A0 80B077E0 3C098016 */  lui     $t1, 0x8016                ## $t1 = 80160000
/* 002A4 80B077E4 0018C840 */  sll     $t9, $t8,  1
/* 002A8 80B077E8 01194021 */  addu    $t0, $t0, $t9
/* 002AC 80B077EC 9508B484 */  lhu     $t0, %lo(D_80B0B484)($t0)
/* 002B0 80B077F0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 002B4 80B077F4 3C0580B0 */  lui     $a1, %hi(func_80B07D7C)    ## $a1 = 80B00000
/* 002B8 80B077F8 A608029C */  sh      $t0, 0x029C($s0)           ## 0000029C
/* 002BC 80B077FC 9529F54E */  lhu     $t1, -0x0AB2($t1)          ## 8015F54E
/* 002C0 80B07800 310AFFFF */  andi    $t2, $t0, 0xFFFF           ## $t2 = 00000000
/* 002C4 80B07804 012A5824 */  and     $t3, $t1, $t2
/* 002C8 80B07808 1160000A */  beq     $t3, $zero, .L80B07834
/* 002CC 80B0780C 00000000 */  nop
/* 002D0 80B07810 860C001C */  lh      $t4, 0x001C($s0)           ## 0000001C
/* 002D4 80B07814 3C0580B1 */  lui     $a1, %hi(D_80B0B46C)       ## $a1 = 80B10000
/* 002D8 80B07818 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 002DC 80B0781C 000C6880 */  sll     $t5, $t4,  2
/* 002E0 80B07820 00AD2821 */  addu    $a1, $a1, $t5
/* 002E4 80B07824 0C2C1D50 */  jal     func_80B07540
/* 002E8 80B07828 8CA5B46C */  lw      $a1, %lo(D_80B0B46C)($a1)
/* 002EC 80B0782C 10000004 */  beq     $zero, $zero, .L80B07840
/* 002F0 80B07830 8FBF002C */  lw      $ra, 0x002C($sp)
.L80B07834:
/* 002F4 80B07834 0C2C1D50 */  jal     func_80B07540
/* 002F8 80B07838 24A57D7C */  addiu   $a1, $a1, %lo(func_80B07D7C) ## $a1 = 00007D7C
/* 002FC 80B0783C 8FBF002C */  lw      $ra, 0x002C($sp)
.L80B07840:
/* 00300 80B07840 8FB00028 */  lw      $s0, 0x0028($sp)
/* 00304 80B07844 27BD0040 */  addiu   $sp, $sp, 0x0040           ## $sp = 00000000
/* 00308 80B07848 03E00008 */  jr      $ra
/* 0030C 80B0784C 00000000 */  nop
