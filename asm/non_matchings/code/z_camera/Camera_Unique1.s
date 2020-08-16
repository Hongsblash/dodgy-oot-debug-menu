.late_rodata
glabel D_8013A1C4
    .float 0.01

glabel D_8013A1C8
    .incbin "baserom.z64", 0xBB1368, 0xC

glabel D_8013A1D4
    .float 0.01

glabel D_8013A1D8
    .float 0.01

glabel D_8013A1DC
    .float 0.01

glabel D_8013A1E0
    .float 0.005

glabel D_8013A1E4
    .float 0.01

glabel D_8013A1E8
    .float 0.01

.text
glabel Camera_Unique1
/* AC8744 800515A4 27BDFF50 */  addiu $sp, $sp, -0xb0
/* AC8748 800515A8 AFB00014 */  sw    $s0, 0x14($sp)
/* AC874C 800515AC 00808025 */  move  $s0, $a0
/* AC8750 800515B0 AFBF001C */  sw    $ra, 0x1c($sp)
/* AC8754 800515B4 AFB10018 */  sw    $s1, 0x18($sp)
/* AC8758 800515B8 0C00B721 */  jal   Player_GetCameraYOffset
/* AC875C 800515BC 8C840090 */   lw    $a0, 0x90($a0)
/* AC8760 800515C0 8602015E */  lh    $v0, 0x15e($s0)
/* AC8764 800515C4 10400008 */  beqz  $v0, .L800515E8
/* AC8768 800515C8 2401000A */   li    $at, 10
/* AC876C 800515CC 10410006 */  beq   $v0, $at, .L800515E8
/* AC8770 800515D0 24010014 */   li    $at, 20
/* AC8774 800515D4 10410004 */  beq   $v0, $at, .L800515E8
/* AC8778 800515D8 3C0E8016 */   lui   $t6, %hi(gGameInfo) # $t6, 0x8016
/* AC877C 800515DC 8DCEFA90 */  lw    $t6, %lo(gGameInfo)($t6)
/* AC8780 800515E0 85C30314 */  lh    $v1, 0x314($t6)
/* AC8784 800515E4 10600054 */  beqz  $v1, .L80051738
.L800515E8:
/* AC8788 800515E8 3C018014 */   lui   $at, %hi(D_8013A1C4)
/* AC878C 800515EC C42EA1C4 */  lwc1  $f14, %lo(D_8013A1C4)($at)
/* AC8790 800515F0 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* AC8794 800515F4 44814000 */  mtc1  $at, $f8
/* AC8798 800515F8 3C014288 */  lui   $at, 0x4288
/* AC879C 800515FC 3C0B8016 */  lui   $t3, %hi(gGameInfo)
/* AC87A0 80051600 8D6BFA90 */  lw    $t3, %lo(gGameInfo)($t3)
/* AC87A4 80051604 44818000 */  mtc1  $at, $f16
/* AC87A8 80051608 860F0142 */  lh    $t7, 0x142($s0)
/* AC87AC 8005160C 856C01F0 */  lh    $t4, 0x1f0($t3)
/* AC87B0 80051610 46008483 */  div.s $f18, $f16, $f0
/* AC87B4 80051614 3C198012 */  lui   $t9, %hi(sCameraSettings+4)
/* AC87B8 80051618 448C2000 */  mtc1  $t4, $f4
/* AC87BC 8005161C 000FC0C0 */  sll   $t8, $t7, 3
/* AC87C0 80051620 86080144 */  lh    $t0, 0x144($s0)
/* AC87C4 80051624 468021A0 */  cvt.s.w $f6, $f4
/* AC87C8 80051628 0338C821 */  addu  $t9, $t9, $t8
/* AC87CC 8005162C 8F39D068 */  lw    $t9, %lo(sCameraSettings+4)($t9)
/* AC87D0 80051630 000848C0 */  sll   $t1, $t0, 3
/* AC87D4 80051634 3C018014 */  lui   $at, %hi(D_8013A1C8)
/* AC87D8 80051638 03295021 */  addu  $t2, $t9, $t1
/* AC87DC 8005163C 8D420004 */  lw    $v0, 4($t2)
/* AC87E0 80051640 3C0C8016 */  lui   $t4, %hi(gGameInfo)
/* AC87E4 80051644 844D0000 */  lh    $t5, ($v0)
/* AC87E8 80051648 24420018 */  addiu $v0, $v0, 0x18
/* AC87EC 8005164C 460E3302 */  mul.s $f12, $f6, $f14
/* AC87F0 80051650 448D3000 */  mtc1  $t5, $f6
/* AC87F4 80051654 460C4280 */  add.s $f10, $f8, $f12
/* AC87F8 80051658 46126102 */  mul.s $f4, $f12, $f18
/* AC87FC 8005165C 46803220 */  cvt.s.w $f8, $f6
/* AC8800 80051660 46045081 */  sub.s $f2, $f10, $f4
/* AC8804 80051664 460E4402 */  mul.s $f16, $f8, $f14
/* AC8808 80051668 00000000 */  nop
/* AC880C 8005166C 46008482 */  mul.s $f18, $f16, $f0
/* AC8810 80051670 00000000 */  nop
/* AC8814 80051674 46029282 */  mul.s $f10, $f18, $f2
/* AC8818 80051678 E60A0000 */  swc1  $f10, ($s0)
/* AC881C 8005167C 844EFFEC */  lh    $t6, -0x14($v0)
/* AC8820 80051680 448E2000 */  mtc1  $t6, $f4
/* AC8824 80051684 00000000 */  nop
/* AC8828 80051688 468021A0 */  cvt.s.w $f6, $f4
/* AC882C 8005168C 460E3202 */  mul.s $f8, $f6, $f14
/* AC8830 80051690 00000000 */  nop
/* AC8834 80051694 46004402 */  mul.s $f16, $f8, $f0
/* AC8838 80051698 00000000 */  nop
/* AC883C 8005169C 46028482 */  mul.s $f18, $f16, $f2
/* AC8840 800516A0 E6120004 */  swc1  $f18, 4($s0)
/* AC8844 800516A4 844FFFF0 */  lh    $t7, -0x10($v0)
/* AC8848 800516A8 448F5000 */  mtc1  $t7, $f10
/* AC884C 800516AC 00000000 */  nop
/* AC8850 800516B0 46805120 */  cvt.s.w $f4, $f10
/* AC8854 800516B4 460E2182 */  mul.s $f6, $f4, $f14
/* AC8858 800516B8 00000000 */  nop
/* AC885C 800516BC 46003202 */  mul.s $f8, $f6, $f0
/* AC8860 800516C0 00000000 */  nop
/* AC8864 800516C4 46024402 */  mul.s $f16, $f8, $f2
/* AC8868 800516C8 E6100008 */  swc1  $f16, 8($s0)
/* AC886C 800516CC 8458FFF4 */  lh    $t8, -0xc($v0)
/* AC8870 800516D0 C424A1C8 */  lwc1  $f4, %lo(D_8013A1C8)($at)
/* AC8874 800516D4 3C013F00 */  li    $at, 0x3F000000 # 0.000000
/* AC8878 800516D8 44989000 */  mtc1  $t8, $f18
/* AC887C 800516DC 44814000 */  mtc1  $at, $f8
/* AC8880 800516E0 468092A0 */  cvt.s.w $f10, $f18
/* AC8884 800516E4 46045182 */  mul.s $f6, $f10, $f4
/* AC8888 800516E8 46083400 */  add.s $f16, $f6, $f8
/* AC888C 800516EC 4600848D */  trunc.w.s $f18, $f16
/* AC8890 800516F0 44199000 */  mfc1  $t9, $f18
/* AC8894 800516F4 00000000 */  nop
/* AC8898 800516F8 A6190018 */  sh    $t9, 0x18($s0)
/* AC889C 800516FC 8449FFF8 */  lh    $t1, -8($v0)
/* AC88A0 80051700 44895000 */  mtc1  $t1, $f10
/* AC88A4 80051704 00000000 */  nop
/* AC88A8 80051708 46805120 */  cvt.s.w $f4, $f10
/* AC88AC 8005170C E6040010 */  swc1  $f4, 0x10($s0)
/* AC88B0 80051710 844AFFFC */  lh    $t2, -4($v0)
/* AC88B4 80051714 448A3000 */  mtc1  $t2, $f6
/* AC88B8 80051718 00000000 */  nop
/* AC88BC 8005171C 46803220 */  cvt.s.w $f8, $f6
/* AC88C0 80051720 460E4402 */  mul.s $f16, $f8, $f14
/* AC88C4 80051724 E6100014 */  swc1  $f16, 0x14($s0)
/* AC88C8 80051728 844B0000 */  lh    $t3, ($v0)
/* AC88CC 8005172C A60B001A */  sh    $t3, 0x1a($s0)
/* AC88D0 80051730 8D8CFA90 */  lw    $t4, %lo(gGameInfo)($t4)
/* AC88D4 80051734 85830314 */  lh    $v1, 0x314($t4)
.L80051738:
/* AC88D8 80051738 50600004 */  beql  $v1, $zero, .L8005174C
/* AC88DC 8005173C 240D0001 */   li    $t5, 1
/* AC88E0 80051740 0C011495 */  jal   Camera_CopyPREGToModeValues
/* AC88E4 80051744 02002025 */   move  $a0, $s0
/* AC88E8 80051748 240D0001 */  li    $t5, 1
.L8005174C:
/* AC88EC 8005174C 3C018012 */  lui   $at, %hi(D_8011D3E8) # $at, 0x8012
/* AC88F0 80051750 26050050 */  addiu $a1, $s0, 0x50
/* AC88F4 80051754 2606005C */  addiu $a2, $s0, 0x5c
/* AC88F8 80051758 AC2DD3E8 */  sw    $t5, %lo(D_8011D3E8)($at)
/* AC88FC 8005175C AFA60038 */  sw    $a2, 0x38($sp)
/* AC8900 80051760 AFA50034 */  sw    $a1, 0x34($sp)
/* AC8904 80051764 0C01F124 */  jal   OLib_Vec3fDiffToVecSphRot90
/* AC8908 80051768 27A4007C */   addiu $a0, $sp, 0x7c
/* AC890C 8005176C 26060074 */  addiu $a2, $s0, 0x74
/* AC8910 80051770 AFA60030 */  sw    $a2, 0x30($sp)
/* AC8914 80051774 27A40074 */  addiu $a0, $sp, 0x74
/* AC8918 80051778 0C01F124 */  jal   OLib_Vec3fDiffToVecSphRot90
/* AC891C 8005177C 8FA50034 */   lw    $a1, 0x34($sp)
/* AC8920 80051780 860E001A */  lh    $t6, 0x1a($s0)
/* AC8924 80051784 3C018012 */  lui   $at, %hi(D_8011D3A0)
/* AC8928 80051788 AC2ED3A0 */  sw    $t6, %lo(D_8011D3A0)($at)
/* AC892C 8005178C 860F015E */  lh    $t7, 0x15e($s0)
/* AC8930 80051790 55E0003E */  bnezl $t7, .L8005188C
/* AC8934 80051794 27A4005C */   addiu $a0, $sp, 0x5c
/* AC8938 80051798 C61200E8 */  lwc1  $f18, 0xe8($s0)
/* AC893C 8005179C C60A00F4 */  lwc1  $f10, 0xf4($s0)
/* AC8940 800517A0 44803000 */  mtc1  $zero, $f6
/* AC8944 800517A4 2611001C */  addiu $s1, $s0, 0x1c
/* AC8948 800517A8 460A9101 */  sub.s $f4, $f18, $f10
/* AC894C 800517AC 27A60098 */  addiu $a2, $sp, 0x98
/* AC8950 800517B0 27A40084 */  addiu $a0, $sp, 0x84
/* AC8954 800517B4 26050094 */  addiu $a1, $s0, 0x94
/* AC8958 800517B8 E60400E8 */  swc1  $f4, 0xe8($s0)
/* AC895C 800517BC 87B8007A */  lh    $t8, 0x7a($sp)
/* AC8960 800517C0 E6260000 */  swc1  $f6, ($s1)
/* AC8964 800517C4 A6380004 */  sh    $t8, 4($s1)
/* AC8968 800517C8 8E080090 */  lw    $t0, 0x90($s0)
/* AC896C 800517CC 8D090908 */  lw    $t1, 0x908($t0)
/* AC8970 800517D0 ACC90000 */  sw    $t1, ($a2)
/* AC8974 800517D4 8D19090C */  lw    $t9, 0x90c($t0)
/* AC8978 800517D8 ACD90004 */  sw    $t9, 4($a2)
/* AC897C 800517DC 8D090910 */  lw    $t1, 0x910($t0)
/* AC8980 800517E0 0C01F124 */  jal   OLib_Vec3fDiffToVecSphRot90
/* AC8984 800517E4 ACC90008 */   sw    $t1, 8($a2)
/* AC8988 800517E8 3C0A8016 */  lui   $t2, %hi(gGameInfo) # $t2, 0x8016
/* AC898C 800517EC 8D4AFA90 */  lw    $t2, %lo(gGameInfo)($t2)
/* AC8990 800517F0 854B01C2 */  lh    $t3, 0x1c2($t2)
/* AC8994 800517F4 A62B0008 */  sh    $t3, 8($s1)
/* AC8998 800517F8 87AD0082 */  lh    $t5, 0x82($sp)
/* AC899C 800517FC 87AC008A */  lh    $t4, 0x8a($sp)
/* AC89A0 80051800 018D1823 */  subu  $v1, $t4, $t5
/* AC89A4 80051804 00031C00 */  sll   $v1, $v1, 0x10
/* AC89A8 80051808 00031C03 */  sra   $v1, $v1, 0x10
/* AC89AC 8005180C 04600003 */  bltz  $v1, .L8005181C
/* AC89B0 80051810 00031023 */   negu  $v0, $v1
/* AC89B4 80051814 10000001 */  b     .L8005181C
/* AC89B8 80051818 00601025 */   move  $v0, $v1
.L8005181C:
/* AC89BC 8005181C 28413A98 */  slti  $at, $v0, 0x3a98
/* AC89C0 80051820 50200004 */  beql  $at, $zero, .L80051834
/* AC89C4 80051824 862E0008 */   lh    $t6, 8($s1)
/* AC89C8 80051828 10000014 */  b     .L8005187C
/* AC89CC 8005182C A6200006 */   sh    $zero, 6($s1)
/* AC89D0 80051830 862E0008 */  lh    $t6, 8($s1)
.L80051834:
/* AC89D4 80051834 006E001A */  div   $zero, $v1, $t6
/* AC89D8 80051838 15C00002 */  bnez  $t6, .L80051844
/* AC89DC 8005183C 00000000 */   nop
/* AC89E0 80051840 0007000D */  break 7
.L80051844:
/* AC89E4 80051844 2401FFFF */  li    $at, -1
/* AC89E8 80051848 15C10004 */  bne   $t6, $at, .L8005185C
/* AC89EC 8005184C 3C018000 */   lui   $at, 0x8000
/* AC89F0 80051850 14610002 */  bne   $v1, $at, .L8005185C
/* AC89F4 80051854 00000000 */   nop
/* AC89F8 80051858 0006000D */  break 6
.L8005185C:
/* AC89FC 8005185C 00007812 */  mflo  $t7
/* AC8A00 80051860 05E10003 */  bgez  $t7, .L80051870
/* AC8A04 80051864 000FC083 */   sra   $t8, $t7, 2
/* AC8A08 80051868 25E10003 */  addiu $at, $t7, 3
/* AC8A0C 8005186C 0001C083 */  sra   $t8, $at, 2
.L80051870:
/* AC8A10 80051870 00184080 */  sll   $t0, $t8, 2
/* AC8A14 80051874 01184023 */  subu  $t0, $t0, $t8
/* AC8A18 80051878 A6280006 */  sh    $t0, 6($s1)
.L8005187C:
/* AC8A1C 8005187C 8619015E */  lh    $t9, 0x15e($s0)
/* AC8A20 80051880 27290001 */  addiu $t1, $t9, 1
/* AC8A24 80051884 A609015E */  sh    $t1, 0x15e($s0)
/* AC8A28 80051888 27A4005C */  addiu $a0, $sp, 0x5c
.L8005188C:
/* AC8A2C 8005188C 8E050090 */  lw    $a1, 0x90($s0)
/* AC8A30 80051890 0C00BBB9 */  jal   func_8002EEE4
/* AC8A34 80051894 2611001C */   addiu $s1, $s0, 0x1c
/* AC8A38 80051898 3C0A8016 */  lui   $t2, %hi(gGameInfo) # $t2, 0x8016
/* AC8A3C 8005189C 8D4AFA90 */  lw    $t2, %lo(gGameInfo)($t2)
/* AC8A40 800518A0 3C0142C8 */  li    $at, 0x42C80000 # 0.000000
/* AC8A44 800518A4 44816000 */  mtc1  $at, $f12
/* AC8A48 800518A8 854B01C6 */  lh    $t3, 0x1c6($t2)
/* AC8A4C 800518AC 3C018014 */  lui   $at, %hi(D_8013A1D4)
/* AC8A50 800518B0 C432A1D4 */  lwc1  $f18, %lo(D_8013A1D4)($at)
/* AC8A54 800518B4 448B4000 */  mtc1  $t3, $f8
/* AC8A58 800518B8 3C073DCC */  lui   $a3, (0x3DCCCCCD >> 16) # lui $a3, 0x3dcc
/* AC8A5C 800518BC 34E7CCCD */  ori   $a3, (0x3DCCCCCD & 0xFFFF) # ori $a3, $a3, 0xcccd
/* AC8A60 800518C0 46804420 */  cvt.s.w $f16, $f8
/* AC8A64 800518C4 C60E00C8 */  lwc1  $f14, 0xc8($s0)
/* AC8A68 800518C8 46128282 */  mul.s $f10, $f16, $f18
/* AC8A6C 800518CC 44065000 */  mfc1  $a2, $f10
/* AC8A70 800518D0 0C010E27 */  jal   func_8004389C
/* AC8A74 800518D4 00000000 */   nop
/* AC8A78 800518D8 E60000C8 */  swc1  $f0, 0xc8($s0)
/* AC8A7C 800518DC 3C0C8016 */  lui   $t4, %hi(gGameInfo) # $t4, 0x8016
/* AC8A80 800518E0 8D8CFA90 */  lw    $t4, %lo(gGameInfo)($t4)
/* AC8A84 800518E4 3C0142C8 */  li    $at, 0x42C80000 # 0.000000
/* AC8A88 800518E8 44816000 */  mtc1  $at, $f12
/* AC8A8C 800518EC 858D01C6 */  lh    $t5, 0x1c6($t4)
/* AC8A90 800518F0 3C018014 */  lui   $at, %hi(D_8013A1D8)
/* AC8A94 800518F4 C428A1D8 */  lwc1  $f8, %lo(D_8013A1D8)($at)
/* AC8A98 800518F8 448D2000 */  mtc1  $t5, $f4
/* AC8A9C 800518FC 3C073DCC */  lui   $a3, (0x3DCCCCCD >> 16) # lui $a3, 0x3dcc
/* AC8AA0 80051900 34E7CCCD */  ori   $a3, (0x3DCCCCCD & 0xFFFF) # ori $a3, $a3, 0xcccd
/* AC8AA4 80051904 468021A0 */  cvt.s.w $f6, $f4
/* AC8AA8 80051908 C60E00C4 */  lwc1  $f14, 0xc4($s0)
/* AC8AAC 8005190C 46083402 */  mul.s $f16, $f6, $f8
/* AC8AB0 80051910 44068000 */  mfc1  $a2, $f16
/* AC8AB4 80051914 0C010E27 */  jal   func_8004389C
/* AC8AB8 80051918 00000000 */   nop
/* AC8ABC 8005191C 3C018014 */  lui   $at, %hi(D_8013A1DC)
/* AC8AC0 80051920 C422A1DC */  lwc1  $f2, %lo(D_8013A1DC)($at)
/* AC8AC4 80051924 E60000C4 */  swc1  $f0, 0xc4($s0)
/* AC8AC8 80051928 3C0E8016 */  lui   $t6, %hi(gGameInfo) # $t6, 0x8016
/* AC8ACC 8005192C 8DCEFA90 */  lw    $t6, %lo(gGameInfo)($t6)
/* AC8AD0 80051930 3C018014 */  lui   $at, %hi(D_8013A1E0)
/* AC8AD4 80051934 44071000 */  mfc1  $a3, $f2
/* AC8AD8 80051938 85CF01C6 */  lh    $t7, 0x1c6($t6)
/* AC8ADC 8005193C C42CA1E0 */  lwc1  $f12, %lo(D_8013A1E0)($at)
/* AC8AE0 80051940 C60E00CC */  lwc1  $f14, 0xcc($s0)
/* AC8AE4 80051944 448F9000 */  mtc1  $t7, $f18
/* AC8AE8 80051948 00000000 */  nop
/* AC8AEC 8005194C 468092A0 */  cvt.s.w $f10, $f18
/* AC8AF0 80051950 46025102 */  mul.s $f4, $f10, $f2
/* AC8AF4 80051954 44062000 */  mfc1  $a2, $f4
/* AC8AF8 80051958 0C010E27 */  jal   func_8004389C
/* AC8AFC 8005195C 00000000 */   nop
/* AC8B00 80051960 3C018014 */  lui   $at, %hi(D_8013A1E4)
/* AC8B04 80051964 C42CA1E4 */  lwc1  $f12, %lo(D_8013A1E4)($at)
/* AC8B08 80051968 E60000CC */  swc1  $f0, 0xcc($s0)
/* AC8B0C 8005196C 3C188016 */  lui   $t8, %hi(gGameInfo) # $t8, 0x8016
/* AC8B10 80051970 8F18FA90 */  lw    $t8, %lo(gGameInfo)($t8)
/* AC8B14 80051974 44076000 */  mfc1  $a3, $f12
/* AC8B18 80051978 C60E00D0 */  lwc1  $f14, 0xd0($s0)
/* AC8B1C 8005197C 870801C8 */  lh    $t0, 0x1c8($t8)
/* AC8B20 80051980 44883000 */  mtc1  $t0, $f6
/* AC8B24 80051984 00000000 */  nop
/* AC8B28 80051988 46803220 */  cvt.s.w $f8, $f6
/* AC8B2C 8005198C 460C4402 */  mul.s $f16, $f8, $f12
/* AC8B30 80051990 44068000 */  mfc1  $a2, $f16
/* AC8B34 80051994 0C010E27 */  jal   func_8004389C
/* AC8B38 80051998 00000000 */   nop
/* AC8B3C 8005199C E60000D0 */  swc1  $f0, 0xd0($s0)
/* AC8B40 800519A0 3C198016 */  lui   $t9, %hi(gGameInfo) # $t9, 0x8016
/* AC8B44 800519A4 8F39FA90 */  lw    $t9, %lo(gGameInfo)($t9)
/* AC8B48 800519A8 3C018014 */  lui   $at, %hi(D_8013A1E8)
/* AC8B4C 800519AC C424A1E8 */  lwc1  $f4, %lo(D_8013A1E8)($at)
/* AC8B50 800519B0 8729019C */  lh    $t1, 0x19c($t9)
/* AC8B54 800519B4 3C063D4C */  lui   $a2, (0x3D4CCCCD >> 16) # lui $a2, 0x3d4c
/* AC8B58 800519B8 3C073DCC */  li    $a3, 0x3DCC0000 # 0.000000
/* AC8B5C 800519BC 44899000 */  mtc1  $t1, $f18
/* AC8B60 800519C0 34E7CCCD */  ori   $a3, (0x3DCCCCCD & 0xFFFF) # ori $a3, $a3, 0xcccd
/* AC8B64 800519C4 34C6CCCD */  ori   $a2, (0x3D4CCCCD & 0xFFFF) # ori $a2, $a2, 0xcccd
/* AC8B68 800519C8 468092A0 */  cvt.s.w $f10, $f18
/* AC8B6C 800519CC C60E00D4 */  lwc1  $f14, 0xd4($s0)
/* AC8B70 800519D0 46045302 */  mul.s $f12, $f10, $f4
/* AC8B74 800519D4 0C010E27 */  jal   func_8004389C
/* AC8B78 800519D8 00000000 */   nop
/* AC8B7C 800519DC E60000D4 */  swc1  $f0, 0xd4($s0)
/* AC8B80 800519E0 02002025 */  move  $a0, $s0
/* AC8B84 800519E4 27A50074 */  addiu $a1, $sp, 0x74
/* AC8B88 800519E8 8E060000 */  lw    $a2, ($s0)
/* AC8B8C 800519EC 0C0115EA */  jal   func_800457A8
/* AC8B90 800519F0 24070001 */   li    $a3, 1
/* AC8B94 800519F4 27A4008C */  addiu $a0, $sp, 0x8c
/* AC8B98 800519F8 8FA50034 */  lw    $a1, 0x34($sp)
/* AC8B9C 800519FC 0C01F124 */  jal   OLib_Vec3fDiffToVecSphRot90
/* AC8BA0 80051A00 8FA60030 */   lw    $a2, 0x30($sp)
/* AC8BA4 80051A04 02002025 */  move  $a0, $s0
/* AC8BA8 80051A08 8FA5008C */  lw    $a1, 0x8c($sp)
/* AC8BAC 80051A0C 8E060004 */  lw    $a2, 4($s0)
/* AC8BB0 80051A10 0C011A33 */  jal   func_800468CC
/* AC8BB4 80051A14 8E070008 */   lw    $a3, 8($s0)
/* AC8BB8 80051A18 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* AC8BBC 80051A1C 44813000 */  mtc1  $at, $f6
/* AC8BC0 80051A20 C60800C4 */  lwc1  $f8, 0xc4($s0)
/* AC8BC4 80051A24 E60000DC */  swc1  $f0, 0xdc($s0)
/* AC8BC8 80051A28 86040018 */  lh    $a0, 0x18($s0)
/* AC8BCC 80051A2C 46083403 */  div.s $f16, $f6, $f8
/* AC8BD0 80051A30 87A50078 */  lh    $a1, 0x78($sp)
/* AC8BD4 80051A34 2407000A */  li    $a3, 10
/* AC8BD8 80051A38 44068000 */  mfc1  $a2, $f16
/* AC8BDC 80051A3C 0C010E47 */  jal   func_8004391C
/* AC8BE0 80051A40 00000000 */   nop
/* AC8BE4 80051A44 3C048016 */  lui   $a0, %hi(gGameInfo) # $a0, 0x8016
/* AC8BE8 80051A48 8C84FA90 */  lw    $a0, %lo(gGameInfo)($a0)
/* AC8BEC 80051A4C A7A20090 */  sh    $v0, 0x90($sp)
/* AC8BF0 80051A50 3C063F00 */  lui   $a2, 0x3f00
/* AC8BF4 80051A54 8483019E */  lh    $v1, 0x19e($a0)
/* AC8BF8 80051A58 24072710 */  li    $a3, 10000
/* AC8BFC 80051A5C 0062082A */  slt   $at, $v1, $v0
/* AC8C00 80051A60 50200004 */  beql  $at, $zero, .L80051A74
/* AC8C04 80051A64 87AA0090 */   lh    $t2, 0x90($sp)
/* AC8C08 80051A68 A7A30090 */  sh    $v1, 0x90($sp)
/* AC8C0C 80051A6C 8483019E */  lh    $v1, 0x19e($a0)
/* AC8C10 80051A70 87AA0090 */  lh    $t2, 0x90($sp)
.L80051A74:
/* AC8C14 80051A74 00031023 */  negu  $v0, $v1
/* AC8C18 80051A78 0142082A */  slt   $at, $t2, $v0
/* AC8C1C 80051A7C 50200003 */  beql  $at, $zero, .L80051A8C
/* AC8C20 80051A80 86220008 */   lh    $v0, 8($s1)
/* AC8C24 80051A84 A7A20090 */  sh    $v0, 0x90($sp)
/* AC8C28 80051A88 86220008 */  lh    $v0, 8($s1)
.L80051A8C:
/* AC8C2C 80051A8C 50400008 */  beql  $v0, $zero, .L80051AB0
/* AC8C30 80051A90 86240004 */   lh    $a0, 4($s1)
/* AC8C34 80051A94 862B0004 */  lh    $t3, 4($s1)
/* AC8C38 80051A98 862C0006 */  lh    $t4, 6($s1)
/* AC8C3C 80051A9C 244EFFFF */  addiu $t6, $v0, -1
/* AC8C40 80051AA0 A62E0008 */  sh    $t6, 8($s1)
/* AC8C44 80051AA4 016C6821 */  addu  $t5, $t3, $t4
/* AC8C48 80051AA8 A62D0004 */  sh    $t5, 4($s1)
/* AC8C4C 80051AAC 86240004 */  lh    $a0, 4($s1)
.L80051AB0:
/* AC8C50 80051AB0 0C010E6B */  jal   func_800439AC
/* AC8C54 80051AB4 87A5007A */   lh    $a1, 0x7a($sp)
/* AC8C58 80051AB8 A7A20092 */  sh    $v0, 0x92($sp)
/* AC8C5C 80051ABC 8FA40030 */  lw    $a0, 0x30($sp)
/* AC8C60 80051AC0 8FA50034 */  lw    $a1, 0x34($sp)
/* AC8C64 80051AC4 0C010F0A */  jal   func_80043C28
/* AC8C68 80051AC8 27A6008C */   addiu $a2, $sp, 0x8c
/* AC8C6C 80051ACC 8FAF0030 */  lw    $t7, 0x30($sp)
/* AC8C70 80051AD0 8FA60038 */  lw    $a2, 0x38($sp)
/* AC8C74 80051AD4 02002025 */  move  $a0, $s0
/* AC8C78 80051AD8 8DE80000 */  lw    $t0, ($t7)
/* AC8C7C 80051ADC ACC80000 */  sw    $t0, ($a2)
/* AC8C80 80051AE0 8DF80004 */  lw    $t8, 4($t7)
/* AC8C84 80051AE4 ACD80004 */  sw    $t8, 4($a2)
/* AC8C88 80051AE8 8DE80008 */  lw    $t0, 8($t7)
/* AC8C8C 80051AEC ACC80008 */  sw    $t0, 8($a2)
/* AC8C90 80051AF0 0C010FCD */  jal   func_80043F34
/* AC8C94 80051AF4 8FA50034 */   lw    $a1, 0x34($sp)
/* AC8C98 80051AF8 C60C0010 */  lwc1  $f12, 0x10($s0)
/* AC8C9C 80051AFC C60E00FC */  lwc1  $f14, 0xfc($s0)
/* AC8CA0 80051B00 8E0600D4 */  lw    $a2, 0xd4($s0)
/* AC8CA4 80051B04 0C010E27 */  jal   func_8004389C
/* AC8CA8 80051B08 3C073F80 */   lui   $a3, 0x3f80
/* AC8CAC 80051B0C E60000FC */  swc1  $f0, 0xfc($s0)
/* AC8CB0 80051B10 A600015A */  sh    $zero, 0x15a($s0)
/* AC8CB4 80051B14 02002025 */  move  $a0, $s0
/* AC8CB8 80051B18 0C011429 */  jal   func_800450A4
/* AC8CBC 80051B1C 8E050014 */   lw    $a1, 0x14($s0)
/* AC8CC0 80051B20 E6000100 */  swc1  $f0, 0x100($s0)
/* AC8CC4 80051B24 8FBF001C */  lw    $ra, 0x1c($sp)
/* AC8CC8 80051B28 8FB10018 */  lw    $s1, 0x18($sp)
/* AC8CCC 80051B2C 8FB00014 */  lw    $s0, 0x14($sp)
/* AC8CD0 80051B30 27BD00B0 */  addiu $sp, $sp, 0xb0
/* AC8CD4 80051B34 03E00008 */  jr    $ra
/* AC8CD8 80051B38 24020001 */   li    $v0, 1
