.late_rodata
glabel D_80139DC0
    .float 0.01

glabel D_80139DC4
    .float 0.01

glabel D_80139DC8
    .float 57.295776

glabel D_80139DCC
    .float 182.041672

glabel D_80139DD0
    .float 57.295776

glabel D_80139DD4
    .float 182.041672

.text
glabel func_80044ADC
/* ABBC7C 80044ADC 27BDFF90 */  addiu $sp, $sp, -0x70
/* ABBC80 80044AE0 AFB00018 */  sw    $s0, 0x18($sp)
/* ABBC84 80044AE4 00808025 */  move  $s0, $a0
/* ABBC88 80044AE8 AFBF001C */  sw    $ra, 0x1c($sp)
/* ABBC8C 80044AEC AFA50074 */  sw    $a1, 0x74($sp)
/* ABBC90 80044AF0 AFA60078 */  sw    $a2, 0x78($sp)
/* ABBC94 80044AF4 0C01DE1C */  jal   Math_Sins
/* ABBC98 80044AF8 87A40076 */   lh    $a0, 0x76($sp)
/* ABBC9C 80044AFC E7A0003C */  swc1  $f0, 0x3c($sp)
/* ABBCA0 80044B00 0C01DE0D */  jal   Math_Coss
/* ABBCA4 80044B04 87A40076 */   lh    $a0, 0x76($sp)
/* ABBCA8 80044B08 8E040090 */  lw    $a0, 0x90($s0)
/* ABBCAC 80044B0C 0C00B721 */  jal   Player_GetCameraYOffset
/* ABBCB0 80044B10 E7A00038 */   swc1  $f0, 0x38($sp)
/* ABBCB4 80044B14 3C028016 */  lui   $v0, %hi(gGameInfo) # $v0, 0x8016
/* ABBCB8 80044B18 8C42FA90 */  lw    $v0, %lo(gGameInfo)($v0)
/* ABBCBC 80044B1C 3C018014 */  lui   $at, %hi(D_80139DC0)
/* ABBCC0 80044B20 C4329DC0 */  lwc1  $f18, %lo(D_80139DC0)($at)
/* ABBCC4 80044B24 844E01BA */  lh    $t6, 0x1ba($v0)
/* ABBCC8 80044B28 844F01B6 */  lh    $t7, 0x1b6($v0)
/* ABBCCC 80044B2C 845801B8 */  lh    $t8, 0x1b8($v0)
/* ABBCD0 80044B30 448E2000 */  mtc1  $t6, $f4
/* ABBCD4 80044B34 448F5000 */  mtc1  $t7, $f10
/* ABBCD8 80044B38 C7AC0038 */  lwc1  $f12, 0x38($sp)
/* ABBCDC 80044B3C 468021A0 */  cvt.s.w $f6, $f4
/* ABBCE0 80044B40 87B9007A */  lh    $t9, 0x7a($sp)
/* ABBCE4 80044B44 3C018016 */  lui   $at, %hi(D_8015CE5C)
/* ABBCE8 80044B48 46805120 */  cvt.s.w $f4, $f10
/* ABBCEC 80044B4C 46123202 */  mul.s $f8, $f6, $f18
/* ABBCF0 80044B50 00000000 */  nop
/* ABBCF4 80044B54 46004082 */  mul.s $f2, $f8, $f0
/* ABBCF8 80044B58 44984000 */  mtc1  $t8, $f8
/* ABBCFC 80044B5C 46122182 */  mul.s $f6, $f4, $f18
/* ABBD00 80044B60 468042A0 */  cvt.s.w $f10, $f8
/* ABBD04 80044B64 46003402 */  mul.s $f16, $f6, $f0
/* ABBD08 80044B68 C6060094 */  lwc1  $f6, 0x94($s0)
/* ABBD0C 80044B6C E7A60064 */  swc1  $f6, 0x64($sp)
/* ABBD10 80044B70 46125102 */  mul.s $f4, $f10, $f18
/* ABBD14 80044B74 C6080104 */  lwc1  $f8, 0x104($s0)
/* ABBD18 80044B78 C7A6003C */  lwc1  $f6, 0x3c($sp)
/* ABBD1C 80044B7C 46024280 */  add.s $f10, $f8, $f2
/* ABBD20 80044B80 46002382 */  mul.s $f14, $f4, $f0
/* ABBD24 80044B84 E7AA0068 */  swc1  $f10, 0x68($sp)
/* ABBD28 80044B88 46068202 */  mul.s $f8, $f16, $f6
/* ABBD2C 80044B8C C604009C */  lwc1  $f4, 0x9c($s0)
/* ABBD30 80044B90 C7AA0064 */  lwc1  $f10, 0x64($sp)
/* ABBD34 80044B94 C7A60068 */  lwc1  $f6, 0x68($sp)
/* ABBD38 80044B98 E7A4006C */  swc1  $f4, 0x6c($sp)
/* ABBD3C 80044B9C E7A6005C */  swc1  $f6, 0x5c($sp)
/* ABBD40 80044BA0 460A4100 */  add.s $f4, $f8, $f10
/* ABBD44 80044BA4 460C8202 */  mul.s $f8, $f16, $f12
/* ABBD48 80044BA8 C7AA006C */  lwc1  $f10, 0x6c($sp)
/* ABBD4C 80044BAC E7A40058 */  swc1  $f4, 0x58($sp)
/* ABBD50 80044BB0 460A4100 */  add.s $f4, $f8, $f10
/* ABBD54 80044BB4 17200009 */  bnez  $t9, .L80044BDC
/* ABBD58 80044BB8 E7A40060 */   swc1  $f4, 0x60($sp)
/* ABBD5C 80044BBC 8E08008C */  lw    $t0, 0x8c($s0)
/* ABBD60 80044BC0 3C058016 */  lui   $a1, %hi(D_8015CE58) # $a1, 0x8016
/* ABBD64 80044BC4 24A5CE58 */  addiu $a1, %lo(D_8015CE58) # addiu $a1, $a1, -0x31a8
/* ABBD68 80044BC8 8D09009C */  lw    $t1, 0x9c($t0)
/* ABBD6C 80044BCC 27A40064 */  addiu $a0, $sp, 0x64
/* ABBD70 80044BD0 312A0001 */  andi  $t2, $t1, 1
/* ABBD74 80044BD4 15400022 */  bnez  $t2, .L80044C60
/* ABBD78 80044BD8 00000000 */   nop
.L80044BDC:
/* ABBD7C 80044BDC C7B2003C */  lwc1  $f18, 0x3c($sp)
/* ABBD80 80044BE0 C7AA0064 */  lwc1  $f10, 0x64($sp)
/* ABBD84 80044BE4 C7A8006C */  lwc1  $f8, 0x6c($sp)
/* ABBD88 80044BE8 46127482 */  mul.s $f18, $f14, $f18
/* ABBD8C 80044BEC C7A60068 */  lwc1  $f6, 0x68($sp)
/* ABBD90 80044BF0 3C068016 */  lui   $a2, %hi(D_8015CE58) # $a2, 0x8016
/* ABBD94 80044BF4 460C7102 */  mul.s $f4, $f14, $f12
/* ABBD98 80044BF8 E426CE5C */  swc1  $f6, %lo(D_8015CE5C)($at)
/* ABBD9C 80044BFC 02002025 */  move  $a0, $s0
/* ABBDA0 80044C00 27A50064 */  addiu $a1, $sp, 0x64
/* ABBDA4 80044C04 24C6CE58 */  addiu $a2, %lo(D_8015CE58) # addiu $a2, $a2, -0x31a8
/* ABBDA8 80044C08 460A9280 */  add.s $f10, $f18, $f10
/* ABBDAC 80044C0C E7AE002C */  swc1  $f14, 0x2c($sp)
/* ABBDB0 80044C10 E7B00030 */  swc1  $f16, 0x30($sp)
/* ABBDB4 80044C14 46082200 */  add.s $f8, $f4, $f8
/* ABBDB8 80044C18 E42ACE58 */  swc1  $f10, %lo(D_8015CE58)($at)
/* ABBDBC 80044C1C 3C018016 */  lui   $at, %hi(D_8015CE60)
/* ABBDC0 80044C20 0C010F46 */  jal   func_80043D18
/* ABBDC4 80044C24 E428CE60 */   swc1  $f8, %lo(D_8015CE60)($at)
/* ABBDC8 80044C28 87AB007A */  lh    $t3, 0x7a($sp)
/* ABBDCC 80044C2C 11600006 */  beqz  $t3, .L80044C48
/* ABBDD0 80044C30 00000000 */   nop
/* ABBDD4 80044C34 C6120104 */  lwc1  $f18, 0x104($s0)
/* ABBDD8 80044C38 3C018016 */  lui   $at, %hi(D_8015CE54)
/* ABBDDC 80044C3C E432CE54 */  swc1  $f18, %lo(D_8015CE54)($at)
/* ABBDE0 80044C40 3C018016 */  lui   $at, %hi(D_8015CE50)
/* ABBDE4 80044C44 E432CE50 */  swc1  $f18, %lo(D_8015CE50)($at)
.L80044C48:
/* ABBDE8 80044C48 3C018016 */  lui   $at, %hi(D_8015CE54)
/* ABBDEC 80044C4C C430CE54 */  lwc1  $f16, %lo(D_8015CE54)($at)
/* ABBDF0 80044C50 3C018016 */  lui   $at, %hi(D_8015CE50)
/* ABBDF4 80044C54 C432CE50 */  lwc1  $f18, %lo(D_8015CE50)($at)
/* ABBDF8 80044C58 10000048 */  b     .L80044D7C
/* ABBDFC 80044C5C C6020104 */   lwc1  $f2, 0x104($s0)
.L80044C60:
/* ABBE00 80044C60 0C01F00A */  jal   OLib_Vec3fDistXZ
/* ABBE04 80044C64 E7B00030 */   swc1  $f16, 0x30($sp)
/* ABBE08 80044C68 3C0140A0 */  li    $at, 0x40A00000 # 0.000000
/* ABBE0C 80044C6C 44811000 */  mtc1  $at, $f2
/* ABBE10 80044C70 3C018016 */  lui   $at, %hi(D_8015CE58)
/* ABBE14 80044C74 C426CE58 */  lwc1  $f6, %lo(D_8015CE58)($at)
/* ABBE18 80044C78 3C018016 */  lui   $at, %hi(D_8015CE64)
/* ABBE1C 80044C7C C428CE64 */  lwc1  $f8, %lo(D_8015CE64)($at)
/* ABBE20 80044C80 3C018016 */  lui   $at, %hi(D_8015CE5C)
/* ABBE24 80044C84 E7A0002C */  swc1  $f0, 0x2c($sp)
/* ABBE28 80044C88 46024282 */  mul.s $f10, $f8, $f2
/* ABBE2C 80044C8C C428CE5C */  lwc1  $f8, %lo(D_8015CE5C)($at)
/* ABBE30 80044C90 27A70034 */  addiu $a3, $sp, 0x34
/* ABBE34 80044C94 02002025 */  move  $a0, $s0
/* ABBE38 80044C98 27A5004C */  addiu $a1, $sp, 0x4c
/* ABBE3C 80044C9C 27A60058 */  addiu $a2, $sp, 0x58
/* ABBE40 80044CA0 460A3100 */  add.s $f4, $f6, $f10
/* ABBE44 80044CA4 E424CE58 */  swc1  $f4, %lo(D_8015CE58)($at)
/* ABBE48 80044CA8 3C018016 */  lui   $at, %hi(D_8015CE68)
/* ABBE4C 80044CAC C426CE68 */  lwc1  $f6, %lo(D_8015CE68)($at)
/* ABBE50 80044CB0 3C018016 */  lui   $at, %hi(D_8015CE5C)
/* ABBE54 80044CB4 46023282 */  mul.s $f10, $f6, $f2
/* ABBE58 80044CB8 460A4100 */  add.s $f4, $f8, $f10
/* ABBE5C 80044CBC E424CE5C */  swc1  $f4, %lo(D_8015CE5C)($at)
/* ABBE60 80044CC0 3C018016 */  lui   $at, %hi(D_8015CE60)
/* ABBE64 80044CC4 C426CE60 */  lwc1  $f6, %lo(D_8015CE60)($at)
/* ABBE68 80044CC8 3C018016 */  lui   $at, %hi(D_8015CE6C)
/* ABBE6C 80044CCC C428CE6C */  lwc1  $f8, %lo(D_8015CE6C)($at)
/* ABBE70 80044CD0 3C018016 */  lui   $at, %hi(D_8015CE60)
/* ABBE74 80044CD4 46024282 */  mul.s $f10, $f8, $f2
/* ABBE78 80044CD8 C7A80030 */  lwc1  $f8, 0x30($sp)
/* ABBE7C 80044CDC 4608003C */  c.lt.s $f0, $f8
/* ABBE80 80044CE0 460A3100 */  add.s $f4, $f6, $f10
/* ABBE84 80044CE4 C7A6002C */  lwc1  $f6, 0x2c($sp)
/* ABBE88 80044CE8 4500000A */  bc1f  .L80044D14
/* ABBE8C 80044CEC E424CE60 */   swc1  $f4, %lo(D_8015CE60)($at)
/* ABBE90 80044CF0 3C068016 */  lui   $a2, %hi(D_8015CE58) # $a2, 0x8016
/* ABBE94 80044CF4 E7A60030 */  swc1  $f6, 0x30($sp)
/* ABBE98 80044CF8 24C6CE58 */  addiu $a2, %lo(D_8015CE58) # addiu $a2, $a2, -0x31a8
/* ABBE9C 80044CFC 02002025 */  move  $a0, $s0
/* ABBEA0 80044D00 0C01115A */  jal   func_80044568
/* ABBEA4 80044D04 27A5004C */   addiu $a1, $sp, 0x4c
/* ABBEA8 80044D08 46000406 */  mov.s $f16, $f0
/* ABBEAC 80044D0C 1000000E */  b     .L80044D48
/* ABBEB0 80044D10 46000486 */   mov.s $f18, $f0
.L80044D14:
/* ABBEB4 80044D14 0C01115A */  jal   func_80044568
/* ABBEB8 80044D18 27A70034 */   addiu $a3, $sp, 0x34
/* ABBEBC 80044D1C 3C068016 */  lui   $a2, %hi(D_8015CE58) # $a2, 0x8016
/* ABBEC0 80044D20 3C018016 */  lui   $at, %hi(D_8015CE50)
/* ABBEC4 80044D24 E420CE50 */  swc1  $f0, %lo(D_8015CE50)($at)
/* ABBEC8 80044D28 24C6CE58 */  addiu $a2, %lo(D_8015CE58) # addiu $a2, $a2, -0x31a8
/* ABBECC 80044D2C 02002025 */  move  $a0, $s0
/* ABBED0 80044D30 27A5004C */  addiu $a1, $sp, 0x4c
/* ABBED4 80044D34 0C01115A */  jal   func_80044568
/* ABBED8 80044D38 27A70034 */   addiu $a3, $sp, 0x34
/* ABBEDC 80044D3C 3C018016 */  lui   $at, %hi(D_8015CE50)
/* ABBEE0 80044D40 C432CE50 */  lwc1  $f18, %lo(D_8015CE50)($at)
/* ABBEE4 80044D44 46000406 */  mov.s $f16, $f0
.L80044D48:
/* ABBEE8 80044D48 3C01C6FA */  li    $at, 0xC6FA0000 # 0.000000
/* ABBEEC 80044D4C 44816000 */  mtc1  $at, $f12
/* ABBEF0 80044D50 00000000 */  nop
/* ABBEF4 80044D54 460C9032 */  c.eq.s $f18, $f12
/* ABBEF8 80044D58 00000000 */  nop
/* ABBEFC 80044D5C 45020003 */  bc1fl .L80044D6C
/* ABBF00 80044D60 460C0032 */   c.eq.s $f0, $f12
/* ABBF04 80044D64 C6120104 */  lwc1  $f18, 0x104($s0)
/* ABBF08 80044D68 460C0032 */  c.eq.s $f0, $f12
.L80044D6C:
/* ABBF0C 80044D6C C6020104 */  lwc1  $f2, 0x104($s0)
/* ABBF10 80044D70 45000002 */  bc1f  .L80044D7C
/* ABBF14 80044D74 00000000 */   nop
/* ABBF18 80044D78 46009406 */  mov.s $f16, $f18
.L80044D7C:
/* ABBF1C 80044D7C 3C0C8016 */  lui   $t4, %hi(gGameInfo) # $t4, 0x8016
/* ABBF20 80044D80 8D8CFA90 */  lw    $t4, %lo(gGameInfo)($t4)
/* ABBF24 80044D84 3C018014 */  lui   $at, %hi(D_80139DC4)
/* ABBF28 80044D88 C4289DC4 */  lwc1  $f8, %lo(D_80139DC4)($at)
/* ABBF2C 80044D8C 858D01BC */  lh    $t5, 0x1bc($t4)
/* ABBF30 80044D90 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* ABBF34 80044D94 46029181 */  sub.s $f6, $f18, $f2
/* ABBF38 80044D98 448D5000 */  mtc1  $t5, $f10
/* ABBF3C 80044D9C C7AE0030 */  lwc1  $f14, 0x30($sp)
/* ABBF40 80044DA0 46805120 */  cvt.s.w $f4, $f10
/* ABBF44 80044DA4 44815000 */  mtc1  $at, $f10
/* ABBF48 80044DA8 3C018016 */  lui   $at, %hi(D_8015CE54)
/* ABBF4C 80044DAC E430CE54 */  swc1  $f16, %lo(D_8015CE54)($at)
/* ABBF50 80044DB0 3C018016 */  lui   $at, %hi(D_8015CE50)
/* ABBF54 80044DB4 E432CE50 */  swc1  $f18, %lo(D_8015CE50)($at)
/* ABBF58 80044DB8 46082002 */  mul.s $f0, $f4, $f8
/* ABBF5C 80044DBC 46028201 */  sub.s $f8, $f16, $f2
/* ABBF60 80044DC0 46005101 */  sub.s $f4, $f10, $f0
/* ABBF64 80044DC4 46060302 */  mul.s $f12, $f0, $f6
/* ABBF68 80044DC8 00000000 */  nop
/* ABBF6C 80044DCC 46082182 */  mul.s $f6, $f4, $f8
/* ABBF70 80044DD0 0C03F494 */  jal   Math_atan2f
/* ABBF74 80044DD4 E7A60040 */   swc1  $f6, 0x40($sp)
/* ABBF78 80044DD8 3C018014 */  lui   $at, %hi(D_80139DC8)
/* ABBF7C 80044DDC C42A9DC8 */  lwc1  $f10, %lo(D_80139DC8)($at)
/* ABBF80 80044DE0 3C018014 */  lui   $at, %hi(D_80139DCC)
/* ABBF84 80044DE4 C4289DCC */  lwc1  $f8, %lo(D_80139DCC)($at)
/* ABBF88 80044DE8 460A0102 */  mul.s $f4, $f0, $f10
/* ABBF8C 80044DEC 3C013F00 */  li    $at, 0x3F000000 # 0.000000
/* ABBF90 80044DF0 44815000 */  mtc1  $at, $f10
/* ABBF94 80044DF4 C7AC0040 */  lwc1  $f12, 0x40($sp)
/* ABBF98 80044DF8 C7AE002C */  lwc1  $f14, 0x2c($sp)
/* ABBF9C 80044DFC 46082182 */  mul.s $f6, $f4, $f8
/* ABBFA0 80044E00 460A3100 */  add.s $f4, $f6, $f10
/* ABBFA4 80044E04 4600220D */  trunc.w.s $f8, $f4
/* ABBFA8 80044E08 44104000 */  mfc1  $s0, $f8
/* ABBFAC 80044E0C 00000000 */  nop
/* ABBFB0 80044E10 00108400 */  sll   $s0, $s0, 0x10
/* ABBFB4 80044E14 0C03F494 */  jal   Math_atan2f
/* ABBFB8 80044E18 00108403 */   sra   $s0, $s0, 0x10
/* ABBFBC 80044E1C 3C018014 */  lui   $at, %hi(D_80139DD0)
/* ABBFC0 80044E20 C4269DD0 */  lwc1  $f6, %lo(D_80139DD0)($at)
/* ABBFC4 80044E24 3C018014 */  lui   $at, %hi(D_80139DD4)
/* ABBFC8 80044E28 C4249DD4 */  lwc1  $f4, %lo(D_80139DD4)($at)
/* ABBFCC 80044E2C 46060282 */  mul.s $f10, $f0, $f6
/* ABBFD0 80044E30 3C013F00 */  li    $at, 0x3F000000 # 0.000000
/* ABBFD4 80044E34 44813000 */  mtc1  $at, $f6
/* ABBFD8 80044E38 8FBF001C */  lw    $ra, 0x1c($sp)
/* ABBFDC 80044E3C 46045202 */  mul.s $f8, $f10, $f4
/* ABBFE0 80044E40 46064280 */  add.s $f10, $f8, $f6
/* ABBFE4 80044E44 4600510D */  trunc.w.s $f4, $f10
/* ABBFE8 80044E48 44082000 */  mfc1  $t0, $f4
/* ABBFEC 80044E4C 00000000 */  nop
/* ABBFF0 80044E50 02081021 */  addu  $v0, $s0, $t0
/* ABBFF4 80044E54 00021400 */  sll   $v0, $v0, 0x10
/* ABBFF8 80044E58 8FB00018 */  lw    $s0, 0x18($sp)
/* ABBFFC 80044E5C 27BD0070 */  addiu $sp, $sp, 0x70
/* ABC000 80044E60 03E00008 */  jr    $ra
/* ABC004 80044E64 00021403 */   sra   $v0, $v0, 0x10

