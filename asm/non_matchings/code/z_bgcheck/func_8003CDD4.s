glabel func_8003CDD4
/* AB3F74 8003CDD4 27BDFF20 */  addiu $sp, $sp, -0xe0
/* AB3F78 8003CDD8 8FAF00FC */  lw    $t7, 0xfc($sp)
/* AB3F7C 8003CDDC AFBF003C */  sw    $ra, 0x3c($sp)
/* AB3F80 8003CDE0 AFB00038 */  sw    $s0, 0x38($sp)
/* AB3F84 8003CDE4 AFA400E0 */  sw    $a0, 0xe0($sp)
/* AB3F88 8003CDE8 AFA500E4 */  sw    $a1, 0xe4($sp)
/* AB3F8C 8003CDEC AFA000D4 */  sw    $zero, 0xd4($sp)
/* AB3F90 8003CDF0 240E0032 */  li    $t6, 50
/* AB3F94 8003CDF4 ADEE0000 */  sw    $t6, ($t7)
/* AB3F98 8003CDF8 8FB800F8 */  lw    $t8, 0xf8($sp)
/* AB3F9C 8003CDFC 00C08025 */  move  $s0, $a2
/* AB3FA0 8003CE00 3C058014 */  lui   $a1, %hi(D_80138980)
/* AB3FA4 8003CE04 AF000000 */  sw    $zero, ($t8)
/* AB3FA8 8003CE08 8FB900E0 */  lw    $t9, 0xe0($sp)
/* AB3FAC 8003CE0C 24A58980 */  addiu $a1, %lo(D_80138980) # addiu $a1, $a1, -0x7680
/* AB3FB0 8003CE10 00E02025 */  move  $a0, $a3
/* AB3FB4 8003CE14 8F290040 */  lw    $t1, 0x40($t9)
/* AB3FB8 8003CE18 AFA900DC */  sw    $t1, 0xdc($sp)
/* AB3FBC 8003CE1C 8CEB0000 */  lw    $t3, ($a3)
/* AB3FC0 8003CE20 ACCB0000 */  sw    $t3, ($a2)
/* AB3FC4 8003CE24 8CEA0004 */  lw    $t2, 4($a3)
/* AB3FC8 8003CE28 ACCA0004 */  sw    $t2, 4($a2)
/* AB3FCC 8003CE2C 8CEB0008 */  lw    $t3, 8($a3)
/* AB3FD0 8003CE30 ACCB0008 */  sw    $t3, 8($a2)
/* AB3FD4 8003CE34 8FAC00F0 */  lw    $t4, 0xf0($sp)
/* AB3FD8 8003CE38 C4E40000 */  lwc1  $f4, ($a3)
/* AB3FDC 8003CE3C 240612DF */  li    $a2, 4831
/* AB3FE0 8003CE40 C5860000 */  lwc1  $f6, ($t4)
/* AB3FE4 8003CE44 46062201 */  sub.s $f8, $f4, $f6
/* AB3FE8 8003CE48 E7A800CC */  swc1  $f8, 0xcc($sp)
/* AB3FEC 8003CE4C C4EA0004 */  lwc1  $f10, 4($a3)
/* AB3FF0 8003CE50 C5920004 */  lwc1  $f18, 4($t4)
/* AB3FF4 8003CE54 46125101 */  sub.s $f4, $f10, $f18
/* AB3FF8 8003CE58 E7A400C8 */  swc1  $f4, 0xc8($sp)
/* AB3FFC 8003CE5C C4E60008 */  lwc1  $f6, 8($a3)
/* AB4000 8003CE60 C5880008 */  lwc1  $f8, 8($t4)
/* AB4004 8003CE64 AFA700EC */  sw    $a3, 0xec($sp)
/* AB4008 8003CE68 46083281 */  sub.s $f10, $f6, $f8
/* AB400C 8003CE6C 0C00E180 */  jal   func_80038600
/* AB4010 8003CE70 E7AA00C4 */   swc1  $f10, 0xc4($sp)
/* AB4014 8003CE74 24010001 */  li    $at, 1
/* AB4018 8003CE78 10410008 */  beq   $v0, $at, .L8003CE9C
/* AB401C 8003CE7C 8FA400F0 */   lw    $a0, 0xf0($sp)
/* AB4020 8003CE80 3C058014 */  lui   $a1, %hi(D_80138990)
/* AB4024 8003CE84 24A58990 */  addiu $a1, %lo(D_80138990) # addiu $a1, $a1, -0x7670
/* AB4028 8003CE88 0C00E180 */  jal   func_80038600
/* AB402C 8003CE8C 240612E0 */   li    $a2, 4832
/* AB4030 8003CE90 24010001 */  li    $at, 1
/* AB4034 8003CE94 54410009 */  bnel  $v0, $at, .L8003CEBC
/* AB4038 8003CE98 44800000 */   mtc1  $zero, $f0
.L8003CE9C:
/* AB403C 8003CE9C 8FA80100 */  lw    $t0, 0x100($sp)
/* AB4040 8003CEA0 3C048014 */  lui   $a0, %hi(D_801389A0)
/* AB4044 8003CEA4 248489A0 */  addiu $a0, %lo(D_801389A0) # addiu $a0, $a0, -0x7660
/* AB4048 8003CEA8 51000004 */  beql  $t0, $zero, .L8003CEBC
/* AB404C 8003CEAC 44800000 */   mtc1  $zero, $f0
/* AB4050 8003CEB0 0C00084C */  jal   osSyncPrintf
/* AB4054 8003CEB4 85050000 */   lh    $a1, ($t0)
/* AB4058 8003CEB8 44800000 */  mtc1  $zero, $f0
.L8003CEBC:
/* AB405C 8003CEBC C7A200CC */  lwc1  $f2, 0xcc($sp)
/* AB4060 8003CEC0 93AD010B */  lbu   $t5, 0x10b($sp)
/* AB4064 8003CEC4 8FA80100 */  lw    $t0, 0x100($sp)
/* AB4068 8003CEC8 46001032 */  c.eq.s $f2, $f0
/* AB406C 8003CECC C7AC00C4 */  lwc1  $f12, 0xc4($sp)
/* AB4070 8003CED0 31AE0001 */  andi  $t6, $t5, 1
/* AB4074 8003CED4 45000005 */  bc1f  .L8003CEEC
/* AB4078 8003CED8 00000000 */   nop   
/* AB407C 8003CEDC 46006032 */  c.eq.s $f12, $f0
/* AB4080 8003CEE0 00000000 */  nop   
/* AB4084 8003CEE4 450300CB */  bc1tl .L8003D214
/* AB4088 8003CEE8 8E0C0000 */   lw    $t4, ($s0)
.L8003CEEC:
/* AB408C 8003CEEC 15C000C8 */  bnez  $t6, .L8003D210
/* AB4090 8003CEF0 C7AC00C4 */   lwc1  $f12, 0xc4($sp)
/* AB4094 8003CEF4 C7A00104 */  lwc1  $f0, 0x104($sp)
/* AB4098 8003CEF8 C7B200C8 */  lwc1  $f18, 0xc8($sp)
/* AB409C 8003CEFC 3C0140A0 */  li    $at, 0x40A00000 # 0.000000
/* AB40A0 8003CF00 44813000 */  mtc1  $at, $f6
/* AB40A4 8003CF04 46120100 */  add.s $f4, $f0, $f18
/* AB40A8 8003CF08 8FA400E0 */  lw    $a0, 0xe0($sp)
/* AB40AC 8003CF0C 97A500E6 */  lhu   $a1, 0xe6($sp)
/* AB40B0 8003CF10 00003025 */  move  $a2, $zero
/* AB40B4 8003CF14 4606203C */  c.lt.s $f4, $f6
/* AB40B8 8003CF18 8FA700F0 */  lw    $a3, 0xf0($sp)
/* AB40BC 8003CF1C 8FAF00EC */  lw    $t7, 0xec($sp)
/* AB40C0 8003CF20 C7A400F4 */  lwc1  $f4, 0xf4($sp)
/* AB40C4 8003CF24 45000058 */  bc1f  .L8003D088
/* AB40C8 8003CF28 8FAB00EC */   lw    $t3, 0xec($sp)
/* AB40CC 8003CF2C 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* AB40D0 8003CF30 44814000 */  mtc1  $at, $f8
/* AB40D4 8003CF34 27B800A8 */  addiu $t8, $sp, 0xa8
/* AB40D8 8003CF38 27B900D0 */  addiu $t9, $sp, 0xd0
/* AB40DC 8003CF3C 27A900A4 */  addiu $t1, $sp, 0xa4
/* AB40E0 8003CF40 240A001B */  li    $t2, 27
/* AB40E4 8003CF44 AFAA0028 */  sw    $t2, 0x28($sp)
/* AB40E8 8003CF48 AFA9001C */  sw    $t1, 0x1c($sp)
/* AB40EC 8003CF4C AFB90018 */  sw    $t9, 0x18($sp)
/* AB40F0 8003CF50 AFB80014 */  sw    $t8, 0x14($sp)
/* AB40F4 8003CF54 AFAF0010 */  sw    $t7, 0x10($sp)
/* AB40F8 8003CF58 AFA80020 */  sw    $t0, 0x20($sp)
/* AB40FC 8003CF5C 0C00F5FC */  jal   func_8003D7F0
/* AB4100 8003CF60 E7A80024 */   swc1  $f8, 0x24($sp)
/* AB4104 8003CF64 104000AA */  beqz  $v0, .L8003D210
/* AB4108 8003CF68 AFA200D4 */   sw    $v0, 0xd4($sp)
/* AB410C 8003CF6C 8FAB00D0 */  lw    $t3, 0xd0($sp)
/* AB4110 8003CF70 3C018014 */  lui   $at, %hi(D_80138F78)
/* AB4114 8003CF74 C4248F78 */  lwc1  $f4, %lo(D_80138F78)($at)
/* AB4118 8003CF78 856C000A */  lh    $t4, 0xa($t3)
/* AB411C 8003CF7C 3C013F00 */  li    $at, 0x3F000000 # 0.000000
/* AB4120 8003CF80 44813000 */  mtc1  $at, $f6
/* AB4124 8003CF84 448C5000 */  mtc1  $t4, $f10
/* AB4128 8003CF88 C7A800A8 */  lwc1  $f8, 0xa8($sp)
/* AB412C 8003CF8C 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* AB4130 8003CF90 468054A0 */  cvt.s.w $f18, $f10
/* AB4134 8003CF94 8FAD00D0 */  lw    $t5, 0xd0($sp)
/* AB4138 8003CF98 46049302 */  mul.s $f12, $f18, $f4
/* AB413C 8003CF9C 460C303C */  c.lt.s $f6, $f12
/* AB4140 8003CFA0 00000000 */  nop   
/* AB4144 8003CFA4 45020016 */  bc1fl .L8003D000
/* AB4148 8003CFA8 85AE0008 */   lh    $t6, 8($t5)
/* AB414C 8003CFAC E6080000 */  swc1  $f8, ($s0)
/* AB4150 8003CFB0 C7B20104 */  lwc1  $f18, 0x104($sp)
/* AB4154 8003CFB4 44815000 */  mtc1  $at, $f10
/* AB4158 8003CFB8 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* AB415C 8003CFBC C7A400AC */  lwc1  $f4, 0xac($sp)
/* AB4160 8003CFC0 4612503C */  c.lt.s $f10, $f18
/* AB4164 8003CFC4 C7B20104 */  lwc1  $f18, 0x104($sp)
/* AB4168 8003CFC8 C7AA00AC */  lwc1  $f10, 0xac($sp)
/* AB416C 8003CFCC 45020007 */  bc1fl .L8003CFEC
/* AB4170 8003CFD0 46125101 */   sub.s $f4, $f10, $f18
/* AB4174 8003CFD4 44813000 */  mtc1  $at, $f6
/* AB4178 8003CFD8 00000000 */  nop   
/* AB417C 8003CFDC 46062201 */  sub.s $f8, $f4, $f6
/* AB4180 8003CFE0 10000003 */  b     .L8003CFF0
/* AB4184 8003CFE4 E6080004 */   swc1  $f8, 4($s0)
/* AB4188 8003CFE8 46125101 */  sub.s $f4, $f10, $f18
.L8003CFEC:
/* AB418C 8003CFEC E6040004 */  swc1  $f4, 4($s0)
.L8003CFF0:
/* AB4190 8003CFF0 C7A600B0 */  lwc1  $f6, 0xb0($sp)
/* AB4194 8003CFF4 1000001D */  b     .L8003D06C
/* AB4198 8003CFF8 E6060008 */   swc1  $f6, 8($s0)
/* AB419C 8003CFFC 85AE0008 */  lh    $t6, 8($t5)
.L8003D000:
/* AB41A0 8003D000 85AF000C */  lh    $t7, 0xc($t5)
/* AB41A4 8003D004 3C018014 */  lui   $at, %hi(D_80138F7C)
/* AB41A8 8003D008 448E4000 */  mtc1  $t6, $f8
/* AB41AC 8003D00C 448F2000 */  mtc1  $t7, $f4
/* AB41B0 8003D010 C4328F7C */  lwc1  $f18, %lo(D_80138F7C)($at)
/* AB41B4 8003D014 468042A0 */  cvt.s.w $f10, $f8
/* AB41B8 8003D018 3C018014 */  lui   $at, %hi(D_80138F80)
/* AB41BC 8003D01C C4288F80 */  lwc1  $f8, %lo(D_80138F80)($at)
/* AB41C0 8003D020 468021A0 */  cvt.s.w $f6, $f4
/* AB41C4 8003D024 46125002 */  mul.s $f0, $f10, $f18
/* AB41C8 8003D028 C7AA00F4 */  lwc1  $f10, 0xf4($sp)
/* AB41CC 8003D02C C7A400A8 */  lwc1  $f4, 0xa8($sp)
/* AB41D0 8003D030 46083082 */  mul.s $f2, $f6, $f8
/* AB41D4 8003D034 00000000 */  nop   
/* AB41D8 8003D038 46005482 */  mul.s $f18, $f10, $f0
/* AB41DC 8003D03C 46049180 */  add.s $f6, $f18, $f4
/* AB41E0 8003D040 E6060000 */  swc1  $f6, ($s0)
/* AB41E4 8003D044 C7A800F4 */  lwc1  $f8, 0xf4($sp)
/* AB41E8 8003D048 C7B200AC */  lwc1  $f18, 0xac($sp)
/* AB41EC 8003D04C 460C4282 */  mul.s $f10, $f8, $f12
/* AB41F0 8003D050 46125100 */  add.s $f4, $f10, $f18
/* AB41F4 8003D054 E6040004 */  swc1  $f4, 4($s0)
/* AB41F8 8003D058 C7A600F4 */  lwc1  $f6, 0xf4($sp)
/* AB41FC 8003D05C C7AA00B0 */  lwc1  $f10, 0xb0($sp)
/* AB4200 8003D060 46023202 */  mul.s $f8, $f6, $f2
/* AB4204 8003D064 460A4480 */  add.s $f18, $f8, $f10
/* AB4208 8003D068 E6120008 */  swc1  $f18, 8($s0)
.L8003D06C:
/* AB420C 8003D06C 8FB800D0 */  lw    $t8, 0xd0($sp)
/* AB4210 8003D070 8FB900F8 */  lw    $t9, 0xf8($sp)
/* AB4214 8003D074 AF380000 */  sw    $t8, ($t9)
/* AB4218 8003D078 8FAA00FC */  lw    $t2, 0xfc($sp)
/* AB421C 8003D07C 8FA900A4 */  lw    $t1, 0xa4($sp)
/* AB4220 8003D080 10000063 */  b     .L8003D210
/* AB4224 8003D084 AD490000 */   sw    $t1, ($t2)
.L8003D088:
/* AB4228 8003D088 46042182 */  mul.s $f6, $f4, $f4
/* AB422C 8003D08C 24020019 */  li    $v0, 25
/* AB4230 8003D090 46021202 */  mul.s $f8, $f2, $f2
/* AB4234 8003D094 00000000 */  nop   
/* AB4238 8003D098 460C6282 */  mul.s $f10, $f12, $f12
/* AB423C 8003D09C 460A4480 */  add.s $f18, $f8, $f10
/* AB4240 8003D0A0 4612303C */  c.lt.s $f6, $f18
/* AB4244 8003D0A4 00000000 */  nop   
/* AB4248 8003D0A8 45000003 */  bc1f  .L8003D0B8
/* AB424C 8003D0AC 00000000 */   nop   
/* AB4250 8003D0B0 10000001 */  b     .L8003D0B8
/* AB4254 8003D0B4 2402001B */   li    $v0, 27
.L8003D0B8:
/* AB4258 8003D0B8 8D6E0000 */  lw    $t6, ($t3)
/* AB425C 8003D0BC 27A30088 */  addiu $v1, $sp, 0x88
/* AB4260 8003D0C0 27A7007C */  addiu $a3, $sp, 0x7c
/* AB4264 8003D0C4 AC6E0000 */  sw    $t6, ($v1)
/* AB4268 8003D0C8 8D6C0004 */  lw    $t4, 4($t3)
/* AB426C 8003D0CC 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* AB4270 8003D0D0 44813000 */  mtc1  $at, $f6
/* AB4274 8003D0D4 AC6C0004 */  sw    $t4, 4($v1)
/* AB4278 8003D0D8 8D6E0008 */  lw    $t6, 8($t3)
/* AB427C 8003D0DC 27B900A8 */  addiu $t9, $sp, 0xa8
/* AB4280 8003D0E0 27A900D0 */  addiu $t1, $sp, 0xd0
/* AB4284 8003D0E4 AC6E0008 */  sw    $t6, 8($v1)
/* AB4288 8003D0E8 C7A4008C */  lwc1  $f4, 0x8c($sp)
/* AB428C 8003D0EC 8FAD00F0 */  lw    $t5, 0xf0($sp)
/* AB4290 8003D0F0 27AA00A4 */  addiu $t2, $sp, 0xa4
/* AB4294 8003D0F4 46002200 */  add.s $f8, $f4, $f0
/* AB4298 8003D0F8 00003025 */  move  $a2, $zero
/* AB429C 8003D0FC E7A8008C */  swc1  $f8, 0x8c($sp)
/* AB42A0 8003D100 8DB80000 */  lw    $t8, ($t5)
/* AB42A4 8003D104 ACF80000 */  sw    $t8, ($a3)
/* AB42A8 8003D108 8DAF0004 */  lw    $t7, 4($t5)
/* AB42AC 8003D10C ACEF0004 */  sw    $t7, 4($a3)
/* AB42B0 8003D110 8DB80008 */  lw    $t8, 8($t5)
/* AB42B4 8003D114 ACF80008 */  sw    $t8, 8($a3)
/* AB42B8 8003D118 C7AA008C */  lwc1  $f10, 0x8c($sp)
/* AB42BC 8003D11C AFA20028 */  sw    $v0, 0x28($sp)
/* AB42C0 8003D120 AFA80020 */  sw    $t0, 0x20($sp)
/* AB42C4 8003D124 AFAA001C */  sw    $t2, 0x1c($sp)
/* AB42C8 8003D128 AFA90018 */  sw    $t1, 0x18($sp)
/* AB42CC 8003D12C AFB90014 */  sw    $t9, 0x14($sp)
/* AB42D0 8003D130 AFA30010 */  sw    $v1, 0x10($sp)
/* AB42D4 8003D134 97A500E6 */  lhu   $a1, 0xe6($sp)
/* AB42D8 8003D138 8FA400E0 */  lw    $a0, 0xe0($sp)
/* AB42DC 8003D13C E7A60024 */  swc1  $f6, 0x24($sp)
/* AB42E0 8003D140 0C00F5FC */  jal   func_8003D7F0
/* AB42E4 8003D144 E7AA0080 */   swc1  $f10, 0x80($sp)
/* AB42E8 8003D148 10400031 */  beqz  $v0, .L8003D210
/* AB42EC 8003D14C AFA200D4 */   sw    $v0, 0xd4($sp)
/* AB42F0 8003D150 8FAB00D0 */  lw    $t3, 0xd0($sp)
/* AB42F4 8003D154 3C018014 */  lui   $at, %hi(D_80138F84)
/* AB42F8 8003D158 C4288F84 */  lwc1  $f8, %lo(D_80138F84)($at)
/* AB42FC 8003D15C 856C0008 */  lh    $t4, 8($t3)
/* AB4300 8003D160 856E000C */  lh    $t6, 0xc($t3)
/* AB4304 8003D164 3C018014 */  lui   $at, %hi(D_80138F88)
/* AB4308 8003D168 448C9000 */  mtc1  $t4, $f18
/* AB430C 8003D16C 448E5000 */  mtc1  $t6, $f10
/* AB4310 8003D170 46809120 */  cvt.s.w $f4, $f18
/* AB4314 8003D174 C4328F88 */  lwc1  $f18, %lo(D_80138F88)($at)
/* AB4318 8003D178 3C018014 */  lui   $at, %hi(D_80138F8C)
/* AB431C 8003D17C 468051A0 */  cvt.s.w $f6, $f10
/* AB4320 8003D180 46082082 */  mul.s $f2, $f4, $f8
/* AB4324 8003D184 C42A8F8C */  lwc1  $f10, %lo(D_80138F8C)($at)
/* AB4328 8003D188 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* AB432C 8003D18C 46123302 */  mul.s $f12, $f6, $f18
/* AB4330 8003D190 00000000 */  nop   
/* AB4334 8003D194 46021102 */  mul.s $f4, $f2, $f2
/* AB4338 8003D198 00000000 */  nop   
/* AB433C 8003D19C 460C6202 */  mul.s $f8, $f12, $f12
/* AB4340 8003D1A0 46082000 */  add.s $f0, $f4, $f8
/* AB4344 8003D1A4 46000384 */  sqrt.s $f14, $f0
/* AB4348 8003D1A8 46007005 */  abs.s $f0, $f14
/* AB434C 8003D1AC 460A003C */  c.lt.s $f0, $f10
/* AB4350 8003D1B0 00000000 */  nop   
/* AB4354 8003D1B4 45010016 */  bc1t  .L8003D210
/* AB4358 8003D1B8 00000000 */   nop   
/* AB435C 8003D1BC 44819000 */  mtc1  $at, $f18
/* AB4360 8003D1C0 C7A600F4 */  lwc1  $f6, 0xf4($sp)
/* AB4364 8003D1C4 C7AA00A8 */  lwc1  $f10, 0xa8($sp)
/* AB4368 8003D1C8 460E9103 */  div.s $f4, $f18, $f14
/* AB436C 8003D1CC 24090001 */  li    $t1, 1
/* AB4370 8003D1D0 46043002 */  mul.s $f0, $f6, $f4
/* AB4374 8003D1D4 00000000 */  nop   
/* AB4378 8003D1D8 46020202 */  mul.s $f8, $f0, $f2
/* AB437C 8003D1DC 460A4480 */  add.s $f18, $f8, $f10
/* AB4380 8003D1E0 460C0182 */  mul.s $f6, $f0, $f12
/* AB4384 8003D1E4 E6120000 */  swc1  $f18, ($s0)
/* AB4388 8003D1E8 C7A400B0 */  lwc1  $f4, 0xb0($sp)
/* AB438C 8003D1EC 46043200 */  add.s $f8, $f6, $f4
/* AB4390 8003D1F0 E6080008 */  swc1  $f8, 8($s0)
/* AB4394 8003D1F4 8FAF00F8 */  lw    $t7, 0xf8($sp)
/* AB4398 8003D1F8 8FAD00D0 */  lw    $t5, 0xd0($sp)
/* AB439C 8003D1FC ADED0000 */  sw    $t5, ($t7)
/* AB43A0 8003D200 8FB900FC */  lw    $t9, 0xfc($sp)
/* AB43A4 8003D204 8FB800A4 */  lw    $t8, 0xa4($sp)
/* AB43A8 8003D208 AF380000 */  sw    $t8, ($t9)
/* AB43AC 8003D20C AFA900D4 */  sw    $t1, 0xd4($sp)
.L8003D210:
/* AB43B0 8003D210 8E0C0000 */  lw    $t4, ($s0)
.L8003D214:
/* AB43B4 8003D214 8FA80100 */  lw    $t0, 0x100($sp)
/* AB43B8 8003D218 27A200B8 */  addiu $v0, $sp, 0xb8
/* AB43BC 8003D21C AC4C0000 */  sw    $t4, ($v0)
/* AB43C0 8003D220 8E0A0004 */  lw    $t2, 4($s0)
/* AB43C4 8003D224 26070008 */  addiu $a3, $s0, 8
/* AB43C8 8003D228 02003025 */  move  $a2, $s0
/* AB43CC 8003D22C AC4A0004 */  sw    $t2, 4($v0)
/* AB43D0 8003D230 8E0C0008 */  lw    $t4, 8($s0)
/* AB43D4 8003D234 AC4C0008 */  sw    $t4, 8($v0)
/* AB43D8 8003D238 C7AA00BC */  lwc1  $f10, 0xbc($sp)
/* AB43DC 8003D23C C7B20104 */  lwc1  $f18, 0x104($sp)
/* AB43E0 8003D240 8FAE00FC */  lw    $t6, 0xfc($sp)
/* AB43E4 8003D244 8FAB00F8 */  lw    $t3, 0xf8($sp)
/* AB43E8 8003D248 46125180 */  add.s $f6, $f10, $f18
/* AB43EC 8003D24C C7A400F4 */  lwc1  $f4, 0xf4($sp)
/* AB43F0 8003D250 AFA000B4 */  sw    $zero, 0xb4($sp)
/* AB43F4 8003D254 AFA70044 */  sw    $a3, 0x44($sp)
/* AB43F8 8003D258 E7A600BC */  swc1  $f6, 0xbc($sp)
/* AB43FC 8003D25C AFA20010 */  sw    $v0, 0x10($sp)
/* AB4400 8003D260 97A500E6 */  lhu   $a1, 0xe6($sp)
/* AB4404 8003D264 8FA400E0 */  lw    $a0, 0xe0($sp)
/* AB4408 8003D268 AFA80020 */  sw    $t0, 0x20($sp)
/* AB440C 8003D26C AFAE001C */  sw    $t6, 0x1c($sp)
/* AB4410 8003D270 AFAB0018 */  sw    $t3, 0x18($sp)
/* AB4414 8003D274 0C01026A */  jal   func_800409A8
/* AB4418 8003D278 E7A40014 */   swc1  $f4, 0x14($sp)
/* AB441C 8003D27C 1040000F */  beqz  $v0, .L8003D2BC
/* AB4420 8003D280 240D0001 */   li    $t5, 1
/* AB4424 8003D284 240F0001 */  li    $t7, 1
/* AB4428 8003D288 AFAD00D4 */  sw    $t5, 0xd4($sp)
/* AB442C 8003D28C AFAF00B4 */  sw    $t7, 0xb4($sp)
/* AB4430 8003D290 8E090000 */  lw    $t1, ($s0)
/* AB4434 8003D294 27B800B8 */  addiu $t8, $sp, 0xb8
/* AB4438 8003D298 AF090000 */  sw    $t1, ($t8)
/* AB443C 8003D29C 8E190004 */  lw    $t9, 4($s0)
/* AB4440 8003D2A0 AF190004 */  sw    $t9, 4($t8)
/* AB4444 8003D2A4 8E090008 */  lw    $t1, 8($s0)
/* AB4448 8003D2A8 AF090008 */  sw    $t1, 8($t8)
/* AB444C 8003D2AC C7A800BC */  lwc1  $f8, 0xbc($sp)
/* AB4450 8003D2B0 C7AA0104 */  lwc1  $f10, 0x104($sp)
/* AB4454 8003D2B4 460A4480 */  add.s $f18, $f8, $f10
/* AB4458 8003D2B8 E7B200BC */  swc1  $f18, 0xbc($sp)
.L8003D2BC:
/* AB445C 8003D2BC 8FA400E0 */  lw    $a0, 0xe0($sp)
/* AB4460 8003D2C0 0C00F157 */  jal   func_8003C55C
/* AB4464 8003D2C4 8FA500EC */   lw    $a1, 0xec($sp)
/* AB4468 8003D2C8 24010001 */  li    $at, 1
/* AB446C 8003D2CC 14410017 */  bne   $v0, $at, .L8003D32C
/* AB4470 8003D2D0 8FA400E0 */   lw    $a0, 0xe0($sp)
/* AB4474 8003D2D4 8FA500DC */  lw    $a1, 0xdc($sp)
/* AB4478 8003D2D8 0C00EB15 */  jal   func_8003AC54
/* AB447C 8003D2DC 02003025 */   move  $a2, $s0
/* AB4480 8003D2E0 8FAA0044 */  lw    $t2, 0x44($sp)
/* AB4484 8003D2E4 C7A600F4 */  lwc1  $f6, 0xf4($sp)
/* AB4488 8003D2E8 8FAB00F8 */  lw    $t3, 0xf8($sp)
/* AB448C 8003D2EC 27AC00B8 */  addiu $t4, $sp, 0xb8
/* AB4490 8003D2F0 AFAC0014 */  sw    $t4, 0x14($sp)
/* AB4494 8003D2F4 00402025 */  move  $a0, $v0
/* AB4498 8003D2F8 8FA500E0 */  lw    $a1, 0xe0($sp)
/* AB449C 8003D2FC 97A600E6 */  lhu   $a2, 0xe6($sp)
/* AB44A0 8003D300 02003825 */  move  $a3, $s0
/* AB44A4 8003D304 AFAA0010 */  sw    $t2, 0x10($sp)
/* AB44A8 8003D308 E7A60018 */  swc1  $f6, 0x18($sp)
/* AB44AC 8003D30C 0C00E6BB */  jal   func_80039AEC
/* AB44B0 8003D310 AFAB001C */   sw    $t3, 0x1c($sp)
/* AB44B4 8003D314 10400005 */  beqz  $v0, .L8003D32C
/* AB44B8 8003D318 8FAD00FC */   lw    $t5, 0xfc($sp)
/* AB44BC 8003D31C 240E0032 */  li    $t6, 50
/* AB44C0 8003D320 ADAE0000 */  sw    $t6, ($t5)
/* AB44C4 8003D324 240F0001 */  li    $t7, 1
/* AB44C8 8003D328 AFAF00D4 */  sw    $t7, 0xd4($sp)
.L8003D32C:
/* AB44CC 8003D32C 8FB800B4 */  lw    $t8, 0xb4($sp)
/* AB44D0 8003D330 24010001 */  li    $at, 1
/* AB44D4 8003D334 8FB900FC */  lw    $t9, 0xfc($sp)
/* AB44D8 8003D338 13010004 */  beq   $t8, $at, .L8003D34C
/* AB44DC 8003D33C 8FA400E0 */   lw    $a0, 0xe0($sp)
/* AB44E0 8003D340 8F290000 */  lw    $t1, ($t9)
/* AB44E4 8003D344 24010032 */  li    $at, 50
/* AB44E8 8003D348 11210041 */  beq   $t1, $at, .L8003D450
.L8003D34C:
/* AB44EC 8003D34C 3C013F80 */   li    $at, 0x3F800000 # 0.000000
/* AB44F0 8003D350 44812000 */  mtc1  $at, $f4
/* AB44F4 8003D354 8FAE0100 */  lw    $t6, 0x100($sp)
/* AB44F8 8003D358 27AA005C */  addiu $t2, $sp, 0x5c
/* AB44FC 8003D35C 27AC00D0 */  addiu $t4, $sp, 0xd0
/* AB4500 8003D360 27AB0058 */  addiu $t3, $sp, 0x58
/* AB4504 8003D364 240D0009 */  li    $t5, 9
/* AB4508 8003D368 AFAD0028 */  sw    $t5, 0x28($sp)
/* AB450C 8003D36C AFAB001C */  sw    $t3, 0x1c($sp)
/* AB4510 8003D370 AFAC0018 */  sw    $t4, 0x18($sp)
/* AB4514 8003D374 AFAA0014 */  sw    $t2, 0x14($sp)
/* AB4518 8003D378 97A500E6 */  lhu   $a1, 0xe6($sp)
/* AB451C 8003D37C 00003025 */  move  $a2, $zero
/* AB4520 8003D380 8FA700F0 */  lw    $a3, 0xf0($sp)
/* AB4524 8003D384 AFB00010 */  sw    $s0, 0x10($sp)
/* AB4528 8003D388 AFAE0020 */  sw    $t6, 0x20($sp)
/* AB452C 8003D38C 0C00F5FC */  jal   func_8003D7F0
/* AB4530 8003D390 E7A40024 */   swc1  $f4, 0x24($sp)
/* AB4534 8003D394 1040002E */  beqz  $v0, .L8003D450
/* AB4538 8003D398 8FAF00D0 */   lw    $t7, 0xd0($sp)
/* AB453C 8003D39C 85F80008 */  lh    $t8, 8($t7)
/* AB4540 8003D3A0 85F9000C */  lh    $t9, 0xc($t7)
/* AB4544 8003D3A4 3C018014 */  lui   $at, %hi(D_80138F90)
/* AB4548 8003D3A8 44984000 */  mtc1  $t8, $f8
/* AB454C 8003D3AC 44999000 */  mtc1  $t9, $f18
/* AB4550 8003D3B0 C4308F90 */  lwc1  $f16, %lo(D_80138F90)($at)
/* AB4554 8003D3B4 468042A0 */  cvt.s.w $f10, $f8
/* AB4558 8003D3B8 3C018014 */  lui   $at, %hi(D_80138F94)
/* AB455C 8003D3BC 468091A0 */  cvt.s.w $f6, $f18
/* AB4560 8003D3C0 46105082 */  mul.s $f2, $f10, $f16
/* AB4564 8003D3C4 C42A8F94 */  lwc1  $f10, %lo(D_80138F94)($at)
/* AB4568 8003D3C8 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* AB456C 8003D3CC 46103302 */  mul.s $f12, $f6, $f16
/* AB4570 8003D3D0 00000000 */  nop   
/* AB4574 8003D3D4 46021102 */  mul.s $f4, $f2, $f2
/* AB4578 8003D3D8 00000000 */  nop   
/* AB457C 8003D3DC 460C6202 */  mul.s $f8, $f12, $f12
/* AB4580 8003D3E0 46082000 */  add.s $f0, $f4, $f8
/* AB4584 8003D3E4 46000384 */  sqrt.s $f14, $f0
/* AB4588 8003D3E8 46007005 */  abs.s $f0, $f14
/* AB458C 8003D3EC 460A003C */  c.lt.s $f0, $f10
/* AB4590 8003D3F0 00000000 */  nop   
/* AB4594 8003D3F4 45010016 */  bc1t  .L8003D450
/* AB4598 8003D3F8 00000000 */   nop   
/* AB459C 8003D3FC 44813000 */  mtc1  $at, $f6
/* AB45A0 8003D400 C7B200F4 */  lwc1  $f18, 0xf4($sp)
/* AB45A4 8003D404 C7AA005C */  lwc1  $f10, 0x5c($sp)
/* AB45A8 8003D408 460E3103 */  div.s $f4, $f6, $f14
/* AB45AC 8003D40C 240E0001 */  li    $t6, 1
/* AB45B0 8003D410 46049002 */  mul.s $f0, $f18, $f4
/* AB45B4 8003D414 00000000 */  nop   
/* AB45B8 8003D418 46020202 */  mul.s $f8, $f0, $f2
/* AB45BC 8003D41C 460A4180 */  add.s $f6, $f8, $f10
/* AB45C0 8003D420 460C0482 */  mul.s $f18, $f0, $f12
/* AB45C4 8003D424 E6060000 */  swc1  $f6, ($s0)
/* AB45C8 8003D428 C7A40064 */  lwc1  $f4, 0x64($sp)
/* AB45CC 8003D42C 46049200 */  add.s $f8, $f18, $f4
/* AB45D0 8003D430 E6080008 */  swc1  $f8, 8($s0)
/* AB45D4 8003D434 8FAA00F8 */  lw    $t2, 0xf8($sp)
/* AB45D8 8003D438 8FA900D0 */  lw    $t1, 0xd0($sp)
/* AB45DC 8003D43C AD490000 */  sw    $t1, ($t2)
/* AB45E0 8003D440 8FAB00FC */  lw    $t3, 0xfc($sp)
/* AB45E4 8003D444 8FAC0058 */  lw    $t4, 0x58($sp)
/* AB45E8 8003D448 AD6C0000 */  sw    $t4, ($t3)
/* AB45EC 8003D44C AFAE00D4 */  sw    $t6, 0xd4($sp)
.L8003D450:
/* AB45F0 8003D450 8FBF003C */  lw    $ra, 0x3c($sp)
/* AB45F4 8003D454 8FA200D4 */  lw    $v0, 0xd4($sp)
/* AB45F8 8003D458 8FB00038 */  lw    $s0, 0x38($sp)
/* AB45FC 8003D45C 03E00008 */  jr    $ra
/* AB4600 8003D460 27BD00E0 */   addiu $sp, $sp, 0xe0

