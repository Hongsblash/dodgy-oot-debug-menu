.late_rodata
glabel D_80140668
    .float 0.33333334

.text
glabel func_800A4A20
/* B1BBC0 800A4A20 3C0E8016 */  lui   $t6, %hi(gGameInfo)
/* B1BBC4 800A4A24 8DCEFA90 */  lw    $t6, %lo(gGameInfo)($t6)
/* B1BBC8 800A4A28 27BDFFD0 */  addiu $sp, $sp, -0x30
/* B1BBCC 800A4A2C AFBF0024 */  sw    $ra, 0x24($sp)
/* B1BBD0 800A4A30 AFB00020 */  sw    $s0, 0x20($sp)
/* B1BBD4 800A4A34 85CF0110 */  lh    $t7, 0x110($t6)
/* B1BBD8 800A4A38 3C018014 */  lui   $at, %hi(D_80140668)
/* B1BBDC 800A4A3C C4280668 */  lwc1  $f8, %lo(D_80140668)($at)
/* B1BBE0 800A4A40 448F2000 */  mtc1  $t7, $f4
/* B1BBE4 800A4A44 C48A002C */  lwc1  $f10, 0x2c($a0)
/* B1BBE8 800A4A48 C4800028 */  lwc1  $f0, 0x28($a0)
/* B1BBEC 800A4A4C 468021A0 */  cvt.s.w $f6, $f4
/* B1BBF0 800A4A50 44802000 */  mtc1  $zero, $f4
/* B1BBF4 800A4A54 46000306 */  mov.s $f12, $f0
/* B1BBF8 800A4A58 00808025 */  move  $s0, $a0
/* B1BBFC 800A4A5C 46083082 */  mul.s $f2, $f6, $f8
/* B1BC00 800A4A60 00000000 */  nop   
/* B1BC04 800A4A64 46025402 */  mul.s $f16, $f10, $f2
/* B1BC08 800A4A68 46100481 */  sub.s $f18, $f0, $f16
/* B1BC0C 800A4A6C E4920028 */  swc1  $f18, 0x28($a0)
/* B1BC10 800A4A70 C4800028 */  lwc1  $f0, 0x28($a0)
/* B1BC14 800A4A74 4604003E */  c.le.s $f0, $f4
/* B1BC18 800A4A78 00000000 */  nop   
/* B1BC1C 800A4A7C 45000007 */  bc1f  .L800A4A9C
/* B1BC20 800A4A80 00000000 */   nop   
/* B1BC24 800A4A84 0C02926C */  jal   func_800A49B0
/* B1BC28 800A4A88 E7AC002C */   swc1  $f12, 0x2c($sp)
/* B1BC2C 800A4A8C 44803000 */  mtc1  $zero, $f6
/* B1BC30 800A4A90 C7AC002C */  lwc1  $f12, 0x2c($sp)
/* B1BC34 800A4A94 E6060028 */  swc1  $f6, 0x28($s0)
/* B1BC38 800A4A98 C6000028 */  lwc1  $f0, 0x28($s0)
.L800A4A9C:
/* B1BC3C 800A4A9C 460C0283 */  div.s $f10, $f0, $f12
/* B1BC40 800A4AA0 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* B1BC44 800A4AA4 44814000 */  mtc1  $at, $f8
/* B1BC48 800A4AA8 8E050020 */  lw    $a1, 0x20($s0)
/* B1BC4C 800A4AAC 92040000 */  lbu   $a0, ($s0)
/* B1BC50 800A4AB0 8E070024 */  lw    $a3, 0x24($s0)
/* B1BC54 800A4AB4 00A03025 */  move  $a2, $a1
/* B1BC58 800A4AB8 460A4401 */  sub.s $f16, $f8, $f10
/* B1BC5C 800A4ABC 0C028B9C */  jal   func_800A2E70
/* B1BC60 800A4AC0 E7B00010 */   swc1  $f16, 0x10($sp)
/* B1BC64 800A4AC4 8FBF0024 */  lw    $ra, 0x24($sp)
/* B1BC68 800A4AC8 8FB00020 */  lw    $s0, 0x20($sp)
/* B1BC6C 800A4ACC 27BD0030 */  addiu $sp, $sp, 0x30
/* B1BC70 800A4AD0 03E00008 */  jr    $ra
/* B1BC74 800A4AD4 00001025 */   move  $v0, $zero

