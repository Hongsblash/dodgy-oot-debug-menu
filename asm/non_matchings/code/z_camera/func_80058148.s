.rdata
glabel D_801398C0
    .asciz "\x1B[34mcamera: personalize ---\x1B[m\n"
    .balign 4

.text
glabel func_80058148
/* ACF2E8 80058148 27BDFF98 */  addiu $sp, $sp, -0x68
/* ACF2EC 8005814C AFB00018 */  sw    $s0, 0x18($sp)
/* ACF2F0 80058150 00808025 */  move  $s0, $a0
/* ACF2F4 80058154 AFBF001C */  sw    $ra, 0x1c($sp)
/* ACF2F8 80058158 AFA5006C */  sw    $a1, 0x6c($sp)
/* ACF2FC 8005815C 0C00BBD1 */  jal   func_8002EF44
/* ACF300 80058160 27A40054 */   addiu $a0, $sp, 0x54
/* ACF304 80058164 0C00B721 */  jal   Player_GetCameraYOffset
/* ACF308 80058168 8FA4006C */   lw    $a0, 0x6c($sp)
/* ACF30C 8005816C 8FAE006C */  lw    $t6, 0x6c($sp)
/* ACF310 80058170 27A20054 */  addiu $v0, $sp, 0x54
/* ACF314 80058174 3C014334 */  li    $at, 0x43340000 # 0.000000
/* ACF318 80058178 AE0E0090 */  sw    $t6, 0x90($s0)
/* ACF31C 8005817C 8C580000 */  lw    $t8, ($v0)
/* ACF320 80058180 44816000 */  mtc1  $at, $f12
/* ACF324 80058184 240A071C */  li    $t2, 1820
/* ACF328 80058188 AE180094 */  sw    $t8, 0x94($s0)
/* ACF32C 8005818C 8C4F0004 */  lw    $t7, 4($v0)
/* ACF330 80058190 240B071C */  li    $t3, 1820
/* ACF334 80058194 44801000 */  mtc1  $zero, $f2
/* ACF338 80058198 AE0F0098 */  sw    $t7, 0x98($s0)
/* ACF33C 8005819C 8C580008 */  lw    $t8, 8($v0)
/* ACF340 800581A0 26050050 */  addiu $a1, $s0, 0x50
/* ACF344 800581A4 44807000 */  mtc1  $zero, $f14
/* ACF348 800581A8 AE18009C */  sw    $t8, 0x9c($s0)
/* ACF34C 800581AC 8C4F000C */  lw    $t7, 0xc($v0)
/* ACF350 800581B0 26040074 */  addiu $a0, $s0, 0x74
/* ACF354 800581B4 27A6004C */  addiu $a2, $sp, 0x4c
/* ACF358 800581B8 AE0F00A0 */  sw    $t7, 0xa0($s0)
/* ACF35C 800581BC 8C580010 */  lw    $t8, 0x10($v0)
/* ACF360 800581C0 AE1800A4 */  sw    $t8, 0xa4($s0)
/* ACF364 800581C4 E7AC004C */  swc1  $f12, 0x4c($sp)
/* ACF368 800581C8 E60C00DC */  swc1  $f12, 0xdc($s0)
/* ACF36C 800581CC 87B90062 */  lh    $t9, 0x62($sp)
/* ACF370 800581D0 A6190136 */  sh    $t9, 0x136($s0)
/* ACF374 800581D4 86080136 */  lh    $t0, 0x136($s0)
/* ACF378 800581D8 A7AA0050 */  sh    $t2, 0x50($sp)
/* ACF37C 800581DC 25098001 */  addiu $t1, $t0, -0x7fff
/* ACF380 800581E0 A7A90052 */  sh    $t1, 0x52($sp)
/* ACF384 800581E4 A60B0134 */  sh    $t3, 0x134($s0)
/* ACF388 800581E8 8A0D0134 */  lwl   $t5, 0x134($s0)
/* ACF38C 800581EC 9A0D0137 */  lwr   $t5, 0x137($s0)
/* ACF390 800581F0 A6000138 */  sh    $zero, 0x138($s0)
/* ACF394 800581F4 A600013E */  sh    $zero, 0x13e($s0)
/* ACF398 800581F8 AA0D013A */  swl   $t5, 0x13a($s0)
/* ACF39C 800581FC E60200D8 */  swc1  $f2, 0xd8($s0)
/* ACF3A0 80058200 E60200F4 */  swc1  $f2, 0xf4($s0)
/* ACF3A4 80058204 BA0D013D */  swr   $t5, 0x13d($s0)
/* ACF3A8 80058208 8C4F0000 */  lw    $t7, ($v0)
/* ACF3AC 8005820C ACAF0000 */  sw    $t7, ($a1)
/* ACF3B0 80058210 8C4E0004 */  lw    $t6, 4($v0)
/* ACF3B4 80058214 ACAE0004 */  sw    $t6, 4($a1)
/* ACF3B8 80058218 8C4F0008 */  lw    $t7, 8($v0)
/* ACF3BC 8005821C ACAF0008 */  sw    $t7, 8($a1)
/* ACF3C0 80058220 C6040054 */  lwc1  $f4, 0x54($s0)
/* ACF3C4 80058224 E60000E8 */  swc1  $f0, 0xe8($s0)
/* ACF3C8 80058228 E60E00E4 */  swc1  $f14, 0xe4($s0)
/* ACF3CC 8005822C 46002180 */  add.s $f6, $f4, $f0
/* ACF3D0 80058230 E60E00EC */  swc1  $f14, 0xec($s0)
/* ACF3D4 80058234 E6060054 */  swc1  $f6, 0x54($s0)
/* ACF3D8 80058238 AFA50024 */  sw    $a1, 0x24($sp)
/* ACF3DC 8005823C 0C010F0A */  jal   func_80043C28
/* ACF3E0 80058240 AFA40020 */   sw    $a0, 0x20($sp)
/* ACF3E4 80058244 8FB80020 */  lw    $t8, 0x20($sp)
/* ACF3E8 80058248 44801000 */  mtc1  $zero, $f2
/* ACF3EC 8005824C 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* ACF3F0 80058250 8F080000 */  lw    $t0, ($t8)
/* ACF3F4 80058254 44814000 */  mtc1  $at, $f8
/* ACF3F8 80058258 02002025 */  move  $a0, $s0
/* ACF3FC 8005825C AE08005C */  sw    $t0, 0x5c($s0)
/* ACF400 80058260 8F190004 */  lw    $t9, 4($t8)
/* ACF404 80058264 27A5003C */  addiu $a1, $sp, 0x3c
/* ACF408 80058268 27A70048 */  addiu $a3, $sp, 0x48
/* ACF40C 8005826C AE190060 */  sw    $t9, 0x60($s0)
/* ACF410 80058270 8F080008 */  lw    $t0, 8($t8)
/* ACF414 80058274 A600015A */  sh    $zero, 0x15a($s0)
/* ACF418 80058278 E6020070 */  swc1  $f2, 0x70($s0)
/* ACF41C 8005827C E6020068 */  swc1  $f2, 0x68($s0)
/* ACF420 80058280 E608006C */  swc1  $f8, 0x6c($s0)
/* ACF424 80058284 AE080064 */  sw    $t0, 0x64($s0)
/* ACF428 80058288 0C01110D */  jal   func_80044434
/* ACF42C 8005828C 8FA60024 */   lw    $a2, 0x24($sp)
/* ACF430 80058290 3C01C6FA */  li    $at, 0xC6FA0000 # 0.000000
/* ACF434 80058294 44815000 */  mtc1  $at, $f10
/* ACF438 80058298 2402FFFF */  li    $v0, -1
/* ACF43C 8005829C 8FA90048 */  lw    $t1, 0x48($sp)
/* ACF440 800582A0 460A0032 */  c.eq.s $f0, $f10
/* ACF444 800582A4 340EB200 */  li    $t6, 45568
/* ACF448 800582A8 3C018012 */  lui   $at, %hi(D_8011D3A0)
/* ACF44C 800582AC 45030003 */  bc1tl .L800582BC
/* ACF450 800582B0 860A014C */   lh    $t2, 0x14c($s0)
/* ACF454 800582B4 A6090146 */  sh    $t1, 0x146($s0)
/* ACF458 800582B8 860A014C */  lh    $t2, 0x14c($s0)
.L800582BC:
/* ACF45C 800582BC 8E0C008C */  lw    $t4, 0x8c($s0)
/* ACF460 800582C0 AE020118 */  sw    $v0, 0x118($s0)
/* ACF464 800582C4 354B0004 */  ori   $t3, $t2, 4
/* ACF468 800582C8 258D01E0 */  addiu $t5, $t4, 0x1e0
/* ACF46C 800582CC AE02011C */  sw    $v0, 0x11c($s0)
/* ACF470 800582D0 160D0004 */  bne   $s0, $t5, .L800582E4
/* ACF474 800582D4 A60B014C */   sh    $t3, 0x14c($s0)
/* ACF478 800582D8 3C018012 */  lui   $at, %hi(D_8011D3A0) # $at, 0x8012
/* ACF47C 800582DC 10000002 */  b     .L800582E8
/* ACF480 800582E0 AC2ED3A0 */   sw    $t6, %lo(D_8011D3A0)($at)
.L800582E4:
/* ACF484 800582E4 AC20D3A0 */  sw    $zero, %lo(D_8011D3A0)($at)
.L800582E8:
/* ACF488 800582E8 0C015FF1 */  jal   func_80057FC4
/* ACF48C 800582EC 02002025 */   move  $a0, $s0
/* ACF490 800582F0 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* ACF494 800582F4 44818000 */  mtc1  $at, $f16
/* ACF498 800582F8 2402FFFF */  li    $v0, -1
/* ACF49C 800582FC A600014A */  sh    $zero, 0x14a($s0)
/* ACF4A0 80058300 A600015C */  sh    $zero, 0x15c($s0)
/* ACF4A4 80058304 A6020156 */  sh    $v0, 0x156($s0)
/* ACF4A8 80058308 02002025 */  move  $a0, $s0
/* ACF4AC 8005830C 86050144 */  lh    $a1, 0x144($s0)
/* ACF4B0 80058310 0C01144A */  jal   func_80045128
/* ACF4B4 80058314 E6100100 */   swc1  $f16, 0x100($s0)
/* ACF4B8 80058318 0C016C11 */  jal   func_8005B044
/* ACF4BC 8005831C 00000000 */   nop   
/* ACF4C0 80058320 3C048014 */  lui   $a0, %hi(D_801398C0) # $a0, 0x8014
/* ACF4C4 80058324 0C00084C */  jal   osSyncPrintf
/* ACF4C8 80058328 248498C0 */   addiu $a0, %lo(D_801398C0) # addiu $a0, $a0, -0x6740
/* ACF4CC 8005832C 860F0164 */  lh    $t7, 0x164($s0)
/* ACF4D0 80058330 55E00004 */  bnezl $t7, .L80058344
/* ACF4D4 80058334 8FBF001C */   lw    $ra, 0x1c($sp)
/* ACF4D8 80058338 0C01622D */  jal   func_800588B4
/* ACF4DC 8005833C 02002025 */   move  $a0, $s0
/* ACF4E0 80058340 8FBF001C */  lw    $ra, 0x1c($sp)
.L80058344:
/* ACF4E4 80058344 8FB00018 */  lw    $s0, 0x18($sp)
/* ACF4E8 80058348 27BD0068 */  addiu $sp, $sp, 0x68
/* ACF4EC 8005834C 03E00008 */  jr    $ra
/* ACF4F0 80058350 00000000 */   nop   

