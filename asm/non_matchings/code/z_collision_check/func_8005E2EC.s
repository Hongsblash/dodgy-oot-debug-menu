glabel func_8005E2EC
/* AD548C 8005E2EC 27BDFFE0 */  addiu $sp, $sp, -0x20
/* AD5490 8005E2F0 AFBF001C */  sw    $ra, 0x1c($sp)
/* AD5494 8005E2F4 AFA60028 */  sw    $a2, 0x28($sp)
/* AD5498 8005E2F8 90A20015 */  lbu   $v0, 0x15($a1)
/* AD549C 8005E2FC 30420018 */  andi  $v0, $v0, 0x18
/* AD54A0 8005E300 14400024 */  bnez  $v0, .L8005E394
/* AD54A4 8005E304 00000000 */   nop   
/* AD54A8 8005E308 90CF0014 */  lbu   $t7, 0x14($a2)
/* AD54AC 8005E30C 24010009 */  li    $at, 9
/* AD54B0 8005E310 00002825 */  move  $a1, $zero
/* AD54B4 8005E314 11E1001F */  beq   $t7, $at, .L8005E394
/* AD54B8 8005E318 00000000 */   nop   
/* AD54BC 8005E31C 0C00A729 */  jal   func_80029CA4
/* AD54C0 8005E320 00E03025 */   move  $a2, $a3
/* AD54C4 8005E324 8FB80028 */  lw    $t8, 0x28($sp)
/* AD54C8 8005E328 3C078013 */  lui   $a3, %hi(D_801333E0) # $a3, 0x8013
/* AD54CC 8005E32C 24E733E0 */  addiu $a3, %lo(D_801333E0) # addiu $a3, $a3, 0x33e0
/* AD54D0 8005E330 8F020000 */  lw    $v0, ($t8)
/* AD54D4 8005E334 24041806 */  li    $a0, 6150
/* AD54D8 8005E338 24060004 */  li    $a2, 4
/* AD54DC 8005E33C 1440000E */  bnez  $v0, .L8005E378
/* AD54E0 8005E340 244500E4 */   addiu $a1, $v0, 0xe4
/* AD54E4 8005E344 3C078013 */  lui   $a3, %hi(D_801333E0) # $a3, 0x8013
/* AD54E8 8005E348 3C198013 */  lui   $t9, %hi(D_801333E8) # $t9, 0x8013
/* AD54EC 8005E34C 24E733E0 */  addiu $a3, %lo(D_801333E0) # addiu $a3, $a3, 0x33e0
/* AD54F0 8005E350 273933E8 */  addiu $t9, %lo(D_801333E8) # addiu $t9, $t9, 0x33e8
/* AD54F4 8005E354 3C058013 */  lui   $a1, %hi(D_801333D4) # $a1, 0x8013
/* AD54F8 8005E358 24A533D4 */  addiu $a1, %lo(D_801333D4) # addiu $a1, $a1, 0x33d4
/* AD54FC 8005E35C AFB90014 */  sw    $t9, 0x14($sp)
/* AD5500 8005E360 AFA70010 */  sw    $a3, 0x10($sp)
/* AD5504 8005E364 24041806 */  li    $a0, 6150
/* AD5508 8005E368 0C03DCE3 */  jal   Audio_PlaySoundGeneral
/* AD550C 8005E36C 24060004 */   li    $a2, 4
/* AD5510 8005E370 1000005E */  b     .L8005E4EC
/* AD5514 8005E374 8FBF001C */   lw    $ra, 0x1c($sp)
.L8005E378:
/* AD5518 8005E378 3C088013 */  lui   $t0, %hi(D_801333E8) # $t0, 0x8013
/* AD551C 8005E37C 250833E8 */  addiu $t0, %lo(D_801333E8) # addiu $t0, $t0, 0x33e8
/* AD5520 8005E380 AFA80014 */  sw    $t0, 0x14($sp)
/* AD5524 8005E384 0C03DCE3 */  jal   Audio_PlaySoundGeneral
/* AD5528 8005E388 AFA70010 */   sw    $a3, 0x10($sp)
/* AD552C 8005E38C 10000057 */  b     .L8005E4EC
/* AD5530 8005E390 8FBF001C */   lw    $ra, 0x1c($sp)
.L8005E394:
/* AD5534 8005E394 14400015 */  bnez  $v0, .L8005E3EC
/* AD5538 8005E398 24010008 */   li    $at, 8
/* AD553C 8005E39C 24050003 */  li    $a1, 3
/* AD5540 8005E3A0 00E03025 */  move  $a2, $a3
/* AD5544 8005E3A4 AFA40020 */  sw    $a0, 0x20($sp)
/* AD5548 8005E3A8 0C00A729 */  jal   func_80029CA4
/* AD554C 8005E3AC AFA7002C */   sw    $a3, 0x2c($sp)
/* AD5550 8005E3B0 8FA90028 */  lw    $t1, 0x28($sp)
/* AD5554 8005E3B4 8FA7002C */  lw    $a3, 0x2c($sp)
/* AD5558 8005E3B8 8FA40020 */  lw    $a0, 0x20($sp)
/* AD555C 8005E3BC 8D220000 */  lw    $v0, ($t1)
/* AD5560 8005E3C0 00E02825 */  move  $a1, $a3
/* AD5564 8005E3C4 14400005 */  bnez  $v0, .L8005E3DC
/* AD5568 8005E3C8 00000000 */   nop   
/* AD556C 8005E3CC 0C018B58 */  jal   func_80062D60
/* AD5570 8005E3D0 00E02825 */   move  $a1, $a3
/* AD5574 8005E3D4 10000045 */  b     .L8005E4EC
/* AD5578 8005E3D8 8FBF001C */   lw    $ra, 0x1c($sp)
.L8005E3DC:
/* AD557C 8005E3DC 0C018B6B */  jal   func_80062DAC
/* AD5580 8005E3E0 244600E4 */   addiu $a2, $v0, 0xe4
/* AD5584 8005E3E4 10000041 */  b     .L8005E4EC
/* AD5588 8005E3E8 8FBF001C */   lw    $ra, 0x1c($sp)
.L8005E3EC:
/* AD558C 8005E3EC 1441001F */  bne   $v0, $at, .L8005E46C
/* AD5590 8005E3F0 00002825 */   move  $a1, $zero
/* AD5594 8005E3F4 0C00A729 */  jal   func_80029CA4
/* AD5598 8005E3F8 00E03025 */   move  $a2, $a3
/* AD559C 8005E3FC 8FAA0028 */  lw    $t2, 0x28($sp)
/* AD55A0 8005E400 3C078013 */  lui   $a3, %hi(D_801333E0) # $a3, 0x8013
/* AD55A4 8005E404 24E733E0 */  addiu $a3, %lo(D_801333E0) # addiu $a3, $a3, 0x33e0
/* AD55A8 8005E408 8D420000 */  lw    $v0, ($t2)
/* AD55AC 8005E40C 24041806 */  li    $a0, 6150
/* AD55B0 8005E410 24060004 */  li    $a2, 4
/* AD55B4 8005E414 1440000E */  bnez  $v0, .L8005E450
/* AD55B8 8005E418 244500E4 */   addiu $a1, $v0, 0xe4
/* AD55BC 8005E41C 3C078013 */  lui   $a3, %hi(D_801333E0) # $a3, 0x8013
/* AD55C0 8005E420 3C0B8013 */  lui   $t3, %hi(D_801333E8) # $t3, 0x8013
/* AD55C4 8005E424 24E733E0 */  addiu $a3, %lo(D_801333E0) # addiu $a3, $a3, 0x33e0
/* AD55C8 8005E428 256B33E8 */  addiu $t3, %lo(D_801333E8) # addiu $t3, $t3, 0x33e8
/* AD55CC 8005E42C 3C058013 */  lui   $a1, %hi(D_801333D4) # $a1, 0x8013
/* AD55D0 8005E430 24A533D4 */  addiu $a1, %lo(D_801333D4) # addiu $a1, $a1, 0x33d4
/* AD55D4 8005E434 AFAB0014 */  sw    $t3, 0x14($sp)
/* AD55D8 8005E438 AFA70010 */  sw    $a3, 0x10($sp)
/* AD55DC 8005E43C 24041806 */  li    $a0, 6150
/* AD55E0 8005E440 0C03DCE3 */  jal   Audio_PlaySoundGeneral
/* AD55E4 8005E444 24060004 */   li    $a2, 4
/* AD55E8 8005E448 10000028 */  b     .L8005E4EC
/* AD55EC 8005E44C 8FBF001C */   lw    $ra, 0x1c($sp)
.L8005E450:
/* AD55F0 8005E450 3C0C8013 */  lui   $t4, %hi(D_801333E8) # $t4, 0x8013
/* AD55F4 8005E454 258C33E8 */  addiu $t4, %lo(D_801333E8) # addiu $t4, $t4, 0x33e8
/* AD55F8 8005E458 AFAC0014 */  sw    $t4, 0x14($sp)
/* AD55FC 8005E45C 0C03DCE3 */  jal   Audio_PlaySoundGeneral
/* AD5600 8005E460 AFA70010 */   sw    $a3, 0x10($sp)
/* AD5604 8005E464 10000021 */  b     .L8005E4EC
/* AD5608 8005E468 8FBF001C */   lw    $ra, 0x1c($sp)
.L8005E46C:
/* AD560C 8005E46C 24010010 */  li    $at, 16
/* AD5610 8005E470 1441001D */  bne   $v0, $at, .L8005E4E8
/* AD5614 8005E474 24050001 */   li    $a1, 1
/* AD5618 8005E478 0C00A729 */  jal   func_80029CA4
/* AD561C 8005E47C 00E03025 */   move  $a2, $a3
/* AD5620 8005E480 8FAD0028 */  lw    $t5, 0x28($sp)
/* AD5624 8005E484 3C078013 */  lui   $a3, %hi(D_801333E0) # $a3, 0x8013
/* AD5628 8005E488 24E733E0 */  addiu $a3, %lo(D_801333E0) # addiu $a3, $a3, 0x33e0
/* AD562C 8005E48C 8DA20000 */  lw    $v0, ($t5)
/* AD5630 8005E490 24041837 */  li    $a0, 6199
/* AD5634 8005E494 24060004 */  li    $a2, 4
/* AD5638 8005E498 1440000E */  bnez  $v0, .L8005E4D4
/* AD563C 8005E49C 244500E4 */   addiu $a1, $v0, 0xe4
/* AD5640 8005E4A0 3C078013 */  lui   $a3, %hi(D_801333E0) # $a3, 0x8013
/* AD5644 8005E4A4 3C0E8013 */  lui   $t6, %hi(D_801333E8) # $t6, 0x8013
/* AD5648 8005E4A8 24E733E0 */  addiu $a3, %lo(D_801333E0) # addiu $a3, $a3, 0x33e0
/* AD564C 8005E4AC 25CE33E8 */  addiu $t6, %lo(D_801333E8) # addiu $t6, $t6, 0x33e8
/* AD5650 8005E4B0 3C058013 */  lui   $a1, %hi(D_801333D4) # $a1, 0x8013
/* AD5654 8005E4B4 24A533D4 */  addiu $a1, %lo(D_801333D4) # addiu $a1, $a1, 0x33d4
/* AD5658 8005E4B8 AFAE0014 */  sw    $t6, 0x14($sp)
/* AD565C 8005E4BC AFA70010 */  sw    $a3, 0x10($sp)
/* AD5660 8005E4C0 24041837 */  li    $a0, 6199
/* AD5664 8005E4C4 0C03DCE3 */  jal   Audio_PlaySoundGeneral
/* AD5668 8005E4C8 24060004 */   li    $a2, 4
/* AD566C 8005E4CC 10000007 */  b     .L8005E4EC
/* AD5670 8005E4D0 8FBF001C */   lw    $ra, 0x1c($sp)
.L8005E4D4:
/* AD5674 8005E4D4 3C0F8013 */  lui   $t7, %hi(D_801333E8) # $t7, 0x8013
/* AD5678 8005E4D8 25EF33E8 */  addiu $t7, %lo(D_801333E8) # addiu $t7, $t7, 0x33e8
/* AD567C 8005E4DC AFAF0014 */  sw    $t7, 0x14($sp)
/* AD5680 8005E4E0 0C03DCE3 */  jal   Audio_PlaySoundGeneral
/* AD5684 8005E4E4 AFA70010 */   sw    $a3, 0x10($sp)
.L8005E4E8:
/* AD5688 8005E4E8 8FBF001C */  lw    $ra, 0x1c($sp)
.L8005E4EC:
/* AD568C 8005E4EC 27BD0020 */  addiu $sp, $sp, 0x20
/* AD5690 8005E4F0 03E00008 */  jr    $ra
/* AD5694 8005E4F4 00000000 */   nop   

