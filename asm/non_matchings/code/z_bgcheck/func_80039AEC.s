glabel func_80039AEC
/* AB0C8C 80039AEC 27BDFEF8 */  addiu $sp, $sp, -0x108
/* AB0C90 80039AF0 AFBF008C */  sw    $ra, 0x8c($sp)
/* AB0C94 80039AF4 AFBE0088 */  sw    $fp, 0x88($sp)
/* AB0C98 80039AF8 AFB70084 */  sw    $s7, 0x84($sp)
/* AB0C9C 80039AFC AFB60080 */  sw    $s6, 0x80($sp)
/* AB0CA0 80039B00 AFB5007C */  sw    $s5, 0x7c($sp)
/* AB0CA4 80039B04 AFB40078 */  sw    $s4, 0x78($sp)
/* AB0CA8 80039B08 AFB30074 */  sw    $s3, 0x74($sp)
/* AB0CAC 80039B0C AFB20070 */  sw    $s2, 0x70($sp)
/* AB0CB0 80039B10 AFB1006C */  sw    $s1, 0x6c($sp)
/* AB0CB4 80039B14 AFB00068 */  sw    $s0, 0x68($sp)
/* AB0CB8 80039B18 F7BE0060 */  sdc1  $f30, 0x60($sp)
/* AB0CBC 80039B1C F7BC0058 */  sdc1  $f28, 0x58($sp)
/* AB0CC0 80039B20 F7BA0050 */  sdc1  $f26, 0x50($sp)
/* AB0CC4 80039B24 F7B80048 */  sdc1  $f24, 0x48($sp)
/* AB0CC8 80039B28 F7B60040 */  sdc1  $f22, 0x40($sp)
/* AB0CCC 80039B2C F7B40038 */  sdc1  $f20, 0x38($sp)
/* AB0CD0 80039B30 AFA40108 */  sw    $a0, 0x108($sp)
/* AB0CD4 80039B34 AFA60110 */  sw    $a2, 0x110($sp)
/* AB0CD8 80039B38 AFA70114 */  sw    $a3, 0x114($sp)
/* AB0CDC 80039B3C AFA000E8 */  sw    $zero, 0xe8($sp)
/* AB0CE0 80039B40 948F0002 */  lhu   $t7, 2($a0)
/* AB0CE4 80039B44 3414FFFF */  li    $s4, 65535
/* AB0CE8 80039B48 00A09825 */  move  $s3, $a1
/* AB0CEC 80039B4C 168F0003 */  bne   $s4, $t7, .L80039B5C
/* AB0CF0 80039B50 8FB6011C */   lw    $s6, 0x11c($sp)
/* AB0CF4 80039B54 10000210 */  b     .L8003A398
/* AB0CF8 80039B58 00001025 */   move  $v0, $zero
.L80039B5C:
/* AB0CFC 80039B5C 8ED90000 */  lw    $t9, ($s6)
/* AB0D00 80039B60 27B700FC */  addiu $s7, $sp, 0xfc
/* AB0D04 80039B64 241E0006 */  li    $fp, 6
/* AB0D08 80039B68 AEF90000 */  sw    $t9, ($s7)
/* AB0D0C 80039B6C 8ED80004 */  lw    $t8, 4($s6)
/* AB0D10 80039B70 24150006 */  li    $s5, 6
/* AB0D14 80039B74 AEF80004 */  sw    $t8, 4($s7)
/* AB0D18 80039B78 8ED90008 */  lw    $t9, 8($s6)
/* AB0D1C 80039B7C AEF90008 */  sw    $t9, 8($s7)
/* AB0D20 80039B80 8E620000 */  lw    $v0, ($s3)
/* AB0D24 80039B84 8FA90108 */  lw    $t1, 0x108($sp)
/* AB0D28 80039B88 C7BC0120 */  lwc1  $f28, 0x120($sp)
/* AB0D2C 80039B8C 8C480018 */  lw    $t0, 0x18($v0)
/* AB0D30 80039B90 AFA800E0 */  sw    $t0, 0xe0($sp)
/* AB0D34 80039B94 952A0002 */  lhu   $t2, 2($t1)
/* AB0D38 80039B98 8E630048 */  lw    $v1, 0x48($s3)
/* AB0D3C 80039B9C 8C520010 */  lw    $s2, 0x10($v0)
/* AB0D40 80039BA0 000A5880 */  sll   $t3, $t2, 2
/* AB0D44 80039BA4 006B8821 */  addu  $s1, $v1, $t3
.L80039BA8:
/* AB0D48 80039BA8 86220000 */  lh    $v0, ($s1)
/* AB0D4C 80039BAC 8FAD00E0 */  lw    $t5, 0xe0($sp)
/* AB0D50 80039BB0 C6C00004 */  lwc1  $f0, 4($s6)
/* AB0D54 80039BB4 00026100 */  sll   $t4, $v0, 4
/* AB0D58 80039BB8 018D8021 */  addu  $s0, $t4, $t5
/* AB0D5C 80039BBC 960E0002 */  lhu   $t6, 2($s0)
/* AB0D60 80039BC0 31CF1FFF */  andi  $t7, $t6, 0x1fff
/* AB0D64 80039BC4 01F50019 */  multu $t7, $s5
/* AB0D68 80039BC8 0000C012 */  mflo  $t8
/* AB0D6C 80039BCC 0258C821 */  addu  $t9, $s2, $t8
/* AB0D70 80039BD0 87280002 */  lh    $t0, 2($t9)
/* AB0D74 80039BD4 44882000 */  mtc1  $t0, $f4
/* AB0D78 80039BD8 00000000 */  nop   
/* AB0D7C 80039BDC 468021A0 */  cvt.s.w $f6, $f4
/* AB0D80 80039BE0 4606003C */  c.lt.s $f0, $f6
/* AB0D84 80039BE4 00000000 */  nop   
/* AB0D88 80039BE8 45020020 */  bc1fl .L80039C6C
/* AB0D8C 80039BEC 860B0008 */   lh    $t3, 8($s0)
/* AB0D90 80039BF0 96090004 */  lhu   $t1, 4($s0)
/* AB0D94 80039BF4 312A1FFF */  andi  $t2, $t1, 0x1fff
/* AB0D98 80039BF8 01550019 */  multu $t2, $s5
/* AB0D9C 80039BFC 00005812 */  mflo  $t3
/* AB0DA0 80039C00 024B6021 */  addu  $t4, $s2, $t3
/* AB0DA4 80039C04 858D0002 */  lh    $t5, 2($t4)
/* AB0DA8 80039C08 448D4000 */  mtc1  $t5, $f8
/* AB0DAC 80039C0C 00000000 */  nop   
/* AB0DB0 80039C10 468042A0 */  cvt.s.w $f10, $f8
/* AB0DB4 80039C14 460A003C */  c.lt.s $f0, $f10
/* AB0DB8 80039C18 00000000 */  nop   
/* AB0DBC 80039C1C 45020013 */  bc1fl .L80039C6C
/* AB0DC0 80039C20 860B0008 */   lh    $t3, 8($s0)
/* AB0DC4 80039C24 960E0006 */  lhu   $t6, 6($s0)
/* AB0DC8 80039C28 01DE0019 */  multu $t6, $fp
/* AB0DCC 80039C2C 00007812 */  mflo  $t7
/* AB0DD0 80039C30 024FC021 */  addu  $t8, $s2, $t7
/* AB0DD4 80039C34 87190002 */  lh    $t9, 2($t8)
/* AB0DD8 80039C38 44992000 */  mtc1  $t9, $f4
/* AB0DDC 80039C3C 00000000 */  nop   
/* AB0DE0 80039C40 468021A0 */  cvt.s.w $f6, $f4
/* AB0DE4 80039C44 4606003C */  c.lt.s $f0, $f6
/* AB0DE8 80039C48 00000000 */  nop   
/* AB0DEC 80039C4C 45020007 */  bc1fl .L80039C6C
/* AB0DF0 80039C50 860B0008 */   lh    $t3, 8($s0)
/* AB0DF4 80039C54 8FA80108 */  lw    $t0, 0x108($sp)
/* AB0DF8 80039C58 95090002 */  lhu   $t1, 2($t0)
/* AB0DFC 80039C5C 00095080 */  sll   $t2, $t1, 2
/* AB0E00 80039C60 100000DC */  b     .L80039FD4
/* AB0E04 80039C64 006A8821 */   addu  $s1, $v1, $t2
/* AB0E08 80039C68 860B0008 */  lh    $t3, 8($s0)
.L80039C6C:
/* AB0E0C 80039C6C 860C000A */  lh    $t4, 0xa($s0)
/* AB0E10 80039C70 860D000C */  lh    $t5, 0xc($s0)
/* AB0E14 80039C74 448B4000 */  mtc1  $t3, $f8
/* AB0E18 80039C78 448C2000 */  mtc1  $t4, $f4
/* AB0E1C 80039C7C 3C018014 */  lui   $at, %hi(D_80138F5C)
/* AB0E20 80039C80 468042A0 */  cvt.s.w $f10, $f8
/* AB0E24 80039C84 C4228F5C */  lwc1  $f2, %lo(D_80138F5C)($at)
/* AB0E28 80039C88 448D4000 */  mtc1  $t5, $f8
/* AB0E2C 80039C8C 860E000E */  lh    $t6, 0xe($s0)
/* AB0E30 80039C90 AFB70010 */  sw    $s7, 0x10($sp)
/* AB0E34 80039C94 468021A0 */  cvt.s.w $f6, $f4
/* AB0E38 80039C98 46025582 */  mul.s $f22, $f10, $f2
/* AB0E3C 80039C9C 468042A0 */  cvt.s.w $f10, $f8
/* AB0E40 80039CA0 46023682 */  mul.s $f26, $f6, $f2
/* AB0E44 80039CA4 448E4000 */  mtc1  $t6, $f8
/* AB0E48 80039CA8 4600B306 */  mov.s $f12, $f22
/* AB0E4C 80039CAC 46025602 */  mul.s $f24, $f10, $f2
/* AB0E50 80039CB0 4600D386 */  mov.s $f14, $f26
/* AB0E54 80039CB4 4616B102 */  mul.s $f4, $f22, $f22
/* AB0E58 80039CB8 00000000 */  nop   
/* AB0E5C 80039CBC 4618C182 */  mul.s $f6, $f24, $f24
/* AB0E60 80039CC0 4406C000 */  mfc1  $a2, $f24
/* AB0E64 80039CC4 46804220 */  cvt.s.w $f8, $f8
/* AB0E68 80039CC8 46062000 */  add.s $f0, $f4, $f6
/* AB0E6C 80039CCC 44074000 */  mfc1  $a3, $f8
/* AB0E70 80039CD0 0C0332C3 */  jal   func_800CCB0C
/* AB0E74 80039CD4 46000504 */   sqrt.s $f20, $f0
/* AB0E78 80039CD8 46000786 */  mov.s $f30, $f0
/* AB0E7C 80039CDC 46000005 */  abs.s $f0, $f0
/* AB0E80 80039CE0 97B80112 */  lhu   $t8, 0x112($sp)
/* AB0E84 80039CE4 4600E03C */  c.lt.s $f28, $f0
/* AB0E88 80039CE8 00000000 */  nop   
/* AB0E8C 80039CEC 45030008 */  bc1tl .L80039D10
/* AB0E90 80039CF0 96220002 */   lhu   $v0, 2($s1)
/* AB0E94 80039CF4 960F0002 */  lhu   $t7, 2($s0)
/* AB0E98 80039CF8 33190007 */  andi  $t9, $t8, 7
/* AB0E9C 80039CFC 00194340 */  sll   $t0, $t9, 0xd
/* AB0EA0 80039D00 01E84824 */  and   $t1, $t7, $t0
/* AB0EA4 80039D04 1120000F */  beqz  $t1, .L80039D44
/* AB0EA8 80039D08 3C018014 */   lui   $at, %hi(D_80138F60)
/* AB0EAC 80039D0C 96220002 */  lhu   $v0, 2($s1)
.L80039D10:
/* AB0EB0 80039D10 56820009 */  bnel  $s4, $v0, .L80039D38
/* AB0EB4 80039D14 8E630048 */   lw    $v1, 0x48($s3)
/* AB0EB8 80039D18 8FAB0108 */  lw    $t3, 0x108($sp)
/* AB0EBC 80039D1C 8E6A0048 */  lw    $t2, 0x48($s3)
/* AB0EC0 80039D20 C6C00004 */  lwc1  $f0, 4($s6)
/* AB0EC4 80039D24 956C0002 */  lhu   $t4, 2($t3)
/* AB0EC8 80039D28 000C6880 */  sll   $t5, $t4, 2
/* AB0ECC 80039D2C 100000A9 */  b     .L80039FD4
/* AB0ED0 80039D30 014D8821 */   addu  $s1, $t2, $t5
/* AB0ED4 80039D34 8E630048 */  lw    $v1, 0x48($s3)
.L80039D38:
/* AB0ED8 80039D38 00027080 */  sll   $t6, $v0, 2
/* AB0EDC 80039D3C 1000FF9A */  b     .L80039BA8
/* AB0EE0 80039D40 006E8821 */   addu  $s1, $v1, $t6
.L80039D44:
/* AB0EE4 80039D44 C42A8F60 */  lwc1  $f10, %lo(D_80138F60)($at)
/* AB0EE8 80039D48 4600A005 */  abs.s $f0, $f20
/* AB0EEC 80039D4C 3C048014 */  lui   $a0, %hi(D_80138784) # $a0, 0x8014
/* AB0EF0 80039D50 460A003C */  c.lt.s $f0, $f10
/* AB0EF4 80039D54 3C058014 */  lui   $a1, %hi(D_80138798) # $a1, 0x8014
/* AB0EF8 80039D58 24A58798 */  addiu $a1, %lo(D_80138798) # addiu $a1, $a1, -0x7868
/* AB0EFC 80039D5C 24848784 */  addiu $a0, %lo(D_80138784) # addiu $a0, $a0, -0x787c
/* AB0F00 80039D60 45020004 */  bc1fl .L80039D74
/* AB0F04 80039D64 3C013F80 */   lui   $at, 0x3f80
/* AB0F08 80039D68 0C0007FC */  jal   __assert
/* AB0F0C 80039D6C 24060B26 */   li    $a2, 2854
/* AB0F10 80039D70 3C013F80 */  li    $at, 0x3F800000 # 0.000000
.L80039D74:
/* AB0F14 80039D74 44812000 */  mtc1  $at, $f4
/* AB0F18 80039D78 4600C005 */  abs.s $f0, $f24
/* AB0F1C 80039D7C 3C018014 */  lui   $at, %hi(D_80138F64)
/* AB0F20 80039D80 46142383 */  div.s $f14, $f4, $f20
/* AB0F24 80039D84 C4268F64 */  lwc1  $f6, %lo(D_80138F64)($at)
/* AB0F28 80039D88 460E0402 */  mul.s $f16, $f0, $f14
/* AB0F2C 80039D8C 4606803C */  c.lt.s $f16, $f6
/* AB0F30 80039D90 00000000 */  nop   
/* AB0F34 80039D94 45020010 */  bc1fl .L80039DD8
/* AB0F38 80039D98 960B0002 */   lhu   $t3, 2($s0)
/* AB0F3C 80039D9C 96220002 */  lhu   $v0, 2($s1)
/* AB0F40 80039DA0 56820009 */  bnel  $s4, $v0, .L80039DC8
/* AB0F44 80039DA4 8E630048 */   lw    $v1, 0x48($s3)
/* AB0F48 80039DA8 8FB90108 */  lw    $t9, 0x108($sp)
/* AB0F4C 80039DAC 8E780048 */  lw    $t8, 0x48($s3)
/* AB0F50 80039DB0 C6C00004 */  lwc1  $f0, 4($s6)
/* AB0F54 80039DB4 972F0002 */  lhu   $t7, 2($t9)
/* AB0F58 80039DB8 000F4080 */  sll   $t0, $t7, 2
/* AB0F5C 80039DBC 10000085 */  b     .L80039FD4
/* AB0F60 80039DC0 03088821 */   addu  $s1, $t8, $t0
/* AB0F64 80039DC4 8E630048 */  lw    $v1, 0x48($s3)
.L80039DC8:
/* AB0F68 80039DC8 00024880 */  sll   $t1, $v0, 2
/* AB0F6C 80039DCC 1000FF76 */  b     .L80039BA8
/* AB0F70 80039DD0 00698821 */   addu  $s1, $v1, $t1
/* AB0F74 80039DD4 960B0002 */  lhu   $t3, 2($s0)
.L80039DD8:
/* AB0F78 80039DD8 96190004 */  lhu   $t9, 4($s0)
/* AB0F7C 80039DDC 316C1FFF */  andi  $t4, $t3, 0x1fff
/* AB0F80 80039DE0 01950019 */  multu $t4, $s5
/* AB0F84 80039DE4 332F1FFF */  andi  $t7, $t9, 0x1fff
/* AB0F88 80039DE8 00005012 */  mflo  $t2
/* AB0F8C 80039DEC 024A6821 */  addu  $t5, $s2, $t2
/* AB0F90 80039DF0 85AE0004 */  lh    $t6, 4($t5)
/* AB0F94 80039DF4 01F50019 */  multu $t7, $s5
/* AB0F98 80039DF8 448E4000 */  mtc1  $t6, $f8
/* AB0F9C 80039DFC 00000000 */  nop   
/* AB0FA0 80039E00 46804320 */  cvt.s.w $f12, $f8
/* AB0FA4 80039E04 0000C012 */  mflo  $t8
/* AB0FA8 80039E08 02584021 */  addu  $t0, $s2, $t8
/* AB0FAC 80039E0C 85090004 */  lh    $t1, 4($t0)
/* AB0FB0 80039E10 46006086 */  mov.s $f2, $f12
/* AB0FB4 80039E14 44895000 */  mtc1  $t1, $f10
/* AB0FB8 80039E18 00000000 */  nop   
/* AB0FBC 80039E1C 46805020 */  cvt.s.w $f0, $f10
/* AB0FC0 80039E20 4602003C */  c.lt.s $f0, $f2
/* AB0FC4 80039E24 00000000 */  nop   
/* AB0FC8 80039E28 45020004 */  bc1fl .L80039E3C
/* AB0FCC 80039E2C 4600603C */   c.lt.s $f12, $f0
/* AB0FD0 80039E30 10000006 */  b     .L80039E4C
/* AB0FD4 80039E34 46000086 */   mov.s $f2, $f0
/* AB0FD8 80039E38 4600603C */  c.lt.s $f12, $f0
.L80039E3C:
/* AB0FDC 80039E3C 00000000 */  nop   
/* AB0FE0 80039E40 45020003 */  bc1fl .L80039E50
/* AB0FE4 80039E44 960B0006 */   lhu   $t3, 6($s0)
/* AB0FE8 80039E48 46000306 */  mov.s $f12, $f0
.L80039E4C:
/* AB0FEC 80039E4C 960B0006 */  lhu   $t3, 6($s0)
.L80039E50:
/* AB0FF0 80039E50 017E0019 */  multu $t3, $fp
/* AB0FF4 80039E54 00006012 */  mflo  $t4
/* AB0FF8 80039E58 024C5021 */  addu  $t2, $s2, $t4
/* AB0FFC 80039E5C 854D0004 */  lh    $t5, 4($t2)
/* AB1000 80039E60 448D2000 */  mtc1  $t5, $f4
/* AB1004 80039E64 00000000 */  nop   
/* AB1008 80039E68 46802020 */  cvt.s.w $f0, $f4
/* AB100C 80039E6C 4602003C */  c.lt.s $f0, $f2
/* AB1010 80039E70 00000000 */  nop   
/* AB1014 80039E74 45020004 */  bc1fl .L80039E88
/* AB1018 80039E78 4600603C */   c.lt.s $f12, $f0
/* AB101C 80039E7C 10000006 */  b     .L80039E98
/* AB1020 80039E80 46000086 */   mov.s $f2, $f0
/* AB1024 80039E84 4600603C */  c.lt.s $f12, $f0
.L80039E88:
/* AB1028 80039E88 00000000 */  nop   
/* AB102C 80039E8C 45020003 */  bc1fl .L80039E9C
/* AB1030 80039E90 461C1081 */   sub.s $f2, $f2, $f28
/* AB1034 80039E94 46000306 */  mov.s $f12, $f0
.L80039E98:
/* AB1038 80039E98 461C1081 */  sub.s $f2, $f2, $f28
.L80039E9C:
/* AB103C 80039E9C C7B20104 */  lwc1  $f18, 0x104($sp)
/* AB1040 80039EA0 461C6300 */  add.s $f12, $f12, $f28
/* AB1044 80039EA4 4602903C */  c.lt.s $f18, $f2
/* AB1048 80039EA8 00000000 */  nop   
/* AB104C 80039EAC 45030008 */  bc1tl .L80039ED0
/* AB1050 80039EB0 96220002 */   lhu   $v0, 2($s1)
/* AB1054 80039EB4 4612603C */  c.lt.s $f12, $f18
/* AB1058 80039EB8 02002025 */  move  $a0, $s0
/* AB105C 80039EBC 02402825 */  move  $a1, $s2
/* AB1060 80039EC0 8FA600FC */  lw    $a2, 0xfc($sp)
/* AB1064 80039EC4 4500000F */  bc1f  .L80039F04
/* AB1068 80039EC8 27A900EC */   addiu $t1, $sp, 0xec
/* AB106C 80039ECC 96220002 */  lhu   $v0, 2($s1)
.L80039ED0:
/* AB1070 80039ED0 56820009 */  bnel  $s4, $v0, .L80039EF8
/* AB1074 80039ED4 8E630048 */   lw    $v1, 0x48($s3)
/* AB1078 80039ED8 8FB90108 */  lw    $t9, 0x108($sp)
/* AB107C 80039EDC 8E6E0048 */  lw    $t6, 0x48($s3)
/* AB1080 80039EE0 C6C00004 */  lwc1  $f0, 4($s6)
/* AB1084 80039EE4 972F0002 */  lhu   $t7, 2($t9)
/* AB1088 80039EE8 000FC080 */  sll   $t8, $t7, 2
/* AB108C 80039EEC 10000039 */  b     .L80039FD4
/* AB1090 80039EF0 01D88821 */   addu  $s1, $t6, $t8
/* AB1094 80039EF4 8E630048 */  lw    $v1, 0x48($s3)
.L80039EF8:
/* AB1098 80039EF8 00024080 */  sll   $t0, $v0, 2
/* AB109C 80039EFC 1000FF2A */  b     .L80039BA8
/* AB10A0 80039F00 00688821 */   addu  $s1, $v1, $t0
.L80039F04:
/* AB10A4 80039F04 8EC70004 */  lw    $a3, 4($s6)
/* AB10A8 80039F08 E7B000B8 */  swc1  $f16, 0xb8($sp)
/* AB10AC 80039F0C E7AE009C */  swc1  $f14, 0x9c($sp)
/* AB10B0 80039F10 0C00E400 */  jal   func_80039000
/* AB10B4 80039F14 AFA90010 */   sw    $t1, 0x10($sp)
/* AB10B8 80039F18 C7AE009C */  lwc1  $f14, 0x9c($sp)
/* AB10BC 80039F1C 1040001F */  beqz  $v0, .L80039F9C
/* AB10C0 80039F20 C7B000B8 */   lwc1  $f16, 0xb8($sp)
/* AB10C4 80039F24 4610E203 */  div.s $f8, $f28, $f16
/* AB10C8 80039F28 C7B20104 */  lwc1  $f18, 0x104($sp)
/* AB10CC 80039F2C C7A600EC */  lwc1  $f6, 0xec($sp)
/* AB10D0 80039F30 46123081 */  sub.s $f2, $f6, $f18
/* AB10D4 80039F34 46001005 */  abs.s $f0, $f2
/* AB10D8 80039F38 4608003E */  c.le.s $f0, $f8
/* AB10DC 80039F3C 00000000 */  nop   
/* AB10E0 80039F40 45020017 */  bc1fl .L80039FA0
/* AB10E4 80039F44 96220002 */   lhu   $v0, 2($s1)
/* AB10E8 80039F48 46181102 */  mul.s $f4, $f2, $f24
/* AB10EC 80039F4C 3C014080 */  li    $at, 0x40800000 # 0.000000
/* AB10F0 80039F50 44815000 */  mtc1  $at, $f10
/* AB10F4 80039F54 02602025 */  move  $a0, $s3
/* AB10F8 80039F58 02002825 */  move  $a1, $s0
/* AB10FC 80039F5C 02E03025 */  move  $a2, $s7
/* AB1100 80039F60 27A70104 */  addiu $a3, $sp, 0x104
/* AB1104 80039F64 460A203E */  c.le.s $f4, $f10
/* AB1108 80039F68 240C0001 */  li    $t4, 1
/* AB110C 80039F6C 4502000C */  bc1fl .L80039FA0
/* AB1110 80039F70 96220002 */   lhu   $v0, 2($s1)
/* AB1114 80039F74 8FAB0124 */  lw    $t3, 0x124($sp)
/* AB1118 80039F78 E7B60010 */  swc1  $f22, 0x10($sp)
/* AB111C 80039F7C E7BA0014 */  swc1  $f26, 0x14($sp)
/* AB1120 80039F80 E7B80018 */  swc1  $f24, 0x18($sp)
/* AB1124 80039F84 E7AE001C */  swc1  $f14, 0x1c($sp)
/* AB1128 80039F88 E7BE0020 */  swc1  $f30, 0x20($sp)
/* AB112C 80039F8C E7BC0024 */  swc1  $f28, 0x24($sp)
/* AB1130 80039F90 AFAC00E8 */  sw    $t4, 0xe8($sp)
/* AB1134 80039F94 0C00E68F */  jal   func_80039A3C
/* AB1138 80039F98 AFAB0028 */   sw    $t3, 0x28($sp)
.L80039F9C:
/* AB113C 80039F9C 96220002 */  lhu   $v0, 2($s1)
.L80039FA0:
/* AB1140 80039FA0 56820009 */  bnel  $s4, $v0, .L80039FC8
/* AB1144 80039FA4 8E630048 */   lw    $v1, 0x48($s3)
/* AB1148 80039FA8 8FAD0108 */  lw    $t5, 0x108($sp)
/* AB114C 80039FAC 8E6A0048 */  lw    $t2, 0x48($s3)
/* AB1150 80039FB0 C6C00004 */  lwc1  $f0, 4($s6)
/* AB1154 80039FB4 95B90002 */  lhu   $t9, 2($t5)
/* AB1158 80039FB8 00197880 */  sll   $t7, $t9, 2
/* AB115C 80039FBC 10000005 */  b     .L80039FD4
/* AB1160 80039FC0 014F8821 */   addu  $s1, $t2, $t7
/* AB1164 80039FC4 8E630048 */  lw    $v1, 0x48($s3)
.L80039FC8:
/* AB1168 80039FC8 00027080 */  sll   $t6, $v0, 2
/* AB116C 80039FCC 1000FEF6 */  b     .L80039BA8
/* AB1170 80039FD0 006E8821 */   addu  $s1, $v1, $t6
.L80039FD4:
/* AB1174 80039FD4 86220000 */  lh    $v0, ($s1)
/* AB1178 80039FD8 8FA800E0 */  lw    $t0, 0xe0($sp)
/* AB117C 80039FDC 0002C100 */  sll   $t8, $v0, 4
/* AB1180 80039FE0 03088021 */  addu  $s0, $t8, $t0
/* AB1184 80039FE4 96090002 */  lhu   $t1, 2($s0)
/* AB1188 80039FE8 312B1FFF */  andi  $t3, $t1, 0x1fff
/* AB118C 80039FEC 01750019 */  multu $t3, $s5
/* AB1190 80039FF0 00006012 */  mflo  $t4
/* AB1194 80039FF4 024C6821 */  addu  $t5, $s2, $t4
/* AB1198 80039FF8 85B90002 */  lh    $t9, 2($t5)
/* AB119C 80039FFC 44993000 */  mtc1  $t9, $f6
/* AB11A0 8003A000 00000000 */  nop   
/* AB11A4 8003A004 46803220 */  cvt.s.w $f8, $f6
/* AB11A8 8003A008 4608003C */  c.lt.s $f0, $f8
/* AB11AC 8003A00C 00000000 */  nop   
/* AB11B0 8003A010 4502001B */  bc1fl .L8003A080
/* AB11B4 8003A014 86190008 */   lh    $t9, 8($s0)
/* AB11B8 8003A018 960A0004 */  lhu   $t2, 4($s0)
/* AB11BC 8003A01C 314F1FFF */  andi  $t7, $t2, 0x1fff
/* AB11C0 8003A020 01F50019 */  multu $t7, $s5
/* AB11C4 8003A024 00007012 */  mflo  $t6
/* AB11C8 8003A028 024EC021 */  addu  $t8, $s2, $t6
/* AB11CC 8003A02C 87080002 */  lh    $t0, 2($t8)
/* AB11D0 8003A030 44885000 */  mtc1  $t0, $f10
/* AB11D4 8003A034 00000000 */  nop   
/* AB11D8 8003A038 46805120 */  cvt.s.w $f4, $f10
/* AB11DC 8003A03C 4604003C */  c.lt.s $f0, $f4
/* AB11E0 8003A040 00000000 */  nop   
/* AB11E4 8003A044 4502000E */  bc1fl .L8003A080
/* AB11E8 8003A048 86190008 */   lh    $t9, 8($s0)
/* AB11EC 8003A04C 96090006 */  lhu   $t1, 6($s0)
/* AB11F0 8003A050 013E0019 */  multu $t1, $fp
/* AB11F4 8003A054 00005812 */  mflo  $t3
/* AB11F8 8003A058 024B6021 */  addu  $t4, $s2, $t3
/* AB11FC 8003A05C 858D0002 */  lh    $t5, 2($t4)
/* AB1200 8003A060 448D3000 */  mtc1  $t5, $f6
/* AB1204 8003A064 00000000 */  nop   
/* AB1208 8003A068 46803220 */  cvt.s.w $f8, $f6
/* AB120C 8003A06C 4608003C */  c.lt.s $f0, $f8
/* AB1210 8003A070 00000000 */  nop   
/* AB1214 8003A074 450300C2 */  bc1tl .L8003A380
/* AB1218 8003A078 C7A800FC */   lwc1  $f8, 0xfc($sp)
/* AB121C 8003A07C 86190008 */  lh    $t9, 8($s0)
.L8003A080:
/* AB1220 8003A080 860A000A */  lh    $t2, 0xa($s0)
/* AB1224 8003A084 860F000C */  lh    $t7, 0xc($s0)
/* AB1228 8003A088 44995000 */  mtc1  $t9, $f10
/* AB122C 8003A08C 448A3000 */  mtc1  $t2, $f6
/* AB1230 8003A090 3C018014 */  lui   $at, %hi(D_80138F68)
/* AB1234 8003A094 46805120 */  cvt.s.w $f4, $f10
/* AB1238 8003A098 C4228F68 */  lwc1  $f2, %lo(D_80138F68)($at)
/* AB123C 8003A09C 448F5000 */  mtc1  $t7, $f10
/* AB1240 8003A0A0 860E000E */  lh    $t6, 0xe($s0)
/* AB1244 8003A0A4 AFB70010 */  sw    $s7, 0x10($sp)
/* AB1248 8003A0A8 46803220 */  cvt.s.w $f8, $f6
/* AB124C 8003A0AC 46022582 */  mul.s $f22, $f4, $f2
/* AB1250 8003A0B0 46805120 */  cvt.s.w $f4, $f10
/* AB1254 8003A0B4 46024682 */  mul.s $f26, $f8, $f2
/* AB1258 8003A0B8 448E5000 */  mtc1  $t6, $f10
/* AB125C 8003A0BC 4600B306 */  mov.s $f12, $f22
/* AB1260 8003A0C0 46022602 */  mul.s $f24, $f4, $f2
/* AB1264 8003A0C4 4600D386 */  mov.s $f14, $f26
/* AB1268 8003A0C8 4616B182 */  mul.s $f6, $f22, $f22
/* AB126C 8003A0CC 00000000 */  nop   
/* AB1270 8003A0D0 4618C202 */  mul.s $f8, $f24, $f24
/* AB1274 8003A0D4 4406C000 */  mfc1  $a2, $f24
/* AB1278 8003A0D8 468052A0 */  cvt.s.w $f10, $f10
/* AB127C 8003A0DC 46083000 */  add.s $f0, $f6, $f8
/* AB1280 8003A0E0 44075000 */  mfc1  $a3, $f10
/* AB1284 8003A0E4 0C0332C3 */  jal   func_800CCB0C
/* AB1288 8003A0E8 46000504 */   sqrt.s $f20, $f0
/* AB128C 8003A0EC 46000786 */  mov.s $f30, $f0
/* AB1290 8003A0F0 46000005 */  abs.s $f0, $f0
/* AB1294 8003A0F4 97A80112 */  lhu   $t0, 0x112($sp)
/* AB1298 8003A0F8 4600E03C */  c.lt.s $f28, $f0
/* AB129C 8003A0FC 00000000 */  nop   
/* AB12A0 8003A100 45030008 */  bc1tl .L8003A124
/* AB12A4 8003A104 96220002 */   lhu   $v0, 2($s1)
/* AB12A8 8003A108 96180002 */  lhu   $t8, 2($s0)
/* AB12AC 8003A10C 31090007 */  andi  $t1, $t0, 7
/* AB12B0 8003A110 00095B40 */  sll   $t3, $t1, 0xd
/* AB12B4 8003A114 030B6024 */  and   $t4, $t8, $t3
/* AB12B8 8003A118 11800008 */  beqz  $t4, .L8003A13C
/* AB12BC 8003A11C 3C018014 */   lui   $at, %hi(D_80138F6C)
/* AB12C0 8003A120 96220002 */  lhu   $v0, 2($s1)
.L8003A124:
/* AB12C4 8003A124 52820096 */  beql  $s4, $v0, .L8003A380
/* AB12C8 8003A128 C7A800FC */   lwc1  $f8, 0xfc($sp)
/* AB12CC 8003A12C 8E6D0048 */  lw    $t5, 0x48($s3)
/* AB12D0 8003A130 0002C880 */  sll   $t9, $v0, 2
/* AB12D4 8003A134 1000008F */  b     .L8003A374
/* AB12D8 8003A138 01B98821 */   addu  $s1, $t5, $t9
.L8003A13C:
/* AB12DC 8003A13C C4248F6C */  lwc1  $f4, %lo(D_80138F6C)($at)
/* AB12E0 8003A140 4600A005 */  abs.s $f0, $f20
/* AB12E4 8003A144 3C048014 */  lui   $a0, %hi(D_801387A8) # $a0, 0x8014
/* AB12E8 8003A148 4604003C */  c.lt.s $f0, $f4
/* AB12EC 8003A14C 3C058014 */  lui   $a1, %hi(D_801387BC) # $a1, 0x8014
/* AB12F0 8003A150 24A587BC */  addiu $a1, %lo(D_801387BC) # addiu $a1, $a1, -0x7844
/* AB12F4 8003A154 248487A8 */  addiu $a0, %lo(D_801387A8) # addiu $a0, $a0, -0x7858
/* AB12F8 8003A158 45020004 */  bc1fl .L8003A16C
/* AB12FC 8003A15C 3C013F80 */   lui   $at, 0x3f80
/* AB1300 8003A160 0C0007FC */  jal   __assert
/* AB1304 8003A164 24060B94 */   li    $a2, 2964
/* AB1308 8003A168 3C013F80 */  li    $at, 0x3F800000 # 0.000000
.L8003A16C:
/* AB130C 8003A16C 44813000 */  mtc1  $at, $f6
/* AB1310 8003A170 4600B005 */  abs.s $f0, $f22
/* AB1314 8003A174 3C018014 */  lui   $at, %hi(D_80138F70)
/* AB1318 8003A178 46143383 */  div.s $f14, $f6, $f20
/* AB131C 8003A17C C4288F70 */  lwc1  $f8, %lo(D_80138F70)($at)
/* AB1320 8003A180 460E0402 */  mul.s $f16, $f0, $f14
/* AB1324 8003A184 4608803C */  c.lt.s $f16, $f8
/* AB1328 8003A188 00000000 */  nop   
/* AB132C 8003A18C 45020009 */  bc1fl .L8003A1B4
/* AB1330 8003A190 960E0002 */   lhu   $t6, 2($s0)
/* AB1334 8003A194 96220002 */  lhu   $v0, 2($s1)
/* AB1338 8003A198 52820079 */  beql  $s4, $v0, .L8003A380
/* AB133C 8003A19C C7A800FC */   lwc1  $f8, 0xfc($sp)
/* AB1340 8003A1A0 8E6A0048 */  lw    $t2, 0x48($s3)
/* AB1344 8003A1A4 00027880 */  sll   $t7, $v0, 2
/* AB1348 8003A1A8 10000072 */  b     .L8003A374
/* AB134C 8003A1AC 014F8821 */   addu  $s1, $t2, $t7
/* AB1350 8003A1B0 960E0002 */  lhu   $t6, 2($s0)
.L8003A1B4:
/* AB1354 8003A1B4 960C0004 */  lhu   $t4, 4($s0)
/* AB1358 8003A1B8 31C81FFF */  andi  $t0, $t6, 0x1fff
/* AB135C 8003A1BC 01150019 */  multu $t0, $s5
/* AB1360 8003A1C0 318D1FFF */  andi  $t5, $t4, 0x1fff
/* AB1364 8003A1C4 00004812 */  mflo  $t1
/* AB1368 8003A1C8 0249C021 */  addu  $t8, $s2, $t1
/* AB136C 8003A1CC 870B0000 */  lh    $t3, ($t8)
/* AB1370 8003A1D0 01B50019 */  multu $t5, $s5
/* AB1374 8003A1D4 448B5000 */  mtc1  $t3, $f10
/* AB1378 8003A1D8 00000000 */  nop   
/* AB137C 8003A1DC 46805320 */  cvt.s.w $f12, $f10
/* AB1380 8003A1E0 0000C812 */  mflo  $t9
/* AB1384 8003A1E4 02595021 */  addu  $t2, $s2, $t9
/* AB1388 8003A1E8 854F0000 */  lh    $t7, ($t2)
/* AB138C 8003A1EC 46006086 */  mov.s $f2, $f12
/* AB1390 8003A1F0 448F2000 */  mtc1  $t7, $f4
/* AB1394 8003A1F4 00000000 */  nop   
/* AB1398 8003A1F8 46802020 */  cvt.s.w $f0, $f4
/* AB139C 8003A1FC 4602003C */  c.lt.s $f0, $f2
/* AB13A0 8003A200 00000000 */  nop   
/* AB13A4 8003A204 45020004 */  bc1fl .L8003A218
/* AB13A8 8003A208 4600603C */   c.lt.s $f12, $f0
/* AB13AC 8003A20C 10000006 */  b     .L8003A228
/* AB13B0 8003A210 46000086 */   mov.s $f2, $f0
/* AB13B4 8003A214 4600603C */  c.lt.s $f12, $f0
.L8003A218:
/* AB13B8 8003A218 00000000 */  nop   
/* AB13BC 8003A21C 45020003 */  bc1fl .L8003A22C
/* AB13C0 8003A220 960E0006 */   lhu   $t6, 6($s0)
/* AB13C4 8003A224 46000306 */  mov.s $f12, $f0
.L8003A228:
/* AB13C8 8003A228 960E0006 */  lhu   $t6, 6($s0)
.L8003A22C:
/* AB13CC 8003A22C 01DE0019 */  multu $t6, $fp
/* AB13D0 8003A230 00004012 */  mflo  $t0
/* AB13D4 8003A234 02484821 */  addu  $t1, $s2, $t0
/* AB13D8 8003A238 85380000 */  lh    $t8, ($t1)
/* AB13DC 8003A23C 44983000 */  mtc1  $t8, $f6
/* AB13E0 8003A240 00000000 */  nop   
/* AB13E4 8003A244 46803020 */  cvt.s.w $f0, $f6
/* AB13E8 8003A248 4602003C */  c.lt.s $f0, $f2
/* AB13EC 8003A24C 00000000 */  nop   
/* AB13F0 8003A250 45020004 */  bc1fl .L8003A264
/* AB13F4 8003A254 4600603C */   c.lt.s $f12, $f0
/* AB13F8 8003A258 10000006 */  b     .L8003A274
/* AB13FC 8003A25C 46000086 */   mov.s $f2, $f0
/* AB1400 8003A260 4600603C */  c.lt.s $f12, $f0
.L8003A264:
/* AB1404 8003A264 00000000 */  nop   
/* AB1408 8003A268 45020003 */  bc1fl .L8003A278
/* AB140C 8003A26C 461C1081 */   sub.s $f2, $f2, $f28
/* AB1410 8003A270 46000306 */  mov.s $f12, $f0
.L8003A274:
/* AB1414 8003A274 461C1081 */  sub.s $f2, $f2, $f28
.L8003A278:
/* AB1418 8003A278 C7B200FC */  lwc1  $f18, 0xfc($sp)
/* AB141C 8003A27C 461C6300 */  add.s $f12, $f12, $f28
/* AB1420 8003A280 4602903C */  c.lt.s $f18, $f2
/* AB1424 8003A284 00000000 */  nop   
/* AB1428 8003A288 45030008 */  bc1tl .L8003A2AC
/* AB142C 8003A28C 96220002 */   lhu   $v0, 2($s1)
/* AB1430 8003A290 4612603C */  c.lt.s $f12, $f18
/* AB1434 8003A294 02002025 */  move  $a0, $s0
/* AB1438 8003A298 02402825 */  move  $a1, $s2
/* AB143C 8003A29C 8FA70104 */  lw    $a3, 0x104($sp)
/* AB1440 8003A2A0 45000008 */  bc1f  .L8003A2C4
/* AB1444 8003A2A4 27AD00EC */   addiu $t5, $sp, 0xec
/* AB1448 8003A2A8 96220002 */  lhu   $v0, 2($s1)
.L8003A2AC:
/* AB144C 8003A2AC 52820034 */  beql  $s4, $v0, .L8003A380
/* AB1450 8003A2B0 C7A800FC */   lwc1  $f8, 0xfc($sp)
/* AB1454 8003A2B4 8E6B0048 */  lw    $t3, 0x48($s3)
/* AB1458 8003A2B8 00026080 */  sll   $t4, $v0, 2
/* AB145C 8003A2BC 1000002D */  b     .L8003A374
/* AB1460 8003A2C0 016C8821 */   addu  $s1, $t3, $t4
.L8003A2C4:
/* AB1464 8003A2C4 8EC60004 */  lw    $a2, 4($s6)
/* AB1468 8003A2C8 E7B000B8 */  swc1  $f16, 0xb8($sp)
/* AB146C 8003A2CC E7AE009C */  swc1  $f14, 0x9c($sp)
/* AB1470 8003A2D0 0C00E3D8 */  jal   func_80038F60
/* AB1474 8003A2D4 AFAD0010 */   sw    $t5, 0x10($sp)
/* AB1478 8003A2D8 C7AE009C */  lwc1  $f14, 0x9c($sp)
/* AB147C 8003A2DC 1040001F */  beqz  $v0, .L8003A35C
/* AB1480 8003A2E0 C7B000B8 */   lwc1  $f16, 0xb8($sp)
/* AB1484 8003A2E4 4610E283 */  div.s $f10, $f28, $f16
/* AB1488 8003A2E8 C7B200FC */  lwc1  $f18, 0xfc($sp)
/* AB148C 8003A2EC C7A800EC */  lwc1  $f8, 0xec($sp)
/* AB1490 8003A2F0 46124081 */  sub.s $f2, $f8, $f18
/* AB1494 8003A2F4 46001005 */  abs.s $f0, $f2
/* AB1498 8003A2F8 460A003E */  c.le.s $f0, $f10
/* AB149C 8003A2FC 00000000 */  nop   
/* AB14A0 8003A300 45020017 */  bc1fl .L8003A360
/* AB14A4 8003A304 96220002 */   lhu   $v0, 2($s1)
/* AB14A8 8003A308 46161182 */  mul.s $f6, $f2, $f22
/* AB14AC 8003A30C 3C014080 */  li    $at, 0x40800000 # 0.000000
/* AB14B0 8003A310 44812000 */  mtc1  $at, $f4
/* AB14B4 8003A314 02602025 */  move  $a0, $s3
/* AB14B8 8003A318 02002825 */  move  $a1, $s0
/* AB14BC 8003A31C 02E03025 */  move  $a2, $s7
/* AB14C0 8003A320 27A70104 */  addiu $a3, $sp, 0x104
/* AB14C4 8003A324 4604303E */  c.le.s $f6, $f4
/* AB14C8 8003A328 240A0001 */  li    $t2, 1
/* AB14CC 8003A32C 4502000C */  bc1fl .L8003A360
/* AB14D0 8003A330 96220002 */   lhu   $v0, 2($s1)
/* AB14D4 8003A334 8FB90124 */  lw    $t9, 0x124($sp)
/* AB14D8 8003A338 E7B60010 */  swc1  $f22, 0x10($sp)
/* AB14DC 8003A33C E7BA0014 */  swc1  $f26, 0x14($sp)
/* AB14E0 8003A340 E7B80018 */  swc1  $f24, 0x18($sp)
/* AB14E4 8003A344 E7AE001C */  swc1  $f14, 0x1c($sp)
/* AB14E8 8003A348 E7BE0020 */  swc1  $f30, 0x20($sp)
/* AB14EC 8003A34C E7BC0024 */  swc1  $f28, 0x24($sp)
/* AB14F0 8003A350 AFAA00E8 */  sw    $t2, 0xe8($sp)
/* AB14F4 8003A354 0C00E68F */  jal   func_80039A3C
/* AB14F8 8003A358 AFB90028 */   sw    $t9, 0x28($sp)
.L8003A35C:
/* AB14FC 8003A35C 96220002 */  lhu   $v0, 2($s1)
.L8003A360:
/* AB1500 8003A360 52820007 */  beql  $s4, $v0, .L8003A380
/* AB1504 8003A364 C7A800FC */   lwc1  $f8, 0xfc($sp)
/* AB1508 8003A368 8E6F0048 */  lw    $t7, 0x48($s3)
/* AB150C 8003A36C 00027080 */  sll   $t6, $v0, 2
/* AB1510 8003A370 01EE8821 */  addu  $s1, $t7, $t6
.L8003A374:
/* AB1514 8003A374 1000FF17 */  b     .L80039FD4
/* AB1518 8003A378 C6C00004 */   lwc1  $f0, 4($s6)
/* AB151C 8003A37C C7A800FC */  lwc1  $f8, 0xfc($sp)
.L8003A380:
/* AB1520 8003A380 8FA80114 */  lw    $t0, 0x114($sp)
/* AB1524 8003A384 E5080000 */  swc1  $f8, ($t0)
/* AB1528 8003A388 8FA90118 */  lw    $t1, 0x118($sp)
/* AB152C 8003A38C C7AA0104 */  lwc1  $f10, 0x104($sp)
/* AB1530 8003A390 E52A0000 */  swc1  $f10, ($t1)
/* AB1534 8003A394 8FA200E8 */  lw    $v0, 0xe8($sp)
.L8003A398:
/* AB1538 8003A398 8FBF008C */  lw    $ra, 0x8c($sp)
/* AB153C 8003A39C D7B40038 */  ldc1  $f20, 0x38($sp)
/* AB1540 8003A3A0 D7B60040 */  ldc1  $f22, 0x40($sp)
/* AB1544 8003A3A4 D7B80048 */  ldc1  $f24, 0x48($sp)
/* AB1548 8003A3A8 D7BA0050 */  ldc1  $f26, 0x50($sp)
/* AB154C 8003A3AC D7BC0058 */  ldc1  $f28, 0x58($sp)
/* AB1550 8003A3B0 D7BE0060 */  ldc1  $f30, 0x60($sp)
/* AB1554 8003A3B4 8FB00068 */  lw    $s0, 0x68($sp)
/* AB1558 8003A3B8 8FB1006C */  lw    $s1, 0x6c($sp)
/* AB155C 8003A3BC 8FB20070 */  lw    $s2, 0x70($sp)
/* AB1560 8003A3C0 8FB30074 */  lw    $s3, 0x74($sp)
/* AB1564 8003A3C4 8FB40078 */  lw    $s4, 0x78($sp)
/* AB1568 8003A3C8 8FB5007C */  lw    $s5, 0x7c($sp)
/* AB156C 8003A3CC 8FB60080 */  lw    $s6, 0x80($sp)
/* AB1570 8003A3D0 8FB70084 */  lw    $s7, 0x84($sp)
/* AB1574 8003A3D4 8FBE0088 */  lw    $fp, 0x88($sp)
/* AB1578 8003A3D8 03E00008 */  jr    $ra
/* AB157C 8003A3DC 27BD0108 */   addiu $sp, $sp, 0x108

