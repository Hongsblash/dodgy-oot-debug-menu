.include "macro.inc"

# assembler directives
.set noat      # allow manual use of $at
.set noreorder # don't insert nops after branches
.set gp=64     # allow use of 64-bit general purposee registers

.section .text

.align 4

glabel func_80103DC0
/* B7AF60 80103DC0 27BDFFE8 */  addiu $sp, $sp, -0x18
/* B7AF64 80103DC4 AFBF0014 */  sw    $ra, 0x14($sp)
/* B7AF68 80103DC8 0C04191C */  jal   __osSpGetStatus
/* B7AF6C 80103DCC AFA40018 */   sw    $a0, 0x18($sp)
/* B7AF70 80103DD0 304E0100 */  andi  $t6, $v0, 0x100
/* B7AF74 80103DD4 11C00003 */  beqz  $t6, .L80103DE4
/* B7AF78 80103DD8 8FA40018 */   lw    $a0, 0x18($sp)
/* B7AF7C 80103DDC 10000002 */  b     .L80103DE8
/* B7AF80 80103DE0 24030001 */   li    $v1, 1
.L80103DE4:
/* B7AF84 80103DE4 00001825 */  move  $v1, $zero
.L80103DE8:
/* B7AF88 80103DE8 304F0080 */  andi  $t7, $v0, 0x80
/* B7AF8C 80103DEC 51E00008 */  beql  $t7, $zero, .L80103E10
/* B7AF90 80103DF0 8FBF0014 */   lw    $ra, 0x14($sp)
/* B7AF94 80103DF4 8C980004 */  lw    $t8, 4($a0)
/* B7AF98 80103DF8 2401FFFD */  li    $at, -3
/* B7AF9C 80103DFC 0303C825 */  or    $t9, $t8, $v1
/* B7AFA0 80103E00 AC990004 */  sw    $t9, 4($a0)
/* B7AFA4 80103E04 03214824 */  and   $t1, $t9, $at
/* B7AFA8 80103E08 AC890004 */  sw    $t1, 4($a0)
/* B7AFAC 80103E0C 8FBF0014 */  lw    $ra, 0x14($sp)
.L80103E10:
/* B7AFB0 80103E10 27BD0018 */  addiu $sp, $sp, 0x18
/* B7AFB4 80103E14 00601025 */  move  $v0, $v1
/* B7AFB8 80103E18 03E00008 */  jr    $ra
/* B7AFBC 80103E1C 00000000 */   nop   

glabel func_80103E20
/* B7AFC0 80103E20 27BDFFB8 */  addiu $sp, $sp, -0x48
/* B7AFC4 80103E24 AFB00018 */  sw    $s0, 0x18($sp)
/* B7AFC8 80103E28 00808025 */  move  $s0, $a0
/* B7AFCC 80103E2C AFBF001C */  sw    $ra, 0x1c($sp)
/* B7AFD0 80103E30 AFA5004C */  sw    $a1, 0x4c($sp)
/* B7AFD4 80103E34 AFA60050 */  sw    $a2, 0x50($sp)
/* B7AFD8 80103E38 AFA70054 */  sw    $a3, 0x54($sp)
/* B7AFDC 80103E3C 27A60058 */  addiu $a2, $sp, 0x58
/* B7AFE0 80103E40 27A50054 */  addiu $a1, $sp, 0x54
/* B7AFE4 80103E44 0C041058 */  jal   func_80104160
/* B7AFE8 80103E48 27A40050 */   addiu $a0, $sp, 0x50
/* B7AFEC 80103E4C 3C018013 */  lui   $at, %hi(D_80134D10)
/* B7AFF0 80103E50 C7AC004C */  lwc1  $f12, 0x4c($sp)
/* B7AFF4 80103E54 C4244D10 */  lwc1  $f4, %lo(D_80134D10)($at)
/* B7AFF8 80103E58 46046302 */  mul.s $f12, $f12, $f4
/* B7AFFC 80103E5C 0C0400A4 */  jal   sinf
/* B7B000 80103E60 E7AC004C */   swc1  $f12, 0x4c($sp)
/* B7B004 80103E64 C7AC004C */  lwc1  $f12, 0x4c($sp)
/* B7B008 80103E68 0C041184 */  jal   cosf
/* B7B00C 80103E6C E7A00044 */   swc1  $f0, 0x44($sp)
/* B7B010 80103E70 C7AC0050 */  lwc1  $f12, 0x50($sp)
/* B7B014 80103E74 C7A80054 */  lwc1  $f8, 0x54($sp)
/* B7B018 80103E78 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* B7B01C 80103E7C 44813000 */  mtc1  $at, $f6
/* B7B020 80103E80 46086282 */  mul.s $f10, $f12, $f8
/* B7B024 80103E84 C7A40058 */  lwc1  $f4, 0x58($sp)
/* B7B028 80103E88 46003081 */  sub.s $f2, $f6, $f0
/* B7B02C 80103E8C 02002025 */  move  $a0, $s0
/* B7B030 80103E90 E7A00040 */  swc1  $f0, 0x40($sp)
/* B7B034 80103E94 46025482 */  mul.s $f18, $f10, $f2
/* B7B038 80103E98 00000000 */  nop   
/* B7B03C 80103E9C 46044182 */  mul.s $f6, $f8, $f4
/* B7B040 80103EA0 E7B2003C */  swc1  $f18, 0x3c($sp)
/* B7B044 80103EA4 46023282 */  mul.s $f10, $f6, $f2
/* B7B048 80103EA8 00000000 */  nop   
/* B7B04C 80103EAC 460C2202 */  mul.s $f8, $f4, $f12
/* B7B050 80103EB0 E7AA0038 */  swc1  $f10, 0x38($sp)
/* B7B054 80103EB4 46024182 */  mul.s $f6, $f8, $f2
/* B7B058 80103EB8 0C0406D0 */  jal   func_80101B40
/* B7B05C 80103EBC E7A60034 */   swc1  $f6, 0x34($sp)
/* B7B060 80103EC0 C7AE0044 */  lwc1  $f14, 0x44($sp)
/* B7B064 80103EC4 C7AA0050 */  lwc1  $f10, 0x50($sp)
/* B7B068 80103EC8 C7A40054 */  lwc1  $f4, 0x54($sp)
/* B7B06C 80103ECC C7A60058 */  lwc1  $f6, 0x58($sp)
/* B7B070 80103ED0 460E5002 */  mul.s $f0, $f10, $f14
/* B7B074 80103ED4 3C013F80 */  li    $at, 0x3F800000 # 0.000000
/* B7B078 80103ED8 C7B00040 */  lwc1  $f16, 0x40($sp)
/* B7B07C 80103EDC 460E2202 */  mul.s $f8, $f4, $f14
/* B7B080 80103EE0 C7B2003C */  lwc1  $f18, 0x3c($sp)
/* B7B084 80103EE4 460E3102 */  mul.s $f4, $f6, $f14
/* B7B088 80103EE8 00000000 */  nop   
/* B7B08C 80103EEC 460A5082 */  mul.s $f2, $f10, $f10
/* B7B090 80103EF0 E7A80028 */  swc1  $f8, 0x28($sp)
/* B7B094 80103EF4 44814000 */  mtc1  $at, $f8
/* B7B098 80103EF8 E7A40024 */  swc1  $f4, 0x24($sp)
/* B7B09C 80103EFC 46024181 */  sub.s $f6, $f8, $f2
/* B7B0A0 80103F00 46103102 */  mul.s $f4, $f6, $f16
/* B7B0A4 80103F04 46022280 */  add.s $f10, $f4, $f2
/* B7B0A8 80103F08 E60A0000 */  swc1  $f10, ($s0)
/* B7B0AC 80103F0C C7A80038 */  lwc1  $f8, 0x38($sp)
/* B7B0B0 80103F10 46004181 */  sub.s $f6, $f8, $f0
/* B7B0B4 80103F14 E6060024 */  swc1  $f6, 0x24($s0)
/* B7B0B8 80103F18 C7A40038 */  lwc1  $f4, 0x38($sp)
/* B7B0BC 80103F1C 44813000 */  mtc1  $at, $f6
/* B7B0C0 80103F20 46002280 */  add.s $f10, $f4, $f0
/* B7B0C4 80103F24 E60A0018 */  swc1  $f10, 0x18($s0)
/* B7B0C8 80103F28 C7A80054 */  lwc1  $f8, 0x54($sp)
/* B7B0CC 80103F2C 46084302 */  mul.s $f12, $f8, $f8
/* B7B0D0 80103F30 460C3101 */  sub.s $f4, $f6, $f12
/* B7B0D4 80103F34 46102282 */  mul.s $f10, $f4, $f16
/* B7B0D8 80103F38 460C5200 */  add.s $f8, $f10, $f12
/* B7B0DC 80103F3C E6080014 */  swc1  $f8, 0x14($s0)
/* B7B0E0 80103F40 C7A20034 */  lwc1  $f2, 0x34($sp)
/* B7B0E4 80103F44 C7A60028 */  lwc1  $f6, 0x28($sp)
/* B7B0E8 80103F48 46061100 */  add.s $f4, $f2, $f6
/* B7B0EC 80103F4C E6040020 */  swc1  $f4, 0x20($s0)
/* B7B0F0 80103F50 C7AA0028 */  lwc1  $f10, 0x28($sp)
/* B7B0F4 80103F54 44812000 */  mtc1  $at, $f4
/* B7B0F8 80103F58 460A1201 */  sub.s $f8, $f2, $f10
/* B7B0FC 80103F5C E6080008 */  swc1  $f8, 8($s0)
/* B7B100 80103F60 C7A60058 */  lwc1  $f6, 0x58($sp)
/* B7B104 80103F64 46063002 */  mul.s $f0, $f6, $f6
/* B7B108 80103F68 46002281 */  sub.s $f10, $f4, $f0
/* B7B10C 80103F6C 46105202 */  mul.s $f8, $f10, $f16
/* B7B110 80103F70 46004180 */  add.s $f6, $f8, $f0
/* B7B114 80103F74 E6060028 */  swc1  $f6, 0x28($s0)
/* B7B118 80103F78 C7A40024 */  lwc1  $f4, 0x24($sp)
/* B7B11C 80103F7C 46049281 */  sub.s $f10, $f18, $f4
/* B7B120 80103F80 E60A0010 */  swc1  $f10, 0x10($s0)
/* B7B124 80103F84 C7A80024 */  lwc1  $f8, 0x24($sp)
/* B7B128 80103F88 46089180 */  add.s $f6, $f18, $f8
/* B7B12C 80103F8C E6060004 */  swc1  $f6, 4($s0)
/* B7B130 80103F90 8FBF001C */  lw    $ra, 0x1c($sp)
/* B7B134 80103F94 8FB00018 */  lw    $s0, 0x18($sp)
/* B7B138 80103F98 27BD0048 */  addiu $sp, $sp, 0x48
/* B7B13C 80103F9C 03E00008 */  jr    $ra
/* B7B140 80103FA0 00000000 */   nop   

glabel func_80103FA4
/* B7B144 80103FA4 27BDFFA0 */  addiu $sp, $sp, -0x60
/* B7B148 80103FA8 44856000 */  mtc1  $a1, $f12
/* B7B14C 80103FAC 44867000 */  mtc1  $a2, $f14
/* B7B150 80103FB0 C7A40070 */  lwc1  $f4, 0x70($sp)
/* B7B154 80103FB4 AFBF001C */  sw    $ra, 0x1c($sp)
/* B7B158 80103FB8 AFA40060 */  sw    $a0, 0x60($sp)
/* B7B15C 80103FBC 44056000 */  mfc1  $a1, $f12
/* B7B160 80103FC0 44067000 */  mfc1  $a2, $f14
/* B7B164 80103FC4 AFA7006C */  sw    $a3, 0x6c($sp)
/* B7B168 80103FC8 27A40020 */  addiu $a0, $sp, 0x20
/* B7B16C 80103FCC 0C040F88 */  jal   func_80103E20
/* B7B170 80103FD0 E7A40010 */   swc1  $f4, 0x10($sp)
/* B7B174 80103FD4 27A40020 */  addiu $a0, $sp, 0x20
/* B7B178 80103FD8 0C041938 */  jal   func_801064E0
/* B7B17C 80103FDC 8FA50060 */   lw    $a1, 0x60($sp)
/* B7B180 80103FE0 8FBF001C */  lw    $ra, 0x1c($sp)
/* B7B184 80103FE4 27BD0060 */  addiu $sp, $sp, 0x60
/* B7B188 80103FE8 03E00008 */  jr    $ra
/* B7B18C 80103FEC 00000000 */   nop   

glabel func_80103FF0
/* B7B190 80103FF0 3C058001 */  lui   $a1, %hi(osViClock)
/* B7B194 80103FF4 24A5ACF8 */  addiu $a1, %lo(osViClock) # addiu $a1, $a1, -0x5308
/* B7B198 80103FF8 8CAE0000 */  lw    $t6, ($a1)
/* B7B19C 80103FFC 44844000 */  mtc1  $a0, $f8
/* B7B1A0 80104000 3C014F80 */  li    $at, 0x4F800000 # 0.000000
/* B7B1A4 80104004 448E2000 */  mtc1  $t6, $f4
/* B7B1A8 80104008 468042A0 */  cvt.s.w $f10, $f8
/* B7B1AC 8010400C 04810004 */  bgez  $a0, .L80104020
/* B7B1B0 80104010 468021A0 */   cvt.s.w $f6, $f4
/* B7B1B4 80104014 44818000 */  mtc1  $at, $f16
/* B7B1B8 80104018 00000000 */  nop   
/* B7B1BC 8010401C 46105280 */  add.s $f10, $f10, $f16
.L80104020:
/* B7B1C0 80104020 460A3483 */  div.s $f18, $f6, $f10
/* B7B1C4 80104024 3C013F00 */  li    $at, 0x3F000000 # 0.000000
/* B7B1C8 80104028 44812000 */  mtc1  $at, $f4
/* B7B1CC 8010402C 24030001 */  li    $v1, 1
/* B7B1D0 80104030 3C014F00 */  lui   $at, 0x4f00
/* B7B1D4 80104034 3C08A450 */  lui   $t0, 0xa450
/* B7B1D8 80104038 3C0AA450 */  li    $t2, 0xA4500000 # 0.000000
/* B7B1DC 8010403C 46049000 */  add.s $f0, $f18, $f4
/* B7B1E0 80104040 444FF800 */  cfc1  $t7, $31
/* B7B1E4 80104044 44C3F800 */  ctc1  $v1, $31
/* B7B1E8 80104048 00000000 */  nop   
/* B7B1EC 8010404C 46000224 */  cvt.w.s $f8, $f0
/* B7B1F0 80104050 4443F800 */  cfc1  $v1, $31
/* B7B1F4 80104054 00000000 */  nop   
/* B7B1F8 80104058 30630078 */  andi  $v1, $v1, 0x78
/* B7B1FC 8010405C 50600013 */  beql  $v1, $zero, .L801040AC
/* B7B200 80104060 44034000 */   mfc1  $v1, $f8
/* B7B204 80104064 44814000 */  mtc1  $at, $f8
/* B7B208 80104068 24030001 */  li    $v1, 1
/* B7B20C 8010406C 46080201 */  sub.s $f8, $f0, $f8
/* B7B210 80104070 44C3F800 */  ctc1  $v1, $31
/* B7B214 80104074 00000000 */  nop   
/* B7B218 80104078 46004224 */  cvt.w.s $f8, $f8
/* B7B21C 8010407C 4443F800 */  cfc1  $v1, $31
/* B7B220 80104080 00000000 */  nop   
/* B7B224 80104084 30630078 */  andi  $v1, $v1, 0x78
/* B7B228 80104088 14600005 */  bnez  $v1, .L801040A0
/* B7B22C 8010408C 00000000 */   nop   
/* B7B230 80104090 44034000 */  mfc1  $v1, $f8
/* B7B234 80104094 3C018000 */  lui   $at, 0x8000
/* B7B238 80104098 10000007 */  b     .L801040B8
/* B7B23C 8010409C 00611825 */   or    $v1, $v1, $at
.L801040A0:
/* B7B240 801040A0 10000005 */  b     .L801040B8
/* B7B244 801040A4 2403FFFF */   li    $v1, -1
/* B7B248 801040A8 44034000 */  mfc1  $v1, $f8
.L801040AC:
/* B7B24C 801040AC 00000000 */  nop   
/* B7B250 801040B0 0460FFFB */  bltz  $v1, .L801040A0
/* B7B254 801040B4 00000000 */   nop   
.L801040B8:
/* B7B258 801040B8 44CFF800 */  ctc1  $t7, $31
/* B7B25C 801040BC 2C610084 */  sltiu $at, $v1, 0x84
/* B7B260 801040C0 10200003 */  beqz  $at, .L801040D0
/* B7B264 801040C4 2479FFFF */   addiu $t9, $v1, -1
/* B7B268 801040C8 03E00008 */  jr    $ra
/* B7B26C 801040CC 2402FFFF */   li    $v0, -1

.L801040D0:
/* B7B270 801040D0 24010042 */  li    $at, 66
/* B7B274 801040D4 0061001B */  divu  $zero, $v1, $at
/* B7B278 801040D8 00001012 */  mflo  $v0
/* B7B27C 801040DC 305800FF */  andi  $t8, $v0, 0xff
/* B7B280 801040E0 2B010011 */  slti  $at, $t8, 0x11
/* B7B284 801040E4 14200002 */  bnez  $at, .L801040F0
/* B7B288 801040E8 304400FF */   andi  $a0, $v0, 0xff
/* B7B28C 801040EC 24040010 */  li    $a0, 16
.L801040F0:
/* B7B290 801040F0 AD190010 */  sw    $t9, 0x10($t0)
/* B7B294 801040F4 2489FFFF */  addiu $t1, $a0, -1
/* B7B298 801040F8 AD490014 */  sw    $t1, 0x14($t2)
/* B7B29C 801040FC 8CAB0000 */  lw    $t3, ($a1)
/* B7B2A0 80104100 0163001A */  div   $zero, $t3, $v1
/* B7B2A4 80104104 00001012 */  mflo  $v0
/* B7B2A8 80104108 14600002 */  bnez  $v1, .L80104114
/* B7B2AC 8010410C 00000000 */   nop   
/* B7B2B0 80104110 0007000D */  break 7
.L80104114:
/* B7B2B4 80104114 2401FFFF */  li    $at, -1
/* B7B2B8 80104118 14610004 */  bne   $v1, $at, .L8010412C
/* B7B2BC 8010411C 3C018000 */   lui   $at, 0x8000
/* B7B2C0 80104120 15610002 */  bne   $t3, $at, .L8010412C
/* B7B2C4 80104124 00000000 */   nop   
/* B7B2C8 80104128 0006000D */  break 6
.L8010412C:
/* B7B2CC 8010412C 03E00008 */  jr    $ra
/* B7B2D0 80104130 00000000 */   nop   
