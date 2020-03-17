.include "macro.inc"

# assembler directives
.set noat      # allow manual use of $at
.set noreorder # don't insert nops after branches
.set gp=64     # allow use of 64-bit general purposee registers

.section .text

.align 4

glabel func_80106170
/* B7D310 80106170 AFA40000 */  sw    $a0, ($sp)
/* B7D314 80106174 308EFFFF */  andi  $t6, $a0, 0xffff
/* B7D318 80106178 01C02025 */  move  $a0, $t6
/* B7D31C 8010617C 00001825 */  move  $v1, $zero
/* B7D320 80106180 24020400 */  li    $v0, 1024
.L80106184:
/* B7D324 80106184 00037840 */  sll   $t7, $v1, 1
/* B7D328 80106188 0082C024 */  and   $t8, $a0, $v0
/* B7D32C 8010618C 13000008 */  beqz  $t8, .L801061B0
/* B7D330 80106190 01E01825 */   move  $v1, $t7
/* B7D334 80106194 31F90020 */  andi  $t9, $t7, 0x20
/* B7D338 80106198 13200003 */  beqz  $t9, .L801061A8
/* B7D33C 8010619C 00000000 */   nop   
/* B7D340 801061A0 10000007 */  b     .L801061C0
/* B7D344 801061A4 39E30014 */   xori  $v1, $t7, 0x14
.L801061A8:
/* B7D348 801061A8 10000005 */  b     .L801061C0
/* B7D34C 801061AC 24630001 */   addiu $v1, $v1, 1
.L801061B0:
/* B7D350 801061B0 30690020 */  andi  $t1, $v1, 0x20
/* B7D354 801061B4 11200002 */  beqz  $t1, .L801061C0
/* B7D358 801061B8 386A0015 */   xori  $t2, $v1, 0x15
/* B7D35C 801061BC 01401825 */  move  $v1, $t2
.L801061C0:
/* B7D360 801061C0 00025842 */  srl   $t3, $v0, 1
/* B7D364 801061C4 1560FFEF */  bnez  $t3, .L80106184
/* B7D368 801061C8 01601025 */   move  $v0, $t3
/* B7D36C 801061CC 00036040 */  sll   $t4, $v1, 1
/* B7D370 801061D0 318D0020 */  andi  $t5, $t4, 0x20
/* B7D374 801061D4 11A00002 */  beqz  $t5, .L801061E0
/* B7D378 801061D8 01801825 */   move  $v1, $t4
/* B7D37C 801061DC 39830015 */  xori  $v1, $t4, 0x15
.L801061E0:
/* B7D380 801061E0 00037840 */  sll   $t7, $v1, 1
/* B7D384 801061E4 31F80020 */  andi  $t8, $t7, 0x20
/* B7D388 801061E8 13000002 */  beqz  $t8, .L801061F4
/* B7D38C 801061EC 01E01825 */   move  $v1, $t7
/* B7D390 801061F0 39E30015 */  xori  $v1, $t7, 0x15
.L801061F4:
/* B7D394 801061F4 00034040 */  sll   $t0, $v1, 1
/* B7D398 801061F8 31090020 */  andi  $t1, $t0, 0x20
/* B7D39C 801061FC 11200002 */  beqz  $t1, .L80106208
/* B7D3A0 80106200 01001825 */   move  $v1, $t0
/* B7D3A4 80106204 39030015 */  xori  $v1, $t0, 0x15
.L80106208:
/* B7D3A8 80106208 00035840 */  sll   $t3, $v1, 1
/* B7D3AC 8010620C 316C0020 */  andi  $t4, $t3, 0x20
/* B7D3B0 80106210 11800002 */  beqz  $t4, .L8010621C
/* B7D3B4 80106214 01601825 */   move  $v1, $t3
/* B7D3B8 80106218 39630015 */  xori  $v1, $t3, 0x15
.L8010621C:
/* B7D3BC 8010621C 00037040 */  sll   $t6, $v1, 1
/* B7D3C0 80106220 31CF0020 */  andi  $t7, $t6, 0x20
/* B7D3C4 80106224 11E00002 */  beqz  $t7, .L80106230
/* B7D3C8 80106228 01C01825 */   move  $v1, $t6
/* B7D3CC 8010622C 39C30015 */  xori  $v1, $t6, 0x15
.L80106230:
/* B7D3D0 80106230 00601025 */  move  $v0, $v1
/* B7D3D4 80106234 3059001F */  andi  $t9, $v0, 0x1f
/* B7D3D8 80106238 03E00008 */  jr    $ra
/* B7D3DC 8010623C 03201025 */   move  $v0, $t9

glabel func_80106240
/* B7D3E0 80106240 00802825 */  move  $a1, $a0
/* B7D3E4 80106244 00001825 */  move  $v1, $zero
/* B7D3E8 80106248 24020020 */  li    $v0, 32
.L8010624C:
/* B7D3EC 8010624C 24040080 */  li    $a0, 128
/* B7D3F0 80106250 90A60000 */  lbu   $a2, ($a1)
.L80106254:
/* B7D3F4 80106254 00037040 */  sll   $t6, $v1, 1
/* B7D3F8 80106258 00C47824 */  and   $t7, $a2, $a0
/* B7D3FC 8010625C 11E00008 */  beqz  $t7, .L80106280
/* B7D400 80106260 01C01825 */   move  $v1, $t6
/* B7D404 80106264 31D80100 */  andi  $t8, $t6, 0x100
/* B7D408 80106268 13000003 */  beqz  $t8, .L80106278
/* B7D40C 8010626C 00000000 */   nop   
/* B7D410 80106270 10000007 */  b     .L80106290
/* B7D414 80106274 39C30084 */   xori  $v1, $t6, 0x84
.L80106278:
/* B7D418 80106278 10000005 */  b     .L80106290
/* B7D41C 8010627C 24630001 */   addiu $v1, $v1, 1
.L80106280:
/* B7D420 80106280 30680100 */  andi  $t0, $v1, 0x100
/* B7D424 80106284 11000002 */  beqz  $t0, .L80106290
/* B7D428 80106288 38690085 */   xori  $t1, $v1, 0x85
/* B7D42C 8010628C 01201825 */  move  $v1, $t1
.L80106290:
/* B7D430 80106290 00045042 */  srl   $t2, $a0, 1
/* B7D434 80106294 1540FFEF */  bnez  $t2, .L80106254
/* B7D438 80106298 01402025 */   move  $a0, $t2
/* B7D43C 8010629C 2442FFFF */  addiu $v0, $v0, -1
/* B7D440 801062A0 1440FFEA */  bnez  $v0, .L8010624C
/* B7D444 801062A4 24A50001 */   addiu $a1, $a1, 1
/* B7D448 801062A8 00035840 */  sll   $t3, $v1, 1
.L801062AC:
/* B7D44C 801062AC 316C0100 */  andi  $t4, $t3, 0x100
/* B7D450 801062B0 11800002 */  beqz  $t4, .L801062BC
/* B7D454 801062B4 01601825 */   move  $v1, $t3
/* B7D458 801062B8 39630085 */  xori  $v1, $t3, 0x85
.L801062BC:
/* B7D45C 801062BC 24420001 */  addiu $v0, $v0, 1
/* B7D460 801062C0 2C410008 */  sltiu $at, $v0, 8
/* B7D464 801062C4 5420FFF9 */  bnezl $at, .L801062AC
/* B7D468 801062C8 00035840 */   sll   $t3, $v1, 1
/* B7D46C 801062CC 03E00008 */  jr    $ra
/* B7D470 801062D0 306200FF */   andi  $v0, $v1, 0xff
