glabel func_8008E9F8
/* B05B98 8008E9F8 3C038012 */  lui   $v1, %hi(D_80125C44)
/* B05B9C 8008E9FC 00651821 */  addu  $v1, $v1, $a1
/* B05BA0 8008EA00 90635C44 */  lbu   $v1, %lo(D_80125C44)($v1)
/* B05BA4 8008EA04 27BDFFE0 */  addiu $sp, $sp, -0x20
/* B05BA8 8008EA08 24010002 */  li    $at, 2
/* B05BAC 8008EA0C 14610007 */  bne   $v1, $at, .L8008EA2C
/* B05BB0 8008EA10 AFBF0014 */   sw    $ra, 0x14($sp)
/* B05BB4 8008EA14 0C023A74 */  jal   func_8008E9D0
/* B05BB8 8008EA18 AFA3001C */   sw    $v1, 0x1c($sp)
/* B05BBC 8008EA1C 10400003 */  beqz  $v0, .L8008EA2C
/* B05BC0 8008EA20 8FA3001C */   lw    $v1, 0x1c($sp)
/* B05BC4 8008EA24 10000002 */  b     .L8008EA30
/* B05BC8 8008EA28 24020001 */   li    $v0, 1
.L8008EA2C:
/* B05BCC 8008EA2C 00601025 */  move  $v0, $v1
.L8008EA30:
/* B05BD0 8008EA30 8FBF0014 */  lw    $ra, 0x14($sp)
/* B05BD4 8008EA34 27BD0020 */  addiu $sp, $sp, 0x20
/* B05BD8 8008EA38 03E00008 */  jr    $ra
/* B05BDC 8008EA3C 00000000 */   nop   

