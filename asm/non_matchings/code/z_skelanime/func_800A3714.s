glabel func_800A3714
/* B1A8B4 800A3714 27BDFFD8 */  addiu $sp, $sp, -0x28
/* B1A8B8 800A3718 AFBF0024 */  sw    $ra, 0x24($sp)
/* B1A8BC 800A371C AFB00020 */  sw    $s0, 0x20($sp)
/* B1A8C0 800A3720 AFA40028 */  sw    $a0, 0x28($sp)
/* B1A8C4 800A3724 3C0F8013 */  lui   $t7, %hi(D_8012A480) # $t7, 0x8013
/* B1A8C8 800A3728 8DEFA480 */  lw    $t7, %lo(D_8012A480)($t7)
/* B1A8CC 800A372C 90AE0000 */  lbu   $t6, ($a1)
/* B1A8D0 800A3730 00A08025 */  move  $s0, $a1
/* B1A8D4 800A3734 01CFC024 */  and   $t8, $t6, $t7
/* B1A8D8 800A3738 57000009 */  bnezl $t8, .L800A3760
/* B1A8DC 800A373C 8FBF0024 */   lw    $ra, 0x24($sp)
/* B1A8E0 800A3740 90A40001 */  lbu   $a0, 1($a1)
/* B1A8E4 800A3744 8CA50004 */  lw    $a1, 4($a1)
/* B1A8E8 800A3748 C604000C */  lwc1  $f4, 0xc($s0)
/* B1A8EC 800A374C 8E070008 */  lw    $a3, 8($s0)
/* B1A8F0 800A3750 00A03025 */  move  $a2, $a1
/* B1A8F4 800A3754 0C028B9C */  jal   func_800A2E70
/* B1A8F8 800A3758 E7A40010 */   swc1  $f4, 0x10($sp)
/* B1A8FC 800A375C 8FBF0024 */  lw    $ra, 0x24($sp)
.L800A3760:
/* B1A900 800A3760 8FB00020 */  lw    $s0, 0x20($sp)
/* B1A904 800A3764 27BD0028 */  addiu $sp, $sp, 0x28
/* B1A908 800A3768 03E00008 */  jr    $ra
/* B1A90C 800A376C 00000000 */   nop   

