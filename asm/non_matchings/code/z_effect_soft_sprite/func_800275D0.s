.rdata
glabel D_801357DC
    .asciz "../z_effect_soft_sprite.c"
    .balign 4

.text
glabel func_800275D0
/* A9E770 800275D0 3C068011 */  lui   $a2, %hi(EffectSS2Info) # $a2, 0x8011
/* A9E774 800275D4 24C658B0 */  addiu $a2, %lo(EffectSS2Info) # addiu $a2, $a2, 0x58b0
/* A9E778 800275D8 8CC20004 */  lw    $v0, 4($a2)
/* A9E77C 800275DC 8CC30008 */  lw    $v1, 8($a2)
/* A9E780 800275E0 27BDFFF8 */  addiu $sp, $sp, -8
/* A9E784 800275E4 AFB00004 */  sw    $s0, 4($sp)
/* A9E788 800275E8 0043082A */  slt   $at, $v0, $v1
/* A9E78C 800275EC 00808025 */  move  $s0, $a0
/* A9E790 800275F0 14200003 */  bnez  $at, .L80027600
/* A9E794 800275F4 AFA5000C */   sw    $a1, 0xc($sp)
/* A9E798 800275F8 ACC00004 */  sw    $zero, 4($a2)
/* A9E79C 800275FC 00001025 */  move  $v0, $zero
.L80027600:
/* A9E7A0 80027600 3C068011 */  lui   $a2, %hi(EffectSS2Info) # $a2, 0x8011
/* A9E7A4 80027604 8CC658B0 */  lw    $a2, %lo(EffectSS2Info)($a2)
/* A9E7A8 80027608 00027080 */  sll   $t6, $v0, 2
/* A9E7AC 8002760C 01C27023 */  subu  $t6, $t6, $v0
/* A9E7B0 80027610 000E7140 */  sll   $t6, $t6, 5
/* A9E7B4 80027614 00402025 */  move  $a0, $v0
/* A9E7B8 80027618 00002825 */  move  $a1, $zero
/* A9E7BC 8002761C 2408FFFF */  li    $t0, -1
/* A9E7C0 80027620 00CE3821 */  addu  $a3, $a2, $t6
.L80027624:
/* A9E7C4 80027624 84EF005C */  lh    $t7, 0x5c($a3)
/* A9E7C8 80027628 550F0004 */  bnel  $t0, $t7, .L8002763C
/* A9E7CC 8002762C 24840001 */   addiu $a0, $a0, 1
/* A9E7D0 80027630 1000000C */  b     .L80027664
/* A9E7D4 80027634 24050001 */   li    $a1, 1
/* A9E7D8 80027638 24840001 */  addiu $a0, $a0, 1
.L8002763C:
/* A9E7DC 8002763C 0083082A */  slt   $at, $a0, $v1
/* A9E7E0 80027640 14200002 */  bnez  $at, .L8002764C
/* A9E7E4 80027644 00000000 */   nop   
/* A9E7E8 80027648 00002025 */  move  $a0, $zero
.L8002764C:
/* A9E7EC 8002764C 10820005 */  beq   $a0, $v0, .L80027664
/* A9E7F0 80027650 0004C080 */   sll   $t8, $a0, 2
/* A9E7F4 80027654 0304C023 */  subu  $t8, $t8, $a0
/* A9E7F8 80027658 0018C140 */  sll   $t8, $t8, 5
/* A9E7FC 8002765C 1000FFF1 */  b     .L80027624
/* A9E800 80027660 00D83821 */   addu  $a3, $a2, $t8
.L80027664:
/* A9E804 80027664 24010001 */  li    $at, 1
/* A9E808 80027668 14A10005 */  bne   $a1, $at, .L80027680
/* A9E80C 8002766C 00024880 */   sll   $t1, $v0, 2
/* A9E810 80027670 8FB9000C */  lw    $t9, 0xc($sp)
/* A9E814 80027674 00001025 */  move  $v0, $zero
/* A9E818 80027678 1000001F */  b     .L800276F8
/* A9E81C 8002767C AF240000 */   sw    $a0, ($t9)
.L80027680:
/* A9E820 80027680 01224823 */  subu  $t1, $t1, $v0
/* A9E824 80027684 00094940 */  sll   $t1, $t1, 5
/* A9E828 80027688 00402025 */  move  $a0, $v0
/* A9E82C 8002768C 00C93821 */  addu  $a3, $a2, $t1
.L80027690:
/* A9E830 80027690 90E5005E */  lbu   $a1, 0x5e($a3)
/* A9E834 80027694 00B0082A */  slt   $at, $a1, $s0
/* A9E838 80027698 54200008 */  bnezl $at, .L800276BC
/* A9E83C 8002769C 24840001 */   addiu $a0, $a0, 1
/* A9E840 800276A0 56050013 */  bnel  $s0, $a1, .L800276F0
/* A9E844 800276A4 8FAD000C */   lw    $t5, 0xc($sp)
/* A9E848 800276A8 94EA005A */  lhu   $t2, 0x5a($a3)
/* A9E84C 800276AC 314B0001 */  andi  $t3, $t2, 1
/* A9E850 800276B0 5160000F */  beql  $t3, $zero, .L800276F0
/* A9E854 800276B4 8FAD000C */   lw    $t5, 0xc($sp)
/* A9E858 800276B8 24840001 */  addiu $a0, $a0, 1
.L800276BC:
/* A9E85C 800276BC 0083082A */  slt   $at, $a0, $v1
/* A9E860 800276C0 14200002 */  bnez  $at, .L800276CC
/* A9E864 800276C4 00000000 */   nop   
/* A9E868 800276C8 00002025 */  move  $a0, $zero
.L800276CC:
/* A9E86C 800276CC 14820003 */  bne   $a0, $v0, .L800276DC
/* A9E870 800276D0 00046080 */   sll   $t4, $a0, 2
/* A9E874 800276D4 10000008 */  b     .L800276F8
/* A9E878 800276D8 24020001 */   li    $v0, 1
.L800276DC:
/* A9E87C 800276DC 01846023 */  subu  $t4, $t4, $a0
/* A9E880 800276E0 000C6140 */  sll   $t4, $t4, 5
/* A9E884 800276E4 1000FFEA */  b     .L80027690
/* A9E888 800276E8 00CC3821 */   addu  $a3, $a2, $t4
/* A9E88C 800276EC 8FAD000C */  lw    $t5, 0xc($sp)
.L800276F0:
/* A9E890 800276F0 00001025 */  move  $v0, $zero
/* A9E894 800276F4 ADA40000 */  sw    $a0, ($t5)
.L800276F8:
/* A9E898 800276F8 8FB00004 */  lw    $s0, 4($sp)
/* A9E89C 800276FC 03E00008 */  jr    $ra
/* A9E8A0 80027700 27BD0008 */   addiu $sp, $sp, 8

