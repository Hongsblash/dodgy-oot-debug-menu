glabel Fault_WaitForInputImpl
/* B4B8A0 800D4700 3C09FFFA */  lui   $t1, (0xFFFA5A5A >> 16) # lui $t1, 0xfffa
/* B4B8A4 800D4704 03A05025 */  move  $t2, $sp
/* B4B8A8 800D4708 27BDFFB8 */  addiu $sp, $sp, -0x48
/* B4B8AC 800D470C 35295A5A */  ori   $t1, (0xFFFA5A5A & 0xFFFF) # ori $t1, $t1, 0x5a5a
.L800D4710:
/* B4B8B0 800D4710 254AFFF8 */  addiu $t2, $t2, -8
/* B4B8B4 800D4714 AD490000 */  sw    $t1, ($t2)
/* B4B8B8 800D4718 155DFFFD */  bne   $t2, $sp, .L800D4710
/* B4B8BC 800D471C AD490004 */   sw    $t1, 4($t2)
/* B4B8C0 800D4720 AFB20020 */  sw    $s2, 0x20($sp)
/* B4B8C4 800D4724 3C128017 */  lui   $s2, %hi(sFaultStructPtr) # $s2, 0x8017
/* B4B8C8 800D4728 2652A800 */  addiu $s2, %lo(sFaultStructPtr) # addiu $s2, $s2, -0x5800
/* B4B8CC 800D472C AFBF003C */  sw    $ra, 0x3c($sp)
/* B4B8D0 800D4730 AFBE0038 */  sw    $fp, 0x38($sp)
/* B4B8D4 800D4734 AFB70034 */  sw    $s7, 0x34($sp)
/* B4B8D8 800D4738 AFB60030 */  sw    $s6, 0x30($sp)
/* B4B8DC 800D473C AFB5002C */  sw    $s5, 0x2c($sp)
/* B4B8E0 800D4740 AFB40028 */  sw    $s4, 0x28($sp)
/* B4B8E4 800D4744 AFB30024 */  sw    $s3, 0x24($sp)
/* B4B8E8 800D4748 AFB1001C */  sw    $s1, 0x1c($sp)
/* B4B8EC 800D474C AFB00018 */  sw    $s0, 0x18($sp)
/* B4B8F0 800D4750 8E4E0000 */  lw    $t6, ($s2)
/* B4B8F4 800D4754 24110258 */  li    $s1, 600
/* B4B8F8 800D4758 241E0400 */  li    $fp, 1024
/* B4B8FC 800D475C 25CF07E4 */  addiu $t7, $t6, 0x7e4
/* B4B900 800D4760 AFAF0044 */  sw    $t7, 0x44($sp)
/* B4B904 800D4764 24170800 */  li    $s7, 2048
/* B4B908 800D4768 24160200 */  li    $s6, 512
/* B4B90C 800D476C 24150100 */  li    $s5, 256
/* B4B910 800D4770 34148000 */  li    $s4, 32768
/* B4B914 800D4774 24130020 */  li    $s3, 32
.L800D4778:
/* B4B918 800D4778 0C03518F */  jal   Fault_Sleep
/* B4B91C 800D477C 24040010 */   li    $a0, 16
/* B4B920 800D4780 0C0351AD */  jal   Fault_UpdatePadImpl
/* B4B924 800D4784 00000000 */   nop   
/* B4B928 800D4788 8FB80044 */  lw    $t8, 0x44($sp)
/* B4B92C 800D478C 9710000C */  lhu   $s0, 0xc($t8)
/* B4B930 800D4790 56130006 */  bnel  $s0, $s3, .L800D47AC
/* B4B934 800D4794 8E590000 */   lw    $t9, ($s2)
/* B4B938 800D4798 8E420000 */  lw    $v0, ($s2)
/* B4B93C 800D479C 904307CF */  lbu   $v1, 0x7cf($v0)
/* B4B940 800D47A0 2C630001 */  sltiu $v1, $v1, 1
/* B4B944 800D47A4 A04307CF */  sb    $v1, 0x7cf($v0)
/* B4B948 800D47A8 8E590000 */  lw    $t9, ($s2)
.L800D47AC:
/* B4B94C 800D47AC 2A230001 */  slti  $v1, $s1, 1
/* B4B950 800D47B0 932807CF */  lbu   $t0, 0x7cf($t9)
/* B4B954 800D47B4 11000005 */  beqz  $t0, .L800D47CC
/* B4B958 800D47B8 00000000 */   nop   
/* B4B95C 800D47BC 1060FFEE */  beqz  $v1, .L800D4778
/* B4B960 800D47C0 2631FFFF */   addiu $s1, $s1, -1
/* B4B964 800D47C4 10000015 */  b     .L800D481C
/* B4B968 800D47C8 00001025 */   move  $v0, $zero
.L800D47CC:
/* B4B96C 800D47CC 12140003 */  beq   $s0, $s4, .L800D47DC
/* B4B970 800D47D0 00000000 */   nop   
/* B4B974 800D47D4 16150003 */  bne   $s0, $s5, .L800D47E4
/* B4B978 800D47D8 00000000 */   nop   
.L800D47DC:
/* B4B97C 800D47DC 1000000F */  b     .L800D481C
/* B4B980 800D47E0 00001025 */   move  $v0, $zero
.L800D47E4:
/* B4B984 800D47E4 16160003 */  bne   $s0, $s6, .L800D47F4
/* B4B988 800D47E8 00000000 */   nop   
/* B4B98C 800D47EC 1000000B */  b     .L800D481C
/* B4B990 800D47F0 24020001 */   li    $v0, 1
.L800D47F4:
/* B4B994 800D47F4 16170003 */  bne   $s0, $s7, .L800D4804
/* B4B998 800D47F8 00000000 */   nop   
/* B4B99C 800D47FC 0C0359DC */  jal   FaultDrawer_SetOsSyncPrintfEnabled
/* B4B9A0 800D4800 24040001 */   li    $a0, 1
.L800D4804:
/* B4B9A4 800D4804 161EFFDC */  bne   $s0, $fp, .L800D4778
/* B4B9A8 800D4808 00000000 */   nop   
/* B4B9AC 800D480C 0C0359DC */  jal   FaultDrawer_SetOsSyncPrintfEnabled
/* B4B9B0 800D4810 00002025 */   move  $a0, $zero
/* B4B9B4 800D4814 1000FFD8 */  b     .L800D4778
/* B4B9B8 800D4818 00000000 */   nop   
.L800D481C:
/* B4B9BC 800D481C 8FBF003C */  lw    $ra, 0x3c($sp)
/* B4B9C0 800D4820 8FB00018 */  lw    $s0, 0x18($sp)
/* B4B9C4 800D4824 8FB1001C */  lw    $s1, 0x1c($sp)
/* B4B9C8 800D4828 8FB20020 */  lw    $s2, 0x20($sp)
/* B4B9CC 800D482C 8FB30024 */  lw    $s3, 0x24($sp)
/* B4B9D0 800D4830 8FB40028 */  lw    $s4, 0x28($sp)
/* B4B9D4 800D4834 8FB5002C */  lw    $s5, 0x2c($sp)
/* B4B9D8 800D4838 8FB60030 */  lw    $s6, 0x30($sp)
/* B4B9DC 800D483C 8FB70034 */  lw    $s7, 0x34($sp)
/* B4B9E0 800D4840 8FBE0038 */  lw    $fp, 0x38($sp)
/* B4B9E4 800D4844 03E00008 */  jr    $ra
/* B4B9E8 800D4848 27BD0048 */   addiu $sp, $sp, 0x48

