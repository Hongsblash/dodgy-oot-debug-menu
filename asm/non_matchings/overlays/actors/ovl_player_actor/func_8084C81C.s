glabel func_8084C81C
/* 1A60C 8084C81C 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 1A610 8084C820 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 1A614 8084C824 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 1A618 8084C828 8C8E0680 */  lw      $t6, 0x0680($a0)           ## 00000680
/* 1A61C 8084C82C 00803025 */  or      $a2, $a0, $zero            ## $a2 = 00000000
/* 1A620 8084C830 24C501B4 */  addiu   $a1, $a2, 0x01B4           ## $a1 = 000001B4
/* 1A624 8084C834 35CF0040 */  ori     $t7, $t6, 0x0040           ## $t7 = 00000040
/* 1A628 8084C838 AC8F0680 */  sw      $t7, 0x0680($a0)           ## 00000680
/* 1A62C 8084C83C AFA60018 */  sw      $a2, 0x0018($sp)           
/* 1A630 8084C840 0C028EF0 */  jal     func_800A3BC0              
/* 1A634 8084C844 8FA4001C */  lw      $a0, 0x001C($sp)           
/* 1A638 8084C848 1040000C */  beq     $v0, $zero, .L8084C87C     
/* 1A63C 8084C84C 8FA60018 */  lw      $a2, 0x0018($sp)           
/* 1A640 8084C850 00C02025 */  or      $a0, $a2, $zero            ## $a0 = 00000000
/* 1A644 8084C854 8FA5001C */  lw      $a1, 0x001C($sp)           
/* 1A648 8084C858 0C20F03A */  jal     func_8083C0E8              
/* 1A64C 8084C85C AFA60018 */  sw      $a2, 0x0018($sp)           
/* 1A650 8084C860 8FA60018 */  lw      $a2, 0x0018($sp)           
/* 1A654 8084C864 3C01FFFB */  lui     $at, 0xFFFB                ## $at = FFFB0000
/* 1A658 8084C868 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = FFFBFFFF
/* 1A65C 8084C86C 8CD80680 */  lw      $t8, 0x0680($a2)           ## 00000680
/* 1A660 8084C870 0301C824 */  and     $t9, $t8, $at              
/* 1A664 8084C874 10000005 */  beq     $zero, $zero, .L8084C88C   
/* 1A668 8084C878 ACD90680 */  sw      $t9, 0x0680($a2)           ## 00000680
.L8084C87C:
/* 1A66C 8084C87C 3C058085 */  lui     $a1, %hi(D_808548D8)       ## $a1 = 80850000
/* 1A670 8084C880 24A548D8 */  addiu   $a1, $a1, %lo(D_808548D8)  ## $a1 = 808548D8
/* 1A674 8084C884 0C20CA49 */  jal     func_80832924              
/* 1A678 8084C888 00C02025 */  or      $a0, $a2, $zero            ## $a0 = 00000000
.L8084C88C:
/* 1A67C 8084C88C 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 1A680 8084C890 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 1A684 8084C894 03E00008 */  jr      $ra                        
/* 1A688 8084C898 00000000 */  nop


