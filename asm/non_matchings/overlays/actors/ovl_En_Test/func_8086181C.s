glabel func_8086181C
/* 021CC 8086181C 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 021D0 80861820 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 021D4 80861824 00803025 */  or      $a2, $a0, $zero            ## $a2 = 00000000
/* 021D8 80861828 3C050601 */  lui     $a1, 0x0601                ## $a1 = 06010000
/* 021DC 8086182C 24A5BE4C */  addiu   $a1, $a1, 0xBE4C           ## $a1 = 0600BE4C
/* 021E0 80861830 AFA60018 */  sw      $a2, 0x0018($sp)           
/* 021E4 80861834 0C02947A */  jal     SkelAnimeChangeAnimationDefaultStop              
/* 021E8 80861838 24840188 */  addiu   $a0, $a0, 0x0188           ## $a0 = 00000188
/* 021EC 8086183C 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 021F0 80861840 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 021F4 80861844 24180011 */  addiu   $t8, $zero, 0x0011         ## $t8 = 00000011
/* 021F8 80861848 908E086C */  lbu     $t6, 0x086C($a0)           ## 0000086C
/* 021FC 8086184C 24190010 */  addiu   $t9, $zero, 0x0010         ## $t9 = 00000010
/* 02200 80861850 3C058086 */  lui     $a1, %hi(func_80861898)    ## $a1 = 80860000
/* 02204 80861854 31CFFFFB */  andi    $t7, $t6, 0xFFFB           ## $t7 = 00000000
/* 02208 80861858 24A51898 */  addiu   $a1, $a1, %lo(func_80861898) ## $a1 = 80861898
/* 0220C 8086185C A08F086C */  sb      $t7, 0x086C($a0)           ## 0000086C
/* 02210 80861860 A09807C8 */  sb      $t8, 0x07C8($a0)           ## 000007C8
/* 02214 80861864 A0990879 */  sb      $t9, 0x0879($a0)           ## 00000879
/* 02218 80861868 0C217D94 */  jal     func_8085F650              
/* 0221C 8086186C E4840068 */  swc1    $f4, 0x0068($a0)           ## 00000068
/* 02220 80861870 8FA60018 */  lw      $a2, 0x0018($sp)           
/* 02224 80861874 24090003 */  addiu   $t1, $zero, 0x0003         ## $t1 = 00000003
/* 02228 80861878 90C807DE */  lbu     $t0, 0x07DE($a2)           ## 000007DE
/* 0222C 8086187C 51000003 */  beql    $t0, $zero, .L8086188C     
/* 02230 80861880 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 02234 80861884 A0C907DE */  sb      $t1, 0x07DE($a2)           ## 000007DE
/* 02238 80861888 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L8086188C:
/* 0223C 8086188C 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 02240 80861890 03E00008 */  jr      $ra                        
/* 02244 80861894 00000000 */  nop


