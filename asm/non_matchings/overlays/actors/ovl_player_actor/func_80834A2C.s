glabel func_80834A2C
/* 0281C 80834A2C 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 02820 80834A30 AFB00014 */  sw      $s0, 0x0014($sp)           
/* 02824 80834A34 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 02828 80834A38 AFB10018 */  sw      $s1, 0x0018($sp)           
/* 0282C 80834A3C 00A08825 */  or      $s1, $a1, $zero            ## $s1 = 00000000
/* 02830 80834A40 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 02834 80834A44 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 02838 80834A48 0C028EF0 */  jal     func_800A3BC0              
/* 0283C 80834A4C 260506C8 */  addiu   $a1, $s0, 0x06C8           ## $a1 = 000006C8
/* 02840 80834A50 54400018 */  bnel    $v0, $zero, .L80834AB4     
/* 02844 80834A54 82180151 */  lb      $t8, 0x0151($s0)           ## 00000151
/* 02848 80834A58 0C20CDC9 */  jal     func_80833724              
/* 0284C 80834A5C 92040152 */  lbu     $a0, 0x0152($s0)           ## 00000152
/* 02850 80834A60 820E0151 */  lb      $t6, 0x0151($s0)           ## 00000151
/* 02854 80834A64 3C038085 */  lui     $v1, %hi(D_80853614)       ## $v1 = 80850000
/* 02858 80834A68 24633614 */  addiu   $v1, $v1, %lo(D_80853614)  ## $v1 = 80853614
/* 0285C 80834A6C 144E0024 */  bne     $v0, $t6, .L80834B00       
/* 02860 80834A70 00000000 */  nop
/* 02864 80834A74 8C620000 */  lw      $v0, 0x0000($v1)           ## 80853614
/* 02868 80834A78 0002102B */  sltu    $v0, $zero, $v0            
/* 0286C 80834A7C 1440000A */  bne     $v0, $zero, .L80834AA8     
/* 02870 80834A80 00000000 */  nop
/* 02874 80834A84 9202015B */  lbu     $v0, 0x015B($s0)           ## 0000015B
/* 02878 80834A88 38420003 */  xori    $v0, $v0, 0x0003           ## $v0 = 00000003
/* 0287C 80834A8C 0002102B */  sltu    $v0, $zero, $v0            
/* 02880 80834A90 10400005 */  beq     $v0, $zero, .L80834AA8     
/* 02884 80834A94 00000000 */  nop
/* 02888 80834A98 3C020001 */  lui     $v0, 0x0001                ## $v0 = 00010000
/* 0288C 80834A9C 00511021 */  addu    $v0, $v0, $s1              
/* 02890 80834AA0 80421E5C */  lb      $v0, 0x1E5C($v0)           ## 00011E5C
/* 02894 80834AA4 2C420001 */  sltiu   $v0, $v0, 0x0001           
.L80834AA8:
/* 02898 80834AA8 10400015 */  beq     $v0, $zero, .L80834B00     
/* 0289C 80834AAC AC620000 */  sw      $v0, 0x0000($v1)           ## 80853614
/* 028A0 80834AB0 82180151 */  lb      $t8, 0x0151($s0)           ## 00000151
.L80834AB4:
/* 028A4 80834AB4 3C058085 */  lui     $a1, %hi(D_80853EDC)       ## $a1 = 80850000
/* 028A8 80834AB8 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 028AC 80834ABC 0018C880 */  sll     $t9, $t8,  2               
/* 028B0 80834AC0 00B92821 */  addu    $a1, $a1, $t9              
/* 028B4 80834AC4 0C20CD8E */  jal     func_80833638              
/* 028B8 80834AC8 8CA53EDC */  lw      $a1, %lo(D_80853EDC)($a1)  
/* 028BC 80834ACC A6000834 */  sh      $zero, 0x0834($s0)         ## 00000834
/* 028C0 80834AD0 A20006AC */  sb      $zero, 0x06AC($s0)         ## 000006AC
/* 028C4 80834AD4 3C088085 */  lui     $t0, %hi(D_80853614)       ## $t0 = 80850000
/* 028C8 80834AD8 8D083614 */  lw      $t0, %lo(D_80853614)($t0)  
/* 028CC 80834ADC 3C018085 */  lui     $at, %hi(D_80853618)       ## $at = 80850000
/* 028D0 80834AE0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 028D4 80834AE4 AC283618 */  sw      $t0, %lo(D_80853618)($at)  
/* 028D8 80834AE8 8E19082C */  lw      $t9, 0x082C($s0)           ## 0000082C
/* 028DC 80834AEC 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 028E0 80834AF0 0320F809 */  jalr    $ra, $t9                   
/* 028E4 80834AF4 00000000 */  nop
/* 028E8 80834AF8 10000014 */  beq     $zero, $zero, .L80834B4C   
/* 028EC 80834AFC 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80834B00:
/* 028F0 80834B00 0C20CCD4 */  jal     func_80833350              
/* 028F4 80834B04 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 028F8 80834B08 1040000C */  beq     $v0, $zero, .L80834B3C     
/* 028FC 80834B0C 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 02900 80834B10 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 02904 80834B14 0C20D23B */  jal     func_808348EC              
/* 02908 80834B18 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0290C 80834B1C 0C20CCCE */  jal     func_80833338              
/* 02910 80834B20 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 02914 80834B24 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 02918 80834B28 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0291C 80834B2C 0C20C899 */  jal     func_80832264              
/* 02920 80834B30 00403025 */  or      $a2, $v0, $zero            ## $a2 = 00000000
/* 02924 80834B34 10000003 */  beq     $zero, $zero, .L80834B44   
/* 02928 80834B38 A20006AC */  sb      $zero, 0x06AC($s0)         ## 000006AC
.L80834B3C:
/* 0292C 80834B3C 0C20D23B */  jal     func_808348EC              
/* 02930 80834B40 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
.L80834B44:
/* 02934 80834B44 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
/* 02938 80834B48 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80834B4C:
/* 0293C 80834B4C 8FB00014 */  lw      $s0, 0x0014($sp)           
/* 02940 80834B50 8FB10018 */  lw      $s1, 0x0018($sp)           
/* 02944 80834B54 03E00008 */  jr      $ra                        
/* 02948 80834B58 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000


