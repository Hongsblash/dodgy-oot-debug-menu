glabel func_80835C58
/* 03A48 80835C58 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 03A4C 80835C5C AFBF001C */  sw      $ra, 0x001C($sp)           
/* 03A50 80835C60 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 03A54 80835C64 AFA40020 */  sw      $a0, 0x0020($sp)           
/* 03A58 80835C68 AFA60028 */  sw      $a2, 0x0028($sp)           
/* 03A5C 80835C6C AFA7002C */  sw      $a3, 0x002C($sp)           
/* 03A60 80835C70 8CA20674 */  lw      $v0, 0x0674($a1)           ## 00000674
/* 03A64 80835C74 3C0F8085 */  lui     $t7, %hi(func_8084E3C4)    ## $t7 = 80850000
/* 03A68 80835C78 00A08025 */  or      $s0, $a1, $zero            ## $s0 = 00000000
/* 03A6C 80835C7C 14C20003 */  bne     $a2, $v0, .L80835C8C       
/* 03A70 80835C80 25EFE3C4 */  addiu   $t7, $t7, %lo(func_8084E3C4) ## $t7 = 8084E3C4
/* 03A74 80835C84 10000044 */  beq     $zero, $zero, .L80835D98   
/* 03A78 80835C88 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L80835C8C:
/* 03A7C 80835C8C 144F0009 */  bne     $v0, $t7, .L80835CB4       
/* 03A80 80835C90 3C088085 */  lui     $t0, %hi(func_808507F4)    ## $t0 = 80850000
/* 03A84 80835C94 0C03B616 */  jal     func_800ED858              
/* 03A88 80835C98 00002025 */  or      $a0, $zero, $zero          ## $a0 = 00000000
/* 03A8C 80835C9C 8E180680 */  lw      $t8, 0x0680($s0)           ## 00000680
/* 03A90 80835CA0 3C01FCFF */  lui     $at, 0xFCFF                ## $at = FCFF0000
/* 03A94 80835CA4 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = FCFFFFFF
/* 03A98 80835CA8 0301C824 */  and     $t9, $t8, $at              
/* 03A9C 80835CAC 10000006 */  beq     $zero, $zero, .L80835CC8   
/* 03AA0 80835CB0 AE190680 */  sw      $t9, 0x0680($s0)           ## 00000680
.L80835CB4:
/* 03AA4 80835CB4 250807F4 */  addiu   $t0, $t0, %lo(func_808507F4) ## $t0 = 000007F4
/* 03AA8 80835CB8 14480003 */  bne     $v0, $t0, .L80835CC8       
/* 03AAC 80835CBC 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 03AB0 80835CC0 0C20C8D0 */  jal     func_80832340              
/* 03AB4 80835CC4 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
.L80835CC8:
/* 03AB8 80835CC8 820A0151 */  lb      $t2, 0x0151($s0)           ## 00000151
/* 03ABC 80835CCC 820B0154 */  lb      $t3, 0x0154($s0)           ## 00000154
/* 03AC0 80835CD0 8FA90028 */  lw      $t1, 0x0028($sp)           
/* 03AC4 80835CD4 114B000B */  beq     $t2, $t3, .L80835D04       
/* 03AC8 80835CD8 AE090674 */  sw      $t1, 0x0674($s0)           ## 00000674
/* 03ACC 80835CDC 8FAC002C */  lw      $t4, 0x002C($sp)           
/* 03AD0 80835CE0 318D0001 */  andi    $t5, $t4, 0x0001           ## $t5 = 00000000
/* 03AD4 80835CE4 11A00005 */  beq     $t5, $zero, .L80835CFC     
/* 03AD8 80835CE8 00000000 */  nop
/* 03ADC 80835CEC 8E0E067C */  lw      $t6, 0x067C($s0)           ## 0000067C
/* 03AE0 80835CF0 000E7A40 */  sll     $t7, $t6,  9               
/* 03AE4 80835CF4 05E20004 */  bltzl   $t7, .L80835D08            
/* 03AE8 80835CF8 8FB8002C */  lw      $t8, 0x002C($sp)           
.L80835CFC:
/* 03AEC 80835CFC 0C023B1C */  jal     func_8008EC70              
/* 03AF0 80835D00 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
.L80835D04:
/* 03AF4 80835D04 8FB8002C */  lw      $t8, 0x002C($sp)           
.L80835D08:
/* 03AF8 80835D08 33190001 */  andi    $t9, $t8, 0x0001           ## $t9 = 00000000
/* 03AFC 80835D0C 1720000D */  bne     $t9, $zero, .L80835D44     
/* 03B00 80835D10 00000000 */  nop
/* 03B04 80835D14 8E08067C */  lw      $t0, 0x067C($s0)           ## 0000067C
/* 03B08 80835D18 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 03B0C 80835D1C 31090800 */  andi    $t1, $t0, 0x0800           ## $t1 = 00000000
/* 03B10 80835D20 15200008 */  bne     $t1, $zero, .L80835D44     
/* 03B14 80835D24 00000000 */  nop
/* 03B18 80835D28 0C20D191 */  jal     func_80834644              
/* 03B1C 80835D2C 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 03B20 80835D30 8E0A067C */  lw      $t2, 0x067C($s0)           ## 0000067C
/* 03B24 80835D34 3C01FFBF */  lui     $at, 0xFFBF                ## $at = FFBF0000
/* 03B28 80835D38 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = FFBFFFFF
/* 03B2C 80835D3C 01415824 */  and     $t3, $t2, $at              
/* 03B30 80835D40 AE0B067C */  sw      $t3, 0x067C($s0)           ## 0000067C
.L80835D44:
/* 03B34 80835D44 0C20CB6F */  jal     func_80832DBC              
/* 03B38 80835D48 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03B3C 80835D4C 8E0C067C */  lw      $t4, 0x067C($s0)           ## 0000067C
/* 03B40 80835D50 3C014BFF */  lui     $at, 0x4BFF                ## $at = 4BFF0000
/* 03B44 80835D54 3421FFBB */  ori     $at, $at, 0xFFBB           ## $at = 4BFFFFBB
/* 03B48 80835D58 8E0E0680 */  lw      $t6, 0x0680($s0)           ## 00000680
/* 03B4C 80835D5C 92180692 */  lbu     $t8, 0x0692($s0)           ## 00000692
/* 03B50 80835D60 01816824 */  and     $t5, $t4, $at              
/* 03B54 80835D64 3C01E7F7 */  lui     $at, 0xE7F7                ## $at = E7F70000
/* 03B58 80835D68 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = E7F7FFFF
/* 03B5C 80835D6C 01C17824 */  and     $t7, $t6, $at              
/* 03B60 80835D70 3319FF75 */  andi    $t9, $t8, 0xFF75           ## $t9 = 00000000
/* 03B64 80835D74 AE0D067C */  sw      $t5, 0x067C($s0)           ## 0000067C
/* 03B68 80835D78 AE0F0680 */  sw      $t7, 0x0680($s0)           ## 00000680
/* 03B6C 80835D7C A2190692 */  sb      $t9, 0x0692($s0)           ## 00000692
/* 03B70 80835D80 A200084F */  sb      $zero, 0x084F($s0)         ## 0000084F
/* 03B74 80835D84 A6000850 */  sh      $zero, 0x0850($s0)         ## 00000850
/* 03B78 80835D88 A20006AC */  sb      $zero, 0x06AC($s0)         ## 000006AC
/* 03B7C 80835D8C 0C20C9BC */  jal     func_808326F0              
/* 03B80 80835D90 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03B84 80835D94 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L80835D98:
/* 03B88 80835D98 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 03B8C 80835D9C 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 03B90 80835DA0 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 03B94 80835DA4 03E00008 */  jr      $ra                        
/* 03B98 80835DA8 00000000 */  nop


