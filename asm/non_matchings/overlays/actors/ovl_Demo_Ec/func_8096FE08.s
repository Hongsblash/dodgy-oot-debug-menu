glabel func_8096FE08
/* 02958 8096FE08 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 0295C 8096FE0C AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 02960 8096FE10 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 02964 8096FE14 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 02968 8096FE18 0C25B743 */  jal     func_8096DD0C              
/* 0296C 8096FE1C AFA5003C */  sw      $a1, 0x003C($sp)           
/* 02970 8096FE20 3C060601 */  lui     $a2, 0x0601                ## $a2 = 06010000
/* 02974 8096FE24 24C6FEF0 */  addiu   $a2, $a2, 0xFEF0           ## $a2 = 0600FEF0
/* 02978 8096FE28 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 0297C 8096FE2C 0C25B5CA */  jal     func_8096D728              
/* 02980 8096FE30 8FA5003C */  lw      $a1, 0x003C($sp)           
/* 02984 8096FE34 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 02988 8096FE38 0C25B76F */  jal     func_8096DDBC              
/* 0298C 8096FE3C 8FA5003C */  lw      $a1, 0x003C($sp)           
/* 02990 8096FE40 8602001C */  lh      $v0, 0x001C($s0)           ## 0000001C
/* 02994 8096FE44 2401001E */  addiu   $at, $zero, 0x001E         ## $at = 0000001E
/* 02998 8096FE48 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 0299C 8096FE4C 14410006 */  bne     $v0, $at, .L8096FE68       
/* 029A0 8096FE50 00003025 */  or      $a2, $zero, $zero          ## $a2 = 00000000
/* 029A4 8096FE54 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 029A8 8096FE58 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 029AC 8096FE5C 44810000 */  mtc1    $at, $f0                   ## $f0 = 1.00
/* 029B0 8096FE60 10000013 */  beq     $zero, $zero, .L8096FEB0   
/* 029B4 8096FE64 24A52FA0 */  addiu   $a1, $a1, 0x2FA0           ## $a1 = 06002FA0
.L8096FE68:
/* 029B8 8096FE68 2401001F */  addiu   $at, $zero, 0x001F         ## $at = 0000001F
/* 029BC 8096FE6C 14410005 */  bne     $v0, $at, .L8096FE84       
/* 029C0 8096FE70 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 029C4 8096FE74 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 029C8 8096FE78 44810000 */  mtc1    $at, $f0                   ## $f0 = 1.00
/* 029CC 8096FE7C 1000000C */  beq     $zero, $zero, .L8096FEB0   
/* 029D0 8096FE80 24A53A98 */  addiu   $a1, $a1, 0x3A98           ## $a1 = 06003A98
.L8096FE84:
/* 029D4 8096FE84 24010020 */  addiu   $at, $zero, 0x0020         ## $at = 00000020
/* 029D8 8096FE88 14410006 */  bne     $v0, $at, .L8096FEA4       
/* 029DC 8096FE8C 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 029E0 8096FE90 3C014170 */  lui     $at, 0x4170                ## $at = 41700000
/* 029E4 8096FE94 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 029E8 8096FE98 44810000 */  mtc1    $at, $f0                   ## $f0 = 15.00
/* 029EC 8096FE9C 10000004 */  beq     $zero, $zero, .L8096FEB0   
/* 029F0 8096FEA0 24A52FA0 */  addiu   $a1, $a1, 0x2FA0           ## $a1 = 06002FA0
.L8096FEA4:
/* 029F4 8096FEA4 3C0140A0 */  lui     $at, 0x40A0                ## $at = 40A00000
/* 029F8 8096FEA8 44810000 */  mtc1    $at, $f0                   ## $f0 = 5.00
/* 029FC 8096FEAC 24A502B8 */  addiu   $a1, $a1, 0x02B8           ## $a1 = 06003258
.L8096FEB0:
/* 02A00 8096FEB0 24070000 */  addiu   $a3, $zero, 0x0000         ## $a3 = 00000000
/* 02A04 8096FEB4 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 02A08 8096FEB8 0C25B5E7 */  jal     func_8096D79C              
/* 02A0C 8096FEBC E7A00028 */  swc1    $f0, 0x0028($sp)           
/* 02A10 8096FEC0 26020050 */  addiu   $v0, $s0, 0x0050           ## $v0 = 00000050
/* 02A14 8096FEC4 C7A00028 */  lwc1    $f0, 0x0028($sp)           
/* 02A18 8096FEC8 C4440000 */  lwc1    $f4, 0x0000($v0)           ## 00000050
/* 02A1C 8096FECC C4480004 */  lwc1    $f8, 0x0004($v0)           ## 00000054
/* 02A20 8096FED0 C4500008 */  lwc1    $f16, 0x0008($v0)          ## 00000058
/* 02A24 8096FED4 46002182 */  mul.s   $f6, $f4, $f0              
/* 02A28 8096FED8 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 02A2C 8096FEDC 46004282 */  mul.s   $f10, $f8, $f0             
/* 02A30 8096FEE0 00000000 */  nop
/* 02A34 8096FEE4 46008482 */  mul.s   $f18, $f16, $f0            
/* 02A38 8096FEE8 E4460000 */  swc1    $f6, 0x0000($v0)           ## 00000050
/* 02A3C 8096FEEC E44A0004 */  swc1    $f10, 0x0004($v0)          ## 00000054
/* 02A40 8096FEF0 E4520008 */  swc1    $f18, 0x0008($v0)          ## 00000058
/* 02A44 8096FEF4 0C25B593 */  jal     func_8096D64C              
/* 02A48 8096FEF8 8FA5003C */  lw      $a1, 0x003C($sp)           
/* 02A4C 8096FEFC 3C068003 */  lui     $a2, 0x8003                ## $a2 = 80030000
/* 02A50 8096FF00 24C6B5EC */  addiu   $a2, $a2, 0xB5EC           ## $a2 = 8002B5EC
/* 02A54 8096FF04 260400B4 */  addiu   $a0, $s0, 0x00B4           ## $a0 = 000000B4
/* 02A58 8096FF08 24050000 */  addiu   $a1, $zero, 0x0000         ## $a1 = 00000000
/* 02A5C 8096FF0C 0C00AC78 */  jal     ActorShape_Init
              
/* 02A60 8096FF10 3C0741F0 */  lui     $a3, 0x41F0                ## $a3 = 41F00000
/* 02A64 8096FF14 240E001B */  addiu   $t6, $zero, 0x001B         ## $t6 = 0000001B
/* 02A68 8096FF18 240F0016 */  addiu   $t7, $zero, 0x0016         ## $t7 = 00000016
/* 02A6C 8096FF1C AE0E0194 */  sw      $t6, 0x0194($s0)           ## 00000194
/* 02A70 8096FF20 AE0F0198 */  sw      $t7, 0x0198($s0)           ## 00000198
/* 02A74 8096FF24 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 02A78 8096FF28 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 02A7C 8096FF2C 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 02A80 8096FF30 03E00008 */  jr      $ra                        
/* 02A84 8096FF34 00000000 */  nop


