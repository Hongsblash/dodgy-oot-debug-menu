glabel func_80853080
/* 20E70 80853080 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 20E74 80853084 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 20E78 80853088 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 20E7C 8085308C AFBF001C */  sw      $ra, 0x001C($sp)           
/* 20E80 80853090 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 20E84 80853094 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 20E88 80853098 3C068084 */  lui     $a2, %hi(func_80840BC8)    ## $a2 = 80840000
/* 20E8C 8085309C 24C60BC8 */  addiu   $a2, $a2, %lo(func_80840BC8) ## $a2 = 80840BC8
/* 20E90 808530A0 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 20E94 808530A4 0C20D716 */  jal     func_80835C58              
/* 20E98 808530A8 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 20E9C 808530AC 0C20CCCE */  jal     func_80833338              
/* 20EA0 808530B0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 20EA4 808530B4 8FA40024 */  lw      $a0, 0x0024($sp)           
/* 20EA8 808530B8 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 20EAC 808530BC 0C20CAC3 */  jal     func_80832B0C              
/* 20EB0 808530C0 00403025 */  or      $a2, $v0, $zero            ## $a2 = 00000000
/* 20EB4 808530C4 860E00B6 */  lh      $t6, 0x00B6($s0)           ## 000000B6
/* 20EB8 808530C8 A60E083C */  sh      $t6, 0x083C($s0)           ## 0000083C
/* 20EBC 808530CC 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 20EC0 808530D0 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 20EC4 808530D4 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 20EC8 808530D8 03E00008 */  jr      $ra                        
/* 20ECC 808530DC 00000000 */  nop


