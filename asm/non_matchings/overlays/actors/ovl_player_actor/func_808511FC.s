glabel func_808511FC
/* 1EFEC 808511FC 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 1EFF0 80851200 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 1EFF4 80851204 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 1EFF8 80851208 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 1EFFC 8085120C AFA60020 */  sw      $a2, 0x0020($sp)           
/* 1F000 80851210 0C028EF0 */  jal     func_800A3BC0              
/* 1F004 80851214 24A501B4 */  addiu   $a1, $a1, 0x01B4           ## $a1 = 000001B4
/* 1F008 80851218 10400007 */  beq     $v0, $zero, .L80851238     
/* 1F00C 8085121C 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 1F010 80851220 8FA5001C */  lw      $a1, 0x001C($sp)           
/* 1F014 80851224 0C2143E7 */  jal     func_80850F9C              
/* 1F018 80851228 8FA60020 */  lw      $a2, 0x0020($sp)           
/* 1F01C 8085122C 8FAF001C */  lw      $t7, 0x001C($sp)           
/* 1F020 80851230 240E0001 */  addiu   $t6, $zero, 0x0001         ## $t6 = 00000001
/* 1F024 80851234 A5EE0850 */  sh      $t6, 0x0850($t7)           ## 00000850
.L80851238:
/* 1F028 80851238 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 1F02C 8085123C 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 1F030 80851240 03E00008 */  jr      $ra                        
/* 1F034 80851244 00000000 */  nop


