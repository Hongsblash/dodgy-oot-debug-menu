glabel func_80A8910C
/* 002FC 80A8910C 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00300 80A89110 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00304 80A89114 0C00BCCD */  jal     func_8002F334              
/* 00308 80A89118 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 0030C 80A8911C 1040000C */  beq     $v0, $zero, .L80A89150     
/* 00310 80A89120 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 00314 80A89124 240E6078 */  addiu   $t6, $zero, 0x6078         ## $t6 = 00006078
/* 00318 80A89128 A48E010E */  sh      $t6, 0x010E($a0)           ## 0000010E
/* 0031C 80A8912C 3C0580A9 */  lui     $a1, %hi(func_80A890C0)    ## $a1 = 80A90000
/* 00320 80A89130 24A590C0 */  addiu   $a1, $a1, %lo(func_80A890C0) ## $a1 = 80A890C0
/* 00324 80A89134 0C2A2384 */  jal     func_80A88E10              
/* 00328 80A89138 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 0032C 80A8913C 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 00330 80A89140 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00334 80A89144 8C8F0004 */  lw      $t7, 0x0004($a0)           ## 00000004
/* 00338 80A89148 01E1C025 */  or      $t8, $t7, $at              ## $t8 = 00010000
/* 0033C 80A8914C AC980004 */  sw      $t8, 0x0004($a0)           ## 00000004
.L80A89150:
/* 00340 80A89150 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 00344 80A89154 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 00348 80A89158 03E00008 */  jr      $ra                        
/* 0034C 80A8915C 00000000 */  nop


