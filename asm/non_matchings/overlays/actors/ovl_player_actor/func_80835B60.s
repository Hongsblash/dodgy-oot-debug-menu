glabel func_80835B60
/* 03950 80835B60 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 03954 80835B64 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 03958 80835B68 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0395C 80835B6C AFBF001C */  sw      $ra, 0x001C($sp)           
/* 03960 80835B70 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 03964 80835B74 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 03968 80835B78 0C20D1D6 */  jal     func_80834758              
/* 0396C 80835B7C 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 03970 80835B80 50400004 */  beql    $v0, $zero, .L80835B94     
/* 03974 80835B84 8E0E067C */  lw      $t6, 0x067C($s0)           ## 0000067C
/* 03978 80835B88 1000001A */  beq     $zero, $zero, .L80835BF4   
/* 0397C 80835B8C 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
/* 03980 80835B90 8E0E067C */  lw      $t6, 0x067C($s0)           ## 0000067C
.L80835B94:
/* 03984 80835B94 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03988 80835B98 3C058083 */  lui     $a1, %hi(func_80835C08)    ## $a1 = 80830000
/* 0398C 80835B9C 000E7980 */  sll     $t7, $t6,  6               
/* 03990 80835BA0 05E00014 */  bltz    $t7, .L80835BF4            
/* 03994 80835BA4 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 03998 80835BA8 0C20CD8E */  jal     func_80833638              
/* 0399C 80835BAC 24A55C08 */  addiu   $a1, $a1, %lo(func_80835C08) ## $a1 = 80835C08
/* 039A0 80835BB0 3C060400 */  lui     $a2, 0x0400                ## $a2 = 04000000
/* 039A4 80835BB4 24C625F8 */  addiu   $a2, $a2, 0x25F8           ## $a2 = 040025F8
/* 039A8 80835BB8 8FA40024 */  lw      $a0, 0x0024($sp)           
/* 039AC 80835BBC 0C02901F */  jal     func_800A407C              
/* 039B0 80835BC0 260506C8 */  addiu   $a1, $s0, 0x06C8           ## $a1 = 000006C8
/* 039B4 80835BC4 3C058012 */  lui     $a1, 0x8012                ## $a1 = 80120000
/* 039B8 80835BC8 24A55EF8 */  addiu   $a1, $a1, 0x5EF8           ## $a1 = 80125EF8
/* 039BC 80835BCC 0C20D5FA */  jal     func_808357E8              
/* 039C0 80835BD0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 039C4 80835BD4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 039C8 80835BD8 0C00BDF7 */  jal     func_8002F7DC              
/* 039CC 80835BDC 24050836 */  addiu   $a1, $zero, 0x0836         ## $a1 = 00000836
/* 039D0 80835BE0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 039D4 80835BE4 0C20C9A6 */  jal     func_80832698              
/* 039D8 80835BE8 24056800 */  addiu   $a1, $zero, 0x6800         ## $a1 = 00006800
/* 039DC 80835BEC 10000001 */  beq     $zero, $zero, .L80835BF4   
/* 039E0 80835BF0 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L80835BF4:
/* 039E4 80835BF4 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 039E8 80835BF8 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 039EC 80835BFC 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 039F0 80835C00 03E00008 */  jr      $ra                        
/* 039F4 80835C04 00000000 */  nop


