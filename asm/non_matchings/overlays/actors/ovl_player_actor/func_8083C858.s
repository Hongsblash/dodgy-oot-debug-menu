glabel func_8083C858
/* 0A648 8083C858 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 0A64C 8083C85C AFBF001C */  sw      $ra, 0x001C($sp)           
/* 0A650 8083C860 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 0A654 8083C864 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0A658 8083C868 0C20CEF3 */  jal     func_80833BCC              
/* 0A65C 8083C86C AFA50024 */  sw      $a1, 0x0024($sp)           
/* 0A660 8083C870 10400004 */  beq     $v0, $zero, .L8083C884     
/* 0A664 8083C874 8FA40024 */  lw      $a0, 0x0024($sp)           
/* 0A668 8083C878 3C068084 */  lui     $a2, %hi(func_8084227C)    ## $a2 = 80840000
/* 0A66C 8083C87C 10000003 */  beq     $zero, $zero, .L8083C88C   
/* 0A670 8083C880 24C6227C */  addiu   $a2, $a2, %lo(func_8084227C) ## $a2 = 8084227C
.L8083C884:
/* 0A674 8083C884 3C068084 */  lui     $a2, %hi(func_80842180)    ## $a2 = 80840000
/* 0A678 8083C888 24C62180 */  addiu   $a2, $a2, %lo(func_80842180) ## $a2 = 80842180
.L8083C88C:
/* 0A67C 8083C88C 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0A680 8083C890 0C20D716 */  jal     func_80835C58              
/* 0A684 8083C894 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 0A688 8083C898 920E015B */  lbu     $t6, 0x015B($s0)           ## 0000015B
/* 0A68C 8083C89C 3C068085 */  lui     $a2, %hi(D_80853944)       ## $a2 = 80850000
/* 0A690 8083C8A0 8FA40024 */  lw      $a0, 0x0024($sp)           
/* 0A694 8083C8A4 000E7880 */  sll     $t7, $t6,  2               
/* 0A698 8083C8A8 00CF3021 */  addu    $a2, $a2, $t7              
/* 0A69C 8083C8AC 8CC63944 */  lw      $a2, %lo(D_80853944)($a2)  
/* 0A6A0 8083C8B0 0C20CAFA */  jal     func_80832BE8              
/* 0A6A4 8083C8B4 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0A6A8 8083C8B8 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 0A6AC 8083C8BC A600089C */  sh      $zero, 0x089C($s0)         ## 0000089C
/* 0A6B0 8083C8C0 E6000868 */  swc1    $f0, 0x0868($s0)           ## 00000868
/* 0A6B4 8083C8C4 E6000864 */  swc1    $f0, 0x0864($s0)           ## 00000864
/* 0A6B8 8083C8C8 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 0A6BC 8083C8CC 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 0A6C0 8083C8D0 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 0A6C4 8083C8D4 03E00008 */  jr      $ra                        
/* 0A6C8 8083C8D8 00000000 */  nop


