glabel func_8087C65C
/* 0086C 8087C65C 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00870 8087C660 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 00874 8087C664 8FAE0018 */  lw      $t6, 0x0018($sp)           
/* 00878 8087C668 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 0087C 8087C66C AFA5001C */  sw      $a1, 0x001C($sp)           
/* 00880 8087C670 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 00884 8087C674 0C00B2D0 */  jal     Flags_GetSwitch
              
/* 00888 8087C678 91C50168 */  lbu     $a1, 0x0168($t6)           ## 00000168
/* 0088C 8087C67C 10400007 */  beq     $v0, $zero, .L8087C69C     
/* 00890 8087C680 8FA4001C */  lw      $a0, 0x001C($sp)           
/* 00894 8087C684 0C020120 */  jal     func_80080480              
/* 00898 8087C688 8FA50018 */  lw      $a1, 0x0018($sp)           
/* 0089C 8087C68C 8FB80018 */  lw      $t8, 0x0018($sp)           
/* 008A0 8087C690 3C0F8088 */  lui     $t7, %hi(func_8087C6AC)    ## $t7 = 80880000
/* 008A4 8087C694 25EFC6AC */  addiu   $t7, $t7, %lo(func_8087C6AC) ## $t7 = 8087C6AC
/* 008A8 8087C698 AF0F0164 */  sw      $t7, 0x0164($t8)           ## 00000164
.L8087C69C:
/* 008AC 8087C69C 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 008B0 8087C6A0 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 008B4 8087C6A4 03E00008 */  jr      $ra                        
/* 008B8 8087C6A8 00000000 */  nop


