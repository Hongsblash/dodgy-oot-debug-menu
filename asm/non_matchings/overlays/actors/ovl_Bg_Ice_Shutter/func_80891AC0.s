glabel func_80891AC0
/* 00000 80891AC0 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 00004 80891AC4 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 00008 80891AC8 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0000C 80891ACC AFBF001C */  sw      $ra, 0x001C($sp)           
/* 00010 80891AD0 0C01DE1C */  jal     Math_Sins
              ## sins?
/* 00014 80891AD4 848400B4 */  lh      $a0, 0x00B4($a0)           ## 000000B4
/* 00018 80891AD8 C6040060 */  lwc1    $f4, 0x0060($s0)           ## 00000060
/* 0001C 80891ADC 46040182 */  mul.s   $f6, $f0, $f4              
/* 00020 80891AE0 E7A60024 */  swc1    $f6, 0x0024($sp)           
/* 00024 80891AE4 0C01DE0D */  jal     Math_Coss
              ## coss?
/* 00028 80891AE8 860400B4 */  lh      $a0, 0x00B4($s0)           ## 000000B4
/* 0002C 80891AEC C6080060 */  lwc1    $f8, 0x0060($s0)           ## 00000060
/* 00030 80891AF0 C610000C */  lwc1    $f16, 0x000C($s0)          ## 0000000C
/* 00034 80891AF4 860400B6 */  lh      $a0, 0x00B6($s0)           ## 000000B6
/* 00038 80891AF8 46080282 */  mul.s   $f10, $f0, $f8             
/* 0003C 80891AFC 46105480 */  add.s   $f18, $f10, $f16           
/* 00040 80891B00 0C01DE1C */  jal     Math_Sins
              ## sins?
/* 00044 80891B04 E6120028 */  swc1    $f18, 0x0028($s0)          ## 00000028
/* 00048 80891B08 C7A40024 */  lwc1    $f4, 0x0024($sp)           
/* 0004C 80891B0C C6080008 */  lwc1    $f8, 0x0008($s0)           ## 00000008
/* 00050 80891B10 860400B6 */  lh      $a0, 0x00B6($s0)           ## 000000B6
/* 00054 80891B14 46040182 */  mul.s   $f6, $f0, $f4              
/* 00058 80891B18 46083280 */  add.s   $f10, $f6, $f8             
/* 0005C 80891B1C 0C01DE0D */  jal     Math_Coss
              ## coss?
/* 00060 80891B20 E60A0024 */  swc1    $f10, 0x0024($s0)          ## 00000024
/* 00064 80891B24 C7B00024 */  lwc1    $f16, 0x0024($sp)          
/* 00068 80891B28 C6040010 */  lwc1    $f4, 0x0010($s0)           ## 00000010
/* 0006C 80891B2C 46100482 */  mul.s   $f18, $f0, $f16            
/* 00070 80891B30 46049180 */  add.s   $f6, $f18, $f4             
/* 00074 80891B34 E606002C */  swc1    $f6, 0x002C($s0)           ## 0000002C
/* 00078 80891B38 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 0007C 80891B3C 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 00080 80891B40 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 00084 80891B44 03E00008 */  jr      $ra                        
/* 00088 80891B48 00000000 */  nop


