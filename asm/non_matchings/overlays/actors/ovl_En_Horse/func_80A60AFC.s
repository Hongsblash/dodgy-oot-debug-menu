glabel func_80A60AFC
/* 0580C 80A60AFC 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 05810 80A60B00 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 05814 80A60B04 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 05818 80A60B08 AFA60020 */  sw      $a2, 0x0020($sp)           
/* 0581C 80A60B0C 8CCF000C */  lw      $t7, 0x000C($a2)           ## 0000000C
/* 05820 80A60B10 00803825 */  or      $a3, $a0, $zero            ## $a3 = 00000000
/* 05824 80A60B14 3C0180A6 */  lui     $at, %hi(D_80A66914)       ## $at = 80A60000
/* 05828 80A60B18 448F2000 */  mtc1    $t7, $f4                   ## $f4 = 0.00
/* 0582C 80A60B1C 248401AC */  addiu   $a0, $a0, 0x01AC           ## $a0 = 000001AC
/* 05830 80A60B20 468021A0 */  cvt.s.w $f6, $f4                   
/* 05834 80A60B24 E486FE78 */  swc1    $f6, -0x0188($a0)          ## 00000024
/* 05838 80A60B28 8FB80020 */  lw      $t8, 0x0020($sp)           
/* 0583C 80A60B2C 8C8BFE78 */  lw      $t3, -0x0188($a0)          ## 00000024
/* 05840 80A60B30 8F190010 */  lw      $t9, 0x0010($t8)           ## 00000010
/* 05844 80A60B34 24180006 */  addiu   $t8, $zero, 0x0006         ## $t8 = 00000006
/* 05848 80A60B38 44994000 */  mtc1    $t9, $f8                   ## $f8 = 0.00
/* 0584C 80A60B3C 24190004 */  addiu   $t9, $zero, 0x0004         ## $t9 = 00000004
/* 05850 80A60B40 468042A0 */  cvt.s.w $f10, $f8                  
/* 05854 80A60B44 E48AFE7C */  swc1    $f10, -0x0184($a0)         ## 00000028
/* 05858 80A60B48 8FA80020 */  lw      $t0, 0x0020($sp)           
/* 0585C 80A60B4C 8C8AFE7C */  lw      $t2, -0x0184($a0)          ## 00000028
/* 05860 80A60B50 8D090014 */  lw      $t1, 0x0014($t0)           ## 00000014
/* 05864 80A60B54 AC8BFF54 */  sw      $t3, -0x00AC($a0)          ## 00000100
/* 05868 80A60B58 AC8AFF58 */  sw      $t2, -0x00A8($a0)          ## 00000104
/* 0586C 80A60B5C 44898000 */  mtc1    $t1, $f16                  ## $f16 = 0.00
/* 05870 80A60B60 3C0A80A6 */  lui     $t2, %hi(D_80A65E58)       ## $t2 = 80A60000
/* 05874 80A60B64 468084A0 */  cvt.s.w $f18, $f16                 
/* 05878 80A60B68 E492FE80 */  swc1    $f18, -0x0180($a0)         ## 0000002C
/* 0587C 80A60B6C 8C8BFE80 */  lw      $t3, -0x0180($a0)          ## 0000002C
/* 05880 80A60B70 AC8BFF5C */  sw      $t3, -0x00A4($a0)          ## 00000108
/* 05884 80A60B74 8FAC0020 */  lw      $t4, 0x0020($sp)           
/* 05888 80A60B78 958D0008 */  lhu     $t5, 0x0008($t4)           ## 00000008
/* 0588C 80A60B7C AC980064 */  sw      $t8, 0x0064($a0)           ## 00000210
/* 05890 80A60B80 AC9901D4 */  sw      $t9, 0x01D4($a0)           ## 00000380
/* 05894 80A60B84 A48DFE86 */  sh      $t5, -0x017A($a0)          ## 00000032
/* 05898 80A60B88 888FFE84 */  lwl     $t7, -0x017C($a0)          ## 00000030
/* 0589C 80A60B8C 988FFE87 */  lwr     $t7, -0x0179($a0)          ## 00000033
/* 058A0 80A60B90 A88FFF08 */  swl     $t7, -0x00F8($a0)          ## 000000B4
/* 058A4 80A60B94 B88FFF0B */  swr     $t7, -0x00F5($a0)          ## 000000B7
/* 058A8 80A60B98 948FFE88 */  lhu     $t7, -0x0178($a0)          ## 00000034
/* 058AC 80A60B9C A48FFF0C */  sh      $t7, -0x00F4($a0)          ## 000000B8
/* 058B0 80A60BA0 C4266914 */  lwc1    $f6, %lo(D_80A66914)($at)  
/* 058B4 80A60BA4 C4E40068 */  lwc1    $f4, 0x0068($a3)           ## 00000068
/* 058B8 80A60BA8 8CE80158 */  lw      $t0, 0x0158($a3)           ## 00000158
/* 058BC 80A60BAC 46062202 */  mul.s   $f8, $f4, $f6              
/* 058C0 80A60BB0 00084880 */  sll     $t1, $t0,  2               
/* 058C4 80A60BB4 01495021 */  addu    $t2, $t2, $t1              
/* 058C8 80A60BB8 8D4A5E58 */  lw      $t2, %lo(D_80A65E58)($t2)  
/* 058CC 80A60BBC 8D450018 */  lw      $a1, 0x0018($t2)           ## 80A60018
/* 058D0 80A60BC0 44064000 */  mfc1    $a2, $f8                   
/* 058D4 80A60BC4 0C0294A7 */  jal     func_800A529C              
/* 058D8 80A60BC8 00000000 */  nop
/* 058DC 80A60BCC 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 058E0 80A60BD0 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 058E4 80A60BD4 03E00008 */  jr      $ra                        
/* 058E8 80A60BD8 00000000 */  nop


