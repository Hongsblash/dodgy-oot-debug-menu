glabel func_80A35A08
/* 006F8 80A35A08 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 006FC 80A35A0C AFBF0034 */  sw      $ra, 0x0034($sp)           
/* 00700 80A35A10 AFB00030 */  sw      $s0, 0x0030($sp)           
/* 00704 80A35A14 AFA5003C */  sw      $a1, 0x003C($sp)           
/* 00708 80A35A18 848E0318 */  lh      $t6, 0x0318($a0)           ## 00000318
/* 0070C 80A35A1C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00710 80A35A20 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 00714 80A35A24 51C00005 */  beql    $t6, $zero, .L80A35A3C     
/* 00718 80A35A28 3C014396 */  lui     $at, 0x4396                ## $at = 43960000
/* 0071C 80A35A2C 0C00B2D0 */  jal     Flags_GetSwitch
              
/* 00720 80A35A30 86050018 */  lh      $a1, 0x0018($s0)           ## 00000018
/* 00724 80A35A34 10400009 */  beq     $v0, $zero, .L80A35A5C     
/* 00728 80A35A38 3C014396 */  lui     $at, 0x4396                ## $at = 43960000
.L80A35A3C:
/* 0072C 80A35A3C 44812000 */  mtc1    $at, $f4                   ## $f4 = 300.00
/* 00730 80A35A40 C6060090 */  lwc1    $f6, 0x0090($s0)           ## 00000090
/* 00734 80A35A44 3C0142B4 */  lui     $at, 0x42B4                ## $at = 42B40000
/* 00738 80A35A48 24040038 */  addiu   $a0, $zero, 0x0038         ## $a0 = 00000038
/* 0073C 80A35A4C 4606203C */  c.lt.s  $f4, $f6                   
/* 00740 80A35A50 00000000 */  nop
/* 00744 80A35A54 4502000B */  bc1fl   .L80A35A84                 
/* 00748 80A35A58 44819000 */  mtc1    $at, $f18                  ## $f18 = 90.00
.L80A35A5C:
/* 0074C 80A35A5C 3C0142F0 */  lui     $at, 0x42F0                ## $at = 42F00000
/* 00750 80A35A60 44815000 */  mtc1    $at, $f10                  ## $f10 = 120.00
/* 00754 80A35A64 C6080080 */  lwc1    $f8, 0x0080($s0)           ## 00000080
/* 00758 80A35A68 8602008A */  lh      $v0, 0x008A($s0)           ## 0000008A
/* 0075C 80A35A6C 460A4400 */  add.s   $f16, $f8, $f10            
/* 00760 80A35A70 A6020032 */  sh      $v0, 0x0032($s0)           ## 00000032
/* 00764 80A35A74 A60200B6 */  sh      $v0, 0x00B6($s0)           ## 000000B6
/* 00768 80A35A78 10000005 */  beq     $zero, $zero, .L80A35A90   
/* 0076C 80A35A7C E6100028 */  swc1    $f16, 0x0028($s0)          ## 00000028
/* 00770 80A35A80 44819000 */  mtc1    $at, $f18                  ## $f18 = 120.00
.L80A35A84:
/* 00774 80A35A84 A6000318 */  sh      $zero, 0x0318($s0)         ## 00000318
/* 00778 80A35A88 0C03D6B3 */  jal     func_800F5ACC              
/* 0077C 80A35A8C E61200C4 */  swc1    $f18, 0x00C4($s0)          ## 000000C4
.L80A35A90:
/* 00780 80A35A90 960F0088 */  lhu     $t7, 0x0088($s0)           ## 00000088
/* 00784 80A35A94 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00788 80A35A98 31F80002 */  andi    $t8, $t7, 0x0002           ## $t8 = 00000000
/* 0078C 80A35A9C 13000030 */  beq     $t8, $zero, .L80A35B60     
/* 00790 80A35AA0 00000000 */  nop
/* 00794 80A35AA4 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 00798 80A35AA8 2405387A */  addiu   $a1, $zero, 0x387A         ## $a1 = 0000387A
/* 0079C 80A35AAC C6060080 */  lwc1    $f6, 0x0080($s0)           ## 00000080
/* 007A0 80A35AB0 8E0A0024 */  lw      $t2, 0x0024($s0)           ## 00000024
/* 007A4 80A35AB4 8E190004 */  lw      $t9, 0x0004($s0)           ## 00000004
/* 007A8 80A35AB8 960B0088 */  lhu     $t3, 0x0088($s0)           ## 00000088
/* 007AC 80A35ABC 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 007B0 80A35AC0 E6060028 */  swc1    $f6, 0x0028($s0)           ## 00000028
/* 007B4 80A35AC4 AE0A0038 */  sw      $t2, 0x0038($s0)           ## 00000038
/* 007B8 80A35AC8 44812000 */  mtc1    $at, $f4                   ## $f4 = 1.00
/* 007BC 80A35ACC 8E0A002C */  lw      $t2, 0x002C($s0)           ## 0000002C
/* 007C0 80A35AD0 8E090028 */  lw      $t1, 0x0028($s0)           ## 00000028
/* 007C4 80A35AD4 44804000 */  mtc1    $zero, $f8                 ## $f8 = 0.00
/* 007C8 80A35AD8 37280001 */  ori     $t0, $t9, 0x0001           ## $t0 = 00000001
/* 007CC 80A35ADC 316CFFFD */  andi    $t4, $t3, 0xFFFD           ## $t4 = 00000000
/* 007D0 80A35AE0 AE080004 */  sw      $t0, 0x0004($s0)           ## 00000004
/* 007D4 80A35AE4 A60C0088 */  sh      $t4, 0x0088($s0)           ## 00000088
/* 007D8 80A35AE8 E60401A4 */  swc1    $f4, 0x01A4($s0)           ## 000001A4
/* 007DC 80A35AEC AE0A0040 */  sw      $t2, 0x0040($s0)           ## 00000040
/* 007E0 80A35AF0 AE09003C */  sw      $t1, 0x003C($s0)           ## 0000003C
/* 007E4 80A35AF4 E6080060 */  swc1    $f8, 0x0060($s0)           ## 00000060
/* 007E8 80A35AF8 3C014000 */  lui     $at, 0x4000                ## $at = 40000000
/* 007EC 80A35AFC 44815000 */  mtc1    $at, $f10                  ## $f10 = 2.00
/* 007F0 80A35B00 240D0002 */  addiu   $t5, $zero, 0x0002         ## $t5 = 00000002
/* 007F4 80A35B04 AFAD0010 */  sw      $t5, 0x0010($sp)           
/* 007F8 80A35B08 AFA00020 */  sw      $zero, 0x0020($sp)         
/* 007FC 80A35B0C AFA0001C */  sw      $zero, 0x001C($sp)         
/* 00800 80A35B10 AFA00018 */  sw      $zero, 0x0018($sp)         
/* 00804 80A35B14 8FA4003C */  lw      $a0, 0x003C($sp)           
/* 00808 80A35B18 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0080C 80A35B1C 260604D0 */  addiu   $a2, $s0, 0x04D0           ## $a2 = 000004D0
/* 00810 80A35B20 3C074040 */  lui     $a3, 0x4040                ## $a3 = 40400000
/* 00814 80A35B24 0C00CC98 */  jal     func_80033260              
/* 00818 80A35B28 E7AA0014 */  swc1    $f10, 0x0014($sp)          
/* 0081C 80A35B2C 3C014000 */  lui     $at, 0x4000                ## $at = 40000000
/* 00820 80A35B30 44818000 */  mtc1    $at, $f16                  ## $f16 = 2.00
/* 00824 80A35B34 240E0002 */  addiu   $t6, $zero, 0x0002         ## $t6 = 00000002
/* 00828 80A35B38 AFAE0010 */  sw      $t6, 0x0010($sp)           
/* 0082C 80A35B3C 8FA4003C */  lw      $a0, 0x003C($sp)           
/* 00830 80A35B40 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 00834 80A35B44 260604C4 */  addiu   $a2, $s0, 0x04C4           ## $a2 = 000004C4
/* 00838 80A35B48 3C074040 */  lui     $a3, 0x4040                ## $a3 = 40400000
/* 0083C 80A35B4C AFA00018 */  sw      $zero, 0x0018($sp)         
/* 00840 80A35B50 AFA0001C */  sw      $zero, 0x001C($sp)         
/* 00844 80A35B54 AFA00020 */  sw      $zero, 0x0020($sp)         
/* 00848 80A35B58 0C00CC98 */  jal     func_80033260              
/* 0084C 80A35B5C E7B00014 */  swc1    $f16, 0x0014($sp)          
.L80A35B60:
/* 00850 80A35B60 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 00854 80A35B64 26040188 */  addiu   $a0, $s0, 0x0188           ## $a0 = 00000188
/* 00858 80A35B68 50400004 */  beql    $v0, $zero, .L80A35B7C     
/* 0085C 80A35B6C 8FBF0034 */  lw      $ra, 0x0034($sp)           
/* 00860 80A35B70 0C28D752 */  jal     func_80A35D48              
/* 00864 80A35B74 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00868 80A35B78 8FBF0034 */  lw      $ra, 0x0034($sp)           
.L80A35B7C:
/* 0086C 80A35B7C 8FB00030 */  lw      $s0, 0x0030($sp)           
/* 00870 80A35B80 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 00874 80A35B84 03E00008 */  jr      $ra                        
/* 00878 80A35B88 00000000 */  nop


