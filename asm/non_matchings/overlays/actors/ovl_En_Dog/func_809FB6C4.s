glabel func_809FB6C4
/* 006F4 809FB6C4 27BDFFB0 */  addiu   $sp, $sp, 0xFFB0           ## $sp = FFFFFFB0
/* 006F8 809FB6C8 3C0F80A0 */  lui     $t7, %hi(D_809FC008)       ## $t7 = 80A00000
/* 006FC 809FB6CC AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 00700 809FB6D0 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 00704 809FB6D4 AFA50054 */  sw      $a1, 0x0054($sp)           
/* 00708 809FB6D8 25EFC008 */  addiu   $t7, $t7, %lo(D_809FC008)  ## $t7 = 809FC008
/* 0070C 809FB6DC 8DF90000 */  lw      $t9, 0x0000($t7)           ## 809FC008
/* 00710 809FB6E0 27AE0044 */  addiu   $t6, $sp, 0x0044           ## $t6 = FFFFFFF4
/* 00714 809FB6E4 8DF80004 */  lw      $t8, 0x0004($t7)           ## 809FC00C
/* 00718 809FB6E8 ADD90000 */  sw      $t9, 0x0000($t6)           ## FFFFFFF4
/* 0071C 809FB6EC 8DF90008 */  lw      $t9, 0x0008($t7)           ## 809FC010
/* 00720 809FB6F0 3C0980A0 */  lui     $t1, %hi(D_809FC014)       ## $t1 = 80A00000
/* 00724 809FB6F4 2529C014 */  addiu   $t1, $t1, %lo(D_809FC014)  ## $t1 = 809FC014
/* 00728 809FB6F8 ADD80004 */  sw      $t8, 0x0004($t6)           ## FFFFFFF8
/* 0072C 809FB6FC ADD90008 */  sw      $t9, 0x0008($t6)           ## FFFFFFFC
/* 00730 809FB700 8D2B0000 */  lw      $t3, 0x0000($t1)           ## 809FC014
/* 00734 809FB704 27A80038 */  addiu   $t0, $sp, 0x0038           ## $t0 = FFFFFFE8
/* 00738 809FB708 8D2A0004 */  lw      $t2, 0x0004($t1)           ## 809FC018
/* 0073C 809FB70C AD0B0000 */  sw      $t3, 0x0000($t0)           ## FFFFFFE8
/* 00740 809FB710 8D2B0008 */  lw      $t3, 0x0008($t1)           ## 809FC01C
/* 00744 809FB714 AD0A0004 */  sw      $t2, 0x0004($t0)           ## FFFFFFEC
/* 00748 809FB718 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0074C 809FB71C AD0B0008 */  sw      $t3, 0x0008($t0)           ## FFFFFFF0
/* 00750 809FB720 0C27ECA7 */  jal     func_809FB29C              
/* 00754 809FB724 8FA50054 */  lw      $a1, 0x0054($sp)           
/* 00758 809FB728 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 0075C 809FB72C 14410003 */  bne     $v0, $at, .L809FB73C       
/* 00760 809FB730 3C0C80A0 */  lui     $t4, %hi(func_809FB940)    ## $t4 = 80A00000
/* 00764 809FB734 258CB940 */  addiu   $t4, $t4, %lo(func_809FB940) ## $t4 = 809FB940
/* 00768 809FB738 AE0C0190 */  sw      $t4, 0x0190($s0)           ## 00000190
.L809FB73C:
/* 0076C 809FB73C 860201EC */  lh      $v0, 0x01EC($s0)           ## 000001EC
/* 00770 809FB740 26040068 */  addiu   $a0, $s0, 0x0068           ## $a0 = 00000068
/* 00774 809FB744 3C063ECC */  lui     $a2, 0x3ECC                ## $a2 = 3ECC0000
/* 00778 809FB748 14400003 */  bne     $v0, $zero, .L809FB758     
/* 0077C 809FB74C 244DFFFF */  addiu   $t5, $v0, 0xFFFF           ## $t5 = FFFFFFFF
/* 00780 809FB750 10000003 */  beq     $zero, $zero, .L809FB760   
/* 00784 809FB754 00001825 */  or      $v1, $zero, $zero          ## $v1 = 00000000
.L809FB758:
/* 00788 809FB758 A60D01EC */  sh      $t5, 0x01EC($s0)           ## 000001EC
/* 0078C 809FB75C 860301EC */  lh      $v1, 0x01EC($s0)           ## 000001EC
.L809FB760:
/* 00790 809FB760 10600028 */  beq     $v1, $zero, .L809FB804     
/* 00794 809FB764 8FAC0054 */  lw      $t4, 0x0054($sp)           
/* 00798 809FB768 860E01F0 */  lh      $t6, 0x01F0($s0)           ## 000001F0
/* 0079C 809FB76C 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 007A0 809FB770 34C6CCCD */  ori     $a2, $a2, 0xCCCD           ## $a2 = 3ECCCCCD
/* 007A4 809FB774 15C00005 */  bne     $t6, $zero, .L809FB78C     
/* 007A8 809FB778 3C073F80 */  lui     $a3, 0x3F80                ## $a3 = 3F800000
/* 007AC 809FB77C 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 007B0 809FB780 44810000 */  mtc1    $at, $f0                   ## $f0 = 1.00
/* 007B4 809FB784 10000005 */  beq     $zero, $zero, .L809FB79C   
/* 007B8 809FB788 44050000 */  mfc1    $a1, $f0                   
.L809FB78C:
/* 007BC 809FB78C 3C014080 */  lui     $at, 0x4080                ## $at = 40800000
/* 007C0 809FB790 44810000 */  mtc1    $at, $f0                   ## $f0 = 4.00
/* 007C4 809FB794 00000000 */  nop
/* 007C8 809FB798 44050000 */  mfc1    $a1, $f0                   
.L809FB79C:
/* 007CC 809FB79C 0C01E0C4 */  jal     Math_SmoothScaleMaxMinF
              
/* 007D0 809FB7A0 E7A40010 */  swc1    $f4, 0x0010($sp)           
/* 007D4 809FB7A4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 007D8 809FB7A8 0C27ECEB */  jal     func_809FB3AC              
/* 007DC 809FB7AC 8FA50054 */  lw      $a1, 0x0054($sp)           
/* 007E0 809FB7B0 8A180030 */  lwl     $t8, 0x0030($s0)           ## 00000030
/* 007E4 809FB7B4 9A180033 */  lwr     $t8, 0x0033($s0)           ## 00000033
/* 007E8 809FB7B8 861901E6 */  lh      $t9, 0x01E6($s0)           ## 000001E6
/* 007EC 809FB7BC 3C028016 */  lui     $v0, 0x8016                ## $v0 = 80160000
/* 007F0 809FB7C0 AA1800B4 */  swl     $t8, 0x00B4($s0)           ## 000000B4
/* 007F4 809FB7C4 BA1800B7 */  swr     $t8, 0x00B7($s0)           ## 000000B7
/* 007F8 809FB7C8 96180034 */  lhu     $t8, 0x0034($s0)           ## 00000034
/* 007FC 809FB7CC 2B210009 */  slti    $at, $t9, 0x0009           
/* 00800 809FB7D0 10200007 */  beq     $at, $zero, .L809FB7F0     
/* 00804 809FB7D4 A61800B8 */  sh      $t8, 0x00B8($s0)           ## 000000B8
/* 00808 809FB7D8 3C028016 */  lui     $v0, 0x8016                ## $v0 = 80160000
/* 0080C 809FB7DC 2442E660 */  addiu   $v0, $v0, 0xE660           ## $v0 = 8015E660
/* 00810 809FB7E0 94481400 */  lhu     $t0, 0x1400($v0)           ## 8015FA60
/* 00814 809FB7E4 35090001 */  ori     $t1, $t0, 0x0001           ## $t1 = 00000001
/* 00818 809FB7E8 10000016 */  beq     $zero, $zero, .L809FB844   
/* 0081C 809FB7EC A4491400 */  sh      $t1, 0x1400($v0)           ## 8015FA60
.L809FB7F0:
/* 00820 809FB7F0 2442E660 */  addiu   $v0, $v0, 0xE660           ## $v0 = 8015CCC0
/* 00824 809FB7F4 944A1400 */  lhu     $t2, 0x1400($v0)           ## 8015E0C0
/* 00828 809FB7F8 314BFFFE */  andi    $t3, $t2, 0xFFFE           ## $t3 = 00000000
/* 0082C 809FB7FC 10000011 */  beq     $zero, $zero, .L809FB844   
/* 00830 809FB800 A44B1400 */  sh      $t3, 0x1400($v0)           ## 8015E0C0
.L809FB804:
/* 00834 809FB804 8D83009C */  lw      $v1, 0x009C($t4)           ## 0000009C
/* 00838 809FB808 24010003 */  addiu   $at, $zero, 0x0003         ## $at = 00000003
/* 0083C 809FB80C 27AE0044 */  addiu   $t6, $sp, 0x0044           ## $t6 = FFFFFFF4
/* 00840 809FB810 0061001B */  divu    $zero, $v1, $at            
/* 00844 809FB814 00001810 */  mfhi    $v1                        
/* 00848 809FB818 00036880 */  sll     $t5, $v1,  2               
/* 0084C 809FB81C 01AE1021 */  addu    $v0, $t5, $t6              
/* 00850 809FB820 8C4F0000 */  lw      $t7, 0x0000($v0)           ## 8015CCC0
/* 00854 809FB824 2404003C */  addiu   $a0, $zero, 0x003C         ## $a0 = 0000003C
/* 00858 809FB828 A60F01F0 */  sh      $t7, 0x01F0($s0)           ## 000001F0
/* 0085C 809FB82C 0C01DF64 */  jal     Math_Rand_S16Offset
              
/* 00860 809FB830 84450002 */  lh      $a1, 0x0002($v0)           ## 8015CCC2
/* 00864 809FB834 3C1880A0 */  lui     $t8, %hi(func_809FB858)    ## $t8 = 80A00000
/* 00868 809FB838 2718B858 */  addiu   $t8, $t8, %lo(func_809FB858) ## $t8 = 809FB858
/* 0086C 809FB83C A60201EC */  sh      $v0, 0x01EC($s0)           ## 000001EC
/* 00870 809FB840 AE180190 */  sw      $t8, 0x0190($s0)           ## 00000190
.L809FB844:
/* 00874 809FB844 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 00878 809FB848 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 0087C 809FB84C 27BD0050 */  addiu   $sp, $sp, 0x0050           ## $sp = 00000000
/* 00880 809FB850 03E00008 */  jr      $ra                        
/* 00884 809FB854 00000000 */  nop


