glabel func_809F49A4
/* 016E4 809F49A4 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 016E8 809F49A8 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 016EC 809F49AC 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 016F0 809F49B0 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 016F4 809F49B4 AFA50034 */  sw      $a1, 0x0034($sp)           
/* 016F8 809F49B8 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 016FC 809F49BC 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 01700 809F49C0 86030258 */  lh      $v1, 0x0258($s0)           ## 00000258
/* 01704 809F49C4 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 01708 809F49C8 14600037 */  bne     $v1, $zero, .L809F4AA8     
/* 0170C 809F49CC 00000000 */  nop
/* 01710 809F49D0 860E0264 */  lh      $t6, 0x0264($s0)           ## 00000264
/* 01714 809F49D4 15C00034 */  bne     $t6, $zero, .L809F4AA8     
/* 01718 809F49D8 00000000 */  nop
/* 0171C 809F49DC C6040288 */  lwc1    $f4, 0x0288($s0)           ## 00000288
/* 01720 809F49E0 C6060024 */  lwc1    $f6, 0x0024($s0)           ## 00000024
/* 01724 809F49E4 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 01728 809F49E8 44811000 */  mtc1    $at, $f2                   ## $f2 = 10.00
/* 0172C 809F49EC 46062301 */  sub.s   $f12, $f4, $f6             
/* 01730 809F49F0 C6080290 */  lwc1    $f8, 0x0290($s0)           ## 00000290
/* 01734 809F49F4 C60A002C */  lwc1    $f10, 0x002C($s0)          ## 0000002C
/* 01738 809F49F8 46006005 */  abs.s   $f0, $f12                  
/* 0173C 809F49FC 460A4381 */  sub.s   $f14, $f8, $f10            
/* 01740 809F4A00 4602003C */  c.lt.s  $f0, $f2                   
/* 01744 809F4A04 00000000 */  nop
/* 01748 809F4A08 45000015 */  bc1f    .L809F4A60                 
/* 0174C 809F4A0C 00000000 */  nop
/* 01750 809F4A10 46007005 */  abs.s   $f0, $f14                  
/* 01754 809F4A14 8FA40034 */  lw      $a0, 0x0034($sp)           
/* 01758 809F4A18 4602003C */  c.lt.s  $f0, $f2                   
/* 0175C 809F4A1C 248420D8 */  addiu   $a0, $a0, 0x20D8           ## $a0 = 000020D8
/* 01760 809F4A20 4500000F */  bc1f    .L809F4A60                 
/* 01764 809F4A24 00000000 */  nop
/* 01768 809F4A28 E7AC002C */  swc1    $f12, 0x002C($sp)          
/* 0176C 809F4A2C 0C042F6F */  jal     func_8010BDBC              
/* 01770 809F4A30 E7AE0028 */  swc1    $f14, 0x0028($sp)          
/* 01774 809F4A34 C7AC002C */  lwc1    $f12, 0x002C($sp)          
/* 01778 809F4A38 10400009 */  beq     $v0, $zero, .L809F4A60     
/* 0177C 809F4A3C C7AE0028 */  lwc1    $f14, 0x0028($sp)          
/* 01780 809F4A40 44808000 */  mtc1    $zero, $f16                ## $f16 = 0.00
/* 01784 809F4A44 3C18809F */  lui     $t8, %hi(func_809F4BA4)    ## $t8 = 809F0000
/* 01788 809F4A48 240F0005 */  addiu   $t7, $zero, 0x0005         ## $t7 = 00000005
/* 0178C 809F4A4C 27184BA4 */  addiu   $t8, $t8, %lo(func_809F4BA4) ## $t8 = 809F4BA4
/* 01790 809F4A50 A60F0274 */  sh      $t7, 0x0274($s0)           ## 00000274
/* 01794 809F4A54 AE180214 */  sw      $t8, 0x0214($s0)           ## 00000214
/* 01798 809F4A58 1000004D */  beq     $zero, $zero, .L809F4B90   
/* 0179C 809F4A5C E6100068 */  swc1    $f16, 0x0068($s0)          ## 00000068
.L809F4A60:
/* 017A0 809F4A60 0C03F494 */  jal     Math_atan2f              
/* 017A4 809F4A64 00000000 */  nop
/* 017A8 809F4A68 3C01809F */  lui     $at, %hi(D_809F6050)       ## $at = 809F0000
/* 017AC 809F4A6C C4326050 */  lwc1    $f18, %lo(D_809F6050)($at) 
/* 017B0 809F4A70 260400B6 */  addiu   $a0, $s0, 0x00B6           ## $a0 = 000000B6
/* 017B4 809F4A74 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 017B8 809F4A78 46120102 */  mul.s   $f4, $f0, $f18             
/* 017BC 809F4A7C 24070BB8 */  addiu   $a3, $zero, 0x0BB8         ## $a3 = 00000BB8
/* 017C0 809F4A80 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 017C4 809F4A84 4600218D */  trunc.w.s $f6, $f4                   
/* 017C8 809F4A88 44053000 */  mfc1    $a1, $f6                   
/* 017CC 809F4A8C 00000000 */  nop
/* 017D0 809F4A90 00052C00 */  sll     $a1, $a1, 16               
/* 017D4 809F4A94 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 017D8 809F4A98 00052C03 */  sra     $a1, $a1, 16               
/* 017DC 809F4A9C 860800B6 */  lh      $t0, 0x00B6($s0)           ## 000000B6
/* 017E0 809F4AA0 10000021 */  beq     $zero, $zero, .L809F4B28   
/* 017E4 809F4AA4 A6080032 */  sh      $t0, 0x0032($s0)           ## 00000032
.L809F4AA8:
/* 017E8 809F4AA8 14610013 */  bne     $v1, $at, .L809F4AF8       
/* 017EC 809F4AAC 3C0141A0 */  lui     $at, 0x41A0                ## $at = 41A00000
/* 017F0 809F4AB0 44816000 */  mtc1    $at, $f12                  ## $f12 = 20.00
/* 017F4 809F4AB4 0C00CFBE */  jal     Math_Rand_ZeroFloat
              
/* 017F8 809F4AB8 00000000 */  nop
/* 017FC 809F4ABC 4600020D */  trunc.w.s $f8, $f0                   
/* 01800 809F4AC0 3C0141A0 */  lui     $at, 0x41A0                ## $at = 41A00000
/* 01804 809F4AC4 44819000 */  mtc1    $at, $f18                  ## $f18 = 20.00
/* 01808 809F4AC8 440A4000 */  mfc1    $t2, $f8                   
/* 0180C 809F4ACC 00000000 */  nop
/* 01810 809F4AD0 000A5C00 */  sll     $t3, $t2, 16               
/* 01814 809F4AD4 000B6403 */  sra     $t4, $t3, 16               
/* 01818 809F4AD8 448C5000 */  mtc1    $t4, $f10                  ## $f10 = 0.00
/* 0181C 809F4ADC 00000000 */  nop
/* 01820 809F4AE0 46805420 */  cvt.s.w $f16, $f10                 
/* 01824 809F4AE4 46128100 */  add.s   $f4, $f16, $f18            
/* 01828 809F4AE8 4600218D */  trunc.w.s $f6, $f4                   
/* 0182C 809F4AEC 440E3000 */  mfc1    $t6, $f6                   
/* 01830 809F4AF0 00000000 */  nop
/* 01834 809F4AF4 A60E0264 */  sh      $t6, 0x0264($s0)           ## 00000264
.L809F4AF8:
/* 01838 809F4AF8 8605008A */  lh      $a1, 0x008A($s0)           ## 0000008A
/* 0183C 809F4AFC AFA00010 */  sw      $zero, 0x0010($sp)         
/* 01840 809F4B00 26040032 */  addiu   $a0, $s0, 0x0032           ## $a0 = 00000032
/* 01844 809F4B04 24060014 */  addiu   $a2, $zero, 0x0014         ## $a2 = 00000014
/* 01848 809F4B08 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 0184C 809F4B0C 24071388 */  addiu   $a3, $zero, 0x1388         ## $a3 = 00001388
/* 01850 809F4B10 8605008A */  lh      $a1, 0x008A($s0)           ## 0000008A
/* 01854 809F4B14 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 01858 809F4B18 260400B6 */  addiu   $a0, $s0, 0x00B6           ## $a0 = 000000B6
/* 0185C 809F4B1C 24060003 */  addiu   $a2, $zero, 0x0003         ## $a2 = 00000003
/* 01860 809F4B20 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 01864 809F4B24 24071388 */  addiu   $a3, $zero, 0x1388         ## $a3 = 00001388
.L809F4B28:
/* 01868 809F4B28 8602025E */  lh      $v0, 0x025E($s0)           ## 0000025E
/* 0186C 809F4B2C 1440000B */  bne     $v0, $zero, .L809F4B5C     
/* 01870 809F4B30 30480003 */  andi    $t0, $v0, 0x0003           ## $t0 = 00000000
/* 01874 809F4B34 8618026A */  lh      $t8, 0x026A($s0)           ## 0000026A
/* 01878 809F4B38 240F0014 */  addiu   $t7, $zero, 0x0014         ## $t7 = 00000014
/* 0187C 809F4B3C A60F025E */  sh      $t7, 0x025E($s0)           ## 0000025E
/* 01880 809F4B40 33190001 */  andi    $t9, $t8, 0x0001           ## $t9 = 00000000
/* 01884 809F4B44 17200009 */  bne     $t9, $zero, .L809F4B6C     
/* 01888 809F4B48 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 0188C 809F4B4C 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 01890 809F4B50 24053880 */  addiu   $a1, $zero, 0x3880         ## $a1 = 00003880
/* 01894 809F4B54 10000006 */  beq     $zero, $zero, .L809F4B70   
/* 01898 809F4B58 96020088 */  lhu     $v0, 0x0088($s0)           ## 00000088
.L809F4B5C:
/* 0189C 809F4B5C 15000003 */  bne     $t0, $zero, .L809F4B6C     
/* 018A0 809F4B60 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 018A4 809F4B64 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 018A8 809F4B68 2405387F */  addiu   $a1, $zero, 0x387F         ## $a1 = 0000387F
.L809F4B6C:
/* 018AC 809F4B6C 96020088 */  lhu     $v0, 0x0088($s0)           ## 00000088
.L809F4B70:
/* 018B0 809F4B70 30490008 */  andi    $t1, $v0, 0x0008           ## $t1 = 00000000
/* 018B4 809F4B74 11200006 */  beq     $t1, $zero, .L809F4B90     
/* 018B8 809F4B78 304A0001 */  andi    $t2, $v0, 0x0001           ## $t2 = 00000000
/* 018BC 809F4B7C 11400004 */  beq     $t2, $zero, .L809F4B90     
/* 018C0 809F4B80 3C0140F0 */  lui     $at, 0x40F0                ## $at = 40F00000
/* 018C4 809F4B84 44814000 */  mtc1    $at, $f8                   ## $f8 = 7.50
/* 018C8 809F4B88 00000000 */  nop
/* 018CC 809F4B8C E6080060 */  swc1    $f8, 0x0060($s0)           ## 00000060
.L809F4B90:
/* 018D0 809F4B90 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 018D4 809F4B94 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 018D8 809F4B98 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000
/* 018DC 809F4B9C 03E00008 */  jr      $ra                        
/* 018E0 809F4BA0 00000000 */  nop


