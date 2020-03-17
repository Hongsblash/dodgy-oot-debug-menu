glabel func_80B16854
/* 02DB4 80B16854 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 02DB8 80B16858 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 02DBC 80B1685C AFB00018 */  sw      $s0, 0x0018($sp)           
/* 02DC0 80B16860 848202E2 */  lh      $v0, 0x02E2($a0)           ## 000002E2
/* 02DC4 80B16864 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 02DC8 80B16868 18400003 */  blez    $v0, .L80B16878            
/* 02DCC 80B1686C 244EFFFF */  addiu   $t6, $v0, 0xFFFF           ## $t6 = FFFFFFFF
/* 02DD0 80B16870 10000029 */  beq     $zero, $zero, .L80B16918   
/* 02DD4 80B16874 A48E02E2 */  sh      $t6, 0x02E2($a0)           ## 000002E2
.L80B16878:
/* 02DD8 80B16878 2604014C */  addiu   $a0, $s0, 0x014C           ## $a0 = 0000014C
/* 02DDC 80B1687C 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 02DE0 80B16880 AFA40024 */  sw      $a0, 0x0024($sp)           
/* 02DE4 80B16884 1040000F */  beq     $v0, $zero, .L80B168C4     
/* 02DE8 80B16888 8FA40024 */  lw      $a0, 0x0024($sp)           
/* 02DEC 80B1688C 0C02947A */  jal     func_800A51E8              
/* 02DF0 80B16890 8E0502E4 */  lw      $a1, 0x02E4($s0)           ## 000002E4
/* 02DF4 80B16894 3C0142C8 */  lui     $at, 0x42C8                ## $at = 42C80000
/* 02DF8 80B16898 44816000 */  mtc1    $at, $f12                  ## $f12 = 100.00
/* 02DFC 80B1689C 0C00CFBE */  jal     Math_Rand_ZeroFloat
              
/* 02E00 80B168A0 00000000 */  nop
/* 02E04 80B168A4 3C0142C8 */  lui     $at, 0x42C8                ## $at = 42C80000
/* 02E08 80B168A8 44812000 */  mtc1    $at, $f4                   ## $f4 = 100.00
/* 02E0C 80B168AC 00000000 */  nop
/* 02E10 80B168B0 46040180 */  add.s   $f6, $f0, $f4              
/* 02E14 80B168B4 4600320D */  trunc.w.s $f8, $f6                   
/* 02E18 80B168B8 44184000 */  mfc1    $t8, $f8                   
/* 02E1C 80B168BC 00000000 */  nop
/* 02E20 80B168C0 A61802E2 */  sh      $t8, 0x02E2($s0)           ## 000002E2
.L80B168C4:
/* 02E24 80B168C4 3C0142C0 */  lui     $at, 0x42C0                ## $at = 42C00000
/* 02E28 80B168C8 44815000 */  mtc1    $at, $f10                  ## $f10 = 96.00
/* 02E2C 80B168CC C6000164 */  lwc1    $f0, 0x0164($s0)           ## 00000164
/* 02E30 80B168D0 3C014254 */  lui     $at, 0x4254                ## $at = 42540000
/* 02E34 80B168D4 24080002 */  addiu   $t0, $zero, 0x0002         ## $t0 = 00000002
/* 02E38 80B168D8 460A003C */  c.lt.s  $f0, $f10                  
/* 02E3C 80B168DC 00000000 */  nop
/* 02E40 80B168E0 4502000A */  bc1fl   .L80B1690C                 
/* 02E44 80B168E4 A60802B4 */  sh      $t0, 0x02B4($s0)           ## 000002B4
/* 02E48 80B168E8 44818000 */  mtc1    $at, $f16                  ## $f16 = 53.00
/* 02E4C 80B168EC 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 02E50 80B168F0 4600803E */  c.le.s  $f16, $f0                  
/* 02E54 80B168F4 00000000 */  nop
/* 02E58 80B168F8 45020004 */  bc1fl   .L80B1690C                 
/* 02E5C 80B168FC A60802B4 */  sh      $t0, 0x02B4($s0)           ## 000002B4
/* 02E60 80B16900 10000002 */  beq     $zero, $zero, .L80B1690C   
/* 02E64 80B16904 A61902B4 */  sh      $t9, 0x02B4($s0)           ## 000002B4
/* 02E68 80B16908 A60802B4 */  sh      $t0, 0x02B4($s0)           ## 000002B4
.L80B1690C:
/* 02E6C 80B1690C 960902E0 */  lhu     $t1, 0x02E0($s0)           ## 000002E0
/* 02E70 80B16910 352A0008 */  ori     $t2, $t1, 0x0008           ## $t2 = 00000008
/* 02E74 80B16914 A60A02E0 */  sh      $t2, 0x02E0($s0)           ## 000002E0
.L80B16918:
/* 02E78 80B16918 960B02E0 */  lhu     $t3, 0x02E0($s0)           ## 000002E0
/* 02E7C 80B1691C 356C0004 */  ori     $t4, $t3, 0x0004           ## $t4 = 00000004
/* 02E80 80B16920 A60C02E0 */  sh      $t4, 0x02E0($s0)           ## 000002E0
/* 02E84 80B16924 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 02E88 80B16928 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 02E8C 80B1692C 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 02E90 80B16930 03E00008 */  jr      $ra                        
/* 02E94 80B16934 00000000 */  nop


