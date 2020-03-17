glabel func_80B1C00C
/* 0070C 80B1C00C 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00710 80B1C010 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00714 80B1C014 84820224 */  lh      $v0, 0x0224($a0)           ## 00000224
/* 00718 80B1C018 00803025 */  or      $a2, $a0, $zero            ## $a2 = 00000000
/* 0071C 80B1C01C 14400003 */  bne     $v0, $zero, .L80B1C02C     
/* 00720 80B1C020 244EFFFF */  addiu   $t6, $v0, 0xFFFF           ## $t6 = FFFFFFFF
/* 00724 80B1C024 10000003 */  beq     $zero, $zero, .L80B1C034   
/* 00728 80B1C028 00001825 */  or      $v1, $zero, $zero          ## $v1 = 00000000
.L80B1C02C:
/* 0072C 80B1C02C A4CE0224 */  sh      $t6, 0x0224($a2)           ## 00000224
/* 00730 80B1C030 84C30224 */  lh      $v1, 0x0224($a2)           ## 00000224
.L80B1C034:
/* 00734 80B1C034 54600023 */  bnel    $v1, $zero, .L80B1C0C4     
/* 00738 80B1C038 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 0073C 80B1C03C 84CF0222 */  lh      $t7, 0x0222($a2)           ## 00000222
/* 00740 80B1C040 25F80001 */  addiu   $t8, $t7, 0x0001           ## $t8 = 00000001
/* 00744 80B1C044 A4D80222 */  sh      $t8, 0x0222($a2)           ## 00000222
/* 00748 80B1C048 84D90222 */  lh      $t9, 0x0222($a2)           ## 00000222
/* 0074C 80B1C04C 2B210003 */  slti    $at, $t9, 0x0003           
/* 00750 80B1C050 5420001C */  bnel    $at, $zero, .L80B1C0C4     
/* 00754 80B1C054 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 00758 80B1C058 84C80218 */  lh      $t0, 0x0218($a2)           ## 00000218
/* 0075C 80B1C05C 2404001E */  addiu   $a0, $zero, 0x001E         ## $a0 = 0000001E
/* 00760 80B1C060 2405001E */  addiu   $a1, $zero, 0x001E         ## $a1 = 0000001E
/* 00764 80B1C064 2509FFFF */  addiu   $t1, $t0, 0xFFFF           ## $t1 = FFFFFFFF
/* 00768 80B1C068 A4C90218 */  sh      $t1, 0x0218($a2)           ## 00000218
/* 0076C 80B1C06C 84CA0218 */  lh      $t2, 0x0218($a2)           ## 00000218
/* 00770 80B1C070 05430013 */  bgezl   $t2, .L80B1C0C0            
/* 00774 80B1C074 A4C00222 */  sh      $zero, 0x0222($a2)         ## 00000222
/* 00778 80B1C078 0C01DF64 */  jal     Math_Rand_S16Offset
              
/* 0077C 80B1C07C AFA60018 */  sw      $a2, 0x0018($sp)           
/* 00780 80B1C080 8FA60018 */  lw      $a2, 0x0018($sp)           
/* 00784 80B1C084 240B0002 */  addiu   $t3, $zero, 0x0002         ## $t3 = 00000002
/* 00788 80B1C088 A4C20224 */  sh      $v0, 0x0224($a2)           ## 00000224
/* 0078C 80B1C08C 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 00790 80B1C090 A4CB0218 */  sh      $t3, 0x0218($a2)           ## 00000218
/* 00794 80B1C094 3C013F00 */  lui     $at, 0x3F00                ## $at = 3F000000
/* 00798 80B1C098 44812000 */  mtc1    $at, $f4                   ## $f4 = 0.50
/* 0079C 80B1C09C 8FA60018 */  lw      $a2, 0x0018($sp)           
/* 007A0 80B1C0A0 4600203C */  c.lt.s  $f4, $f0                   
/* 007A4 80B1C0A4 00000000 */  nop
/* 007A8 80B1C0A8 45020005 */  bc1fl   .L80B1C0C0                 
/* 007AC 80B1C0AC A4C00222 */  sh      $zero, 0x0222($a2)         ## 00000222
/* 007B0 80B1C0B0 84CC0218 */  lh      $t4, 0x0218($a2)           ## 00000218
/* 007B4 80B1C0B4 258D0001 */  addiu   $t5, $t4, 0x0001           ## $t5 = 00000001
/* 007B8 80B1C0B8 A4CD0218 */  sh      $t5, 0x0218($a2)           ## 00000218
/* 007BC 80B1C0BC A4C00222 */  sh      $zero, 0x0222($a2)         ## 00000222
.L80B1C0C0:
/* 007C0 80B1C0C0 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80B1C0C4:
/* 007C4 80B1C0C4 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 007C8 80B1C0C8 03E00008 */  jr      $ra                        
/* 007CC 80B1C0CC 00000000 */  nop


