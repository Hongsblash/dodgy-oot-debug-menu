glabel func_8084CBF4
/* 1A9E4 8084CBF4 44808000 */  mtc1    $zero, $f16                ## $f16 = 0.00
/* 1A9E8 8084CBF8 C4820878 */  lwc1    $f2, 0x0878($a0)           ## 00000878
/* 1A9EC 8084CBFC 44857000 */  mtc1    $a1, $f14                  ## $f14 = 0.00
/* 1A9F0 8084CC00 44866000 */  mtc1    $a2, $f12                  ## $f12 = 0.00
/* 1A9F4 8084CC04 46028032 */  c.eq.s  $f16, $f2                  
/* 1A9F8 8084CC08 00000000 */  nop
/* 1A9FC 8084CC0C 45010020 */  bc1t    .L8084CC90                 
/* 1AA00 8084CC10 00000000 */  nop
/* 1AA04 8084CC14 C48401CC */  lwc1    $f4, 0x01CC($a0)           ## 000001CC
/* 1AA08 8084CC18 4604603E */  c.le.s  $f12, $f4                  
/* 1AA0C 8084CC1C 00000000 */  nop
/* 1AA10 8084CC20 4500001B */  bc1f    .L8084CC90                 
/* 1AA14 8084CC24 00000000 */  nop
/* 1AA18 8084CC28 46001005 */  abs.s   $f0, $f2                   
/* 1AA1C 8084CC2C 4600703C */  c.lt.s  $f14, $f0                  
/* 1AA20 8084CC30 00000000 */  nop
/* 1AA24 8084CC34 45020011 */  bc1fl   .L8084CC7C                 
/* 1AA28 8084CC38 46001306 */  mov.s   $f12, $f2                  
/* 1AA2C 8084CC3C 4602803E */  c.le.s  $f16, $f2                  
/* 1AA30 8084CC40 3C01BF80 */  lui     $at, 0xBF80                ## $at = BF800000
/* 1AA34 8084CC44 45020008 */  bc1fl   .L8084CC68                 
/* 1AA38 8084CC48 44810000 */  mtc1    $at, $f0                   ## $f0 = -1.00
/* 1AA3C 8084CC4C 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 1AA40 8084CC50 44810000 */  mtc1    $at, $f0                   ## $f0 = 1.00
/* 1AA44 8084CC54 00000000 */  nop
/* 1AA48 8084CC58 460E0302 */  mul.s   $f12, $f0, $f14            
/* 1AA4C 8084CC5C 10000008 */  beq     $zero, $zero, .L8084CC80   
/* 1AA50 8084CC60 C4860028 */  lwc1    $f6, 0x0028($a0)           ## 00000028
/* 1AA54 8084CC64 44810000 */  mtc1    $at, $f0                   ## $f0 = 1.00
.L8084CC68:
/* 1AA58 8084CC68 00000000 */  nop
/* 1AA5C 8084CC6C 460E0302 */  mul.s   $f12, $f0, $f14            
/* 1AA60 8084CC70 10000003 */  beq     $zero, $zero, .L8084CC80   
/* 1AA64 8084CC74 C4860028 */  lwc1    $f6, 0x0028($a0)           ## 00000028
/* 1AA68 8084CC78 46001306 */  mov.s   $f12, $f2                  
.L8084CC7C:
/* 1AA6C 8084CC7C C4860028 */  lwc1    $f6, 0x0028($a0)           ## 00000028
.L8084CC80:
/* 1AA70 8084CC80 460C1281 */  sub.s   $f10, $f2, $f12            
/* 1AA74 8084CC84 460C3200 */  add.s   $f8, $f6, $f12             
/* 1AA78 8084CC88 E48A0878 */  swc1    $f10, 0x0878($a0)          ## 00000878
/* 1AA7C 8084CC8C E4880028 */  swc1    $f8, 0x0028($a0)           ## 00000028
.L8084CC90:
/* 1AA80 8084CC90 03E00008 */  jr      $ra                        
/* 1AA84 8084CC94 00000000 */  nop


