glabel func_809869F8
/* 01E18 809869F8 3C0E8016 */  lui     $t6, 0x8016                ## $t6 = 80160000
/* 01E1C 809869FC 8DCEFA90 */  lw      $t6, -0x0570($t6)          ## 8015FA90
/* 01E20 80986A00 3C0141F0 */  lui     $at, 0x41F0                ## $at = 41F00000
/* 01E24 80986A04 44814000 */  mtc1    $at, $f8                   ## $f8 = 30.00
/* 01E28 80986A08 85CF1474 */  lh      $t7, 0x1474($t6)           ## 80161474
/* 01E2C 80986A0C C4820024 */  lwc1    $f2, 0x0024($a0)           ## 00000024
/* 01E30 80986A10 8CA21C44 */  lw      $v0, 0x1C44($a1)           ## 00001C44
/* 01E34 80986A14 448F2000 */  mtc1    $t7, $f4                   ## $f4 = 0.00
/* 01E38 80986A18 C4400024 */  lwc1    $f0, 0x0024($v0)           ## 00000024
/* 01E3C 80986A1C 468021A0 */  cvt.s.w $f6, $f4                   
/* 01E40 80986A20 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 01E44 80986A24 46083280 */  add.s   $f10, $f6, $f8             
/* 01E48 80986A28 460A1401 */  sub.s   $f16, $f2, $f10            
/* 01E4C 80986A2C 4610003C */  c.lt.s  $f0, $f16                  
/* 01E50 80986A30 00000000 */  nop
/* 01E54 80986A34 45000007 */  bc1f    .L80986A54                 
/* 01E58 80986A38 00000000 */  nop
/* 01E5C 80986A3C 8C980004 */  lw      $t8, 0x0004($a0)           ## 00000004
/* 01E60 80986A40 33190040 */  andi    $t9, $t8, 0x0040           ## $t9 = 00000000
/* 01E64 80986A44 17200003 */  bne     $t9, $zero, .L80986A54     
/* 01E68 80986A48 00000000 */  nop
/* 01E6C 80986A4C 03E00008 */  jr      $ra                        
/* 01E70 80986A50 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L80986A54:
/* 01E74 80986A54 03E00008 */  jr      $ra                        
/* 01E78 80986A58 00000000 */  nop


