glabel func_808B39E8
/* 00088 808B39E8 44866000 */  mtc1    $a2, $f12                  ## $f12 = 0.00
/* 0008C 808B39EC C4A40008 */  lwc1    $f4, 0x0008($a1)           ## 00000008
/* 00090 808B39F0 44877000 */  mtc1    $a3, $f14                  ## $f14 = 0.00
/* 00094 808B39F4 C4A80000 */  lwc1    $f8, 0x0000($a1)           ## 00000000
/* 00098 808B39F8 460C2182 */  mul.s   $f6, $f4, $f12             
/* 0009C 808B39FC 00000000 */  nop
/* 000A0 808B3A00 460E4282 */  mul.s   $f10, $f8, $f14            
/* 000A4 808B3A04 460A3400 */  add.s   $f16, $f6, $f10            
/* 000A8 808B3A08 E4900000 */  swc1    $f16, 0x0000($a0)          ## 00000000
/* 000AC 808B3A0C C4B20004 */  lwc1    $f18, 0x0004($a1)          ## 00000004
/* 000B0 808B3A10 E4920004 */  swc1    $f18, 0x0004($a0)          ## 00000004
/* 000B4 808B3A14 C4A40008 */  lwc1    $f4, 0x0008($a1)           ## 00000008
/* 000B8 808B3A18 C4A60000 */  lwc1    $f6, 0x0000($a1)           ## 00000000
/* 000BC 808B3A1C 460E2202 */  mul.s   $f8, $f4, $f14             
/* 000C0 808B3A20 00000000 */  nop
/* 000C4 808B3A24 460C3282 */  mul.s   $f10, $f6, $f12            
/* 000C8 808B3A28 460A4401 */  sub.s   $f16, $f8, $f10            
/* 000CC 808B3A2C 03E00008 */  jr      $ra                        
/* 000D0 808B3A30 E4900008 */  swc1    $f16, 0x0008($a0)          ## 00000008


