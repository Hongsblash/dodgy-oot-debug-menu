glabel func_809A9DC0
/* 00740 809A9DC0 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 00744 809A9DC4 00000000 */  nop
/* 00748 809A9DC8 E4800020 */  swc1    $f0, 0x0020($a0)           ## 00000020
/* 0074C 809A9DCC E480001C */  swc1    $f0, 0x001C($a0)           ## 0000001C
/* 00750 809A9DD0 03E00008 */  jr      $ra                        
/* 00754 809A9DD4 E4800018 */  swc1    $f0, 0x0018($a0)           ## 00000018


