glabel func_80B5BFE4
/* 00834 80B5BFE4 C482015C */  lwc1    $f2, 0x015C($a0)           ## 0000015C
/* 00838 80B5BFE8 3C01BF80 */  lui     $at, 0xBF80                ## $at = BF800000
/* 0083C 80B5BFEC C4800158 */  lwc1    $f0, 0x0158($a0)           ## 00000158
/* 00840 80B5BFF0 44812000 */  mtc1    $at, $f4                   ## $f4 = -1.00
/* 00844 80B5BFF4 E4820158 */  swc1    $f2, 0x0158($a0)           ## 00000158
/* 00848 80B5BFF8 E4820164 */  swc1    $f2, 0x0164($a0)           ## 00000164
/* 0084C 80B5BFFC E480015C */  swc1    $f0, 0x015C($a0)           ## 0000015C
/* 00850 80B5C000 03E00008 */  jr      $ra                        
/* 00854 80B5C004 E4840168 */  swc1    $f4, 0x0168($a0)           ## 00000168
