glabel func_80BA0DC0
/* 00060 80BA0DC0 C4840060 */  lwc1    $f4, 0x0060($a0)           ## 00000060
/* 00064 80BA0DC4 C486006C */  lwc1    $f6, 0x006C($a0)           ## 0000006C
/* 00068 80BA0DC8 C4800070 */  lwc1    $f0, 0x0070($a0)           ## 00000070
/* 0006C 80BA0DCC 46062200 */  add.s   $f8, $f4, $f6              
/* 00070 80BA0DD0 E4880060 */  swc1    $f8, 0x0060($a0)           ## 00000060
/* 00074 80BA0DD4 C48A0060 */  lwc1    $f10, 0x0060($a0)          ## 00000060
/* 00078 80BA0DD8 4600503C */  c.lt.s  $f10, $f0                  
/* 0007C 80BA0DDC 00000000 */  nop
/* 00080 80BA0DE0 45000002 */  bc1f    .L80BA0DEC                 
/* 00084 80BA0DE4 00000000 */  nop
/* 00088 80BA0DE8 E4800060 */  swc1    $f0, 0x0060($a0)           ## 00000060
.L80BA0DEC:
/* 0008C 80BA0DEC 03E00008 */  jr      $ra                        
/* 00090 80BA0DF0 00000000 */  nop


