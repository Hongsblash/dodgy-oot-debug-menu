glabel func_80AF2690
/* 00140 80AF2690 8483001C */  lh      $v1, 0x001C($a0)           ## 0000001C
/* 00144 80AF2694 00031A03 */  sra     $v1, $v1,  8               
/* 00148 80AF2698 03E00008 */  jr      $ra                        
/* 0014C 80AF269C 306200FF */  andi    $v0, $v1, 0x00FF           ## $v0 = 00000000


