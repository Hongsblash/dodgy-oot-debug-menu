glabel func_80A68660
/* 00000 80A68660 000570C0 */  sll     $t6, $a1,  3               
/* 00004 80A68664 008E1021 */  addu    $v0, $a0, $t6              
/* 00008 80A68668 844F0000 */  lh      $t7, 0x0000($v0)           ## 00000000
/* 0000C 80A6866C 448F2000 */  mtc1    $t7, $f4                   ## $f4 = 0.00
/* 00010 80A68670 00000000 */  nop
/* 00014 80A68674 468021A0 */  cvt.s.w $f6, $f4                   
/* 00018 80A68678 E4C60000 */  swc1    $f6, 0x0000($a2)           ## 00000000
/* 0001C 80A6867C 84580002 */  lh      $t8, 0x0002($v0)           ## 00000002
/* 00020 80A68680 44984000 */  mtc1    $t8, $f8                   ## $f8 = 0.00
/* 00024 80A68684 00000000 */  nop
/* 00028 80A68688 468042A0 */  cvt.s.w $f10, $f8                  
/* 0002C 80A6868C E4CA0004 */  swc1    $f10, 0x0004($a2)          ## 00000004
/* 00030 80A68690 84590004 */  lh      $t9, 0x0004($v0)           ## 00000004
/* 00034 80A68694 44998000 */  mtc1    $t9, $f16                  ## $f16 = 0.00
/* 00038 80A68698 00000000 */  nop
/* 0003C 80A6869C 468084A0 */  cvt.s.w $f18, $f16                 
/* 00040 80A686A0 03E00008 */  jr      $ra                        
/* 00044 80A686A4 E4D20008 */  swc1    $f18, 0x0008($a2)          ## 00000008


