glabel func_8098F8F8
/* 01498 8098F8F8 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 0149C 8098F8FC 44813000 */  mtc1    $at, $f6                   ## $f6 = 1.00
/* 014A0 8098F900 C48401A0 */  lwc1    $f4, 0x01A0($a0)           ## 000001A0
/* 014A4 8098F904 3C0E8016 */  lui     $t6, 0x8016                ## $t6 = 80160000
/* 014A8 8098F908 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 014AC 8098F90C 46062200 */  add.s   $f8, $f4, $f6              
/* 014B0 8098F910 44819000 */  mtc1    $at, $f18                  ## $f18 = 10.00
/* 014B4 8098F914 240200FF */  addiu   $v0, $zero, 0x00FF         ## $v0 = 000000FF
/* 014B8 8098F918 E48801A0 */  swc1    $f8, 0x01A0($a0)           ## 000001A0
/* 014BC 8098F91C 8DCEFA90 */  lw      $t6, -0x0570($t6)          ## 8015FA90
/* 014C0 8098F920 C48401A0 */  lwc1    $f4, 0x01A0($a0)           ## 000001A0
/* 014C4 8098F924 85CF1476 */  lh      $t7, 0x1476($t6)           ## 80161476
/* 014C8 8098F928 448F5000 */  mtc1    $t7, $f10                  ## $f10 = 0.00
/* 014CC 8098F92C 00000000 */  nop
/* 014D0 8098F930 46805420 */  cvt.s.w $f16, $f10                 
/* 014D4 8098F934 46128000 */  add.s   $f0, $f16, $f18            
/* 014D8 8098F938 4604003E */  c.le.s  $f0, $f4                   
/* 014DC 8098F93C 00000000 */  nop
/* 014E0 8098F940 45020005 */  bc1fl   .L8098F958                 
/* 014E4 8098F944 C48601A0 */  lwc1    $f6, 0x01A0($a0)           ## 000001A0
/* 014E8 8098F948 AC8201A4 */  sw      $v0, 0x01A4($a0)           ## 000001A4
/* 014EC 8098F94C 03E00008 */  jr      $ra                        
/* 014F0 8098F950 A08200C8 */  sb      $v0, 0x00C8($a0)           ## 000000C8
.L8098F954:
/* 014F4 8098F954 C48601A0 */  lwc1    $f6, 0x01A0($a0)           ## 000001A0
.L8098F958:
/* 014F8 8098F958 3C01437F */  lui     $at, 0x437F                ## $at = 437F0000
/* 014FC 8098F95C 44815000 */  mtc1    $at, $f10                  ## $f10 = 255.00
/* 01500 8098F960 46003203 */  div.s   $f8, $f6, $f0              
/* 01504 8098F964 460A4402 */  mul.s   $f16, $f8, $f10            
/* 01508 8098F968 4600848D */  trunc.w.s $f18, $f16                 
/* 0150C 8098F96C 44029000 */  mfc1    $v0, $f18                  
/* 01510 8098F970 00000000 */  nop
/* 01514 8098F974 AC8201A4 */  sw      $v0, 0x01A4($a0)           ## 000001A4
/* 01518 8098F978 A08200C8 */  sb      $v0, 0x00C8($a0)           ## 000000C8
/* 0151C 8098F97C 03E00008 */  jr      $ra                        
/* 01520 8098F980 00000000 */  nop
