glabel func_8096BA98
/* 01FC8 8096BA98 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 01FCC 8096BA9C 44813000 */  mtc1    $at, $f6                   ## $f6 = 1.00
/* 01FD0 8096BAA0 C48401A4 */  lwc1    $f4, 0x01A4($a0)           ## 000001A4
/* 01FD4 8096BAA4 3C0E8016 */  lui     $t6, 0x8016                ## $t6 = 80160000
/* 01FD8 8096BAA8 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 01FDC 8096BAAC 46062200 */  add.s   $f8, $f4, $f6              
/* 01FE0 8096BAB0 44819000 */  mtc1    $at, $f18                  ## $f18 = 10.00
/* 01FE4 8096BAB4 240200FF */  addiu   $v0, $zero, 0x00FF         ## $v0 = 000000FF
/* 01FE8 8096BAB8 E48801A4 */  swc1    $f8, 0x01A4($a0)           ## 000001A4
/* 01FEC 8096BABC 8DCEFA90 */  lw      $t6, -0x0570($t6)          ## 8015FA90
/* 01FF0 8096BAC0 C48401A4 */  lwc1    $f4, 0x01A4($a0)           ## 000001A4
/* 01FF4 8096BAC4 85CF1476 */  lh      $t7, 0x1476($t6)           ## 80161476
/* 01FF8 8096BAC8 448F5000 */  mtc1    $t7, $f10                  ## $f10 = 0.00
/* 01FFC 8096BACC 00000000 */  nop
/* 02000 8096BAD0 46805420 */  cvt.s.w $f16, $f10                 
/* 02004 8096BAD4 46128000 */  add.s   $f0, $f16, $f18            
/* 02008 8096BAD8 4604003E */  c.le.s  $f0, $f4                   
/* 0200C 8096BADC 00000000 */  nop
/* 02010 8096BAE0 45020005 */  bc1fl   .L8096BAF8                 
/* 02014 8096BAE4 C48601A4 */  lwc1    $f6, 0x01A4($a0)           ## 000001A4
/* 02018 8096BAE8 AC8201A8 */  sw      $v0, 0x01A8($a0)           ## 000001A8
/* 0201C 8096BAEC 03E00008 */  jr      $ra                        
/* 02020 8096BAF0 A08200C8 */  sb      $v0, 0x00C8($a0)           ## 000000C8
.L8096BAF4:
/* 02024 8096BAF4 C48601A4 */  lwc1    $f6, 0x01A4($a0)           ## 000001A4
.L8096BAF8:
/* 02028 8096BAF8 3C01437F */  lui     $at, 0x437F                ## $at = 437F0000
/* 0202C 8096BAFC 44815000 */  mtc1    $at, $f10                  ## $f10 = 255.00
/* 02030 8096BB00 46003203 */  div.s   $f8, $f6, $f0              
/* 02034 8096BB04 460A4402 */  mul.s   $f16, $f8, $f10            
/* 02038 8096BB08 4600848D */  trunc.w.s $f18, $f16                 
/* 0203C 8096BB0C 44029000 */  mfc1    $v0, $f18                  
/* 02040 8096BB10 00000000 */  nop
/* 02044 8096BB14 AC8201A8 */  sw      $v0, 0x01A8($a0)           ## 000001A8
/* 02048 8096BB18 A08200C8 */  sb      $v0, 0x00C8($a0)           ## 000000C8
/* 0204C 8096BB1C 03E00008 */  jr      $ra                        
/* 02050 8096BB20 00000000 */  nop


