glabel func_8080BEF8
/* 081B8 8080BEF8 3C0E8016 */  lui     $t6, 0x8016                ## $t6 = 80160000
/* 081BC 8080BEFC 8DCEFA90 */  lw      $t6, -0x0570($t6)          ## 8015FA90
/* 081C0 8080BF00 3C01439D */  lui     $at, 0x439D                ## $at = 439D0000
/* 081C4 8080BF04 44810000 */  mtc1    $at, $f0                   ## $f0 = 314.00
/* 081C8 8080BF08 85CF0F34 */  lh      $t7, 0x0F34($t6)           ## 80160F34
/* 081CC 8080BF0C 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 081D0 8080BF10 34218000 */  ori     $at, $at, 0x8000           ## $at = 00018000
/* 081D4 8080BF14 448F3000 */  mtc1    $t7, $f6                   ## $f6 = 0.00
/* 081D8 8080BF18 00811021 */  addu    $v0, $a0, $at              
/* 081DC 8080BF1C C4444AC4 */  lwc1    $f4, 0x4AC4($v0)           ## 00004AC4
/* 081E0 8080BF20 46803220 */  cvt.s.w $f8, $f6                   
/* 081E4 8080BF24 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 081E8 8080BF28 00240821 */  addu    $at, $at, $a0              
/* 081EC 8080BF2C 24180026 */  addiu   $t8, $zero, 0x0026         ## $t8 = 00000026
/* 081F0 8080BF30 46082280 */  add.s   $f10, $f4, $f8             
/* 081F4 8080BF34 E42ACAC4 */  swc1    $f10, -0x353C($at)         ## 0001CAC4
/* 081F8 8080BF38 C4504AC4 */  lwc1    $f16, 0x4AC4($v0)          ## 00004AC4
/* 081FC 8080BF3C 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 08200 8080BF40 00240821 */  addu    $at, $at, $a0              
/* 08204 8080BF44 4610003E */  c.le.s  $f0, $f16                  
/* 08208 8080BF48 00000000 */  nop
/* 0820C 8080BF4C 45000005 */  bc1f    .L8080BF64                 
/* 08210 8080BF50 00000000 */  nop
/* 08214 8080BF54 E420CAC4 */  swc1    $f0, -0x353C($at)          ## 0001CAC4
/* 08218 8080BF58 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 0821C 8080BF5C 00240821 */  addu    $at, $at, $a0              
/* 08220 8080BF60 A438CA3E */  sh      $t8, -0x35C2($at)          ## 0001CA3E
.L8080BF64:
/* 08224 8080BF64 03E00008 */  jr      $ra                        
/* 08228 8080BF68 00000000 */  nop
