glabel func_80808F84
/* 05244 80808F84 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 05248 80808F88 34218000 */  ori     $at, $at, 0x8000           ## $at = 00018000
/* 0524C 80808F8C 00811021 */  addu    $v0, $a0, $at              
/* 05250 80808F90 844E4A9C */  lh      $t6, 0x4A9C($v0)           ## 00004A9C
/* 05254 80808F94 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 05258 80808F98 00240821 */  addu    $at, $at, $a0              
/* 0525C 80808F9C 25CF0019 */  addiu   $t7, $t6, 0x0019           ## $t7 = 00000019
/* 05260 80808FA0 A42FCA9C */  sh      $t7, -0x3564($at)          ## 0001CA9C
/* 05264 80808FA4 84584A9C */  lh      $t8, 0x4A9C($v0)           ## 00004A9C
/* 05268 80808FA8 241900FF */  addiu   $t9, $zero, 0x00FF         ## $t9 = 000000FF
/* 0526C 80808FAC 240B00FF */  addiu   $t3, $zero, 0x00FF         ## $t3 = 000000FF
/* 05270 80808FB0 2B0100FF */  slti    $at, $t8, 0x00FF           
/* 05274 80808FB4 14200004 */  bne     $at, $zero, .L80808FC8     
/* 05278 80808FB8 240C0063 */  addiu   $t4, $zero, 0x0063         ## $t4 = 00000063
/* 0527C 80808FBC 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 05280 80808FC0 00240821 */  addu    $at, $at, $a0              
/* 05284 80808FC4 A439CA9C */  sh      $t9, -0x3564($at)          ## 0001CA9C
.L80808FC8:
/* 05288 80808FC8 84484ABE */  lh      $t0, 0x4ABE($v0)           ## 00004ABE
/* 0528C 80808FCC 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 05290 80808FD0 00240821 */  addu    $at, $at, $a0              
/* 05294 80808FD4 2509FFE2 */  addiu   $t1, $t0, 0xFFE2           ## $t1 = FFFFFFE2
/* 05298 80808FD8 A429CABE */  sh      $t1, -0x3542($at)          ## 0001CABE
/* 0529C 80808FDC 844A4ABE */  lh      $t2, 0x4ABE($v0)           ## 00004ABE
/* 052A0 80808FE0 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 052A4 80808FE4 00240821 */  addu    $at, $at, $a0              
/* 052A8 80808FE8 1D400011 */  bgtz    $t2, .L80809030            
/* 052AC 80808FEC 240D0021 */  addiu   $t5, $zero, 0x0021         ## $t5 = 00000021
/* 052B0 80808FF0 A420CABE */  sh      $zero, -0x3542($at)        ## 0001CABE
/* 052B4 80808FF4 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 052B8 80808FF8 00240821 */  addu    $at, $at, $a0              
/* 052BC 80808FFC A42BCA9C */  sh      $t3, -0x3564($at)          ## 0001CA9C
/* 052C0 80809000 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 052C4 80809004 00240821 */  addu    $at, $at, $a0              
/* 052C8 80809008 A420CAD0 */  sh      $zero, -0x3530($at)        ## 0001CAD0
/* 052CC 8080900C 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 052D0 80809010 00240821 */  addu    $at, $at, $a0              
/* 052D4 80809014 A420CAD2 */  sh      $zero, -0x352E($at)        ## 0001CAD2
/* 052D8 80809018 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 052DC 8080901C 00240821 */  addu    $at, $at, $a0              
/* 052E0 80809020 A42CCAC8 */  sh      $t4, -0x3538($at)          ## 0001CAC8
/* 052E4 80809024 3C010002 */  lui     $at, 0x0002                ## $at = 00020000
/* 052E8 80809028 00240821 */  addu    $at, $at, $a0              
/* 052EC 8080902C A42DCA3E */  sh      $t5, -0x35C2($at)          ## 0001CA3E
.L80809030:
/* 052F0 80809030 03E00008 */  jr      $ra                        
/* 052F4 80809034 00000000 */  nop


