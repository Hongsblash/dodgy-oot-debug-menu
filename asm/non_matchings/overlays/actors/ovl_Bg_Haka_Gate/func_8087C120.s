glabel func_8087C120
/* 00330 8087C120 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 00334 8087C124 C4840150 */  lwc1    $f4, 0x0150($a0)           ## 00000150
/* 00338 8087C128 8CA21C44 */  lw      $v0, 0x1C44($a1)           ## 00001C44
/* 0033C 8087C12C 46040032 */  c.eq.s  $f0, $f4                   
/* 00340 8087C130 00000000 */  nop
/* 00344 8087C134 45010006 */  bc1t    .L8087C150                 
/* 00348 8087C138 00000000 */  nop
/* 0034C 8087C13C 8C4E0680 */  lw      $t6, 0x0680($v0)           ## 00000680
/* 00350 8087C140 2401FFEF */  addiu   $at, $zero, 0xFFEF         ## $at = FFFFFFEF
/* 00354 8087C144 01C17824 */  and     $t7, $t6, $at              
/* 00358 8087C148 AC4F0680 */  sw      $t7, 0x0680($v0)           ## 00000680
/* 0035C 8087C14C E4800150 */  swc1    $f0, 0x0150($a0)           ## 00000150
.L8087C150:
/* 00360 8087C150 03E00008 */  jr      $ra                        
/* 00364 8087C154 00000000 */  nop


