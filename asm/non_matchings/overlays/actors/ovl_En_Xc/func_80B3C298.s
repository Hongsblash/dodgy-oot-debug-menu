glabel func_80B3C298
/* 000B8 80B3C298 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 000BC 80B3C29C AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 000C0 80B3C2A0 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 000C4 80B3C2A4 8CA21C44 */  lw      $v0, 0x1C44($a1)           ## 00001C44
/* 000C8 80B3C2A8 3C038016 */  lui     $v1, 0x8016                ## $v1 = 80160000
/* 000CC 80B3C2AC 2463FA90 */  addiu   $v1, $v1, 0xFA90           ## $v1 = 8015FA90
/* 000D0 80B3C2B0 8C580024 */  lw      $t8, 0x0024($v0)           ## 00000024
/* 000D4 80B3C2B4 3C014040 */  lui     $at, 0x4040                ## $at = 40400000
/* 000D8 80B3C2B8 44814000 */  mtc1    $at, $f8                   ## $f8 = 3.00
/* 000DC 80B3C2BC AC98032C */  sw      $t8, 0x032C($a0)           ## 0000032C
/* 000E0 80B3C2C0 8C4F0028 */  lw      $t7, 0x0028($v0)           ## 00000028
/* 000E4 80B3C2C4 24850314 */  addiu   $a1, $a0, 0x0314           ## $a1 = 00000314
/* 000E8 80B3C2C8 24070002 */  addiu   $a3, $zero, 0x0002         ## $a3 = 00000002
/* 000EC 80B3C2CC AC8F0330 */  sw      $t7, 0x0330($a0)           ## 00000330
/* 000F0 80B3C2D0 8C58002C */  lw      $t8, 0x002C($v0)           ## 0000002C
/* 000F4 80B3C2D4 AC980334 */  sw      $t8, 0x0334($a0)           ## 00000334
/* 000F8 80B3C2D8 8C790000 */  lw      $t9, 0x0000($v1)           ## 8015FA90
/* 000FC 80B3C2DC 87281474 */  lh      $t0, 0x1474($t9)           ## 00001474
/* 00100 80B3C2E0 44882000 */  mtc1    $t0, $f4                   ## $f4 = 0.00
/* 00104 80B3C2E4 00000000 */  nop
/* 00108 80B3C2E8 468021A0 */  cvt.s.w $f6, $f4                   
/* 0010C 80B3C2EC 46083281 */  sub.s   $f10, $f6, $f8             
/* 00110 80B3C2F0 E48A0328 */  swc1    $f10, 0x0328($a0)          ## 00000328
/* 00114 80B3C2F4 8C690000 */  lw      $t1, 0x0000($v1)           ## 8015FA90
/* 00118 80B3C2F8 85261476 */  lh      $a2, 0x1476($t1)           ## 00001476
/* 0011C 80B3C2FC 24C6000C */  addiu   $a2, $a2, 0x000C           ## $a2 = 0000000C
/* 00120 80B3C300 00063400 */  sll     $a2, $a2, 16               
/* 00124 80B3C304 0C00D285 */  jal     func_80034A14              
/* 00128 80B3C308 00063403 */  sra     $a2, $a2, 16               
/* 0012C 80B3C30C 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 00130 80B3C310 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 00134 80B3C314 03E00008 */  jr      $ra                        
/* 00138 80B3C318 00000000 */  nop


