glabel func_8083264C
/* 0043C 8083264C 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00440 80832650 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00444 80832654 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 00448 80832658 AFA60020 */  sw      $a2, 0x0020($sp)           
/* 0044C 8083265C AFA70024 */  sw      $a3, 0x0024($sp)           
/* 00450 80832660 908E0002 */  lbu     $t6, 0x0002($a0)           ## 00000002
/* 00454 80832664 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
/* 00458 80832668 93A5001F */  lbu     $a1, 0x001F($sp)           
/* 0045C 8083266C 15C10006 */  bne     $t6, $at, .L80832688       
/* 00460 80832670 8FAF0028 */  lw      $t7, 0x0028($sp)           
/* 00464 80832674 448F2000 */  mtc1    $t7, $f4                   ## $f4 = 0.00
/* 00468 80832678 93A60023 */  lbu     $a2, 0x0023($sp)           
/* 0046C 8083267C 93A70027 */  lbu     $a3, 0x0027($sp)           
/* 00470 80832680 0C02A800 */  jal     func_800AA000              
/* 00474 80832684 46802320 */  cvt.s.w $f12, $f4                  
.L80832688:
/* 00478 80832688 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 0047C 8083268C 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 00480 80832690 03E00008 */  jr      $ra                        
/* 00484 80832694 00000000 */  nop
