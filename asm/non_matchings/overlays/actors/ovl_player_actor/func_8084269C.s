glabel func_8084269C
/* 1048C 8084269C 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 10490 808426A0 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 10494 808426A4 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 10498 808426A8 AFA40038 */  sw      $a0, 0x0038($sp)           
/* 1049C 808426AC 94A2089E */  lhu     $v0, 0x089E($a1)           ## 0000089E
/* 104A0 808426B0 00A08025 */  or      $s0, $a1, $zero            ## $s0 = 00000000
/* 104A4 808426B4 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 104A8 808426B8 50400004 */  beql    $v0, $zero, .L808426CC     
/* 104AC 808426BC C6040080 */  lwc1    $f4, 0x0080($s0)           ## 00000080
/* 104B0 808426C0 54410031 */  bnel    $v0, $at, .L80842788       
/* 104B4 808426C4 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 104B8 808426C8 C6040080 */  lwc1    $f4, 0x0080($s0)           ## 00000080
.L808426CC:
/* 104BC 808426CC C60600D0 */  lwc1    $f6, 0x00D0($s0)           ## 000000D0
/* 104C0 808426D0 3C0140A0 */  lui     $at, 0x40A0                ## $at = 40A00000
/* 104C4 808426D4 44815000 */  mtc1    $at, $f10                  ## $f10 = 5.00
/* 104C8 808426D8 46062201 */  sub.s   $f8, $f4, $f6              
/* 104CC 808426DC 260400CC */  addiu   $a0, $s0, 0x00CC           ## $a0 = 000000CC
/* 104D0 808426E0 27A5002C */  addiu   $a1, $sp, 0x002C           ## $a1 = FFFFFFF4
/* 104D4 808426E4 3C0740E0 */  lui     $a3, 0x40E0                ## $a3 = 40E00000
/* 104D8 808426E8 44064000 */  mfc1    $a2, $f8                   
/* 104DC 808426EC 0C210983 */  jal     func_8084260C              
/* 104E0 808426F0 E7AA0010 */  swc1    $f10, 0x0010($sp)          
/* 104E4 808426F4 3C068085 */  lui     $a2, %hi(D_808545B4)       ## $a2 = 80850000
/* 104E8 808426F8 3C078085 */  lui     $a3, %hi(D_808545C0)       ## $a3 = 80850000
/* 104EC 808426FC 240E0032 */  addiu   $t6, $zero, 0x0032         ## $t6 = 00000032
/* 104F0 80842700 240F001E */  addiu   $t7, $zero, 0x001E         ## $t7 = 0000001E
/* 104F4 80842704 AFAF0014 */  sw      $t7, 0x0014($sp)           
/* 104F8 80842708 AFAE0010 */  sw      $t6, 0x0010($sp)           
/* 104FC 8084270C 24E745C0 */  addiu   $a3, $a3, %lo(D_808545C0)  ## $a3 = 808545C0
/* 10500 80842710 24C645B4 */  addiu   $a2, $a2, %lo(D_808545B4)  ## $a2 = 808545B4
/* 10504 80842714 8FA40038 */  lw      $a0, 0x0038($sp)           
/* 10508 80842718 0C00A1B3 */  jal     func_800286CC              
/* 1050C 8084271C 27A5002C */  addiu   $a1, $sp, 0x002C           ## $a1 = FFFFFFF4
/* 10510 80842720 C6100080 */  lwc1    $f16, 0x0080($s0)          ## 00000080
/* 10514 80842724 C61200DC */  lwc1    $f18, 0x00DC($s0)          ## 000000DC
/* 10518 80842728 3C0140A0 */  lui     $at, 0x40A0                ## $at = 40A00000
/* 1051C 8084272C 44813000 */  mtc1    $at, $f6                   ## $f6 = 5.00
/* 10520 80842730 46128101 */  sub.s   $f4, $f16, $f18            
/* 10524 80842734 260400D8 */  addiu   $a0, $s0, 0x00D8           ## $a0 = 000000D8
/* 10528 80842738 AFA40028 */  sw      $a0, 0x0028($sp)           
/* 1052C 8084273C 27A5002C */  addiu   $a1, $sp, 0x002C           ## $a1 = FFFFFFF4
/* 10530 80842740 44062000 */  mfc1    $a2, $f4                   
/* 10534 80842744 3C0740E0 */  lui     $a3, 0x40E0                ## $a3 = 40E00000
/* 10538 80842748 0C210983 */  jal     func_8084260C              
/* 1053C 8084274C E7A60010 */  swc1    $f6, 0x0010($sp)           
/* 10540 80842750 3C068085 */  lui     $a2, %hi(D_808545B4)       ## $a2 = 80850000
/* 10544 80842754 3C078085 */  lui     $a3, %hi(D_808545C0)       ## $a3 = 80850000
/* 10548 80842758 24180032 */  addiu   $t8, $zero, 0x0032         ## $t8 = 00000032
/* 1054C 8084275C 2419001E */  addiu   $t9, $zero, 0x001E         ## $t9 = 0000001E
/* 10550 80842760 AFB90014 */  sw      $t9, 0x0014($sp)           
/* 10554 80842764 AFB80010 */  sw      $t8, 0x0010($sp)           
/* 10558 80842768 24E745C0 */  addiu   $a3, $a3, %lo(D_808545C0)  ## $a3 = 808545C0
/* 1055C 8084276C 24C645B4 */  addiu   $a2, $a2, %lo(D_808545B4)  ## $a2 = 808545B4
/* 10560 80842770 8FA40038 */  lw      $a0, 0x0038($sp)           
/* 10564 80842774 0C00A1B3 */  jal     func_800286CC              
/* 10568 80842778 8FA50028 */  lw      $a1, 0x0028($sp)           
/* 1056C 8084277C 10000002 */  beq     $zero, $zero, .L80842788   
/* 10570 80842780 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
/* 10574 80842784 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L80842788:
/* 10578 80842788 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 1057C 8084278C 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 10580 80842790 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 10584 80842794 03E00008 */  jr      $ra                        
/* 10588 80842798 00000000 */  nop
