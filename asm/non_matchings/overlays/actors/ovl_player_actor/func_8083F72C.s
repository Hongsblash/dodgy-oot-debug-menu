glabel func_8083F72C
/* 0D51C 8083F72C 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 0D520 8083F730 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 0D524 8083F734 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0D528 8083F738 AFA60028 */  sw      $a2, 0x0028($sp)           
/* 0D52C 8083F73C 00C02025 */  or      $a0, $a2, $zero            ## $a0 = 00000000
/* 0D530 8083F740 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 0D534 8083F744 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 0D538 8083F748 3C068084 */  lui     $a2, %hi(func_8083A388)    ## $a2 = 80840000
/* 0D53C 8083F74C 24C6A388 */  addiu   $a2, $a2, %lo(func_8083A388) ## $a2 = 8083A388
/* 0D540 8083F750 0C20DA26 */  jal     func_80836898              
/* 0D544 8083F754 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0D548 8083F758 14400006 */  bne     $v0, $zero, .L8083F774     
/* 0D54C 8083F75C 8FA40028 */  lw      $a0, 0x0028($sp)           
/* 0D550 8083F760 3C068085 */  lui     $a2, %hi(func_8084B78C)    ## $a2 = 80850000
/* 0D554 8083F764 24C6B78C */  addiu   $a2, $a2, %lo(func_8084B78C) ## $a2 = 8084B78C
/* 0D558 8083F768 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0D55C 8083F76C 0C20D716 */  jal     func_80835C58              
/* 0D560 8083F770 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
.L8083F774:
/* 0D564 8083F774 8FA40028 */  lw      $a0, 0x0028($sp)           
/* 0D568 8083F778 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0D56C 8083F77C 0C20C899 */  jal     func_80832264              
/* 0D570 8083F780 8FA60024 */  lw      $a2, 0x0024($sp)           
/* 0D574 8083F784 0C20C889 */  jal     func_80832224              
/* 0D578 8083F788 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 0D57C 8083F78C 8603007E */  lh      $v1, 0x007E($s0)           ## 0000007E
/* 0D580 8083F790 34018000 */  ori     $at, $zero, 0x8000         ## $at = 00008000
/* 0D584 8083F794 00611821 */  addu    $v1, $v1, $at              
/* 0D588 8083F798 00031C00 */  sll     $v1, $v1, 16               
/* 0D58C 8083F79C 00031C03 */  sra     $v1, $v1, 16               
/* 0D590 8083F7A0 A603083C */  sh      $v1, 0x083C($s0)           ## 0000083C
/* 0D594 8083F7A4 A60300B6 */  sh      $v1, 0x00B6($s0)           ## 000000B6
/* 0D598 8083F7A8 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 0D59C 8083F7AC 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 0D5A0 8083F7B0 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 0D5A4 8083F7B4 03E00008 */  jr      $ra                        
/* 0D5A8 8083F7B8 00000000 */  nop


