glabel func_8083BA90
/* 09880 8083BA90 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 09884 8083BA94 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 09888 8083BA98 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 0988C 8083BA9C 00A08025 */  or      $s0, $a1, $zero            ## $s0 = 00000000
/* 09890 8083BAA0 AFA40020 */  sw      $a0, 0x0020($sp)           
/* 09894 8083BAA4 0C20DE52 */  jal     func_80837948              
/* 09898 8083BAA8 AFA7002C */  sw      $a3, 0x002C($sp)           
/* 0989C 8083BAAC 3C068084 */  lui     $a2, %hi(func_80844AF4)    ## $a2 = 80840000
/* 098A0 8083BAB0 24C64AF4 */  addiu   $a2, $a2, %lo(func_80844AF4) ## $a2 = 80844AF4
/* 098A4 8083BAB4 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 098A8 8083BAB8 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 098AC 8083BABC 0C20D716 */  jal     func_80835C58              
/* 098B0 8083BAC0 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 098B4 8083BAC4 920E0692 */  lbu     $t6, 0x0692($s0)           ## 00000692
/* 098B8 8083BAC8 861800B6 */  lh      $t8, 0x00B6($s0)           ## 000000B6
/* 098BC 8083BACC 96190088 */  lhu     $t9, 0x0088($s0)           ## 00000088
/* 098C0 8083BAD0 35CF0002 */  ori     $t7, $t6, 0x0002           ## $t7 = 00000002
/* 098C4 8083BAD4 A20F0692 */  sb      $t7, 0x0692($s0)           ## 00000692
/* 098C8 8083BAD8 A618083C */  sh      $t8, 0x083C($s0)           ## 0000083C
/* 098CC 8083BADC C7A4002C */  lwc1    $f4, 0x002C($sp)           
/* 098D0 8083BAE0 3328FFFE */  andi    $t0, $t9, 0xFFFE           ## $t0 = 00000000
/* 098D4 8083BAE4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 098D8 8083BAE8 E6040838 */  swc1    $f4, 0x0838($s0)           ## 00000838
/* 098DC 8083BAEC C7A60030 */  lwc1    $f6, 0x0030($sp)           
/* 098E0 8083BAF0 A6080088 */  sh      $t0, 0x0088($s0)           ## 00000088
/* 098E4 8083BAF4 A2000893 */  sb      $zero, 0x0893($s0)         ## 00000893
/* 098E8 8083BAF8 0C20CA15 */  jal     func_80832854              
/* 098EC 8083BAFC E6060060 */  swc1    $f6, 0x0060($s0)           ## 00000060
/* 098F0 8083BB00 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 098F4 8083BB04 0C20C9A6 */  jal     func_80832698              
/* 098F8 8083BB08 24056801 */  addiu   $a1, $zero, 0x6801         ## $a1 = 00006801
/* 098FC 8083BB0C 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 09900 8083BB10 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 09904 8083BB14 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 09908 8083BB18 03E00008 */  jr      $ra                        
/* 0990C 8083BB1C 00000000 */  nop


