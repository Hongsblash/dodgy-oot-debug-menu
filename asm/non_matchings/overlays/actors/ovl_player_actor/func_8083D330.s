glabel func_8083D330
/* 0B120 8083D330 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 0B124 8083D334 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 0B128 8083D338 3C060400 */  lui     $a2, 0x0400                ## $a2 = 04000000
/* 0B12C 8083D33C 24C632F0 */  addiu   $a2, $a2, 0x32F0           ## $a2 = 040032F0
/* 0B130 8083D340 0C20C8A1 */  jal     func_80832284              
/* 0B134 8083D344 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 0B138 8083D348 8FA5001C */  lw      $a1, 0x001C($sp)           
/* 0B13C 8083D34C 240E3E80 */  addiu   $t6, $zero, 0x3E80         ## $t6 = 00003E80
/* 0B140 8083D350 240F0001 */  addiu   $t7, $zero, 0x0001         ## $t7 = 00000001
/* 0B144 8083D354 A4AE06C2 */  sh      $t6, 0x06C2($a1)           ## 000006C2
/* 0B148 8083D358 A4AF0850 */  sh      $t7, 0x0850($a1)           ## 00000850
/* 0B14C 8083D35C 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 0B150 8083D360 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 0B154 8083D364 03E00008 */  jr      $ra                        
/* 0B158 8083D368 00000000 */  nop


