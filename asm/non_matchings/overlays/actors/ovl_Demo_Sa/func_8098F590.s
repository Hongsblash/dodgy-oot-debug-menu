glabel func_8098F590
/* 01130 8098F590 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 01134 8098F594 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 01138 8098F598 3C050601 */  lui     $a1, 0x0601                ## $a1 = 06010000
/* 0113C 8098F59C AFA40020 */  sw      $a0, 0x0020($sp)           
/* 01140 8098F5A0 24A5F580 */  addiu   $a1, $a1, 0xF580           ## $a1 = 0600F580
/* 01144 8098F5A4 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 01148 8098F5A8 24060002 */  addiu   $a2, $zero, 0x0002         ## $a2 = 00000002
/* 0114C 8098F5AC 0C2639DB */  jal     func_8098E76C              
/* 01150 8098F5B0 3C07C100 */  lui     $a3, 0xC100                ## $a3 = C1000000
/* 01154 8098F5B4 8FAF0020 */  lw      $t7, 0x0020($sp)           
/* 01158 8098F5B8 240E000E */  addiu   $t6, $zero, 0x000E         ## $t6 = 0000000E
/* 0115C 8098F5BC ADEE0198 */  sw      $t6, 0x0198($t7)           ## 00000198
/* 01160 8098F5C0 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 01164 8098F5C4 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 01168 8098F5C8 03E00008 */  jr      $ra                        
/* 0116C 8098F5CC 00000000 */  nop


