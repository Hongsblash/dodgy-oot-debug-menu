glabel func_8098652C
/* 0194C 8098652C 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 01950 80986530 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 01954 80986534 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 01958 80986538 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 0195C 8098653C AFA40020 */  sw      $a0, 0x0020($sp)           
/* 01960 80986540 24A51868 */  addiu   $a1, $a1, 0x1868           ## $a1 = 06001868
/* 01964 80986544 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 01968 80986548 00003025 */  or      $a2, $zero, $zero          ## $a2 = 00000000
/* 0196C 8098654C 0C2614A0 */  jal     func_80985280              
/* 01970 80986550 24070000 */  addiu   $a3, $zero, 0x0000         ## $a3 = 00000000
/* 01974 80986554 8FAF0020 */  lw      $t7, 0x0020($sp)           
/* 01978 80986558 240E000F */  addiu   $t6, $zero, 0x000F         ## $t6 = 0000000F
/* 0197C 8098655C ADEE0260 */  sw      $t6, 0x0260($t7)           ## 00000260
/* 01980 80986560 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 01984 80986564 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 01988 80986568 03E00008 */  jr      $ra                        
/* 0198C 8098656C 00000000 */  nop


