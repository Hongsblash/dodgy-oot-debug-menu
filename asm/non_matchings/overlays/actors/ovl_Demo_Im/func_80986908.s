glabel func_80986908
/* 01D28 80986908 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 01D2C 8098690C AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 01D30 80986910 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 01D34 80986914 0C261406 */  jal     func_80985018              
/* 01D38 80986918 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 01D3C 8098691C 0C261418 */  jal     func_80985060              
/* 01D40 80986920 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 01D44 80986924 0C2612F8 */  jal     func_80984BE0              
/* 01D48 80986928 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 01D4C 8098692C 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 01D50 80986930 0C261A03 */  jal     func_8098680C              
/* 01D54 80986934 8FA5001C */  lw      $a1, 0x001C($sp)           
/* 01D58 80986938 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 01D5C 8098693C 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 01D60 80986940 03E00008 */  jr      $ra                        
/* 01D64 80986944 00000000 */  nop


