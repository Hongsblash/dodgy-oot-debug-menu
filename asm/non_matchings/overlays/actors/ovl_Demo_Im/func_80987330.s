glabel func_80987330
/* 02750 80987330 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 02754 80987334 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 02758 80987338 0C261406 */  jal     func_80985018              
/* 0275C 8098733C AFA40020 */  sw      $a0, 0x0020($sp)           
/* 02760 80987340 0C261418 */  jal     func_80985060              
/* 02764 80987344 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 02768 80987348 AFA2001C */  sw      $v0, 0x001C($sp)           
/* 0276C 8098734C 0C2612F8 */  jal     func_80984BE0              
/* 02770 80987350 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 02774 80987354 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 02778 80987358 0C261C6D */  jal     func_809871B4              
/* 0277C 8098735C 8FA5001C */  lw      $a1, 0x001C($sp)           
/* 02780 80987360 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 02784 80987364 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 02788 80987368 03E00008 */  jr      $ra                        
/* 0278C 8098736C 00000000 */  nop


