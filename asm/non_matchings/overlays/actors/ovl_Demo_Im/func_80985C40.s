glabel func_80985C40
/* 01060 80985C40 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 01064 80985C44 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 01068 80985C48 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 0106C 80985C4C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 01070 80985C50 0C261406 */  jal     func_80985018              
/* 01074 80985C54 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 01078 80985C58 0C261418 */  jal     func_80985060              
/* 0107C 80985C5C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01080 80985C60 0C2612F8 */  jal     func_80984BE0              
/* 01084 80985C64 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01088 80985C68 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 0108C 80985C6C 0C261678 */  jal     func_809859E0              
/* 01090 80985C70 8FA50024 */  lw      $a1, 0x0024($sp)           
/* 01094 80985C74 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01098 80985C78 0C261323 */  jal     func_80984C8C              
/* 0109C 80985C7C 8FA50024 */  lw      $a1, 0x0024($sp)           
/* 010A0 80985C80 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 010A4 80985C84 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 010A8 80985C88 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 010AC 80985C8C 03E00008 */  jr      $ra                        
/* 010B0 80985C90 00000000 */  nop


