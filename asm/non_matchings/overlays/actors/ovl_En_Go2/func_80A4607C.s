glabel func_80A4607C
/* 0334C 80A4607C 848E001C */  lh      $t6, 0x001C($a0)           ## 0000001C
/* 03350 80A46080 24010003 */  addiu   $at, $zero, 0x0003         ## $at = 00000003
/* 03354 80A46084 31CF001F */  andi    $t7, $t6, 0x001F           ## $t7 = 00000000
/* 03358 80A46088 15E10006 */  bne     $t7, $at, .L80A460A4       
/* 0335C 80A4608C 00000000 */  nop
/* 03360 80A46090 84980194 */  lh      $t8, 0x0194($a0)           ## 00000194
/* 03364 80A46094 3C1980A4 */  lui     $t9, %hi(func_80A47578)    ## $t9 = 80A40000
/* 03368 80A46098 27397578 */  addiu   $t9, $t9, %lo(func_80A47578) ## $t9 = 80A47578
/* 0336C 80A4609C 17000003 */  bne     $t8, $zero, .L80A460AC     
/* 03370 80A460A0 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L80A460A4:
/* 03374 80A460A4 03E00008 */  jr      $ra                        
/* 03378 80A460A8 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L80A460AC:
/* 0337C 80A460AC AC990190 */  sw      $t9, 0x0190($a0)           ## 00000190
/* 03380 80A460B0 03E00008 */  jr      $ra                        
/* 03384 80A460B4 00000000 */  nop


