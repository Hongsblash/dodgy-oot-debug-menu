glabel func_80B32434
/* 001A4 80B32434 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 001A8 80B32438 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 001AC 80B3243C 00803825 */  or      $a3, $a0, $zero            ## $a3 = 00000000
/* 001B0 80B32440 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 001B4 80B32444 24A50FC0 */  addiu   $a1, $a1, 0x0FC0           ## $a1 = 06000FC0
/* 001B8 80B32448 AFA70018 */  sw      $a3, 0x0018($sp)           
/* 001BC 80B3244C 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 001C0 80B32450 0C0294D3 */  jal     func_800A534C              
/* 001C4 80B32454 3C06C0A0 */  lui     $a2, 0xC0A0                ## $a2 = C0A00000
/* 001C8 80B32458 8FA70018 */  lw      $a3, 0x0018($sp)           
/* 001CC 80B3245C 3C0140A0 */  lui     $at, 0x40A0                ## $at = 40A00000
/* 001D0 80B32460 44812000 */  mtc1    $at, $f4                   ## $f4 = 5.00
/* 001D4 80B32464 90EE0294 */  lbu     $t6, 0x0294($a3)           ## 00000294
/* 001D8 80B32468 3C1880B3 */  lui     $t8, %hi(func_80B32C2C)    ## $t8 = 80B30000
/* 001DC 80B3246C 27182C2C */  addiu   $t8, $t8, %lo(func_80B32C2C) ## $t8 = 80B32C2C
/* 001E0 80B32470 35CF0001 */  ori     $t7, $t6, 0x0001           ## $t7 = 00000001
/* 001E4 80B32474 A0EF0294 */  sb      $t7, 0x0294($a3)           ## 00000294
/* 001E8 80B32478 A4E00194 */  sh      $zero, 0x0194($a3)         ## 00000194
/* 001EC 80B3247C ACF80190 */  sw      $t8, 0x0190($a3)           ## 00000190
/* 001F0 80B32480 E4E40068 */  swc1    $f4, 0x0068($a3)           ## 00000068
/* 001F4 80B32484 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 001F8 80B32488 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 001FC 80B3248C 03E00008 */  jr      $ra                        
/* 00200 80B32490 00000000 */  nop


