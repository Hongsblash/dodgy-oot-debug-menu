glabel func_809991C8
/* 00A48 809991C8 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00A4C 809991CC AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00A50 809991D0 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 00A54 809991D4 94820192 */  lhu     $v0, 0x0192($a0)           ## 00000192
/* 00A58 809991D8 3C05809A */  lui     $a1, %hi(func_80999214)    ## $a1 = 809A0000
/* 00A5C 809991DC 24A59214 */  addiu   $a1, $a1, %lo(func_80999214) ## $a1 = 80999214
/* 00A60 809991E0 10400003 */  beq     $v0, $zero, .L809991F0     
/* 00A64 809991E4 244EFFFF */  addiu   $t6, $v0, 0xFFFF           ## $t6 = FFFFFFFF
/* 00A68 809991E8 10000004 */  beq     $zero, $zero, .L809991FC   
/* 00A6C 809991EC A48E0192 */  sh      $t6, 0x0192($a0)           ## 00000192
.L809991F0:
/* 00A70 809991F0 0C2661E0 */  jal     func_80998780              
/* 00A74 809991F4 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 00A78 809991F8 8FA40018 */  lw      $a0, 0x0018($sp)           
.L809991FC:
/* 00A7C 809991FC 0C266465 */  jal     func_80999194              
/* 00A80 80999200 8FA5001C */  lw      $a1, 0x001C($sp)           
/* 00A84 80999204 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 00A88 80999208 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 00A8C 8099920C 03E00008 */  jr      $ra                        
/* 00A90 80999210 00000000 */  nop
