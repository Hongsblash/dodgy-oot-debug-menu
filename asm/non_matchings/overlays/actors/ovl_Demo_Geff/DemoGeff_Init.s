.rdata


glabel D_809785E0
    .asciz "[31mDemo_Geff_Actor_ct:arg_dataがおかしい!!!!!!!!!!!!\n[m"
    .balign 4

.text

glabel DemoGeff_Init
/* 0000C 80977E4C 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00010 80977E50 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00014 80977E54 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 00018 80977E58 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 0001C 80977E5C 8482001C */  lh      $v0, 0x001C($a0)           ## 0000001C
/* 00020 80977E60 00802825 */  or      $a1, $a0, $zero            ## $a1 = 00000000
/* 00024 80977E64 3C048098 */  lui     $a0, %hi(D_809785E0)       ## $a0 = 80980000
/* 00028 80977E68 04400002 */  bltz    $v0, .L80977E74            
/* 0002C 80977E6C 28410009 */  slti    $at, $v0, 0x0009           
/* 00030 80977E70 14200007 */  bne     $at, $zero, .L80977E90     
.L80977E74:
/* 00034 80977E74 248485E0 */  addiu   $a0, $a0, %lo(D_809785E0)  ## $a0 = 809785E0
/* 00038 80977E78 0C00084C */  jal     osSyncPrintf
              
/* 0003C 80977E7C AFA50018 */  sw      $a1, 0x0018($sp)           
/* 00040 80977E80 0C00B55C */  jal     Actor_Kill
              
/* 00044 80977E84 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 00048 80977E88 10000004 */  beq     $zero, $zero, .L80977E9C   
/* 0004C 80977E8C 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80977E90:
/* 00050 80977E90 ACA0014C */  sw      $zero, 0x014C($a1)         ## 0000014C
/* 00054 80977E94 ACA00150 */  sw      $zero, 0x0150($a1)         ## 00000150
/* 00058 80977E98 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80977E9C:
/* 0005C 80977E9C 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 00060 80977EA0 03E00008 */  jr      $ra                        
/* 00064 80977EA4 00000000 */  nop


