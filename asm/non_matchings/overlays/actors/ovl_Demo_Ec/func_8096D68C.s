glabel func_8096D68C
/* 001DC 8096D68C 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 001E0 8096D690 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 001E4 8096D694 848E0192 */  lh      $t6, 0x0192($a0)           ## 00000192
/* 001E8 8096D698 00803025 */  or      $a2, $a0, $zero            ## $a2 = 00000000
/* 001EC 8096D69C 24C30192 */  addiu   $v1, $a2, 0x0192           ## $v1 = 00000192
/* 001F0 8096D6A0 15C00003 */  bne     $t6, $zero, .L8096D6B0     
/* 001F4 8096D6A4 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 001F8 8096D6A8 10000005 */  beq     $zero, $zero, .L8096D6C0   
/* 001FC 8096D6AC 24830192 */  addiu   $v1, $a0, 0x0192           ## $v1 = 00000192
.L8096D6B0:
/* 00200 8096D6B0 846F0000 */  lh      $t7, 0x0000($v1)           ## 00000192
/* 00204 8096D6B4 25F8FFFF */  addiu   $t8, $t7, 0xFFFF           ## $t8 = FFFFFFFF
/* 00208 8096D6B8 A4780000 */  sh      $t8, 0x0000($v1)           ## 00000192
/* 0020C 8096D6BC 84620000 */  lh      $v0, 0x0000($v1)           ## 00000192
.L8096D6C0:
/* 00210 8096D6C0 14400008 */  bne     $v0, $zero, .L8096D6E4     
/* 00214 8096D6C4 2404003C */  addiu   $a0, $zero, 0x003C         ## $a0 = 0000003C
/* 00218 8096D6C8 2405003C */  addiu   $a1, $zero, 0x003C         ## $a1 = 0000003C
/* 0021C 8096D6CC AFA3001C */  sw      $v1, 0x001C($sp)           
/* 00220 8096D6D0 0C01DF64 */  jal     Math_Rand_S16Offset
              
/* 00224 8096D6D4 AFA60038 */  sw      $a2, 0x0038($sp)           
/* 00228 8096D6D8 8FA3001C */  lw      $v1, 0x001C($sp)           
/* 0022C 8096D6DC 8FA60038 */  lw      $a2, 0x0038($sp)           
/* 00230 8096D6E0 A4620000 */  sh      $v0, 0x0000($v1)           ## 00000000
.L8096D6E4:
/* 00234 8096D6E4 84790000 */  lh      $t9, 0x0000($v1)           ## 00000000
/* 00238 8096D6E8 24C20190 */  addiu   $v0, $a2, 0x0190           ## $v0 = 00000190
/* 0023C 8096D6EC A4590000 */  sh      $t9, 0x0000($v0)           ## 00000190
/* 00240 8096D6F0 84480000 */  lh      $t0, 0x0000($v0)           ## 00000190
/* 00244 8096D6F4 29010003 */  slti    $at, $t0, 0x0003           
/* 00248 8096D6F8 54200003 */  bnel    $at, $zero, .L8096D708     
/* 0024C 8096D6FC 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 00250 8096D700 A4400000 */  sh      $zero, 0x0000($v0)         ## 00000190
/* 00254 8096D704 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L8096D708:
/* 00258 8096D708 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 0025C 8096D70C 03E00008 */  jr      $ra                        
/* 00260 8096D710 00000000 */  nop


