glabel BgSpot08Bakudankabe_Draw
/* 004DC 808B07AC 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 004E0 808B07B0 AFA40020 */  sw      $a0, 0x0020($sp)           
/* 004E4 808B07B4 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 004E8 808B07B8 8FA50020 */  lw      $a1, 0x0020($sp)           
/* 004EC 808B07BC AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 004F0 808B07C0 00002025 */  or      $a0, $zero, $zero          ## $a0 = 00000000
/* 004F4 808B07C4 24A50164 */  addiu   $a1, $a1, 0x0164           ## $a1 = 00000164
/* 004F8 808B07C8 0C018A29 */  jal     func_800628A4              
/* 004FC 808B07CC AFA50018 */  sw      $a1, 0x0018($sp)           
/* 00500 808B07D0 24040001 */  addiu   $a0, $zero, 0x0001         ## $a0 = 00000001
/* 00504 808B07D4 0C018A29 */  jal     func_800628A4              
/* 00508 808B07D8 8FA50018 */  lw      $a1, 0x0018($sp)           
/* 0050C 808B07DC 24040002 */  addiu   $a0, $zero, 0x0002         ## $a0 = 00000002
/* 00510 808B07E0 0C018A29 */  jal     func_800628A4              
/* 00514 808B07E4 8FA50018 */  lw      $a1, 0x0018($sp)           
/* 00518 808B07E8 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 0051C 808B07EC 24A53898 */  addiu   $a1, $a1, 0x3898           ## $a1 = 06003898
/* 00520 808B07F0 0C00D498 */  jal     Gfx_DrawDListOpa
              
/* 00524 808B07F4 8FA40024 */  lw      $a0, 0x0024($sp)           
/* 00528 808B07F8 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 0052C 808B07FC 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 00530 808B0800 03E00008 */  jr      $ra                        
/* 00534 808B0804 00000000 */  nop
/* 00538 808B0808 00000000 */  nop
/* 0053C 808B080C 00000000 */  nop

