glabel func_80B96560
/* 00490 80B96560 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 00494 80B96564 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 00498 80B96568 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 0049C 80B9656C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 004A0 80B96570 0C010D5B */  jal     func_8004356C              
/* 004A4 80B96574 AFA5003C */  sw      $a1, 0x003C($sp)           
/* 004A8 80B96578 10400027 */  beq     $v0, $zero, .L80B96618     
/* 004AC 80B9657C 8FA3003C */  lw      $v1, 0x003C($sp)           
/* 004B0 80B96580 860E016E */  lh      $t6, 0x016E($s0)           ## 0000016E
/* 004B4 80B96584 5DC0002D */  bgtzl   $t6, .L80B9663C            
/* 004B8 80B96588 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 004BC 80B9658C 860F001C */  lh      $t7, 0x001C($s0)           ## 0000001C
/* 004C0 80B96590 24010007 */  addiu   $at, $zero, 0x0007         ## $at = 00000007
/* 004C4 80B96594 000FC203 */  sra     $t8, $t7,  8               
/* 004C8 80B96598 33190007 */  andi    $t9, $t8, 0x0007           ## $t9 = 00000000
/* 004CC 80B9659C 57210006 */  bnel    $t9, $at, .L80B965B8       
/* 004D0 80B965A0 846807A0 */  lh      $t0, 0x07A0($v1)           ## 000007A0
/* 004D4 80B965A4 0C2E59F0 */  jal     func_80B967C0              
/* 004D8 80B965A8 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 004DC 80B965AC 10000023 */  beq     $zero, $zero, .L80B9663C   
/* 004E0 80B965B0 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 004E4 80B965B4 846807A0 */  lh      $t0, 0x07A0($v1)           ## 000007A0
.L80B965B8:
/* 004E8 80B965B8 24050001 */  addiu   $a1, $zero, 0x0001         ## $a1 = 00000001
/* 004EC 80B965BC 00084880 */  sll     $t1, $t0,  2               
/* 004F0 80B965C0 00695021 */  addu    $t2, $v1, $t1              
/* 004F4 80B965C4 0C024BE2 */  jal     Quake_Add              
/* 004F8 80B965C8 8D440790 */  lw      $a0, 0x0790($t2)           ## 00000790
/* 004FC 80B965CC 00022400 */  sll     $a0, $v0, 16               
/* 00500 80B965D0 AFA20030 */  sw      $v0, 0x0030($sp)           
/* 00504 80B965D4 00042403 */  sra     $a0, $a0, 16               
/* 00508 80B965D8 0C024B6B */  jal     Quake_SetSpeed              
/* 0050C 80B965DC 24052710 */  addiu   $a1, $zero, 0x2710         ## $a1 = 00002710
/* 00510 80B965E0 87A40032 */  lh      $a0, 0x0032($sp)           
/* 00514 80B965E4 24050002 */  addiu   $a1, $zero, 0x0002         ## $a1 = 00000002
/* 00518 80B965E8 00003025 */  or      $a2, $zero, $zero          ## $a2 = 00000000
/* 0051C 80B965EC 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 00520 80B965F0 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 00524 80B965F4 0C024B9C */  jal     Quake_SetQuakeValues              
/* 00528 80B965F8 AFA4002C */  sw      $a0, 0x002C($sp)           
/* 0052C 80B965FC 8FA4002C */  lw      $a0, 0x002C($sp)           
/* 00530 80B96600 0C024B7C */  jal     Quake_SetCountdown              
/* 00534 80B96604 24050014 */  addiu   $a1, $zero, 0x0014         ## $a1 = 00000014
/* 00538 80B96608 0C2E5993 */  jal     func_80B9664C              
/* 0053C 80B9660C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00540 80B96610 1000000A */  beq     $zero, $zero, .L80B9663C   
/* 00544 80B96614 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80B96618:
/* 00548 80B96618 860B001C */  lh      $t3, 0x001C($s0)           ## 0000001C
/* 0054C 80B9661C 3C0F80B9 */  lui     $t7, %hi(D_80B969C0)       ## $t7 = 80B90000
/* 00550 80B96620 000B6203 */  sra     $t4, $t3,  8               
/* 00554 80B96624 318D0007 */  andi    $t5, $t4, 0x0007           ## $t5 = 00000000
/* 00558 80B96628 000D7040 */  sll     $t6, $t5,  1               
/* 0055C 80B9662C 01EE7821 */  addu    $t7, $t7, $t6              
/* 00560 80B96630 85EF69C0 */  lh      $t7, %lo(D_80B969C0)($t7)  
/* 00564 80B96634 A60F016E */  sh      $t7, 0x016E($s0)           ## 0000016E
/* 00568 80B96638 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80B9663C:
/* 0056C 80B9663C 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 00570 80B96640 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 00574 80B96644 03E00008 */  jr      $ra                        
/* 00578 80B96648 00000000 */  nop
