glabel EnGe2_Draw
/* 01610 80A341E0 27BDFFA8 */  addiu   $sp, $sp, 0xFFA8           ## $sp = FFFFFFA8
/* 01614 80A341E4 AFB10028 */  sw      $s1, 0x0028($sp)           
/* 01618 80A341E8 00A08825 */  or      $s1, $a1, $zero            ## $s1 = 00000000
/* 0161C 80A341EC AFBF002C */  sw      $ra, 0x002C($sp)           
/* 01620 80A341F0 AFB00024 */  sw      $s0, 0x0024($sp)           
/* 01624 80A341F4 8CA50000 */  lw      $a1, 0x0000($a1)           ## 00000000
/* 01628 80A341F8 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0162C 80A341FC 3C0680A3 */  lui     $a2, %hi(D_80A3441C)       ## $a2 = 80A30000
/* 01630 80A34200 24C6441C */  addiu   $a2, $a2, %lo(D_80A3441C)  ## $a2 = 80A3441C
/* 01634 80A34204 27A4003C */  addiu   $a0, $sp, 0x003C           ## $a0 = FFFFFFE4
/* 01638 80A34208 240704FA */  addiu   $a3, $zero, 0x04FA         ## $a3 = 000004FA
/* 0163C 80A3420C 0C031AB1 */  jal     func_800C6AC4              
/* 01640 80A34210 AFA5004C */  sw      $a1, 0x004C($sp)           
/* 01644 80A34214 0C0250F2 */  jal     func_800943C8              
/* 01648 80A34218 8E240000 */  lw      $a0, 0x0000($s1)           ## 00000000
/* 0164C 80A3421C 8FA5004C */  lw      $a1, 0x004C($sp)           
/* 01650 80A34220 3C0FDB06 */  lui     $t7, 0xDB06                ## $t7 = DB060000
/* 01654 80A34224 35EF0020 */  ori     $t7, $t7, 0x0020           ## $t7 = DB060020
/* 01658 80A34228 8CA302C0 */  lw      $v1, 0x02C0($a1)           ## 000002C0
/* 0165C 80A3422C 3C0480A3 */  lui     $a0, %hi(D_80A343BC)       ## $a0 = 80A30000
/* 01660 80A34230 3C0C8016 */  lui     $t4, 0x8016                ## $t4 = 80160000
/* 01664 80A34234 246E0008 */  addiu   $t6, $v1, 0x0008           ## $t6 = 00000008
/* 01668 80A34238 ACAE02C0 */  sw      $t6, 0x02C0($a1)           ## 000002C0
/* 0166C 80A3423C AC6F0000 */  sw      $t7, 0x0000($v1)           ## 00000000
/* 01670 80A34240 861802E4 */  lh      $t8, 0x02E4($s0)           ## 000002E4
/* 01674 80A34244 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 01678 80A34248 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 0167C 80A3424C 0018C880 */  sll     $t9, $t8,  2               
/* 01680 80A34250 00992021 */  addu    $a0, $a0, $t9              
/* 01684 80A34254 8C8443BC */  lw      $a0, %lo(D_80A343BC)($a0)  
/* 01688 80A34258 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 0168C 80A3425C 00003025 */  or      $a2, $zero, $zero          ## $a2 = 00000000
/* 01690 80A34260 00044900 */  sll     $t1, $a0,  4               
/* 01694 80A34264 00095702 */  srl     $t2, $t1, 28               
/* 01698 80A34268 000A5880 */  sll     $t3, $t2,  2               
/* 0169C 80A3426C 018B6021 */  addu    $t4, $t4, $t3              
/* 016A0 80A34270 8D8C6FA8 */  lw      $t4, 0x6FA8($t4)           ## 80166FA8
/* 016A4 80A34274 00814024 */  and     $t0, $a0, $at              
/* 016A8 80A34278 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 016AC 80A3427C 010C6821 */  addu    $t5, $t0, $t4              
/* 016B0 80A34280 01A17021 */  addu    $t6, $t5, $at              
/* 016B4 80A34284 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 016B8 80A34288 0C00BAF3 */  jal     func_8002EBCC              
/* 016BC 80A3428C AC6E0004 */  sw      $t6, 0x0004($v1)           ## 00000004
/* 016C0 80A34290 8E05019C */  lw      $a1, 0x019C($s0)           ## 0000019C
/* 016C4 80A34294 8E0601B8 */  lw      $a2, 0x01B8($s0)           ## 000001B8
/* 016C8 80A34298 9207019A */  lbu     $a3, 0x019A($s0)           ## 0000019A
/* 016CC 80A3429C 3C0F80A3 */  lui     $t7, %hi(func_80A3415C)    ## $t7 = 80A30000
/* 016D0 80A342A0 3C1880A3 */  lui     $t8, %hi(func_80A341A0)    ## $t8 = 80A30000
/* 016D4 80A342A4 271841A0 */  addiu   $t8, $t8, %lo(func_80A341A0) ## $t8 = 80A341A0
/* 016D8 80A342A8 25EF415C */  addiu   $t7, $t7, %lo(func_80A3415C) ## $t7 = 80A3415C
/* 016DC 80A342AC AFAF0010 */  sw      $t7, 0x0010($sp)           
/* 016E0 80A342B0 AFB80014 */  sw      $t8, 0x0014($sp)           
/* 016E4 80A342B4 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 016E8 80A342B8 0C0286B2 */  jal     func_800A1AC8              
/* 016EC 80A342BC 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 016F0 80A342C0 3C0680A3 */  lui     $a2, %hi(D_80A3442C)       ## $a2 = 80A30000
/* 016F4 80A342C4 24C6442C */  addiu   $a2, $a2, %lo(D_80A3442C)  ## $a2 = 80A3442C
/* 016F8 80A342C8 27A4003C */  addiu   $a0, $sp, 0x003C           ## $a0 = FFFFFFE4
/* 016FC 80A342CC 8E250000 */  lw      $a1, 0x0000($s1)           ## 00000000
/* 01700 80A342D0 0C031AD5 */  jal     func_800C6B54              
/* 01704 80A342D4 2407050B */  addiu   $a3, $zero, 0x050B         ## $a3 = 0000050B
/* 01708 80A342D8 8FBF002C */  lw      $ra, 0x002C($sp)           
/* 0170C 80A342DC 8FB00024 */  lw      $s0, 0x0024($sp)           
/* 01710 80A342E0 8FB10028 */  lw      $s1, 0x0028($sp)           
/* 01714 80A342E4 03E00008 */  jr      $ra                        
/* 01718 80A342E8 27BD0058 */  addiu   $sp, $sp, 0x0058           ## $sp = 00000000
/* 0171C 80A342EC 00000000 */  nop

