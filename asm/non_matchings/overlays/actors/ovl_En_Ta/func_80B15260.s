glabel func_80B15260
/* 017C0 80B15260 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 017C4 80B15264 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 017C8 80B15268 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 017CC 80B1526C 0C00BC65 */  jal     func_8002F194              
/* 017D0 80B15270 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 017D4 80B15274 1040000A */  beq     $v0, $zero, .L80B152A0     
/* 017D8 80B15278 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 017DC 80B1527C 8C8F0004 */  lw      $t7, 0x0004($a0)           ## 00000004
/* 017E0 80B15280 3C01FFFE */  lui     $at, 0xFFFE                ## $at = FFFE0000
/* 017E4 80B15284 3C0E80B1 */  lui     $t6, %hi(func_80B15100)    ## $t6 = 80B10000
/* 017E8 80B15288 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = FFFEFFFF
/* 017EC 80B1528C 25CE5100 */  addiu   $t6, $t6, %lo(func_80B15100) ## $t6 = 80B15100
/* 017F0 80B15290 01E1C024 */  and     $t8, $t7, $at              
/* 017F4 80B15294 AC8E025C */  sw      $t6, 0x025C($a0)           ## 0000025C
/* 017F8 80B15298 10000006 */  beq     $zero, $zero, .L80B152B4   
/* 017FC 80B1529C AC980004 */  sw      $t8, 0x0004($a0)           ## 00000004
.L80B152A0:
/* 01800 80B152A0 8FA5001C */  lw      $a1, 0x001C($sp)           
/* 01804 80B152A4 3C06447A */  lui     $a2, 0x447A                ## $a2 = 447A0000
/* 01808 80B152A8 0C00BCB3 */  jal     func_8002F2CC              
/* 0180C 80B152AC AFA40018 */  sw      $a0, 0x0018($sp)           
/* 01810 80B152B0 8FA40018 */  lw      $a0, 0x0018($sp)           
.L80B152B4:
/* 01814 80B152B4 949902E0 */  lhu     $t9, 0x02E0($a0)           ## 000002E0
/* 01818 80B152B8 37280001 */  ori     $t0, $t9, 0x0001           ## $t0 = 00000001
/* 0181C 80B152BC A48802E0 */  sh      $t0, 0x02E0($a0)           ## 000002E0
/* 01820 80B152C0 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 01824 80B152C4 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 01828 80B152C8 03E00008 */  jr      $ra                        
/* 0182C 80B152CC 00000000 */  nop


