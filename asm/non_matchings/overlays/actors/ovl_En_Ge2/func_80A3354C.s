glabel func_80A3354C
/* 0097C 80A3354C 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 00980 80A33550 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 00984 80A33554 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 00988 80A33558 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 0098C 80A3355C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00990 80A33560 E4840068 */  swc1    $f4, 0x0068($a0)           ## 00000068
/* 00994 80A33564 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 00998 80A33568 0C28CBB3 */  jal     func_80A32ECC              
/* 0099C 80A3356C 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 009A0 80A33570 1040000B */  beq     $v0, $zero, .L80A335A0     
/* 009A4 80A33574 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 009A8 80A33578 24050002 */  addiu   $a1, $zero, 0x0002         ## $a1 = 00000002
/* 009AC 80A3357C 0C28CAF4 */  jal     func_80A32BD0              
/* 009B0 80A33580 AFA2002C */  sw      $v0, 0x002C($sp)           
/* 009B4 80A33584 8FA3002C */  lw      $v1, 0x002C($sp)           
/* 009B8 80A33588 860F008A */  lh      $t7, 0x008A($s0)           ## 0000008A
/* 009BC 80A3358C 240E0064 */  addiu   $t6, $zero, 0x0064         ## $t6 = 00000064
/* 009C0 80A33590 A20E0305 */  sb      $t6, 0x0305($s0)           ## 00000305
/* 009C4 80A33594 A2030306 */  sb      $v1, 0x0306($s0)           ## 00000306
/* 009C8 80A33598 1000000D */  beq     $zero, $zero, .L80A335D0   
/* 009CC 80A3359C A60F02F8 */  sh      $t7, 0x02F8($s0)           ## 000002F8
.L80A335A0:
/* 009D0 80A335A0 961802F4 */  lhu     $t8, 0x02F4($s0)           ## 000002F4
/* 009D4 80A335A4 26040032 */  addiu   $a0, $s0, 0x0032           ## $a0 = 00000032
/* 009D8 80A335A8 24060002 */  addiu   $a2, $zero, 0x0002         ## $a2 = 00000002
/* 009DC 80A335AC 33190002 */  andi    $t9, $t8, 0x0002           ## $t9 = 00000000
/* 009E0 80A335B0 13200007 */  beq     $t9, $zero, .L80A335D0     
/* 009E4 80A335B4 24070400 */  addiu   $a3, $zero, 0x0400         ## $a3 = 00000400
/* 009E8 80A335B8 860502F6 */  lh      $a1, 0x02F6($s0)           ## 000002F6
/* 009EC 80A335BC 24080200 */  addiu   $t0, $zero, 0x0200         ## $t0 = 00000200
/* 009F0 80A335C0 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 009F4 80A335C4 AFA80010 */  sw      $t0, 0x0010($sp)           
/* 009F8 80A335C8 86090032 */  lh      $t1, 0x0032($s0)           ## 00000032
/* 009FC 80A335CC A60900B6 */  sh      $t1, 0x00B6($s0)           ## 000000B6
.L80A335D0:
/* 00A00 80A335D0 860A02F6 */  lh      $t2, 0x02F6($s0)           ## 000002F6
/* 00A04 80A335D4 860B00B6 */  lh      $t3, 0x00B6($s0)           ## 000000B6
/* 00A08 80A335D8 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00A0C 80A335DC 554B0004 */  bnel    $t2, $t3, .L80A335F0       
/* 00A10 80A335E0 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 00A14 80A335E4 0C28CAF4 */  jal     func_80A32BD0              
/* 00A18 80A335E8 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 00A1C 80A335EC 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80A335F0:
/* 00A20 80A335F0 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 00A24 80A335F4 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000
/* 00A28 80A335F8 03E00008 */  jr      $ra                        
/* 00A2C 80A335FC 00000000 */  nop


