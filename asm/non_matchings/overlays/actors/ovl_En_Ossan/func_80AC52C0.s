glabel func_80AC52C0
/* 02620 80AC52C0 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 02624 80AC52C4 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 02628 80AC52C8 AFB20020 */  sw      $s2, 0x0020($sp)           
/* 0262C 80AC52CC AFB1001C */  sw      $s1, 0x001C($sp)           
/* 02630 80AC52D0 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 02634 80AC52D4 90AE0252 */  lbu     $t6, 0x0252($a1)           ## 00000252
/* 02638 80AC52D8 00A08825 */  or      $s1, $a1, $zero            ## $s1 = 00000000
/* 0263C 80AC52DC 00809025 */  or      $s2, $a0, $zero            ## $s2 = 00000000
/* 02640 80AC52E0 000E7880 */  sll     $t7, $t6,  2               
/* 02644 80AC52E4 00AFC021 */  addu    $t8, $a1, $t7              
/* 02648 80AC52E8 8F100200 */  lw      $s0, 0x0200($t8)           ## 00000200
/* 0264C 80AC52EC 8E1901BC */  lw      $t9, 0x01BC($s0)           ## 000001BC
/* 02650 80AC52F0 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 02654 80AC52F4 0320F809 */  jalr    $ra, $t9                   
/* 02658 80AC52F8 00000000 */  nop
/* 0265C 80AC52FC 1040000A */  beq     $v0, $zero, .L80AC5328     
/* 02660 80AC5300 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 02664 80AC5304 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 02668 80AC5308 10410013 */  beq     $v0, $at, .L80AC5358       
/* 0266C 80AC530C 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
/* 02670 80AC5310 10410024 */  beq     $v0, $at, .L80AC53A4       
/* 02674 80AC5314 24010004 */  addiu   $at, $zero, 0x0004         ## $at = 00000004
/* 02678 80AC5318 1041002A */  beq     $v0, $at, .L80AC53C4       
/* 0267C 80AC531C 00000000 */  nop
/* 02680 80AC5320 1000002F */  beq     $zero, $zero, .L80AC53E0   
/* 02684 80AC5324 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80AC5328:
/* 02688 80AC5328 0C2B13B2 */  jal     func_80AC4EC8              
/* 0268C 80AC532C 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 02690 80AC5330 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 02694 80AC5334 A2200251 */  sb      $zero, 0x0251($s1)         ## 00000251
/* 02698 80AC5338 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 0269C 80AC533C E62402D0 */  swc1    $f4, 0x02D0($s1)           ## 000002D0
/* 026A0 80AC5340 8E1901AC */  lw      $t9, 0x01AC($s0)           ## 000001AC
/* 026A4 80AC5344 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 026A8 80AC5348 0320F809 */  jalr    $ra, $t9                   
/* 026AC 80AC534C 00000000 */  nop
/* 026B0 80AC5350 10000023 */  beq     $zero, $zero, .L80AC53E0   
/* 026B4 80AC5354 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80AC5358:
/* 026B8 80AC5358 8E1901C0 */  lw      $t9, 0x01C0($s0)           ## 000001C0
/* 026BC 80AC535C 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 026C0 80AC5360 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 026C4 80AC5364 0320F809 */  jalr    $ra, $t9                   
/* 026C8 80AC5368 00000000 */  nop
/* 026CC 80AC536C 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 026D0 80AC5370 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 026D4 80AC5374 0C2B13F8 */  jal     func_80AC4FE0              
/* 026D8 80AC5378 2406009A */  addiu   $a2, $zero, 0x009A         ## $a2 = 0000009A
/* 026DC 80AC537C 44803000 */  mtc1    $zero, $f6                 ## $f6 = 0.00
/* 026E0 80AC5380 A2200251 */  sb      $zero, 0x0251($s1)         ## 00000251
/* 026E4 80AC5384 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 026E8 80AC5388 E62602D0 */  swc1    $f6, 0x02D0($s1)           ## 000002D0
/* 026EC 80AC538C 8E1901AC */  lw      $t9, 0x01AC($s0)           ## 000001AC
/* 026F0 80AC5390 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 026F4 80AC5394 0320F809 */  jalr    $ra, $t9                   
/* 026F8 80AC5398 00000000 */  nop
/* 026FC 80AC539C 10000010 */  beq     $zero, $zero, .L80AC53E0   
/* 02700 80AC53A0 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80AC53A4:
/* 02704 80AC53A4 0C01E221 */  jal     func_80078884              
/* 02708 80AC53A8 24044806 */  addiu   $a0, $zero, 0x4806         ## $a0 = 00004806
/* 0270C 80AC53AC 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 02710 80AC53B0 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 02714 80AC53B4 0C2B13EB */  jal     func_80AC4FAC              
/* 02718 80AC53B8 2406009D */  addiu   $a2, $zero, 0x009D         ## $a2 = 0000009D
/* 0271C 80AC53BC 10000008 */  beq     $zero, $zero, .L80AC53E0   
/* 02720 80AC53C0 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80AC53C4:
/* 02724 80AC53C4 0C01E221 */  jal     func_80078884              
/* 02728 80AC53C8 24044806 */  addiu   $a0, $zero, 0x4806         ## $a0 = 00004806
/* 0272C 80AC53CC 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 02730 80AC53D0 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 02734 80AC53D4 0C2B13EB */  jal     func_80AC4FAC              
/* 02738 80AC53D8 24060085 */  addiu   $a2, $zero, 0x0085         ## $a2 = 00000085
/* 0273C 80AC53DC 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80AC53E0:
/* 02740 80AC53E0 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 02744 80AC53E4 8FB1001C */  lw      $s1, 0x001C($sp)           
/* 02748 80AC53E8 8FB20020 */  lw      $s2, 0x0020($sp)           
/* 0274C 80AC53EC 03E00008 */  jr      $ra                        
/* 02750 80AC53F0 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000


