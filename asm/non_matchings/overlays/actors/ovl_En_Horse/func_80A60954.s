glabel func_80A60954
/* 05664 80A60954 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 05668 80A60958 AFBF002C */  sw      $ra, 0x002C($sp)           
/* 0566C 80A6095C AFB00028 */  sw      $s0, 0x0028($sp)           
/* 05670 80A60960 AFA5003C */  sw      $a1, 0x003C($sp)           
/* 05674 80A60964 AFA60040 */  sw      $a2, 0x0040($sp)           
/* 05678 80A60968 3C0141C8 */  lui     $at, 0x41C8                ## $at = 41C80000
/* 0567C 80A6096C 44813000 */  mtc1    $at, $f6                   ## $f6 = 25.00
/* 05680 80A60970 C4880214 */  lwc1    $f8, 0x0214($a0)           ## 00000214
/* 05684 80A60974 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 05688 80A60978 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0568C 80A6097C 4608303C */  c.lt.s  $f6, $f8                   
/* 05690 80A60980 E4840068 */  swc1    $f4, 0x0068($a0)           ## 00000068
/* 05694 80A60984 45020011 */  bc1fl   .L80A609CC                 
/* 05698 80A60988 260401AC */  addiu   $a0, $s0, 0x01AC           ## $a0 = 000001AC
/* 0569C 80A6098C 8C8201F0 */  lw      $v0, 0x01F0($a0)           ## 0000039C
/* 056A0 80A60990 3C078013 */  lui     $a3, 0x8013                ## $a3 = 80130000
/* 056A4 80A60994 24E733E0 */  addiu   $a3, $a3, 0x33E0           ## $a3 = 801333E0
/* 056A8 80A60998 304E0800 */  andi    $t6, $v0, 0x0800           ## $t6 = 00000000
/* 056AC 80A6099C 15C0000A */  bne     $t6, $zero, .L80A609C8     
/* 056B0 80A609A0 344F0800 */  ori     $t7, $v0, 0x0800           ## $t7 = 00000800
/* 056B4 80A609A4 AC8F01F0 */  sw      $t7, 0x01F0($a0)           ## 0000039C
/* 056B8 80A609A8 3C188013 */  lui     $t8, 0x8013                ## $t8 = 80130000
/* 056BC 80A609AC 271833E8 */  addiu   $t8, $t8, 0x33E8           ## $t8 = 801333E8
/* 056C0 80A609B0 AFB80014 */  sw      $t8, 0x0014($sp)           
/* 056C4 80A609B4 AFA70010 */  sw      $a3, 0x0010($sp)           
/* 056C8 80A609B8 2404282B */  addiu   $a0, $zero, 0x282B         ## $a0 = 0000282B
/* 056CC 80A609BC 260500E4 */  addiu   $a1, $s0, 0x00E4           ## $a1 = 000000E4
/* 056D0 80A609C0 0C03DCE3 */  jal     Audio_PlaySoundGeneral
              
/* 056D4 80A609C4 24060004 */  addiu   $a2, $zero, 0x0004         ## $a2 = 00000004
.L80A609C8:
/* 056D8 80A609C8 260401AC */  addiu   $a0, $s0, 0x01AC           ## $a0 = 000001AC
.L80A609CC:
/* 056DC 80A609CC 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 056E0 80A609D0 AFA40030 */  sw      $a0, 0x0030($sp)           
/* 056E4 80A609D4 50400045 */  beql    $v0, $zero, .L80A60AEC     
/* 056E8 80A609D8 8FBF002C */  lw      $ra, 0x002C($sp)           
/* 056EC 80A609DC 96020384 */  lhu     $v0, 0x0384($s0)           ## 00000384
/* 056F0 80A609E0 AE000210 */  sw      $zero, 0x0210($s0)         ## 00000210
/* 056F4 80A609E4 30590004 */  andi    $t9, $v0, 0x0004           ## $t9 = 00000000
/* 056F8 80A609E8 17200022 */  bne     $t9, $zero, .L80A60A74     
/* 056FC 80A609EC 34480004 */  ori     $t0, $v0, 0x0004           ## $t0 = 00000004
/* 05700 80A609F0 8E090158 */  lw      $t1, 0x0158($s0)           ## 00000158
/* 05704 80A609F4 3C0B80A6 */  lui     $t3, %hi(D_80A65E58)       ## $t3 = 80A60000
/* 05708 80A609F8 A6080384 */  sh      $t0, 0x0384($s0)           ## 00000384
/* 0570C 80A609FC 00095080 */  sll     $t2, $t1,  2               
/* 05710 80A60A00 016A5821 */  addu    $t3, $t3, $t2              
/* 05714 80A60A04 8D6B5E58 */  lw      $t3, %lo(D_80A65E58)($t3)  
/* 05718 80A60A08 00006880 */  sll     $t5, $zero,  2             
/* 0571C 80A60A0C 016D7021 */  addu    $t6, $t3, $t5              
/* 05720 80A60A10 0C028800 */  jal     SkelAnime_GetFrameCount
              
/* 05724 80A60A14 8DC40000 */  lw      $a0, 0x0000($t6)           ## 00000000
/* 05728 80A60A18 8E0F0158 */  lw      $t7, 0x0158($s0)           ## 00000158
/* 0572C 80A60A1C 44825000 */  mtc1    $v0, $f10                  ## $f10 = 0.00
/* 05730 80A60A20 3C1980A6 */  lui     $t9, %hi(D_80A65E58)       ## $t9 = 80A60000
/* 05734 80A60A24 000FC080 */  sll     $t8, $t7,  2               
/* 05738 80A60A28 8E080210 */  lw      $t0, 0x0210($s0)           ## 00000210
/* 0573C 80A60A2C 0338C821 */  addu    $t9, $t9, $t8              
/* 05740 80A60A30 8F395E58 */  lw      $t9, %lo(D_80A65E58)($t9)  
/* 05744 80A60A34 46805420 */  cvt.s.w $f16, $f10                 
/* 05748 80A60A38 00084880 */  sll     $t1, $t0,  2               
/* 0574C 80A60A3C 3C01C040 */  lui     $at, 0xC040                ## $at = C0400000
/* 05750 80A60A40 03295021 */  addu    $t2, $t9, $t1              
/* 05754 80A60A44 44819000 */  mtc1    $at, $f18                  ## $f18 = -3.00
/* 05758 80A60A48 8D450000 */  lw      $a1, 0x0000($t2)           ## 00000000
/* 0575C 80A60A4C 240C0002 */  addiu   $t4, $zero, 0x0002         ## $t4 = 00000002
/* 05760 80A60A50 AFAC0014 */  sw      $t4, 0x0014($sp)           
/* 05764 80A60A54 E7B00010 */  swc1    $f16, 0x0010($sp)          
/* 05768 80A60A58 8FA40030 */  lw      $a0, 0x0030($sp)           
/* 0576C 80A60A5C 3C063F80 */  lui     $a2, 0x3F80                ## $a2 = 3F800000
/* 05770 80A60A60 24070000 */  addiu   $a3, $zero, 0x0000         ## $a3 = 00000000
/* 05774 80A60A64 0C029468 */  jal     SkelAnime_ChangeAnimation
              
/* 05778 80A60A68 E7B20018 */  swc1    $f18, 0x0018($sp)          
/* 0577C 80A60A6C 1000001F */  beq     $zero, $zero, .L80A60AEC   
/* 05780 80A60A70 8FBF002C */  lw      $ra, 0x002C($sp)           
.L80A60A74:
/* 05784 80A60A74 8E0B0158 */  lw      $t3, 0x0158($s0)           ## 00000158
/* 05788 80A60A78 3C0E80A6 */  lui     $t6, %hi(D_80A65E58)       ## $t6 = 80A60000
/* 0578C 80A60A7C 8E0F0210 */  lw      $t7, 0x0210($s0)           ## 00000210
/* 05790 80A60A80 000B6880 */  sll     $t5, $t3,  2               
/* 05794 80A60A84 01CD7021 */  addu    $t6, $t6, $t5              
/* 05798 80A60A88 8DCE5E58 */  lw      $t6, %lo(D_80A65E58)($t6)  
/* 0579C 80A60A8C 000FC080 */  sll     $t8, $t7,  2               
/* 057A0 80A60A90 01D84021 */  addu    $t0, $t6, $t8              
/* 057A4 80A60A94 0C028800 */  jal     SkelAnime_GetFrameCount
              
/* 057A8 80A60A98 8D040000 */  lw      $a0, 0x0000($t0)           ## 00000000
/* 057AC 80A60A9C 8E190158 */  lw      $t9, 0x0158($s0)           ## 00000158
/* 057B0 80A60AA0 44822000 */  mtc1    $v0, $f4                   ## $f4 = 0.00
/* 057B4 80A60AA4 3C0A80A6 */  lui     $t2, %hi(D_80A65E58)       ## $t2 = 80A60000
/* 057B8 80A60AA8 00194880 */  sll     $t1, $t9,  2               
/* 057BC 80A60AAC 8E0C0210 */  lw      $t4, 0x0210($s0)           ## 00000210
/* 057C0 80A60AB0 01495021 */  addu    $t2, $t2, $t1              
/* 057C4 80A60AB4 8D4A5E58 */  lw      $t2, %lo(D_80A65E58)($t2)  
/* 057C8 80A60AB8 468021A0 */  cvt.s.w $f6, $f4                   
/* 057CC 80A60ABC 000C5880 */  sll     $t3, $t4,  2               
/* 057D0 80A60AC0 014B6821 */  addu    $t5, $t2, $t3              
/* 057D4 80A60AC4 44804000 */  mtc1    $zero, $f8                 ## $f8 = 0.00
/* 057D8 80A60AC8 8DA50000 */  lw      $a1, 0x0000($t5)           ## 00000000
/* 057DC 80A60ACC AFA00014 */  sw      $zero, 0x0014($sp)         
/* 057E0 80A60AD0 E7A60010 */  swc1    $f6, 0x0010($sp)           
/* 057E4 80A60AD4 8FA40030 */  lw      $a0, 0x0030($sp)           
/* 057E8 80A60AD8 3C063F80 */  lui     $a2, 0x3F80                ## $a2 = 3F800000
/* 057EC 80A60ADC 24070000 */  addiu   $a3, $zero, 0x0000         ## $a3 = 00000000
/* 057F0 80A60AE0 0C029468 */  jal     SkelAnime_ChangeAnimation
              
/* 057F4 80A60AE4 E7A80018 */  swc1    $f8, 0x0018($sp)           
/* 057F8 80A60AE8 8FBF002C */  lw      $ra, 0x002C($sp)           
.L80A60AEC:
/* 057FC 80A60AEC 8FB00028 */  lw      $s0, 0x0028($sp)           
/* 05800 80A60AF0 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 05804 80A60AF4 03E00008 */  jr      $ra                        
/* 05808 80A60AF8 00000000 */  nop


