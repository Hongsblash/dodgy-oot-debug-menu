glabel EnDekubaba_Draw
/* 03478 809E8C48 27BDFFA8 */  addiu   $sp, $sp, 0xFFA8           ## $sp = FFFFFFA8
/* 0347C 809E8C4C AFB10020 */  sw      $s1, 0x0020($sp)           
/* 03480 809E8C50 00A08825 */  or      $s1, $a1, $zero            ## $s1 = 00000000
/* 03484 809E8C54 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 03488 809E8C58 AFB0001C */  sw      $s0, 0x001C($sp)           
/* 0348C 809E8C5C 8CA50000 */  lw      $a1, 0x0000($a1)           ## 00000000
/* 03490 809E8C60 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 03494 809E8C64 3C06809F */  lui     $a2, %hi(D_809E9130)       ## $a2 = 809F0000
/* 03498 809E8C68 24C69130 */  addiu   $a2, $a2, %lo(D_809E9130)  ## $a2 = 809E9130
/* 0349C 809E8C6C 27A4003C */  addiu   $a0, $sp, 0x003C           ## $a0 = FFFFFFE4
/* 034A0 809E8C70 24070AC0 */  addiu   $a3, $zero, 0x0AC0         ## $a3 = 00000AC0
/* 034A4 809E8C74 0C031AB1 */  jal     Graph_OpenDisps              
/* 034A8 809E8C78 AFA5004C */  sw      $a1, 0x004C($sp)           
/* 034AC 809E8C7C 0C024F46 */  jal     func_80093D18              
/* 034B0 809E8C80 8E240000 */  lw      $a0, 0x0000($s1)           ## 00000000
/* 034B4 809E8C84 8E0F01C0 */  lw      $t7, 0x01C0($s0)           ## 000001C0
/* 034B8 809E8C88 3C0E809F */  lui     $t6, %hi(func_809E80D8)    ## $t6 = 809F0000
/* 034BC 809E8C8C 25CE80D8 */  addiu   $t6, $t6, %lo(func_809E80D8) ## $t6 = 809E80D8
/* 034C0 809E8C90 11CF0055 */  beq     $t6, $t7, .L809E8DE8       
/* 034C4 809E8C94 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 034C8 809E8C98 3C18809F */  lui     $t8, %hi(func_809E8C0C)    ## $t8 = 809F0000
/* 034CC 809E8C9C 27188C0C */  addiu   $t8, $t8, %lo(func_809E8C0C) ## $t8 = 809E8C0C
/* 034D0 809E8CA0 8E050180 */  lw      $a1, 0x0180($s0)           ## 00000180
/* 034D4 809E8CA4 8E06019C */  lw      $a2, 0x019C($s0)           ## 0000019C
/* 034D8 809E8CA8 AFB00014 */  sw      $s0, 0x0014($sp)           
/* 034DC 809E8CAC AFB80010 */  sw      $t8, 0x0010($sp)           
/* 034E0 809E8CB0 0C028572 */  jal     SkelAnime_Draw
              
/* 034E4 809E8CB4 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 034E8 809E8CB8 8E0801C0 */  lw      $t0, 0x01C0($s0)           ## 000001C0
/* 034EC 809E8CBC 3C19809E */  lui     $t9, %hi(func_809E64F4)    ## $t9 = 809E0000
/* 034F0 809E8CC0 273964F4 */  addiu   $t9, $t9, %lo(func_809E64F4) ## $t9 = 809E64F4
/* 034F4 809E8CC4 17280006 */  bne     $t9, $t0, .L809E8CE0       
/* 034F8 809E8CC8 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 034FC 809E8CCC 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03500 809E8CD0 0C27A163 */  jal     func_809E858C              
/* 03504 809E8CD4 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 03508 809E8CD8 10000003 */  beq     $zero, $zero, .L809E8CE8   
/* 0350C 809E8CDC 00000000 */  nop
.L809E8CE0:
/* 03510 809E8CE0 0C27A1AE */  jal     func_809E86B8              
/* 03514 809E8CE4 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
.L809E8CE8:
/* 03518 809E8CE8 3C01809F */  lui     $at, %hi(D_809E91D4)       ## $at = 809F0000
/* 0351C 809E8CEC C42691D4 */  lwc1    $f6, %lo(D_809E91D4)($at)  
/* 03520 809E8CF0 C6040230 */  lwc1    $f4, 0x0230($s0)           ## 00000230
/* 03524 809E8CF4 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 03528 809E8CF8 46062202 */  mul.s   $f8, $f4, $f6              
/* 0352C 809E8CFC E7A80050 */  swc1    $f8, 0x0050($sp)           
/* 03530 809E8D00 8E060010 */  lw      $a2, 0x0010($s0)           ## 00000010
/* 03534 809E8D04 C60E000C */  lwc1    $f14, 0x000C($s0)          ## 0000000C
/* 03538 809E8D08 0C034261 */  jal     Matrix_Translate              
/* 0353C 809E8D0C C60C0008 */  lwc1    $f12, 0x0008($s0)          ## 00000008
/* 03540 809E8D10 86090016 */  lh      $t1, 0x0016($s0)           ## 00000016
/* 03544 809E8D14 3C01809F */  lui     $at, %hi(D_809E91D8)       ## $at = 809F0000
/* 03548 809E8D18 C43291D8 */  lwc1    $f18, %lo(D_809E91D8)($at) 
/* 0354C 809E8D1C 44895000 */  mtc1    $t1, $f10                  ## $f10 = 0.00
/* 03550 809E8D20 24050001 */  addiu   $a1, $zero, 0x0001         ## $a1 = 00000001
/* 03554 809E8D24 46805420 */  cvt.s.w $f16, $f10                 
/* 03558 809E8D28 46128302 */  mul.s   $f12, $f16, $f18           
/* 0355C 809E8D2C 0C034348 */  jal     Matrix_RotateY              
/* 03560 809E8D30 00000000 */  nop
/* 03564 809E8D34 C7AC0050 */  lwc1    $f12, 0x0050($sp)          
/* 03568 809E8D38 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 0356C 809E8D3C 44066000 */  mfc1    $a2, $f12                  
/* 03570 809E8D40 0C0342A3 */  jal     Matrix_Scale              
/* 03574 809E8D44 46006386 */  mov.s   $f14, $f12                 
/* 03578 809E8D48 8FA7004C */  lw      $a3, 0x004C($sp)           
/* 0357C 809E8D4C 3C0BDA38 */  lui     $t3, 0xDA38                ## $t3 = DA380000
/* 03580 809E8D50 356B0003 */  ori     $t3, $t3, 0x0003           ## $t3 = DA380003
/* 03584 809E8D54 8CE202C0 */  lw      $v0, 0x02C0($a3)           ## 000002C0
/* 03588 809E8D58 3C05809F */  lui     $a1, %hi(D_809E9144)       ## $a1 = 809F0000
/* 0358C 809E8D5C 24A59144 */  addiu   $a1, $a1, %lo(D_809E9144)  ## $a1 = 809E9144
/* 03590 809E8D60 244A0008 */  addiu   $t2, $v0, 0x0008           ## $t2 = 00000008
/* 03594 809E8D64 ACEA02C0 */  sw      $t2, 0x02C0($a3)           ## 000002C0
/* 03598 809E8D68 AC4B0000 */  sw      $t3, 0x0000($v0)           ## 00000000
/* 0359C 809E8D6C 8E240000 */  lw      $a0, 0x0000($s1)           ## 00000000
/* 035A0 809E8D70 24060ADC */  addiu   $a2, $zero, 0x0ADC         ## $a2 = 00000ADC
/* 035A4 809E8D74 0C0346A2 */  jal     Matrix_NewMtx              
/* 035A8 809E8D78 AFA20038 */  sw      $v0, 0x0038($sp)           
/* 035AC 809E8D7C 8FA30038 */  lw      $v1, 0x0038($sp)           
/* 035B0 809E8D80 3C0E0600 */  lui     $t6, 0x0600                ## $t6 = 06000000
/* 035B4 809E8D84 25CE10F0 */  addiu   $t6, $t6, 0x10F0           ## $t6 = 060010F0
/* 035B8 809E8D88 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 035BC 809E8D8C 8FA4004C */  lw      $a0, 0x004C($sp)           
/* 035C0 809E8D90 3C0DDE00 */  lui     $t5, 0xDE00                ## $t5 = DE000000
/* 035C4 809E8D94 3C0F809E */  lui     $t7, %hi(func_809E7BB0)    ## $t7 = 809E0000
/* 035C8 809E8D98 8C8202C0 */  lw      $v0, 0x02C0($a0)           ## 000002C0
/* 035CC 809E8D9C 25EF7BB0 */  addiu   $t7, $t7, %lo(func_809E7BB0) ## $t7 = 809E7BB0
/* 035D0 809E8DA0 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 00000008
/* 035D4 809E8DA4 AC8C02C0 */  sw      $t4, 0x02C0($a0)           ## 000002C0
/* 035D8 809E8DA8 AC4E0004 */  sw      $t6, 0x0004($v0)           ## 00000004
/* 035DC 809E8DAC AC4D0000 */  sw      $t5, 0x0000($v0)           ## 00000000
/* 035E0 809E8DB0 8E1801C0 */  lw      $t8, 0x01C0($s0)           ## 000001C0
/* 035E4 809E8DB4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 035E8 809E8DB8 55F80004 */  bnel    $t7, $t8, .L809E8DCC       
/* 035EC 809E8DBC 8E190234 */  lw      $t9, 0x0234($s0)           ## 00000234
/* 035F0 809E8DC0 0C27A279 */  jal     func_809E89E4              
/* 035F4 809E8DC4 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 035F8 809E8DC8 8E190234 */  lw      $t9, 0x0234($s0)           ## 00000234
.L809E8DCC:
/* 035FC 809E8DCC 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03600 809E8DD0 13200027 */  beq     $t9, $zero, .L809E8E70     
/* 03604 809E8DD4 00000000 */  nop
/* 03608 809E8DD8 0C27A2B6 */  jal     func_809E8AD8              
/* 0360C 809E8DDC 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 03610 809E8DE0 10000023 */  beq     $zero, $zero, .L809E8E70   
/* 03614 809E8DE4 00000000 */  nop
.L809E8DE8:
/* 03618 809E8DE8 860201C6 */  lh      $v0, 0x01C6($s0)           ## 000001C6
/* 0361C 809E8DEC 3C064348 */  lui     $a2, 0x4348                ## $a2 = 43480000
/* 03620 809E8DF0 28410029 */  slti    $at, $v0, 0x0029           
/* 03624 809E8DF4 10200003 */  beq     $at, $zero, .L809E8E04     
/* 03628 809E8DF8 30480001 */  andi    $t0, $v0, 0x0001           ## $t0 = 00000000
/* 0362C 809E8DFC 1100001C */  beq     $t0, $zero, .L809E8E70     
/* 03630 809E8E00 00000000 */  nop
.L809E8E04:
/* 03634 809E8E04 44806000 */  mtc1    $zero, $f12                ## $f12 = 0.00
/* 03638 809E8E08 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 0363C 809E8E0C 0C034261 */  jal     Matrix_Translate              
/* 03640 809E8E10 46006386 */  mov.s   $f14, $f12                 
/* 03644 809E8E14 8FA3004C */  lw      $v1, 0x004C($sp)           
/* 03648 809E8E18 3C0ADA38 */  lui     $t2, 0xDA38                ## $t2 = DA380000
/* 0364C 809E8E1C 354A0003 */  ori     $t2, $t2, 0x0003           ## $t2 = DA380003
/* 03650 809E8E20 8C6202C0 */  lw      $v0, 0x02C0($v1)           ## 000002C0
/* 03654 809E8E24 3C05809F */  lui     $a1, %hi(D_809E9158)       ## $a1 = 809F0000
/* 03658 809E8E28 24A59158 */  addiu   $a1, $a1, %lo(D_809E9158)  ## $a1 = 809E9158
/* 0365C 809E8E2C 24490008 */  addiu   $t1, $v0, 0x0008           ## $t1 = 00000008
/* 03660 809E8E30 AC6902C0 */  sw      $t1, 0x02C0($v1)           ## 000002C0
/* 03664 809E8E34 AC4A0000 */  sw      $t2, 0x0000($v0)           ## 00000000
/* 03668 809E8E38 8E240000 */  lw      $a0, 0x0000($s1)           ## 00000000
/* 0366C 809E8E3C 24060AED */  addiu   $a2, $zero, 0x0AED         ## $a2 = 00000AED
/* 03670 809E8E40 0C0346A2 */  jal     Matrix_NewMtx              
/* 03674 809E8E44 00408025 */  or      $s0, $v0, $zero            ## $s0 = 00000000
/* 03678 809E8E48 AE020004 */  sw      $v0, 0x0004($s0)           ## 00000004
/* 0367C 809E8E4C 8FAB004C */  lw      $t3, 0x004C($sp)           
/* 03680 809E8E50 3C0E0600 */  lui     $t6, 0x0600                ## $t6 = 06000000
/* 03684 809E8E54 25CE3070 */  addiu   $t6, $t6, 0x3070           ## $t6 = 06003070
/* 03688 809E8E58 8D6202C0 */  lw      $v0, 0x02C0($t3)           ## 000002C0
/* 0368C 809E8E5C 3C0DDE00 */  lui     $t5, 0xDE00                ## $t5 = DE000000
/* 03690 809E8E60 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 00000008
/* 03694 809E8E64 AD6C02C0 */  sw      $t4, 0x02C0($t3)           ## 000002C0
/* 03698 809E8E68 AC4E0004 */  sw      $t6, 0x0004($v0)           ## 00000004
/* 0369C 809E8E6C AC4D0000 */  sw      $t5, 0x0000($v0)           ## 00000000
.L809E8E70:
/* 036A0 809E8E70 3C06809F */  lui     $a2, %hi(D_809E916C)       ## $a2 = 809F0000
/* 036A4 809E8E74 24C6916C */  addiu   $a2, $a2, %lo(D_809E916C)  ## $a2 = 809E916C
/* 036A8 809E8E78 27A4003C */  addiu   $a0, $sp, 0x003C           ## $a0 = FFFFFFE4
/* 036AC 809E8E7C 8E250000 */  lw      $a1, 0x0000($s1)           ## 00000000
/* 036B0 809E8E80 0C031AD5 */  jal     Graph_CloseDisps              
/* 036B4 809E8E84 24070AF4 */  addiu   $a3, $zero, 0x0AF4         ## $a3 = 00000AF4
/* 036B8 809E8E88 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 036BC 809E8E8C 8FB0001C */  lw      $s0, 0x001C($sp)           
/* 036C0 809E8E90 8FB10020 */  lw      $s1, 0x0020($sp)           
/* 036C4 809E8E94 03E00008 */  jr      $ra                        
/* 036C8 809E8E98 27BD0058 */  addiu   $sp, $sp, 0x0058           ## $sp = 00000000
/* 036CC 809E8E9C 00000000 */  nop

