glabel func_809527A4
/* 034E4 809527A4 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 034E8 809527A8 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 034EC 809527AC AFA5001C */  sw      $a1, 0x001C($sp)           
/* 034F0 809527B0 8C8E0004 */  lw      $t6, 0x0004($a0)           ## 00000004
/* 034F4 809527B4 A48001AC */  sh      $zero, 0x01AC($a0)         ## 000001AC
/* 034F8 809527B8 849801AC */  lh      $t8, 0x01AC($a0)           ## 000001AC
/* 034FC 809527BC 8499008A */  lh      $t9, 0x008A($a0)           ## 0000008A
/* 03500 809527C0 35CF0001 */  ori     $t7, $t6, 0x0001           ## $t7 = 00000001
/* 03504 809527C4 AC8F0004 */  sw      $t7, 0x0004($a0)           ## 00000004
/* 03508 809527C8 A49801F8 */  sh      $t8, 0x01F8($a0)           ## 000001F8
/* 0350C 809527CC A4990032 */  sh      $t9, 0x0032($a0)           ## 00000032
/* 03510 809527D0 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 03514 809527D4 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 03518 809527D8 3C014316 */  lui     $at, 0x4316                ## $at = 43160000
/* 0351C 809527DC 44812000 */  mtc1    $at, $f4                   ## $f4 = 150.00
/* 03520 809527E0 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 03524 809527E4 3C018096 */  lui     $at, %hi(D_809668D0)       ## $at = 80960000
/* 03528 809527E8 46040182 */  mul.s   $f6, $f0, $f4              
/* 0352C 809527EC 240D0001 */  addiu   $t5, $zero, 0x0001         ## $t5 = 00000001
/* 03530 809527F0 240E0004 */  addiu   $t6, $zero, 0x0004         ## $t6 = 00000004
/* 03534 809527F4 44805000 */  mtc1    $zero, $f10                ## $f10 = 0.00
/* 03538 809527F8 3C058095 */  lui     $a1, %hi(func_80952858)    ## $a1 = 80950000
/* 0353C 809527FC 240FFFE2 */  addiu   $t7, $zero, 0xFFE2         ## $t7 = FFFFFFE2
/* 03540 80952800 24180037 */  addiu   $t8, $zero, 0x0037         ## $t8 = 00000037
/* 03544 80952804 4600320D */  trunc.w.s $f8, $f6                   
/* 03548 80952808 24A52858 */  addiu   $a1, $a1, %lo(func_80952858) ## $a1 = 80952858
/* 0354C 8095280C 440B4000 */  mfc1    $t3, $f8                   
/* 03550 80952810 00000000 */  nop
/* 03554 80952814 256C012C */  addiu   $t4, $t3, 0x012C           ## $t4 = 0000012C
/* 03558 80952818 A48C019C */  sh      $t4, 0x019C($a0)           ## 0000019C
/* 0355C 8095281C A02D68D0 */  sb      $t5, %lo(D_809668D0)($at)  
/* 03560 80952820 3C018096 */  lui     $at, %hi(D_80966940)       ## $at = 80960000
/* 03564 80952824 A02E6940 */  sb      $t6, %lo(D_80966940)($at)  
/* 03568 80952828 C49000BC */  lwc1    $f16, 0x00BC($a0)          ## 000000BC
/* 0356C 8095282C 46105032 */  c.eq.s  $f10, $f16                 
/* 03570 80952830 00000000 */  nop
/* 03574 80952834 45010002 */  bc1t    .L80952840                 
/* 03578 80952838 00000000 */  nop
/* 0357C 8095283C AC8F0198 */  sw      $t7, 0x0198($a0)           ## 00000198
.L80952840:
/* 03580 80952840 0C253CB0 */  jal     func_8094F2C0              
/* 03584 80952844 A49802CC */  sh      $t8, 0x02CC($a0)           ## 000002CC
/* 03588 80952848 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 0358C 8095284C 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 03590 80952850 03E00008 */  jr      $ra                        
/* 03594 80952854 00000000 */  nop


