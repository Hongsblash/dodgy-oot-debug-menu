glabel func_80973524
/* 025D4 80973524 44856000 */  mtc1    $a1, $f12                  ## $f12 = 0.00
/* 025D8 80973528 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 025DC 8097352C AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 025E0 80973530 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 025E4 80973534 0C25CD09 */  jal     func_80973424              
/* 025E8 80973538 E7AC001C */  swc1    $f12, 0x001C($sp)          
/* 025EC 8097353C 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 025F0 80973540 3C01437F */  lui     $at, 0x437F                ## $at = 437F0000
/* 025F4 80973544 44812000 */  mtc1    $at, $f4                   ## $f4 = 255.00
/* 025F8 80973548 908E0178 */  lbu     $t6, 0x0178($a0)           ## 00000178
/* 025FC 8097354C 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 02600 80973550 C7AC001C */  lwc1    $f12, 0x001C($sp)          
/* 02604 80973554 44813000 */  mtc1    $at, $f6                   ## $f6 = 1.00
/* 02608 80973558 448E5000 */  mtc1    $t6, $f10                  ## $f10 = 0.00
/* 0260C 8097355C 24180001 */  addiu   $t8, $zero, 0x0001         ## $t8 = 00000001
/* 02610 80973560 460C3201 */  sub.s   $f8, $f6, $f12             
/* 02614 80973564 46805420 */  cvt.s.w $f16, $f10                 
/* 02618 80973568 46082002 */  mul.s   $f0, $f4, $f8              
/* 0261C 8097356C 00000000 */  nop
/* 02620 80973570 460C8482 */  mul.s   $f18, $f16, $f12           
/* 02624 80973574 46009180 */  add.s   $f6, $f18, $f0             
/* 02628 80973578 444FF800 */  cfc1    $t7, $31
/* 0262C 8097357C 44D8F800 */  ctc1    $t8, $31
/* 02630 80973580 00000000 */  nop
/* 02634 80973584 46003124 */  cvt.w.s $f4, $f6                   
/* 02638 80973588 4458F800 */  cfc1    $t8, $31
/* 0263C 8097358C 00000000 */  nop
/* 02640 80973590 33180078 */  andi    $t8, $t8, 0x0078           ## $t8 = 00000000
/* 02644 80973594 13000013 */  beq     $t8, $zero, .L809735E4     
/* 02648 80973598 00000000 */  nop
/* 0264C 8097359C 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 02650 809735A0 44812000 */  mtc1    $at, $f4                   ## $f4 = 2147483648.00
/* 02654 809735A4 24180001 */  addiu   $t8, $zero, 0x0001         ## $t8 = 00000001
/* 02658 809735A8 46043101 */  sub.s   $f4, $f6, $f4              
/* 0265C 809735AC 44D8F800 */  ctc1    $t8, $31
/* 02660 809735B0 00000000 */  nop
/* 02664 809735B4 46002124 */  cvt.w.s $f4, $f4                   
/* 02668 809735B8 4458F800 */  cfc1    $t8, $31
/* 0266C 809735BC 00000000 */  nop
/* 02670 809735C0 33180078 */  andi    $t8, $t8, 0x0078           ## $t8 = 00000000
/* 02674 809735C4 17000005 */  bne     $t8, $zero, .L809735DC     
/* 02678 809735C8 00000000 */  nop
/* 0267C 809735CC 44182000 */  mfc1    $t8, $f4                   
/* 02680 809735D0 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 02684 809735D4 10000007 */  beq     $zero, $zero, .L809735F4   
/* 02688 809735D8 0301C025 */  or      $t8, $t8, $at              ## $t8 = 80000000
.L809735DC:
/* 0268C 809735DC 10000005 */  beq     $zero, $zero, .L809735F4   
/* 02690 809735E0 2418FFFF */  addiu   $t8, $zero, 0xFFFF         ## $t8 = FFFFFFFF
.L809735E4:
/* 02694 809735E4 44182000 */  mfc1    $t8, $f4                   
/* 02698 809735E8 00000000 */  nop
/* 0269C 809735EC 0700FFFB */  bltz    $t8, .L809735DC            
/* 026A0 809735F0 00000000 */  nop
.L809735F4:
/* 026A4 809735F4 90990179 */  lbu     $t9, 0x0179($a0)           ## 00000179
/* 026A8 809735F8 44CFF800 */  ctc1    $t7, $31
/* 026AC 809735FC 24090001 */  addiu   $t1, $zero, 0x0001         ## $t1 = 00000001
/* 026B0 80973600 44994000 */  mtc1    $t9, $f8                   ## $f8 = 0.00
/* 026B4 80973604 A0980178 */  sb      $t8, 0x0178($a0)           ## 00000178
/* 026B8 80973608 468042A0 */  cvt.s.w $f10, $f8                  
/* 026BC 8097360C 460C5402 */  mul.s   $f16, $f10, $f12           
/* 026C0 80973610 46008480 */  add.s   $f18, $f16, $f0            
/* 026C4 80973614 4448F800 */  cfc1    $t0, $31
/* 026C8 80973618 44C9F800 */  ctc1    $t1, $31
/* 026CC 8097361C 00000000 */  nop
/* 026D0 80973620 460091A4 */  cvt.w.s $f6, $f18                  
/* 026D4 80973624 4449F800 */  cfc1    $t1, $31
/* 026D8 80973628 00000000 */  nop
/* 026DC 8097362C 31290078 */  andi    $t1, $t1, 0x0078           ## $t1 = 00000000
/* 026E0 80973630 11200012 */  beq     $t1, $zero, .L8097367C     
/* 026E4 80973634 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 026E8 80973638 44813000 */  mtc1    $at, $f6                   ## $f6 = 2147483648.00
/* 026EC 8097363C 24090001 */  addiu   $t1, $zero, 0x0001         ## $t1 = 00000001
/* 026F0 80973640 46069181 */  sub.s   $f6, $f18, $f6             
/* 026F4 80973644 44C9F800 */  ctc1    $t1, $31
/* 026F8 80973648 00000000 */  nop
/* 026FC 8097364C 460031A4 */  cvt.w.s $f6, $f6                   
/* 02700 80973650 4449F800 */  cfc1    $t1, $31
/* 02704 80973654 00000000 */  nop
/* 02708 80973658 31290078 */  andi    $t1, $t1, 0x0078           ## $t1 = 00000000
/* 0270C 8097365C 15200005 */  bne     $t1, $zero, .L80973674     
/* 02710 80973660 00000000 */  nop
/* 02714 80973664 44093000 */  mfc1    $t1, $f6                   
/* 02718 80973668 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 0271C 8097366C 10000007 */  beq     $zero, $zero, .L8097368C   
/* 02720 80973670 01214825 */  or      $t1, $t1, $at              ## $t1 = 80000000
.L80973674:
/* 02724 80973674 10000005 */  beq     $zero, $zero, .L8097368C   
/* 02728 80973678 2409FFFF */  addiu   $t1, $zero, 0xFFFF         ## $t1 = FFFFFFFF
.L8097367C:
/* 0272C 8097367C 44093000 */  mfc1    $t1, $f6                   
/* 02730 80973680 00000000 */  nop
/* 02734 80973684 0520FFFB */  bltz    $t1, .L80973674            
/* 02738 80973688 00000000 */  nop
.L8097368C:
/* 0273C 8097368C 908A017A */  lbu     $t2, 0x017A($a0)           ## 0000017A
/* 02740 80973690 44C8F800 */  ctc1    $t0, $31
/* 02744 80973694 240C0001 */  addiu   $t4, $zero, 0x0001         ## $t4 = 00000001
/* 02748 80973698 448A2000 */  mtc1    $t2, $f4                   ## $f4 = 0.00
/* 0274C 8097369C A0890179 */  sb      $t1, 0x0179($a0)           ## 00000179
/* 02750 809736A0 46802220 */  cvt.s.w $f8, $f4                   
/* 02754 809736A4 460C4282 */  mul.s   $f10, $f8, $f12            
/* 02758 809736A8 46005400 */  add.s   $f16, $f10, $f0            
/* 0275C 809736AC 444BF800 */  cfc1    $t3, $31
/* 02760 809736B0 44CCF800 */  ctc1    $t4, $31
/* 02764 809736B4 00000000 */  nop
/* 02768 809736B8 460084A4 */  cvt.w.s $f18, $f16                 
/* 0276C 809736BC 444CF800 */  cfc1    $t4, $31
/* 02770 809736C0 00000000 */  nop
/* 02774 809736C4 318C0078 */  andi    $t4, $t4, 0x0078           ## $t4 = 00000000
/* 02778 809736C8 11800012 */  beq     $t4, $zero, .L80973714     
/* 0277C 809736CC 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 02780 809736D0 44819000 */  mtc1    $at, $f18                  ## $f18 = 2147483648.00
/* 02784 809736D4 240C0001 */  addiu   $t4, $zero, 0x0001         ## $t4 = 00000001
/* 02788 809736D8 46128481 */  sub.s   $f18, $f16, $f18           
/* 0278C 809736DC 44CCF800 */  ctc1    $t4, $31
/* 02790 809736E0 00000000 */  nop
/* 02794 809736E4 460094A4 */  cvt.w.s $f18, $f18                 
/* 02798 809736E8 444CF800 */  cfc1    $t4, $31
/* 0279C 809736EC 00000000 */  nop
/* 027A0 809736F0 318C0078 */  andi    $t4, $t4, 0x0078           ## $t4 = 00000000
/* 027A4 809736F4 15800005 */  bne     $t4, $zero, .L8097370C     
/* 027A8 809736F8 00000000 */  nop
/* 027AC 809736FC 440C9000 */  mfc1    $t4, $f18                  
/* 027B0 80973700 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 027B4 80973704 10000007 */  beq     $zero, $zero, .L80973724   
/* 027B8 80973708 01816025 */  or      $t4, $t4, $at              ## $t4 = 80000000
.L8097370C:
/* 027BC 8097370C 10000005 */  beq     $zero, $zero, .L80973724   
/* 027C0 80973710 240CFFFF */  addiu   $t4, $zero, 0xFFFF         ## $t4 = FFFFFFFF
.L80973714:
/* 027C4 80973714 440C9000 */  mfc1    $t4, $f18                  
/* 027C8 80973718 00000000 */  nop
/* 027CC 8097371C 0580FFFB */  bltz    $t4, .L8097370C            
/* 027D0 80973720 00000000 */  nop
.L80973724:
/* 027D4 80973724 908D017E */  lbu     $t5, 0x017E($a0)           ## 0000017E
/* 027D8 80973728 44CBF800 */  ctc1    $t3, $31
/* 027DC 8097372C 240F0001 */  addiu   $t7, $zero, 0x0001         ## $t7 = 00000001
/* 027E0 80973730 448D3000 */  mtc1    $t5, $f6                   ## $f6 = 0.00
/* 027E4 80973734 A08C017A */  sb      $t4, 0x017A($a0)           ## 0000017A
/* 027E8 80973738 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 027EC 8097373C 46803120 */  cvt.s.w $f4, $f6                   
/* 027F0 80973740 460C2202 */  mul.s   $f8, $f4, $f12             
/* 027F4 80973744 46004280 */  add.s   $f10, $f8, $f0             
/* 027F8 80973748 444EF800 */  cfc1    $t6, $31
/* 027FC 8097374C 44CFF800 */  ctc1    $t7, $31
/* 02800 80973750 00000000 */  nop
/* 02804 80973754 46005424 */  cvt.w.s $f16, $f10                 
/* 02808 80973758 444FF800 */  cfc1    $t7, $31
/* 0280C 8097375C 00000000 */  nop
/* 02810 80973760 31EF0078 */  andi    $t7, $t7, 0x0078           ## $t7 = 00000000
/* 02814 80973764 51E00013 */  beql    $t7, $zero, .L809737B4     
/* 02818 80973768 440F8000 */  mfc1    $t7, $f16                  
/* 0281C 8097376C 44818000 */  mtc1    $at, $f16                  ## $f16 = 2147483648.00
/* 02820 80973770 240F0001 */  addiu   $t7, $zero, 0x0001         ## $t7 = 00000001
/* 02824 80973774 46105401 */  sub.s   $f16, $f10, $f16           
/* 02828 80973778 44CFF800 */  ctc1    $t7, $31
/* 0282C 8097377C 00000000 */  nop
/* 02830 80973780 46008424 */  cvt.w.s $f16, $f16                 
/* 02834 80973784 444FF800 */  cfc1    $t7, $31
/* 02838 80973788 00000000 */  nop
/* 0283C 8097378C 31EF0078 */  andi    $t7, $t7, 0x0078           ## $t7 = 00000000
/* 02840 80973790 15E00005 */  bne     $t7, $zero, .L809737A8     
/* 02844 80973794 00000000 */  nop
/* 02848 80973798 440F8000 */  mfc1    $t7, $f16                  
/* 0284C 8097379C 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 02850 809737A0 10000007 */  beq     $zero, $zero, .L809737C0   
/* 02854 809737A4 01E17825 */  or      $t7, $t7, $at              ## $t7 = 80000000
.L809737A8:
/* 02858 809737A8 10000005 */  beq     $zero, $zero, .L809737C0   
/* 0285C 809737AC 240FFFFF */  addiu   $t7, $zero, 0xFFFF         ## $t7 = FFFFFFFF
/* 02860 809737B0 440F8000 */  mfc1    $t7, $f16                  
.L809737B4:
/* 02864 809737B4 00000000 */  nop
/* 02868 809737B8 05E0FFFB */  bltz    $t7, .L809737A8            
/* 0286C 809737BC 00000000 */  nop
.L809737C0:
/* 02870 809737C0 9098017F */  lbu     $t8, 0x017F($a0)           ## 0000017F
/* 02874 809737C4 44CEF800 */  ctc1    $t6, $31
/* 02878 809737C8 24080001 */  addiu   $t0, $zero, 0x0001         ## $t0 = 00000001
/* 0287C 809737CC 44989000 */  mtc1    $t8, $f18                  ## $f18 = NaN
/* 02880 809737D0 A08F017E */  sb      $t7, 0x017E($a0)           ## 0000017E
/* 02884 809737D4 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 02888 809737D8 468091A0 */  cvt.s.w $f6, $f18                  
/* 0288C 809737DC 460C3102 */  mul.s   $f4, $f6, $f12             
/* 02890 809737E0 46002200 */  add.s   $f8, $f4, $f0              
/* 02894 809737E4 4459F800 */  cfc1    $t9, $31
/* 02898 809737E8 44C8F800 */  ctc1    $t0, $31
/* 0289C 809737EC 00000000 */  nop
/* 028A0 809737F0 460042A4 */  cvt.w.s $f10, $f8                  
/* 028A4 809737F4 4448F800 */  cfc1    $t0, $31
/* 028A8 809737F8 00000000 */  nop
/* 028AC 809737FC 31080078 */  andi    $t0, $t0, 0x0078           ## $t0 = 00000000
/* 028B0 80973800 51000013 */  beql    $t0, $zero, .L80973850     
/* 028B4 80973804 44085000 */  mfc1    $t0, $f10                  
/* 028B8 80973808 44815000 */  mtc1    $at, $f10                  ## $f10 = 2147483648.00
/* 028BC 8097380C 24080001 */  addiu   $t0, $zero, 0x0001         ## $t0 = 00000001
/* 028C0 80973810 460A4281 */  sub.s   $f10, $f8, $f10            
/* 028C4 80973814 44C8F800 */  ctc1    $t0, $31
/* 028C8 80973818 00000000 */  nop
/* 028CC 8097381C 460052A4 */  cvt.w.s $f10, $f10                 
/* 028D0 80973820 4448F800 */  cfc1    $t0, $31
/* 028D4 80973824 00000000 */  nop
/* 028D8 80973828 31080078 */  andi    $t0, $t0, 0x0078           ## $t0 = 00000000
/* 028DC 8097382C 15000005 */  bne     $t0, $zero, .L80973844     
/* 028E0 80973830 00000000 */  nop
/* 028E4 80973834 44085000 */  mfc1    $t0, $f10                  
/* 028E8 80973838 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 028EC 8097383C 10000007 */  beq     $zero, $zero, .L8097385C   
/* 028F0 80973840 01014025 */  or      $t0, $t0, $at              ## $t0 = 80000000
.L80973844:
/* 028F4 80973844 10000005 */  beq     $zero, $zero, .L8097385C   
/* 028F8 80973848 2408FFFF */  addiu   $t0, $zero, 0xFFFF         ## $t0 = FFFFFFFF
/* 028FC 8097384C 44085000 */  mfc1    $t0, $f10                  
.L80973850:
/* 02900 80973850 00000000 */  nop
/* 02904 80973854 0500FFFB */  bltz    $t0, .L80973844            
/* 02908 80973858 00000000 */  nop
.L8097385C:
/* 0290C 8097385C 90890180 */  lbu     $t1, 0x0180($a0)           ## 00000180
/* 02910 80973860 44D9F800 */  ctc1    $t9, $31
/* 02914 80973864 240B0001 */  addiu   $t3, $zero, 0x0001         ## $t3 = 00000001
/* 02918 80973868 44898000 */  mtc1    $t1, $f16                  ## $f16 = NaN
/* 0291C 8097386C A088017F */  sb      $t0, 0x017F($a0)           ## 0000017F
/* 02920 80973870 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 02924 80973874 468084A0 */  cvt.s.w $f18, $f16                 
/* 02928 80973878 460C9182 */  mul.s   $f6, $f18, $f12            
/* 0292C 8097387C 46003100 */  add.s   $f4, $f6, $f0              
/* 02930 80973880 444AF800 */  cfc1    $t2, $31
/* 02934 80973884 44CBF800 */  ctc1    $t3, $31
/* 02938 80973888 00000000 */  nop
/* 0293C 8097388C 46002224 */  cvt.w.s $f8, $f4                   
/* 02940 80973890 444BF800 */  cfc1    $t3, $31
/* 02944 80973894 00000000 */  nop
/* 02948 80973898 316B0078 */  andi    $t3, $t3, 0x0078           ## $t3 = 00000000
/* 0294C 8097389C 51600013 */  beql    $t3, $zero, .L809738EC     
/* 02950 809738A0 440B4000 */  mfc1    $t3, $f8                   
/* 02954 809738A4 44814000 */  mtc1    $at, $f8                   ## $f8 = 2147483648.00
/* 02958 809738A8 240B0001 */  addiu   $t3, $zero, 0x0001         ## $t3 = 00000001
/* 0295C 809738AC 46082201 */  sub.s   $f8, $f4, $f8              
/* 02960 809738B0 44CBF800 */  ctc1    $t3, $31
/* 02964 809738B4 00000000 */  nop
/* 02968 809738B8 46004224 */  cvt.w.s $f8, $f8                   
/* 0296C 809738BC 444BF800 */  cfc1    $t3, $31
/* 02970 809738C0 00000000 */  nop
/* 02974 809738C4 316B0078 */  andi    $t3, $t3, 0x0078           ## $t3 = 00000000
/* 02978 809738C8 15600005 */  bne     $t3, $zero, .L809738E0     
/* 0297C 809738CC 00000000 */  nop
/* 02980 809738D0 440B4000 */  mfc1    $t3, $f8                   
/* 02984 809738D4 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 02988 809738D8 10000007 */  beq     $zero, $zero, .L809738F8   
/* 0298C 809738DC 01615825 */  or      $t3, $t3, $at              ## $t3 = 80000000
.L809738E0:
/* 02990 809738E0 10000005 */  beq     $zero, $zero, .L809738F8   
/* 02994 809738E4 240BFFFF */  addiu   $t3, $zero, 0xFFFF         ## $t3 = FFFFFFFF
/* 02998 809738E8 440B4000 */  mfc1    $t3, $f8                   
.L809738EC:
/* 0299C 809738EC 00000000 */  nop
/* 029A0 809738F0 0560FFFB */  bltz    $t3, .L809738E0            
/* 029A4 809738F4 00000000 */  nop
.L809738F8:
/* 029A8 809738F8 908C017B */  lbu     $t4, 0x017B($a0)           ## 0000017B
/* 029AC 809738FC 44CAF800 */  ctc1    $t2, $31
/* 029B0 80973900 240E0001 */  addiu   $t6, $zero, 0x0001         ## $t6 = 00000001
/* 029B4 80973904 448C5000 */  mtc1    $t4, $f10                  ## $f10 = NaN
/* 029B8 80973908 A08B0180 */  sb      $t3, 0x0180($a0)           ## 00000180
/* 029BC 8097390C 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 029C0 80973910 46805420 */  cvt.s.w $f16, $f10                 
/* 029C4 80973914 460C8482 */  mul.s   $f18, $f16, $f12           
/* 029C8 80973918 444DF800 */  cfc1    $t5, $31
/* 029CC 8097391C 44CEF800 */  ctc1    $t6, $31
/* 029D0 80973920 00000000 */  nop
/* 029D4 80973924 460091A4 */  cvt.w.s $f6, $f18                  
/* 029D8 80973928 444EF800 */  cfc1    $t6, $31
/* 029DC 8097392C 00000000 */  nop
/* 029E0 80973930 31CE0078 */  andi    $t6, $t6, 0x0078           ## $t6 = 00000000
/* 029E4 80973934 51C00013 */  beql    $t6, $zero, .L80973984     
/* 029E8 80973938 440E3000 */  mfc1    $t6, $f6                   
/* 029EC 8097393C 44813000 */  mtc1    $at, $f6                   ## $f6 = 2147483648.00
/* 029F0 80973940 240E0001 */  addiu   $t6, $zero, 0x0001         ## $t6 = 00000001
/* 029F4 80973944 46069181 */  sub.s   $f6, $f18, $f6             
/* 029F8 80973948 44CEF800 */  ctc1    $t6, $31
/* 029FC 8097394C 00000000 */  nop
/* 02A00 80973950 460031A4 */  cvt.w.s $f6, $f6                   
/* 02A04 80973954 444EF800 */  cfc1    $t6, $31
/* 02A08 80973958 00000000 */  nop
/* 02A0C 8097395C 31CE0078 */  andi    $t6, $t6, 0x0078           ## $t6 = 00000000
/* 02A10 80973960 15C00005 */  bne     $t6, $zero, .L80973978     
/* 02A14 80973964 00000000 */  nop
/* 02A18 80973968 440E3000 */  mfc1    $t6, $f6                   
/* 02A1C 8097396C 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 02A20 80973970 10000007 */  beq     $zero, $zero, .L80973990   
/* 02A24 80973974 01C17025 */  or      $t6, $t6, $at              ## $t6 = 80000000
.L80973978:
/* 02A28 80973978 10000005 */  beq     $zero, $zero, .L80973990   
/* 02A2C 8097397C 240EFFFF */  addiu   $t6, $zero, 0xFFFF         ## $t6 = FFFFFFFF
/* 02A30 80973980 440E3000 */  mfc1    $t6, $f6                   
.L80973984:
/* 02A34 80973984 00000000 */  nop
/* 02A38 80973988 05C0FFFB */  bltz    $t6, .L80973978            
/* 02A3C 8097398C 00000000 */  nop
.L80973990:
/* 02A40 80973990 908F017C */  lbu     $t7, 0x017C($a0)           ## 0000017C
/* 02A44 80973994 44CDF800 */  ctc1    $t5, $31
/* 02A48 80973998 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 02A4C 8097399C 448F2000 */  mtc1    $t7, $f4                   ## $f4 = NaN
/* 02A50 809739A0 A08E017B */  sb      $t6, 0x017B($a0)           ## 0000017B
/* 02A54 809739A4 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 02A58 809739A8 46802220 */  cvt.s.w $f8, $f4                   
/* 02A5C 809739AC 460C4282 */  mul.s   $f10, $f8, $f12            
/* 02A60 809739B0 4458F800 */  cfc1    $t8, $31
/* 02A64 809739B4 44D9F800 */  ctc1    $t9, $31
/* 02A68 809739B8 00000000 */  nop
/* 02A6C 809739BC 46005424 */  cvt.w.s $f16, $f10                 
/* 02A70 809739C0 4459F800 */  cfc1    $t9, $31
/* 02A74 809739C4 00000000 */  nop
/* 02A78 809739C8 33390078 */  andi    $t9, $t9, 0x0078           ## $t9 = 00000000
/* 02A7C 809739CC 53200013 */  beql    $t9, $zero, .L80973A1C     
/* 02A80 809739D0 44198000 */  mfc1    $t9, $f16                  
/* 02A84 809739D4 44818000 */  mtc1    $at, $f16                  ## $f16 = 2147483648.00
/* 02A88 809739D8 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 02A8C 809739DC 46105401 */  sub.s   $f16, $f10, $f16           
/* 02A90 809739E0 44D9F800 */  ctc1    $t9, $31
/* 02A94 809739E4 00000000 */  nop
/* 02A98 809739E8 46008424 */  cvt.w.s $f16, $f16                 
/* 02A9C 809739EC 4459F800 */  cfc1    $t9, $31
/* 02AA0 809739F0 00000000 */  nop
/* 02AA4 809739F4 33390078 */  andi    $t9, $t9, 0x0078           ## $t9 = 00000000
/* 02AA8 809739F8 17200005 */  bne     $t9, $zero, .L80973A10     
/* 02AAC 809739FC 00000000 */  nop
/* 02AB0 80973A00 44198000 */  mfc1    $t9, $f16                  
/* 02AB4 80973A04 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 02AB8 80973A08 10000007 */  beq     $zero, $zero, .L80973A28   
/* 02ABC 80973A0C 0321C825 */  or      $t9, $t9, $at              ## $t9 = 80000000
.L80973A10:
/* 02AC0 80973A10 10000005 */  beq     $zero, $zero, .L80973A28   
/* 02AC4 80973A14 2419FFFF */  addiu   $t9, $zero, 0xFFFF         ## $t9 = FFFFFFFF
/* 02AC8 80973A18 44198000 */  mfc1    $t9, $f16                  
.L80973A1C:
/* 02ACC 80973A1C 00000000 */  nop
/* 02AD0 80973A20 0720FFFB */  bltz    $t9, .L80973A10            
/* 02AD4 80973A24 00000000 */  nop
.L80973A28:
/* 02AD8 80973A28 9088017D */  lbu     $t0, 0x017D($a0)           ## 0000017D
/* 02ADC 80973A2C 44D8F800 */  ctc1    $t8, $31
/* 02AE0 80973A30 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 02AE4 80973A34 44889000 */  mtc1    $t0, $f18                  ## $f18 = NaN
/* 02AE8 80973A38 A099017C */  sb      $t9, 0x017C($a0)           ## 0000017C
/* 02AEC 80973A3C 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 02AF0 80973A40 468091A0 */  cvt.s.w $f6, $f18                  
/* 02AF4 80973A44 460C3102 */  mul.s   $f4, $f6, $f12             
/* 02AF8 80973A48 4449F800 */  cfc1    $t1, $31
/* 02AFC 80973A4C 44CAF800 */  ctc1    $t2, $31
/* 02B00 80973A50 00000000 */  nop
/* 02B04 80973A54 46002224 */  cvt.w.s $f8, $f4                   
/* 02B08 80973A58 444AF800 */  cfc1    $t2, $31
/* 02B0C 80973A5C 00000000 */  nop
/* 02B10 80973A60 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 02B14 80973A64 51400013 */  beql    $t2, $zero, .L80973AB4     
/* 02B18 80973A68 440A4000 */  mfc1    $t2, $f8                   
/* 02B1C 80973A6C 44814000 */  mtc1    $at, $f8                   ## $f8 = 2147483648.00
/* 02B20 80973A70 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 02B24 80973A74 46082201 */  sub.s   $f8, $f4, $f8              
/* 02B28 80973A78 44CAF800 */  ctc1    $t2, $31
/* 02B2C 80973A7C 00000000 */  nop
/* 02B30 80973A80 46004224 */  cvt.w.s $f8, $f8                   
/* 02B34 80973A84 444AF800 */  cfc1    $t2, $31
/* 02B38 80973A88 00000000 */  nop
/* 02B3C 80973A8C 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 02B40 80973A90 15400005 */  bne     $t2, $zero, .L80973AA8     
/* 02B44 80973A94 00000000 */  nop
/* 02B48 80973A98 440A4000 */  mfc1    $t2, $f8                   
/* 02B4C 80973A9C 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 02B50 80973AA0 10000007 */  beq     $zero, $zero, .L80973AC0   
/* 02B54 80973AA4 01415025 */  or      $t2, $t2, $at              ## $t2 = 80000000
.L80973AA8:
/* 02B58 80973AA8 10000005 */  beq     $zero, $zero, .L80973AC0   
/* 02B5C 80973AAC 240AFFFF */  addiu   $t2, $zero, 0xFFFF         ## $t2 = FFFFFFFF
/* 02B60 80973AB0 440A4000 */  mfc1    $t2, $f8                   
.L80973AB4:
/* 02B64 80973AB4 00000000 */  nop
/* 02B68 80973AB8 0540FFFB */  bltz    $t2, .L80973AA8            
/* 02B6C 80973ABC 00000000 */  nop
.L80973AC0:
/* 02B70 80973AC0 908B0181 */  lbu     $t3, 0x0181($a0)           ## 00000181
/* 02B74 80973AC4 44C9F800 */  ctc1    $t1, $31
/* 02B78 80973AC8 240D0001 */  addiu   $t5, $zero, 0x0001         ## $t5 = 00000001
/* 02B7C 80973ACC 448B5000 */  mtc1    $t3, $f10                  ## $f10 = NaN
/* 02B80 80973AD0 A08A017D */  sb      $t2, 0x017D($a0)           ## 0000017D
/* 02B84 80973AD4 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 02B88 80973AD8 46805420 */  cvt.s.w $f16, $f10                 
/* 02B8C 80973ADC 460C8482 */  mul.s   $f18, $f16, $f12           
/* 02B90 80973AE0 444CF800 */  cfc1    $t4, $31
/* 02B94 80973AE4 44CDF800 */  ctc1    $t5, $31
/* 02B98 80973AE8 00000000 */  nop
/* 02B9C 80973AEC 460091A4 */  cvt.w.s $f6, $f18                  
/* 02BA0 80973AF0 444DF800 */  cfc1    $t5, $31
/* 02BA4 80973AF4 00000000 */  nop
/* 02BA8 80973AF8 31AD0078 */  andi    $t5, $t5, 0x0078           ## $t5 = 00000000
/* 02BAC 80973AFC 51A00013 */  beql    $t5, $zero, .L80973B4C     
/* 02BB0 80973B00 440D3000 */  mfc1    $t5, $f6                   
/* 02BB4 80973B04 44813000 */  mtc1    $at, $f6                   ## $f6 = 2147483648.00
/* 02BB8 80973B08 240D0001 */  addiu   $t5, $zero, 0x0001         ## $t5 = 00000001
/* 02BBC 80973B0C 46069181 */  sub.s   $f6, $f18, $f6             
/* 02BC0 80973B10 44CDF800 */  ctc1    $t5, $31
/* 02BC4 80973B14 00000000 */  nop
/* 02BC8 80973B18 460031A4 */  cvt.w.s $f6, $f6                   
/* 02BCC 80973B1C 444DF800 */  cfc1    $t5, $31
/* 02BD0 80973B20 00000000 */  nop
/* 02BD4 80973B24 31AD0078 */  andi    $t5, $t5, 0x0078           ## $t5 = 00000000
/* 02BD8 80973B28 15A00005 */  bne     $t5, $zero, .L80973B40     
/* 02BDC 80973B2C 00000000 */  nop
/* 02BE0 80973B30 440D3000 */  mfc1    $t5, $f6                   
/* 02BE4 80973B34 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 02BE8 80973B38 10000007 */  beq     $zero, $zero, .L80973B58   
/* 02BEC 80973B3C 01A16825 */  or      $t5, $t5, $at              ## $t5 = 80000000
.L80973B40:
/* 02BF0 80973B40 10000005 */  beq     $zero, $zero, .L80973B58   
/* 02BF4 80973B44 240DFFFF */  addiu   $t5, $zero, 0xFFFF         ## $t5 = FFFFFFFF
/* 02BF8 80973B48 440D3000 */  mfc1    $t5, $f6                   
.L80973B4C:
/* 02BFC 80973B4C 00000000 */  nop
/* 02C00 80973B50 05A0FFFB */  bltz    $t5, .L80973B40            
/* 02C04 80973B54 00000000 */  nop
.L80973B58:
/* 02C08 80973B58 908E0182 */  lbu     $t6, 0x0182($a0)           ## 00000182
/* 02C0C 80973B5C 44CCF800 */  ctc1    $t4, $31
/* 02C10 80973B60 24180001 */  addiu   $t8, $zero, 0x0001         ## $t8 = 00000001
/* 02C14 80973B64 448E2000 */  mtc1    $t6, $f4                   ## $f4 = NaN
/* 02C18 80973B68 A08D0181 */  sb      $t5, 0x0181($a0)           ## 00000181
/* 02C1C 80973B6C 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 02C20 80973B70 46802220 */  cvt.s.w $f8, $f4                   
/* 02C24 80973B74 460C4282 */  mul.s   $f10, $f8, $f12            
/* 02C28 80973B78 444FF800 */  cfc1    $t7, $31
/* 02C2C 80973B7C 44D8F800 */  ctc1    $t8, $31
/* 02C30 80973B80 00000000 */  nop
/* 02C34 80973B84 46005424 */  cvt.w.s $f16, $f10                 
/* 02C38 80973B88 4458F800 */  cfc1    $t8, $31
/* 02C3C 80973B8C 00000000 */  nop
/* 02C40 80973B90 33180078 */  andi    $t8, $t8, 0x0078           ## $t8 = 00000000
/* 02C44 80973B94 53000013 */  beql    $t8, $zero, .L80973BE4     
/* 02C48 80973B98 44188000 */  mfc1    $t8, $f16                  
/* 02C4C 80973B9C 44818000 */  mtc1    $at, $f16                  ## $f16 = 2147483648.00
/* 02C50 80973BA0 24180001 */  addiu   $t8, $zero, 0x0001         ## $t8 = 00000001
/* 02C54 80973BA4 46105401 */  sub.s   $f16, $f10, $f16           
/* 02C58 80973BA8 44D8F800 */  ctc1    $t8, $31
/* 02C5C 80973BAC 00000000 */  nop
/* 02C60 80973BB0 46008424 */  cvt.w.s $f16, $f16                 
/* 02C64 80973BB4 4458F800 */  cfc1    $t8, $31
/* 02C68 80973BB8 00000000 */  nop
/* 02C6C 80973BBC 33180078 */  andi    $t8, $t8, 0x0078           ## $t8 = 00000000
/* 02C70 80973BC0 17000005 */  bne     $t8, $zero, .L80973BD8     
/* 02C74 80973BC4 00000000 */  nop
/* 02C78 80973BC8 44188000 */  mfc1    $t8, $f16                  
/* 02C7C 80973BCC 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 02C80 80973BD0 10000007 */  beq     $zero, $zero, .L80973BF0   
/* 02C84 80973BD4 0301C025 */  or      $t8, $t8, $at              ## $t8 = 80000000
.L80973BD8:
/* 02C88 80973BD8 10000005 */  beq     $zero, $zero, .L80973BF0   
/* 02C8C 80973BDC 2418FFFF */  addiu   $t8, $zero, 0xFFFF         ## $t8 = FFFFFFFF
/* 02C90 80973BE0 44188000 */  mfc1    $t8, $f16                  
.L80973BE4:
/* 02C94 80973BE4 00000000 */  nop
/* 02C98 80973BE8 0700FFFB */  bltz    $t8, .L80973BD8            
/* 02C9C 80973BEC 00000000 */  nop
.L80973BF0:
/* 02CA0 80973BF0 90990183 */  lbu     $t9, 0x0183($a0)           ## 00000183
/* 02CA4 80973BF4 44CFF800 */  ctc1    $t7, $31
/* 02CA8 80973BF8 24090001 */  addiu   $t1, $zero, 0x0001         ## $t1 = 00000001
/* 02CAC 80973BFC 44999000 */  mtc1    $t9, $f18                  ## $f18 = NaN
/* 02CB0 80973C00 A0980182 */  sb      $t8, 0x0182($a0)           ## 00000182
/* 02CB4 80973C04 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 02CB8 80973C08 468091A0 */  cvt.s.w $f6, $f18                  
/* 02CBC 80973C0C 460C3102 */  mul.s   $f4, $f6, $f12             
/* 02CC0 80973C10 4448F800 */  cfc1    $t0, $31
/* 02CC4 80973C14 44C9F800 */  ctc1    $t1, $31
/* 02CC8 80973C18 00000000 */  nop
/* 02CCC 80973C1C 46002224 */  cvt.w.s $f8, $f4                   
/* 02CD0 80973C20 4449F800 */  cfc1    $t1, $31
/* 02CD4 80973C24 00000000 */  nop
/* 02CD8 80973C28 31290078 */  andi    $t1, $t1, 0x0078           ## $t1 = 00000000
/* 02CDC 80973C2C 51200013 */  beql    $t1, $zero, .L80973C7C     
/* 02CE0 80973C30 44094000 */  mfc1    $t1, $f8                   
/* 02CE4 80973C34 44814000 */  mtc1    $at, $f8                   ## $f8 = 2147483648.00
/* 02CE8 80973C38 24090001 */  addiu   $t1, $zero, 0x0001         ## $t1 = 00000001
/* 02CEC 80973C3C 46082201 */  sub.s   $f8, $f4, $f8              
/* 02CF0 80973C40 44C9F800 */  ctc1    $t1, $31
/* 02CF4 80973C44 00000000 */  nop
/* 02CF8 80973C48 46004224 */  cvt.w.s $f8, $f8                   
/* 02CFC 80973C4C 4449F800 */  cfc1    $t1, $31
/* 02D00 80973C50 00000000 */  nop
/* 02D04 80973C54 31290078 */  andi    $t1, $t1, 0x0078           ## $t1 = 00000000
/* 02D08 80973C58 15200005 */  bne     $t1, $zero, .L80973C70     
/* 02D0C 80973C5C 00000000 */  nop
/* 02D10 80973C60 44094000 */  mfc1    $t1, $f8                   
/* 02D14 80973C64 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 02D18 80973C68 10000007 */  beq     $zero, $zero, .L80973C88   
/* 02D1C 80973C6C 01214825 */  or      $t1, $t1, $at              ## $t1 = 80000000
.L80973C70:
/* 02D20 80973C70 10000005 */  beq     $zero, $zero, .L80973C88   
/* 02D24 80973C74 2409FFFF */  addiu   $t1, $zero, 0xFFFF         ## $t1 = FFFFFFFF
/* 02D28 80973C78 44094000 */  mfc1    $t1, $f8                   
.L80973C7C:
/* 02D2C 80973C7C 00000000 */  nop
/* 02D30 80973C80 0520FFFB */  bltz    $t1, .L80973C70            
/* 02D34 80973C84 00000000 */  nop
.L80973C88:
/* 02D38 80973C88 A0890183 */  sb      $t1, 0x0183($a0)           ## 00000183
/* 02D3C 80973C8C 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 02D40 80973C90 44C8F800 */  ctc1    $t0, $31
/* 02D44 80973C94 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 02D48 80973C98 03E00008 */  jr      $ra                        
/* 02D4C 80973C9C 00000000 */  nop


