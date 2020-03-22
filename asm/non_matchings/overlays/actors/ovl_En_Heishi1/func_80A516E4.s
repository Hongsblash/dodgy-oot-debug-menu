glabel func_80A516E4
/* 00414 80A516E4 27BDFFB8 */  addiu   $sp, $sp, 0xFFB8           ## $sp = FFFFFFB8
/* 00418 80A516E8 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 0041C 80A516EC 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00420 80A516F0 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 00424 80A516F4 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 00428 80A516F8 AFA5004C */  sw      $a1, 0x004C($sp)           
/* 0042C 80A516FC 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 00430 80A51700 AFA4002C */  sw      $a0, 0x002C($sp)           
/* 00434 80A51704 8FA4002C */  lw      $a0, 0x002C($sp)           
/* 00438 80A51708 0C0295B2 */  jal     func_800A56C8              
/* 0043C 80A5170C 3C053F80 */  lui     $a1, 0x3F80                ## $a1 = 3F800000
/* 00440 80A51710 14400004 */  bne     $v0, $zero, .L80A51724     
/* 00444 80A51714 8FA4002C */  lw      $a0, 0x002C($sp)           
/* 00448 80A51718 0C0295B2 */  jal     func_800A56C8              
/* 0044C 80A5171C 3C054188 */  lui     $a1, 0x4188                ## $a1 = 41880000
/* 00450 80A51720 10400003 */  beq     $v0, $zero, .L80A51730     
.L80A51724:
/* 00454 80A51724 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00458 80A51728 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 0045C 80A5172C 240528F0 */  addiu   $a1, $zero, 0x28F0         ## $a1 = 000028F0
.L80A51730:
/* 00460 80A51730 3C0E80A5 */  lui     $t6, %hi(D_80A527A0)       ## $t6 = 80A50000
/* 00464 80A51734 8DCE27A0 */  lw      $t6, %lo(D_80A527A0)($t6)  
/* 00468 80A51738 8FAF004C */  lw      $t7, 0x004C($sp)           
/* 0046C 80A5173C 3C180001 */  lui     $t8, 0x0001                ## $t8 = 00010000
/* 00470 80A51740 15C000D0 */  bne     $t6, $zero, .L80A51A84     
/* 00474 80A51744 030FC021 */  addu    $t8, $t8, $t7              
/* 00478 80A51748 8619026C */  lh      $t9, 0x026C($s0)           ## 0000026C
/* 0047C 80A5174C 8F181E08 */  lw      $t8, 0x1E08($t8)           ## 00011E08
/* 00480 80A51750 3C0D8016 */  lui     $t5, 0x8016                ## $t5 = 80160000
/* 00484 80A51754 001948C0 */  sll     $t1, $t9,  3               
/* 00488 80A51758 03094021 */  addu    $t0, $t8, $t1              
/* 0048C 80A5175C 8D030004 */  lw      $v1, 0x0004($t0)           ## 00000004
/* 00490 80A51760 860F02AA */  lh      $t7, 0x02AA($s0)           ## 000002AA
/* 00494 80A51764 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 00498 80A51768 00035100 */  sll     $t2, $v1,  4               
/* 0049C 80A5176C 000A5F02 */  srl     $t3, $t2, 28               
/* 004A0 80A51770 000B6080 */  sll     $t4, $t3,  2               
/* 004A4 80A51774 01AC6821 */  addu    $t5, $t5, $t4              
/* 004A8 80A51778 8DAD6FA8 */  lw      $t5, 0x6FA8($t5)           ## 80166FA8
/* 004AC 80A5177C 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 004B0 80A51780 00617024 */  and     $t6, $v1, $at              
/* 004B4 80A51784 000FC880 */  sll     $t9, $t7,  2               
/* 004B8 80A51788 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 004BC 80A5178C 032FC823 */  subu    $t9, $t9, $t7              
/* 004C0 80A51790 01AE1021 */  addu    $v0, $t5, $t6              
/* 004C4 80A51794 00411021 */  addu    $v0, $v0, $at              
/* 004C8 80A51798 0019C840 */  sll     $t9, $t9,  1               
/* 004CC 80A5179C 00591021 */  addu    $v0, $v0, $t9              
/* 004D0 80A517A0 84580000 */  lh      $t8, 0x0000($v0)           ## 00000000
/* 004D4 80A517A4 8E070274 */  lw      $a3, 0x0274($s0)           ## 00000274
/* 004D8 80A517A8 AFA80044 */  sw      $t0, 0x0044($sp)           
/* 004DC 80A517AC 44982000 */  mtc1    $t8, $f4                   ## $f4 = 0.00
/* 004E0 80A517B0 AFA20040 */  sw      $v0, 0x0040($sp)           
/* 004E4 80A517B4 26040024 */  addiu   $a0, $s0, 0x0024           ## $a0 = 00000024
/* 004E8 80A517B8 46802120 */  cvt.s.w $f4, $f4                   
/* 004EC 80A517BC 3C063F80 */  lui     $a2, 0x3F80                ## $a2 = 3F800000
/* 004F0 80A517C0 44052000 */  mfc1    $a1, $f4                   
/* 004F4 80A517C4 0C01E107 */  jal     Math_SmoothScaleMaxF
              
/* 004F8 80A517C8 00000000 */  nop
/* 004FC 80A517CC 8FA20040 */  lw      $v0, 0x0040($sp)           
/* 00500 80A517D0 2604002C */  addiu   $a0, $s0, 0x002C           ## $a0 = 0000002C
/* 00504 80A517D4 3C063F80 */  lui     $a2, 0x3F80                ## $a2 = 3F800000
/* 00508 80A517D8 84490004 */  lh      $t1, 0x0004($v0)           ## 00000004
/* 0050C 80A517DC 8E070274 */  lw      $a3, 0x0274($s0)           ## 00000274
/* 00510 80A517E0 44893000 */  mtc1    $t1, $f6                   ## $f6 = 0.00
/* 00514 80A517E4 00000000 */  nop
/* 00518 80A517E8 468031A0 */  cvt.s.w $f6, $f6                   
/* 0051C 80A517EC 44053000 */  mfc1    $a1, $f6                   
/* 00520 80A517F0 0C01E107 */  jal     Math_SmoothScaleMaxF
              
/* 00524 80A517F4 00000000 */  nop
/* 00528 80A517F8 26040274 */  addiu   $a0, $s0, 0x0274           ## $a0 = 00000274
/* 0052C 80A517FC 8E05028C */  lw      $a1, 0x028C($s0)           ## 0000028C
/* 00530 80A51800 3C063F80 */  lui     $a2, 0x3F80                ## $a2 = 3F800000
/* 00534 80A51804 0C01E107 */  jal     Math_SmoothScaleMaxF
              
/* 00538 80A51808 8E070290 */  lw      $a3, 0x0290($s0)           ## 00000290
/* 0053C 80A5180C 8FA20040 */  lw      $v0, 0x0040($sp)           
/* 00540 80A51810 C6100024 */  lwc1    $f16, 0x0024($s0)          ## 00000024
/* 00544 80A51814 C606002C */  lwc1    $f6, 0x002C($s0)           ## 0000002C
/* 00548 80A51818 844A0000 */  lh      $t2, 0x0000($v0)           ## 00000000
/* 0054C 80A5181C 844B0004 */  lh      $t3, 0x0004($v0)           ## 00000004
/* 00550 80A51820 448A4000 */  mtc1    $t2, $f8                   ## $f8 = 0.00
/* 00554 80A51824 448B9000 */  mtc1    $t3, $f18                  ## $f18 = 0.00
/* 00558 80A51828 468042A0 */  cvt.s.w $f10, $f8                  
/* 0055C 80A5182C 46809120 */  cvt.s.w $f4, $f18                  
/* 00560 80A51830 46105301 */  sub.s   $f12, $f10, $f16           
/* 00564 80A51834 46062381 */  sub.s   $f14, $f4, $f6             
/* 00568 80A51838 E7AC003C */  swc1    $f12, 0x003C($sp)          
/* 0056C 80A5183C 0C03F494 */  jal     Math_atan2f              
/* 00570 80A51840 E7AE0038 */  swc1    $f14, 0x0038($sp)          
/* 00574 80A51844 3C0180A5 */  lui     $at, %hi(D_80A52C2C)       ## $at = 80A50000
/* 00578 80A51848 C4282C2C */  lwc1    $f8, %lo(D_80A52C2C)($at)  
/* 0057C 80A5184C C6120278 */  lwc1    $f18, 0x0278($s0)          ## 00000278
/* 00580 80A51850 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 00584 80A51854 46080282 */  mul.s   $f10, $f0, $f8             
/* 00588 80A51858 4600910D */  trunc.w.s $f4, $f18                  
/* 0058C 80A5185C 260400B6 */  addiu   $a0, $s0, 0x00B6           ## $a0 = 000000B6
/* 00590 80A51860 24060003 */  addiu   $a2, $zero, 0x0003         ## $a2 = 00000003
/* 00594 80A51864 44072000 */  mfc1    $a3, $f4                   
/* 00598 80A51868 4600540D */  trunc.w.s $f16, $f10                 
/* 0059C 80A5186C 00073C00 */  sll     $a3, $a3, 16               
/* 005A0 80A51870 00073C03 */  sra     $a3, $a3, 16               
/* 005A4 80A51874 44058000 */  mfc1    $a1, $f16                  
/* 005A8 80A51878 00000000 */  nop
/* 005AC 80A5187C 00052C00 */  sll     $a1, $a1, 16               
/* 005B0 80A51880 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 005B4 80A51884 00052C03 */  sra     $a1, $a1, 16               
/* 005B8 80A51888 26040278 */  addiu   $a0, $s0, 0x0278           ## $a0 = 00000278
/* 005BC 80A5188C 8E050294 */  lw      $a1, 0x0294($s0)           ## 00000294
/* 005C0 80A51890 3C063F80 */  lui     $a2, 0x3F80                ## $a2 = 3F800000
/* 005C4 80A51894 0C01E107 */  jal     Math_SmoothScaleMaxF
              
/* 005C8 80A51898 8E070298 */  lw      $a3, 0x0298($s0)           ## 00000298
/* 005CC 80A5189C 860E02A4 */  lh      $t6, 0x02A4($s0)           ## 000002A4
/* 005D0 80A518A0 55C00020 */  bnel    $t6, $zero, .L80A51924     
/* 005D4 80A518A4 2604027C */  addiu   $a0, $s0, 0x027C           ## $a0 = 0000027C
/* 005D8 80A518A8 860F0262 */  lh      $t7, 0x0262($s0)           ## 00000262
/* 005DC 80A518AC 3C014600 */  lui     $at, 0x4600                ## $at = 46000000
/* 005E0 80A518B0 44813000 */  mtc1    $at, $f6                   ## $f6 = 8192.00
/* 005E4 80A518B4 25F90001 */  addiu   $t9, $t7, 0x0001           ## $t9 = 00000001
/* 005E8 80A518B8 A6190262 */  sh      $t9, 0x0262($s0)           ## 00000262
/* 005EC 80A518BC 86180262 */  lh      $t8, 0x0262($s0)           ## 00000262
/* 005F0 80A518C0 3C01BF80 */  lui     $at, 0xBF80                ## $at = BF800000
/* 005F4 80A518C4 E6060280 */  swc1    $f6, 0x0280($s0)           ## 00000280
/* 005F8 80A518C8 33090001 */  andi    $t1, $t8, 0x0001           ## $t1 = 00000000
/* 005FC 80A518CC 51200007 */  beql    $t1, $zero, .L80A518EC     
/* 00600 80A518D0 3C0141F0 */  lui     $at, 0x41F0                ## $at = 41F00000
/* 00604 80A518D4 C6080280 */  lwc1    $f8, 0x0280($s0)           ## 00000280
/* 00608 80A518D8 44815000 */  mtc1    $at, $f10                  ## $f10 = 30.00
/* 0060C 80A518DC 00000000 */  nop
/* 00610 80A518E0 460A4402 */  mul.s   $f16, $f8, $f10            
/* 00614 80A518E4 E6100280 */  swc1    $f16, 0x0280($s0)          ## 00000280
/* 00618 80A518E8 3C0141F0 */  lui     $at, 0x41F0                ## $at = 41F00000
.L80A518EC:
/* 0061C 80A518EC 44816000 */  mtc1    $at, $f12                  ## $f12 = 30.00
/* 00620 80A518F0 0C00CFBE */  jal     Math_Rand_ZeroFloat
              
/* 00624 80A518F4 00000000 */  nop
/* 00628 80A518F8 860A026A */  lh      $t2, 0x026A($s0)           ## 0000026A
/* 0062C 80A518FC 4600048D */  trunc.w.s $f18, $f0                  
/* 00630 80A51900 3C0C80A5 */  lui     $t4, %hi(D_80A52844)       ## $t4 = 80A50000
/* 00634 80A51904 000A5840 */  sll     $t3, $t2,  1               
/* 00638 80A51908 018B6021 */  addu    $t4, $t4, $t3              
/* 0063C 80A5190C 858C2844 */  lh      $t4, %lo(D_80A52844)($t4)  
/* 00640 80A51910 44199000 */  mfc1    $t9, $f18                  
/* 00644 80A51914 00000000 */  nop
/* 00648 80A51918 0199C021 */  addu    $t8, $t4, $t9              
/* 0064C 80A5191C A61802A4 */  sh      $t8, 0x02A4($s0)           ## 000002A4
/* 00650 80A51920 2604027C */  addiu   $a0, $s0, 0x027C           ## $a0 = 0000027C
.L80A51924:
/* 00654 80A51924 8E050280 */  lw      $a1, 0x0280($s0)           ## 00000280
/* 00658 80A51928 8E06029C */  lw      $a2, 0x029C($s0)           ## 0000029C
/* 0065C 80A5192C 0C01E107 */  jal     Math_SmoothScaleMaxF
              
/* 00660 80A51930 8E0702A0 */  lw      $a3, 0x02A0($s0)           ## 000002A0
/* 00664 80A51934 3C028016 */  lui     $v0, 0x8016                ## $v0 = 80160000
/* 00668 80A51938 8C42FA90 */  lw      $v0, -0x0570($v0)          ## 8015FA90
/* 0066C 80A5193C 8605026C */  lh      $a1, 0x026C($s0)           ## 0000026C
/* 00670 80A51940 844912D6 */  lh      $t1, 0x12D6($v0)           ## 801612D6
/* 00674 80A51944 54A90020 */  bnel    $a1, $t1, .L80A519C8       
/* 00678 80A51948 C7A0003C */  lwc1    $f0, 0x003C($sp)           
/* 0067C 80A5194C 844A12D4 */  lh      $t2, 0x12D4($v0)           ## 801612D4
/* 00680 80A51950 3C0480A5 */  lui     $a0, %hi(D_80A52B1C)       ## $a0 = 80A50000
/* 00684 80A51954 5140001C */  beql    $t2, $zero, .L80A519C8     
/* 00688 80A51958 C7A0003C */  lwc1    $f0, 0x003C($sp)           
/* 0068C 80A5195C 0C00084C */  jal     osSyncPrintf
              
/* 00690 80A51960 24842B1C */  addiu   $a0, $a0, %lo(D_80A52B1C)  ## $a0 = 80A52B1C
/* 00694 80A51964 3C0480A5 */  lui     $a0, %hi(D_80A52B30)       ## $a0 = 80A50000
/* 00698 80A51968 24842B30 */  addiu   $a0, $a0, %lo(D_80A52B30)  ## $a0 = 80A52B30
/* 0069C 80A5196C 0C00084C */  jal     osSyncPrintf
              
/* 006A0 80A51970 860502AA */  lh      $a1, 0x02AA($s0)           ## 000002AA
/* 006A4 80A51974 C6040278 */  lwc1    $f4, 0x0278($s0)           ## 00000278
/* 006A8 80A51978 3C0480A5 */  lui     $a0, %hi(D_80A52B44)       ## $a0 = 80A50000
/* 006AC 80A5197C 24842B44 */  addiu   $a0, $a0, %lo(D_80A52B44)  ## $a0 = 80A52B44
/* 006B0 80A51980 460021A1 */  cvt.d.s $f6, $f4                   
/* 006B4 80A51984 44073000 */  mfc1    $a3, $f6                   
/* 006B8 80A51988 44063800 */  mfc1    $a2, $f7                   
/* 006BC 80A5198C 0C00084C */  jal     osSyncPrintf
              
/* 006C0 80A51990 00000000 */  nop
/* 006C4 80A51994 3C0480A5 */  lui     $a0, %hi(D_80A52B58)       ## $a0 = 80A50000
/* 006C8 80A51998 24842B58 */  addiu   $a0, $a0, %lo(D_80A52B58)  ## $a0 = 80A52B58
/* 006CC 80A5199C 0C00084C */  jal     osSyncPrintf
              
/* 006D0 80A519A0 86050270 */  lh      $a1, 0x0270($s0)           ## 00000270
/* 006D4 80A519A4 8FAB0044 */  lw      $t3, 0x0044($sp)           
/* 006D8 80A519A8 3C0480A5 */  lui     $a0, %hi(D_80A52B6C)       ## $a0 = 80A50000
/* 006DC 80A519AC 24842B6C */  addiu   $a0, $a0, %lo(D_80A52B6C)  ## $a0 = 80A52B6C
/* 006E0 80A519B0 0C00084C */  jal     osSyncPrintf
              
/* 006E4 80A519B4 91650000 */  lbu     $a1, 0x0000($t3)           ## 00000000
/* 006E8 80A519B8 3C0480A5 */  lui     $a0, %hi(D_80A52B80)       ## $a0 = 80A50000
/* 006EC 80A519BC 0C00084C */  jal     osSyncPrintf
              
/* 006F0 80A519C0 24842B80 */  addiu   $a0, $a0, %lo(D_80A52B80)  ## $a0 = 80A52B80
/* 006F4 80A519C4 C7A0003C */  lwc1    $f0, 0x003C($sp)           
.L80A519C8:
/* 006F8 80A519C8 3C0141A0 */  lui     $at, 0x41A0                ## $at = 41A00000
/* 006FC 80A519CC 44811000 */  mtc1    $at, $f2                   ## $f2 = 20.00
/* 00700 80A519D0 46000005 */  abs.s   $f0, $f0                   
/* 00704 80A519D4 4602003C */  c.lt.s  $f0, $f2                   
/* 00708 80A519D8 C7A00038 */  lwc1    $f0, 0x0038($sp)           
/* 0070C 80A519DC 4502002A */  bc1fl   .L80A51A88                 
/* 00710 80A519E0 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 00714 80A519E4 46000005 */  abs.s   $f0, $f0                   
/* 00718 80A519E8 4602003C */  c.lt.s  $f0, $f2                   
/* 0071C 80A519EC 00000000 */  nop
/* 00720 80A519F0 45020025 */  bc1fl   .L80A51A88                 
/* 00724 80A519F4 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 00728 80A519F8 860D0270 */  lh      $t5, 0x0270($s0)           ## 00000270
/* 0072C 80A519FC 55A00022 */  bnel    $t5, $zero, .L80A51A88     
/* 00730 80A51A00 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 00734 80A51A04 860E026A */  lh      $t6, 0x026A($s0)           ## 0000026A
/* 00738 80A51A08 29C10002 */  slti    $at, $t6, 0x0002           
/* 0073C 80A51A0C 1420001A */  bne     $at, $zero, .L80A51A78     
/* 00740 80A51A10 00000000 */  nop
/* 00744 80A51A14 860F02AA */  lh      $t7, 0x02AA($s0)           ## 000002AA
/* 00748 80A51A18 29E10004 */  slti    $at, $t7, 0x0004           
/* 0074C 80A51A1C 14200016 */  bne     $at, $zero, .L80A51A78     
/* 00750 80A51A20 3C0180A5 */  lui     $at, %hi(D_80A52C30)       ## $at = 80A50000
/* 00754 80A51A24 0C00CFBE */  jal     Math_Rand_ZeroFloat
              
/* 00758 80A51A28 C42C2C30 */  lwc1    $f12, %lo(D_80A52C30)($at) 
/* 0075C 80A51A2C 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 00760 80A51A30 44814000 */  mtc1    $at, $f8                   ## $f8 = 1.00
/* 00764 80A51A34 24190005 */  addiu   $t9, $zero, 0x0005         ## $t9 = 00000005
/* 00768 80A51A38 4600403C */  c.lt.s  $f8, $f0                   
/* 0076C 80A51A3C 00000000 */  nop
/* 00770 80A51A40 4500000D */  bc1f    .L80A51A78                 
/* 00774 80A51A44 00000000 */  nop
/* 00778 80A51A48 860202AA */  lh      $v0, 0x02AA($s0)           ## 000002AA
/* 0077C 80A51A4C 24010007 */  addiu   $at, $zero, 0x0007         ## $at = 00000007
/* 00780 80A51A50 54410004 */  bnel    $v0, $at, .L80A51A64       
/* 00784 80A51A54 28410004 */  slti    $at, $v0, 0x0004           
/* 00788 80A51A58 A60002AA */  sh      $zero, 0x02AA($s0)         ## 000002AA
/* 0078C 80A51A5C 860202AA */  lh      $v0, 0x02AA($s0)           ## 000002AA
/* 00790 80A51A60 28410004 */  slti    $at, $v0, 0x0004           
.L80A51A64:
/* 00794 80A51A64 14200002 */  bne     $at, $zero, .L80A51A70     
/* 00798 80A51A68 244CFFFD */  addiu   $t4, $v0, 0xFFFD           ## $t4 = FFFFFFFD
/* 0079C 80A51A6C A60C02AA */  sh      $t4, 0x02AA($s0)           ## 000002AA
.L80A51A70:
/* 007A0 80A51A70 10000004 */  beq     $zero, $zero, .L80A51A84   
/* 007A4 80A51A74 A6190270 */  sh      $t9, 0x0270($s0)           ## 00000270
.L80A51A78:
/* 007A8 80A51A78 3C1880A5 */  lui     $t8, %hi(func_80A51C4C)    ## $t8 = 80A50000
/* 007AC 80A51A7C 27181C4C */  addiu   $t8, $t8, %lo(func_80A51C4C) ## $t8 = 80A51C4C
/* 007B0 80A51A80 AE18025C */  sw      $t8, 0x025C($s0)           ## 0000025C
.L80A51A84:
/* 007B4 80A51A84 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80A51A88:
/* 007B8 80A51A88 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 007BC 80A51A8C 27BD0048 */  addiu   $sp, $sp, 0x0048           ## $sp = 00000000
/* 007C0 80A51A90 03E00008 */  jr      $ra                        
/* 007C4 80A51A94 00000000 */  nop


