.late_rodata
glabel jtbl_8098CA9C
.word L809896CC
.word L809896CC
.word L809892E8
.word L809892E8
.word L809892E8
.word L809892E8
.word L809892E8
.word L80989300
.word L809896CC
.word L809896CC
.word L809896CC
.word L809896CC
.word L809896CC
.word L80989318
.word L809896CC
.word L80989348
.word L8098964C
.word L809896CC

.text
glabel func_809892A4
/* 00424 809892A4 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00428 809892A8 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 0042C 809892AC 808F001E */  lb      $t7, 0x001E($a0)           ## 0000001E
/* 00430 809892B0 908E014C */  lbu     $t6, 0x014C($a0)           ## 0000014C
/* 00434 809892B4 00803825 */  or      $a3, $a0, $zero            ## $a3 = 00000000
/* 00438 809892B8 00A03025 */  or      $a2, $a1, $zero            ## $a2 = 00000000
/* 0043C 809892BC 15CF0103 */  bne     $t6, $t7, .L809896CC       
/* 00440 809892C0 8CA81C44 */  lw      $t0, 0x1C44($a1)           ## 00001C44
/* 00444 809892C4 9498001C */  lhu     $t8, 0x001C($a0)           ## 0000001C
/* 00448 809892C8 2F010012 */  sltiu   $at, $t8, 0x0012           
/* 0044C 809892CC 102000FF */  beq     $at, $zero, .L809896CC     
/* 00450 809892D0 0018C080 */  sll     $t8, $t8,  2               
/* 00454 809892D4 3C018099 */  lui     $at, %hi(jtbl_8098CA9C)       ## $at = 80990000
/* 00458 809892D8 00380821 */  addu    $at, $at, $t8              
/* 0045C 809892DC 8C38CA9C */  lw      $t8, %lo(jtbl_8098CA9C)($at)  
/* 00460 809892E0 03000008 */  jr      $t8                        
/* 00464 809892E4 00000000 */  nop
glabel L809892E8
/* 00468 809892E8 3C058099 */  lui     $a1, %hi(func_80989800)    ## $a1 = 80990000
/* 0046C 809892EC 24A59800 */  addiu   $a1, $a1, %lo(func_80989800) ## $a1 = 80989800
/* 00470 809892F0 0C2623A0 */  jal     func_80988E80              
/* 00474 809892F4 00E02025 */  or      $a0, $a3, $zero            ## $a0 = 00000000
/* 00478 809892F8 100000F5 */  beq     $zero, $zero, .L809896D0   
/* 0047C 809892FC 8FBF0014 */  lw      $ra, 0x0014($sp)           
glabel L80989300
/* 00480 80989300 3C058099 */  lui     $a1, %hi(func_8098987C)    ## $a1 = 80990000
/* 00484 80989304 24A5987C */  addiu   $a1, $a1, %lo(func_8098987C) ## $a1 = 8098987C
/* 00488 80989308 0C2623A0 */  jal     func_80988E80              
/* 0048C 8098930C 00E02025 */  or      $a0, $a3, $zero            ## $a0 = 00000000
/* 00490 80989310 100000EF */  beq     $zero, $zero, .L809896D0   
/* 00494 80989314 8FBF0014 */  lw      $ra, 0x0014($sp)           
glabel L80989318
/* 00498 80989318 00C02025 */  or      $a0, $a2, $zero            ## $a0 = 00000000
/* 0049C 8098931C 24050002 */  addiu   $a1, $zero, 0x0002         ## $a1 = 00000002
/* 004A0 80989320 0C01B129 */  jal     func_8006C4A4              
/* 004A4 80989324 AFA70018 */  sw      $a3, 0x0018($sp)           
/* 004A8 80989328 104000E8 */  beq     $v0, $zero, .L809896CC     
/* 004AC 8098932C 8FA70018 */  lw      $a3, 0x0018($sp)           
/* 004B0 80989330 3C058099 */  lui     $a1, %hi(func_809898C8)    ## $a1 = 80990000
/* 004B4 80989334 24A598C8 */  addiu   $a1, $a1, %lo(func_809898C8) ## $a1 = 809898C8
/* 004B8 80989338 0C2623A0 */  jal     func_80988E80              
/* 004BC 8098933C 00E02025 */  or      $a0, $a3, $zero            ## $a0 = 00000000
/* 004C0 80989340 100000E3 */  beq     $zero, $zero, .L809896D0   
/* 004C4 80989344 8FBF0014 */  lw      $ra, 0x0014($sp)           
glabel L80989348
/* 004C8 80989348 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 004CC 8098934C 00C12021 */  addu    $a0, $a2, $at              
/* 004D0 80989350 240500FF */  addiu   $a1, $zero, 0x00FF         ## $a1 = 000000FF
/* 004D4 80989354 A0850B06 */  sb      $a1, 0x0B06($a0)           ## 00000B06
/* 004D8 80989358 A0850B07 */  sb      $a1, 0x0B07($a0)           ## 00000B07
/* 004DC 8098935C A0850B08 */  sb      $a1, 0x0B08($a0)           ## 00000B08
/* 004E0 80989360 00260821 */  addu    $at, $at, $a2              
/* 004E4 80989364 A0200B05 */  sb      $zero, 0x0B05($at)         ## 00010B05
/* 004E8 80989368 90E2014E */  lbu     $v0, 0x014E($a3)           ## 0000014E
/* 004EC 8098936C 28410015 */  slti    $at, $v0, 0x0015           
/* 004F0 80989370 1020003E */  beq     $at, $zero, .L8098946C     
/* 004F4 80989374 00401825 */  or      $v1, $v0, $zero            ## $v1 = 00000000
/* 004F8 80989378 2861000F */  slti    $at, $v1, 0x000F           
/* 004FC 8098937C 5420003C */  bnel    $at, $zero, .L80989470     
/* 00500 80989380 3C01437F */  lui     $at, 0x437F                ## $at = 437F0000
/* 00504 80989384 44822000 */  mtc1    $v0, $f4                   ## $f4 = 0.00
/* 00508 80989388 3C01437F */  lui     $at, 0x437F                ## $at = 437F0000
/* 0050C 8098938C 44811000 */  mtc1    $at, $f2                   ## $f2 = 255.00
/* 00510 80989390 04410005 */  bgez    $v0, .L809893A8            
/* 00514 80989394 468021A0 */  cvt.s.w $f6, $f4                   
/* 00518 80989398 3C014F80 */  lui     $at, 0x4F80                ## $at = 4F800000
/* 0051C 8098939C 44814000 */  mtc1    $at, $f8                   ## $f8 = 4294967296.00
/* 00520 809893A0 00000000 */  nop
/* 00524 809893A4 46083180 */  add.s   $f6, $f6, $f8              
.L809893A8:
/* 00528 809893A8 3C014170 */  lui     $at, 0x4170                ## $at = 41700000
/* 0052C 809893AC 44815000 */  mtc1    $at, $f10                  ## $f10 = 15.00
/* 00530 809893B0 3C0140A0 */  lui     $at, 0x40A0                ## $at = 40A00000
/* 00534 809893B4 44819000 */  mtc1    $at, $f18                  ## $f18 = 5.00
/* 00538 809893B8 460A3401 */  sub.s   $f16, $f6, $f10            
/* 0053C 809893BC 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 00540 809893C0 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00544 809893C4 00260821 */  addu    $at, $at, $a2              
/* 00548 809893C8 46128003 */  div.s   $f0, $f16, $f18            
/* 0054C 809893CC 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 00550 809893D0 A0390B05 */  sb      $t9, 0x0B05($at)           ## 00010B05
/* 00554 809893D4 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 00558 809893D8 46001102 */  mul.s   $f4, $f2, $f0              
/* 0055C 809893DC 46041201 */  sub.s   $f8, $f2, $f4              
/* 00560 809893E0 4449F800 */  cfc1    $t1, $31
/* 00564 809893E4 44CAF800 */  ctc1    $t2, $31
/* 00568 809893E8 00000000 */  nop
/* 0056C 809893EC 460041A4 */  cvt.w.s $f6, $f8                   
/* 00570 809893F0 444AF800 */  cfc1    $t2, $31
/* 00574 809893F4 00000000 */  nop
/* 00578 809893F8 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 0057C 809893FC 51400013 */  beql    $t2, $zero, .L8098944C     
/* 00580 80989400 440A3000 */  mfc1    $t2, $f6                   
/* 00584 80989404 44813000 */  mtc1    $at, $f6                   ## $f6 = 2147483648.00
/* 00588 80989408 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 0058C 8098940C 46064181 */  sub.s   $f6, $f8, $f6              
/* 00590 80989410 44CAF800 */  ctc1    $t2, $31
/* 00594 80989414 00000000 */  nop
/* 00598 80989418 460031A4 */  cvt.w.s $f6, $f6                   
/* 0059C 8098941C 444AF800 */  cfc1    $t2, $31
/* 005A0 80989420 00000000 */  nop
/* 005A4 80989424 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 005A8 80989428 15400005 */  bne     $t2, $zero, .L80989440     
/* 005AC 8098942C 00000000 */  nop
/* 005B0 80989430 440A3000 */  mfc1    $t2, $f6                   
/* 005B4 80989434 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 005B8 80989438 10000007 */  beq     $zero, $zero, .L80989458   
/* 005BC 8098943C 01415025 */  or      $t2, $t2, $at              ## $t2 = 80000000
.L80989440:
/* 005C0 80989440 10000005 */  beq     $zero, $zero, .L80989458   
/* 005C4 80989444 240AFFFF */  addiu   $t2, $zero, 0xFFFF         ## $t2 = FFFFFFFF
/* 005C8 80989448 440A3000 */  mfc1    $t2, $f6                   
.L8098944C:
/* 005CC 8098944C 00000000 */  nop
/* 005D0 80989450 0540FFFB */  bltz    $t2, .L80989440            
/* 005D4 80989454 00000000 */  nop
.L80989458:
/* 005D8 80989458 A08A0B09 */  sb      $t2, 0x0B09($a0)           ## 00000B09
/* 005DC 8098945C 90E2014E */  lbu     $v0, 0x014E($a3)           ## 0000014E
/* 005E0 80989460 44C9F800 */  ctc1    $t1, $31
/* 005E4 80989464 00401825 */  or      $v1, $v0, $zero            ## $v1 = 00000000
/* 005E8 80989468 00000000 */  nop
.L8098946C:
/* 005EC 8098946C 3C01437F */  lui     $at, 0x437F                ## $at = 437F0000
.L80989470:
/* 005F0 80989470 44811000 */  mtc1    $at, $f2                   ## $f2 = 255.00
/* 005F4 80989474 2861000F */  slti    $at, $v1, 0x000F           
/* 005F8 80989478 1020003A */  beq     $at, $zero, .L80989564     
/* 005FC 8098947C 28610004 */  slti    $at, $v1, 0x0004           
/* 00600 80989480 54200039 */  bnel    $at, $zero, .L80989568     
/* 00604 80989484 2401000F */  addiu   $at, $zero, 0x000F         ## $at = 0000000F
/* 00608 80989488 44825000 */  mtc1    $v0, $f10                  ## $f10 = 0.00
/* 0060C 8098948C 3C014F80 */  lui     $at, 0x4F80                ## $at = 4F800000
/* 00610 80989490 04410004 */  bgez    $v0, .L809894A4            
/* 00614 80989494 46805420 */  cvt.s.w $f16, $f10                 
/* 00618 80989498 44819000 */  mtc1    $at, $f18                  ## $f18 = 4294967296.00
/* 0061C 8098949C 00000000 */  nop
/* 00620 809894A0 46128400 */  add.s   $f16, $f16, $f18           
.L809894A4:
/* 00624 809894A4 3C014080 */  lui     $at, 0x4080                ## $at = 40800000
/* 00628 809894A8 44812000 */  mtc1    $at, $f4                   ## $f4 = 4.00
/* 0062C 809894AC 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 00630 809894B0 44813000 */  mtc1    $at, $f6                   ## $f6 = 10.00
/* 00634 809894B4 46048201 */  sub.s   $f8, $f16, $f4             
/* 00638 809894B8 240D0001 */  addiu   $t5, $zero, 0x0001         ## $t5 = 00000001
/* 0063C 809894BC 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00640 809894C0 00260821 */  addu    $at, $at, $a2              
/* 00644 809894C4 46064003 */  div.s   $f0, $f8, $f6              
/* 00648 809894C8 240B0001 */  addiu   $t3, $zero, 0x0001         ## $t3 = 00000001
/* 0064C 809894CC A02B0B05 */  sb      $t3, 0x0B05($at)           ## 00010B05
/* 00650 809894D0 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 00654 809894D4 46001282 */  mul.s   $f10, $f2, $f0             
/* 00658 809894D8 444CF800 */  cfc1    $t4, $31
/* 0065C 809894DC 44CDF800 */  ctc1    $t5, $31
/* 00660 809894E0 00000000 */  nop
/* 00664 809894E4 460054A4 */  cvt.w.s $f18, $f10                 
/* 00668 809894E8 444DF800 */  cfc1    $t5, $31
/* 0066C 809894EC 00000000 */  nop
/* 00670 809894F0 31AD0078 */  andi    $t5, $t5, 0x0078           ## $t5 = 00000000
/* 00674 809894F4 51A00013 */  beql    $t5, $zero, .L80989544     
/* 00678 809894F8 440D9000 */  mfc1    $t5, $f18                  
/* 0067C 809894FC 44819000 */  mtc1    $at, $f18                  ## $f18 = 2147483648.00
/* 00680 80989500 240D0001 */  addiu   $t5, $zero, 0x0001         ## $t5 = 00000001
/* 00684 80989504 46125481 */  sub.s   $f18, $f10, $f18           
/* 00688 80989508 44CDF800 */  ctc1    $t5, $31
/* 0068C 8098950C 00000000 */  nop
/* 00690 80989510 460094A4 */  cvt.w.s $f18, $f18                 
/* 00694 80989514 444DF800 */  cfc1    $t5, $31
/* 00698 80989518 00000000 */  nop
/* 0069C 8098951C 31AD0078 */  andi    $t5, $t5, 0x0078           ## $t5 = 00000000
/* 006A0 80989520 15A00005 */  bne     $t5, $zero, .L80989538     
/* 006A4 80989524 00000000 */  nop
/* 006A8 80989528 440D9000 */  mfc1    $t5, $f18                  
/* 006AC 8098952C 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 006B0 80989530 10000007 */  beq     $zero, $zero, .L80989550   
/* 006B4 80989534 01A16825 */  or      $t5, $t5, $at              ## $t5 = 80000000
.L80989538:
/* 006B8 80989538 10000005 */  beq     $zero, $zero, .L80989550   
/* 006BC 8098953C 240DFFFF */  addiu   $t5, $zero, 0xFFFF         ## $t5 = FFFFFFFF
/* 006C0 80989540 440D9000 */  mfc1    $t5, $f18                  
.L80989544:
/* 006C4 80989544 00000000 */  nop
/* 006C8 80989548 05A0FFFB */  bltz    $t5, .L80989538            
/* 006CC 8098954C 00000000 */  nop
.L80989550:
/* 006D0 80989550 A08D0B09 */  sb      $t5, 0x0B09($a0)           ## 00000B09
/* 006D4 80989554 90E2014E */  lbu     $v0, 0x014E($a3)           ## 0000014E
/* 006D8 80989558 44CCF800 */  ctc1    $t4, $31
/* 006DC 8098955C 00401825 */  or      $v1, $v0, $zero            ## $v1 = 00000000
/* 006E0 80989560 00000000 */  nop
.L80989564:
/* 006E4 80989564 2401000F */  addiu   $at, $zero, 0x000F         ## $at = 0000000F
.L80989568:
/* 006E8 80989568 14610003 */  bne     $v1, $at, .L80989578       
/* 006EC 8098956C 00C02025 */  or      $a0, $a2, $zero            ## $a0 = 00000000
/* 006F0 80989570 AD000134 */  sw      $zero, 0x0134($t0)         ## 00000134
/* 006F4 80989574 90E2014E */  lbu     $v0, 0x014E($a3)           ## 0000014E
.L80989578:
/* 006F8 80989578 10400004 */  beq     $v0, $zero, .L8098958C     
/* 006FC 8098957C 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 00700 80989580 244EFFFF */  addiu   $t6, $v0, 0xFFFF           ## $t6 = FFFFFFFF
/* 00704 80989584 A0EE014E */  sb      $t6, 0x014E($a3)           ## 0000014E
/* 00708 80989588 31C200FF */  andi    $v0, $t6, 0x00FF           ## $v0 = 000000FF
.L8098958C:
/* 0070C 8098958C 54410050 */  bnel    $v0, $at, .L809896D0       
/* 00710 80989590 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 00714 80989594 84CF00A4 */  lh      $t7, 0x00A4($a2)           ## 000000A4
/* 00718 80989598 24010043 */  addiu   $at, $zero, 0x0043         ## $at = 00000043
/* 0071C 8098959C 3C038016 */  lui     $v1, %hi(gSaveContext)
/* 00720 809895A0 15E10010 */  bne     $t7, $at, .L809895E4       
/* 00724 809895A4 2463E660 */  addiu   $v1, %lo(gSaveContext)
/* 00728 809895A8 3C038016 */  lui     $v1, %hi(gSaveContext)
/* 0072C 809895AC 24180019 */  addiu   $t8, $zero, 0x0019         ## $t8 = 00000019
/* 00730 809895B0 3C018099 */  lui     $at, %hi(D_8098CF84)       ## $at = 80990000
/* 00734 809895B4 2463E660 */  addiu   $v1, %lo(gSaveContext)
/* 00738 809895B8 A438CF84 */  sh      $t8, %lo(D_8098CF84)($at)  
/* 0073C 809895BC 8C790004 */  lw      $t9, 0x0004($v1)           ## 8015E664
/* 00740 809895C0 3C0A8099 */  lui     $t2, %hi(D_8098BBA0)       ## $t2 = 80990000
/* 00744 809895C4 3C098099 */  lui     $t1, %hi(D_8098C080)       ## $t1 = 80990000
/* 00748 809895C8 13200004 */  beq     $t9, $zero, .L809895DC     
/* 0074C 809895CC 254ABBA0 */  addiu   $t2, $t2, %lo(D_8098BBA0)  ## $t2 = 8098BBA0
/* 00750 809895D0 2529C080 */  addiu   $t1, $t1, %lo(D_8098C080)  ## $t1 = 8098C080
/* 00754 809895D4 1000000F */  beq     $zero, $zero, .L80989614   
/* 00758 809895D8 ACC91D68 */  sw      $t1, 0x1D68($a2)           ## 00001D68
.L809895DC:
/* 0075C 809895DC 1000000D */  beq     $zero, $zero, .L80989614   
/* 00760 809895E0 ACCA1D68 */  sw      $t2, 0x1D68($a2)           ## 00001D68
.L809895E4:
/* 00764 809895E4 240B0020 */  addiu   $t3, $zero, 0x0020         ## $t3 = 00000020
/* 00768 809895E8 3C018099 */  lui     $at, %hi(D_8098CF84)       ## $at = 80990000
/* 0076C 809895EC A42BCF84 */  sh      $t3, %lo(D_8098CF84)($at)  
/* 00770 809895F0 8C6C0004 */  lw      $t4, 0x0004($v1)           ## 8015E664
/* 00774 809895F4 3C0E8099 */  lui     $t6, %hi(D_8098B910)       ## $t6 = 80990000
/* 00778 809895F8 3C0D8099 */  lui     $t5, %hi(D_8098BDD0)       ## $t5 = 80990000
/* 0077C 809895FC 11800004 */  beq     $t4, $zero, .L80989610     
/* 00780 80989600 25CEB910 */  addiu   $t6, $t6, %lo(D_8098B910)  ## $t6 = 8098B910
/* 00784 80989604 25ADBDD0 */  addiu   $t5, $t5, %lo(D_8098BDD0)  ## $t5 = 8098BDD0
/* 00788 80989608 10000002 */  beq     $zero, $zero, .L80989614   
/* 0078C 8098960C ACCD1D68 */  sw      $t5, 0x1D68($a2)           ## 00001D68
.L80989610:
/* 00790 80989610 ACCE1D68 */  sw      $t6, 0x1D68($a2)           ## 00001D68
.L80989614:
/* 00794 80989614 0C03032E */  jal     func_800C0CB8              
/* 00798 80989618 AFA70018 */  sw      $a3, 0x0018($sp)           
/* 0079C 8098961C 3C038016 */  lui     $v1, %hi(gSaveContext)
/* 007A0 80989620 2463E660 */  addiu   $v1, %lo(gSaveContext)
/* 007A4 80989624 10400003 */  beq     $v0, $zero, .L80989634     
/* 007A8 80989628 8FA70018 */  lw      $a3, 0x0018($sp)           
/* 007AC 8098962C 240F0001 */  addiu   $t7, $zero, 0x0001         ## $t7 = 00000001
/* 007B0 80989630 A06F1414 */  sb      $t7, 0x1414($v1)           ## 8015FA74
.L80989634:
/* 007B4 80989634 3C058099 */  lui     $a1, %hi(func_809896DC)    ## $a1 = 80990000
/* 007B8 80989638 24A596DC */  addiu   $a1, $a1, %lo(func_809896DC) ## $a1 = 809896DC
/* 007BC 8098963C 0C2623A0 */  jal     func_80988E80              
/* 007C0 80989640 00E02025 */  or      $a0, $a3, $zero            ## $a0 = 00000000
/* 007C4 80989644 10000022 */  beq     $zero, $zero, .L809896D0   
/* 007C8 80989648 8FBF0014 */  lw      $ra, 0x0014($sp)           
glabel L8098964C
/* 007CC 8098964C 84D800A4 */  lh      $t8, 0x00A4($a2)           ## 000000A4
/* 007D0 80989650 24010043 */  addiu   $at, $zero, 0x0043         ## $at = 00000043
/* 007D4 80989654 240E0001 */  addiu   $t6, $zero, 0x0001         ## $t6 = 00000001
/* 007D8 80989658 1701000D */  bne     $t8, $at, .L80989690       
/* 007DC 8098965C 00E02025 */  or      $a0, $a3, $zero            ## $a0 = 00000000
/* 007E0 80989660 3C038016 */  lui     $v1, %hi(gSaveContext)
/* 007E4 80989664 2463E660 */  addiu   $v1, %lo(gSaveContext)
/* 007E8 80989668 8C790004 */  lw      $t9, 0x0004($v1)           ## 8015E664
/* 007EC 8098966C 3C0A8099 */  lui     $t2, %hi(D_8098BCB0)       ## $t2 = 80990000
/* 007F0 80989670 3C098099 */  lui     $t1, %hi(D_8098C1B0)       ## $t1 = 80990000
/* 007F4 80989674 13200004 */  beq     $t9, $zero, .L80989688     
/* 007F8 80989678 254ABCB0 */  addiu   $t2, $t2, %lo(D_8098BCB0)  ## $t2 = 8098BCB0
/* 007FC 8098967C 2529C1B0 */  addiu   $t1, $t1, %lo(D_8098C1B0)  ## $t1 = 8098C1B0
/* 00800 80989680 1000000E */  beq     $zero, $zero, .L809896BC   
/* 00804 80989684 ACC91D68 */  sw      $t1, 0x1D68($a2)           ## 00001D68
.L80989688:
/* 00808 80989688 1000000C */  beq     $zero, $zero, .L809896BC   
/* 0080C 8098968C ACCA1D68 */  sw      $t2, 0x1D68($a2)           ## 00001D68
.L80989690:
/* 00810 80989690 3C038016 */  lui     $v1, %hi(gSaveContext)
/* 00814 80989694 2463E660 */  addiu   $v1, %lo(gSaveContext)
/* 00818 80989698 8C6B0004 */  lw      $t3, 0x0004($v1)           ## 8015E664
/* 0081C 8098969C 3C0D8099 */  lui     $t5, %hi(D_8098BA20)       ## $t5 = 80990000
/* 00820 809896A0 3C0C8099 */  lui     $t4, %hi(D_8098BF00)       ## $t4 = 80990000
/* 00824 809896A4 11600004 */  beq     $t3, $zero, .L809896B8     
/* 00828 809896A8 25ADBA20 */  addiu   $t5, $t5, %lo(D_8098BA20)  ## $t5 = 8098BA20
/* 0082C 809896AC 258CBF00 */  addiu   $t4, $t4, %lo(D_8098BF00)  ## $t4 = 8098BF00
/* 00830 809896B0 10000002 */  beq     $zero, $zero, .L809896BC   
/* 00834 809896B4 ACCC1D68 */  sw      $t4, 0x1D68($a2)           ## 00001D68
.L809896B8:
/* 00838 809896B8 ACCD1D68 */  sw      $t5, 0x1D68($a2)           ## 00001D68
.L809896BC:
/* 0083C 809896BC 3C058099 */  lui     $a1, %hi(func_809896E8)    ## $a1 = 80990000
/* 00840 809896C0 A06E1414 */  sb      $t6, 0x1414($v1)           ## 8015FA74
/* 00844 809896C4 0C2623A0 */  jal     func_80988E80              
/* 00848 809896C8 24A596E8 */  addiu   $a1, $a1, %lo(func_809896E8) ## $a1 = 809896E8
glabel L809896CC
.L809896CC:
/* 0084C 809896CC 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L809896D0:
/* 00850 809896D0 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 00854 809896D4 03E00008 */  jr      $ra                        
/* 00858 809896D8 00000000 */  nop
