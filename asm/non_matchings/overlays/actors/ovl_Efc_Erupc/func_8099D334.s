glabel func_8099D334
/* 00684 8099D334 27BDFF70 */  addiu   $sp, $sp, 0xFF70           ## $sp = FFFFFF70
/* 00688 8099D338 AFBF0044 */  sw      $ra, 0x0044($sp)           
/* 0068C 8099D33C AFBE0040 */  sw      $s8, 0x0040($sp)           
/* 00690 8099D340 AFB7003C */  sw      $s7, 0x003C($sp)           
/* 00694 8099D344 AFB60038 */  sw      $s6, 0x0038($sp)           
/* 00698 8099D348 AFB50034 */  sw      $s5, 0x0034($sp)           
/* 0069C 8099D34C AFB40030 */  sw      $s4, 0x0030($sp)           
/* 006A0 8099D350 AFB3002C */  sw      $s3, 0x002C($sp)           
/* 006A4 8099D354 AFB20028 */  sw      $s2, 0x0028($sp)           
/* 006A8 8099D358 AFB10024 */  sw      $s1, 0x0024($sp)           
/* 006AC 8099D35C AFB00020 */  sw      $s0, 0x0020($sp)           
/* 006B0 8099D360 F7B40018 */  sdc1    $f20, 0x0018($sp)          
/* 006B4 8099D364 8CB00000 */  lw      $s0, 0x0000($a1)           ## 00000000
/* 006B8 8099D368 00808825 */  or      $s1, $a0, $zero            ## $s1 = 00000000
/* 006BC 8099D36C 00A0B025 */  or      $s6, $a1, $zero            ## $s6 = 00000000
/* 006C0 8099D370 3C06809A */  lui     $a2, %hi(D_8099D7E0)       ## $a2 = 809A0000
/* 006C4 8099D374 24C6D7E0 */  addiu   $a2, $a2, %lo(D_8099D7E0)  ## $a2 = 8099D7E0
/* 006C8 8099D378 27A40070 */  addiu   $a0, $sp, 0x0070           ## $a0 = FFFFFFE0
/* 006CC 8099D37C 24070170 */  addiu   $a3, $zero, 0x0170         ## $a3 = 00000170
/* 006D0 8099D380 0C031AB1 */  jal     Graph_OpenDisps              
/* 006D4 8099D384 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 006D8 8099D388 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 006DC 8099D38C 4481A000 */  mtc1    $at, $f20                  ## $f20 = 1.00
/* 006E0 8099D390 00009825 */  or      $s3, $zero, $zero          ## $s3 = 00000000
/* 006E4 8099D394 3C1EFA00 */  lui     $s8, 0xFA00                ## $s8 = FA000000
/* 006E8 8099D398 3C17DE00 */  lui     $s7, 0xDE00                ## $s7 = DE000000
.L8099D39C:
/* 006EC 8099D39C 922E0024 */  lbu     $t6, 0x0024($s1)           ## 00000024
/* 006F0 8099D3A0 3C120600 */  lui     $s2, 0x0600                ## $s2 = 06000000
/* 006F4 8099D3A4 26522760 */  addiu   $s2, $s2, 0x2760           ## $s2 = 06002760
/* 006F8 8099D3A8 11C00047 */  beq     $t6, $zero, .L8099D4C8     
/* 006FC 8099D3AC 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00700 8099D3B0 34211DA0 */  ori     $at, $at, 0x1DA0           ## $at = 00011DA0
/* 00704 8099D3B4 3C150600 */  lui     $s5, 0x0600                ## $s5 = 06000000
/* 00708 8099D3B8 26B527D8 */  addiu   $s5, $s5, 0x27D8           ## $s5 = 060027D8
/* 0070C 8099D3BC 02C1A021 */  addu    $s4, $s6, $at              
/* 00710 8099D3C0 0C024F61 */  jal     func_80093D84              
/* 00714 8099D3C4 8EC40000 */  lw      $a0, 0x0000($s6)           ## 00000000
/* 00718 8099D3C8 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 0071C 8099D3CC 3C0E9600 */  lui     $t6, 0x9600                ## $t6 = 96000000
/* 00720 8099D3D0 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 00724 8099D3D4 244F0008 */  addiu   $t7, $v0, 0x0008           ## $t7 = 00000008
/* 00728 8099D3D8 AE0F02D0 */  sw      $t7, 0x02D0($s0)           ## 000002D0
/* 0072C 8099D3DC AC520004 */  sw      $s2, 0x0004($v0)           ## 00000004
/* 00730 8099D3E0 AC570000 */  sw      $s7, 0x0000($v0)           ## 00000000
/* 00734 8099D3E4 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 00738 8099D3E8 24580008 */  addiu   $t8, $v0, 0x0008           ## $t8 = 00000008
/* 0073C 8099D3EC AE1802D0 */  sw      $t8, 0x02D0($s0)           ## 000002D0
/* 00740 8099D3F0 AC5E0000 */  sw      $s8, 0x0000($v0)           ## 00000000
/* 00744 8099D3F4 92280026 */  lbu     $t0, 0x0026($s1)           ## 00000026
/* 00748 8099D3F8 922B0027 */  lbu     $t3, 0x0027($s1)           ## 00000027
/* 0074C 8099D3FC 922F0028 */  lbu     $t7, 0x0028($s1)           ## 00000028
/* 00750 8099D400 00084E00 */  sll     $t1, $t0, 24               
/* 00754 8099D404 8628002A */  lh      $t0, 0x002A($s1)           ## 0000002A
/* 00758 8099D408 000B6400 */  sll     $t4, $t3, 16               
/* 0075C 8099D40C 012C6825 */  or      $t5, $t1, $t4              ## $t5 = 00000000
/* 00760 8099D410 000FC200 */  sll     $t8, $t7,  8               
/* 00764 8099D414 01B8C825 */  or      $t9, $t5, $t8              ## $t9 = 00000008
/* 00768 8099D418 310A00FF */  andi    $t2, $t0, 0x00FF           ## $t2 = 00000000
/* 0076C 8099D41C 032A5825 */  or      $t3, $t9, $t2              ## $t3 = 00000008
/* 00770 8099D420 AC4B0004 */  sw      $t3, 0x0004($v0)           ## 00000004
/* 00774 8099D424 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 00778 8099D428 3C0CFB00 */  lui     $t4, 0xFB00                ## $t4 = FB000000
/* 0077C 8099D42C 3C0DE700 */  lui     $t5, 0xE700                ## $t5 = E7000000
/* 00780 8099D430 24490008 */  addiu   $t1, $v0, 0x0008           ## $t1 = 00000008
/* 00784 8099D434 AE0902D0 */  sw      $t1, 0x02D0($s0)           ## 000002D0
/* 00788 8099D438 AC4E0004 */  sw      $t6, 0x0004($v0)           ## 00000004
/* 0078C 8099D43C AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
/* 00790 8099D440 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 00794 8099D444 244F0008 */  addiu   $t7, $v0, 0x0008           ## $t7 = 00000008
/* 00798 8099D448 AE0F02D0 */  sw      $t7, 0x02D0($s0)           ## 000002D0
/* 0079C 8099D44C AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 007A0 8099D450 AC4D0000 */  sw      $t5, 0x0000($v0)           ## 00000000
/* 007A4 8099D454 8E260008 */  lw      $a2, 0x0008($s1)           ## 00000008
/* 007A8 8099D458 C62E0004 */  lwc1    $f14, 0x0004($s1)          ## 00000004
/* 007AC 8099D45C 0C034261 */  jal     Matrix_Translate              
/* 007B0 8099D460 C62C0000 */  lwc1    $f12, 0x0000($s1)          ## 00000000
/* 007B4 8099D464 0C0347F5 */  jal     func_800D1FD4              
/* 007B8 8099D468 02802025 */  or      $a0, $s4, $zero            ## $a0 = 00000000
/* 007BC 8099D46C C62C0030 */  lwc1    $f12, 0x0030($s1)          ## 00000030
/* 007C0 8099D470 4406A000 */  mfc1    $a2, $f20                  
/* 007C4 8099D474 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 007C8 8099D478 0C0342A3 */  jal     Matrix_Scale              
/* 007CC 8099D47C 46006386 */  mov.s   $f14, $f12                 
/* 007D0 8099D480 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 007D4 8099D484 3C08DA38 */  lui     $t0, 0xDA38                ## $t0 = DA380000
/* 007D8 8099D488 35080003 */  ori     $t0, $t0, 0x0003           ## $t0 = DA380003
/* 007DC 8099D48C 24580008 */  addiu   $t8, $v0, 0x0008           ## $t8 = 00000008
/* 007E0 8099D490 AE1802D0 */  sw      $t8, 0x02D0($s0)           ## 000002D0
/* 007E4 8099D494 3C05809A */  lui     $a1, %hi(D_8099D7F4)       ## $a1 = 809A0000
/* 007E8 8099D498 24A5D7F4 */  addiu   $a1, $a1, %lo(D_8099D7F4)  ## $a1 = 8099D7F4
/* 007EC 8099D49C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 007F0 8099D4A0 24060189 */  addiu   $a2, $zero, 0x0189         ## $a2 = 00000189
/* 007F4 8099D4A4 AC480000 */  sw      $t0, 0x0000($v0)           ## 00000000
/* 007F8 8099D4A8 0C0346A2 */  jal     Matrix_NewMtx              
/* 007FC 8099D4AC 00409025 */  or      $s2, $v0, $zero            ## $s2 = 00000000
/* 00800 8099D4B0 AE420004 */  sw      $v0, 0x0004($s2)           ## 00000004
/* 00804 8099D4B4 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 00808 8099D4B8 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 0080C 8099D4BC AE1902D0 */  sw      $t9, 0x02D0($s0)           ## 000002D0
/* 00810 8099D4C0 AC550004 */  sw      $s5, 0x0004($v0)           ## 00000004
/* 00814 8099D4C4 AC570000 */  sw      $s7, 0x0000($v0)           ## 00000000
.L8099D4C8:
/* 00818 8099D4C8 26730001 */  addiu   $s3, $s3, 0x0001           ## $s3 = 00000001
/* 0081C 8099D4CC 00139C00 */  sll     $s3, $s3, 16               
/* 00820 8099D4D0 00139C03 */  sra     $s3, $s3, 16               
/* 00824 8099D4D4 2A610064 */  slti    $at, $s3, 0x0064           
/* 00828 8099D4D8 1420FFB0 */  bne     $at, $zero, .L8099D39C     
/* 0082C 8099D4DC 2631003C */  addiu   $s1, $s1, 0x003C           ## $s1 = 0000003C
/* 00830 8099D4E0 3C06809A */  lui     $a2, %hi(D_8099D808)       ## $a2 = 809A0000
/* 00834 8099D4E4 24C6D808 */  addiu   $a2, $a2, %lo(D_8099D808)  ## $a2 = 8099D808
/* 00838 8099D4E8 27A40070 */  addiu   $a0, $sp, 0x0070           ## $a0 = FFFFFFE0
/* 0083C 8099D4EC 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 00840 8099D4F0 0C031AD5 */  jal     Graph_CloseDisps              
/* 00844 8099D4F4 2407018F */  addiu   $a3, $zero, 0x018F         ## $a3 = 0000018F
/* 00848 8099D4F8 8FBF0044 */  lw      $ra, 0x0044($sp)           
/* 0084C 8099D4FC D7B40018 */  ldc1    $f20, 0x0018($sp)          
/* 00850 8099D500 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 00854 8099D504 8FB10024 */  lw      $s1, 0x0024($sp)           
/* 00858 8099D508 8FB20028 */  lw      $s2, 0x0028($sp)           
/* 0085C 8099D50C 8FB3002C */  lw      $s3, 0x002C($sp)           
/* 00860 8099D510 8FB40030 */  lw      $s4, 0x0030($sp)           
/* 00864 8099D514 8FB50034 */  lw      $s5, 0x0034($sp)           
/* 00868 8099D518 8FB60038 */  lw      $s6, 0x0038($sp)           
/* 0086C 8099D51C 8FB7003C */  lw      $s7, 0x003C($sp)           
/* 00870 8099D520 8FBE0040 */  lw      $s8, 0x0040($sp)           
/* 00874 8099D524 03E00008 */  jr      $ra                        
/* 00878 8099D528 27BD0090 */  addiu   $sp, $sp, 0x0090           ## $sp = 00000000


