glabel EnElf_Draw
/* 03C28 80A05858 27BDFF80 */  addiu   $sp, $sp, 0xFF80           ## $sp = FFFFFF80
/* 03C2C 80A0585C AFBF002C */  sw      $ra, 0x002C($sp)
/* 03C30 80A05860 AFB00028 */  sw      $s0, 0x0028($sp)
/* 03C34 80A05864 848E02A8 */  lh      $t6, 0x02A8($a0)           ## 000002A8
/* 03C38 80A05868 24010008 */  addiu   $at, $zero, 0x0008         ## $at = 00000008
/* 03C3C 80A0586C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 03C40 80A05870 00A03025 */  or      $a2, $a1, $zero            ## $a2 = 00000000
/* 03C44 80A05874 11C101A1 */  beq     $t6, $at, .L80A05EFC
/* 03C48 80A05878 8CA21C44 */  lw      $v0, 0x1C44($a1)           ## 00001C44
/* 03C4C 80A0587C 948F02C4 */  lhu     $t7, 0x02C4($a0)           ## 000002C4
/* 03C50 80A05880 31F80008 */  andi    $t8, $t7, 0x0008           ## $t8 = 00000000
/* 03C54 80A05884 5700019E */  bnel    $t8, $zero, .L80A05F00
/* 03C58 80A05888 8FBF002C */  lw      $ra, 0x002C($sp)
/* 03C5C 80A0588C 8C59067C */  lw      $t9, 0x067C($v0)           ## 0000067C
/* 03C60 80A05890 3C0B8016 */  lui     $t3, 0x8016                ## $t3 = 80160000
/* 03C64 80A05894 24050020 */  addiu   $a1, $zero, 0x0020         ## $a1 = 00000020
/* 03C68 80A05898 001952C0 */  sll     $t2, $t9, 11
/* 03C6C 80A0589C 0543000C */  bgezl   $t2, .L80A058D0
/* 03C70 80A058A0 8CC40000 */  lw      $a0, 0x0000($a2)           ## 00000000
/* 03C74 80A058A4 8D6BFA90 */  lw      $t3, -0x0570($t3)          ## 8015FA90
/* 03C78 80A058A8 C48800EC */  lwc1    $f8, 0x00EC($a0)           ## 000000EC
/* 03C7C 80A058AC 856C1508 */  lh      $t4, 0x1508($t3)           ## 80161508
/* 03C80 80A058B0 448C2000 */  mtc1    $t4, $f4                   ## $f4 = 0.00
/* 03C84 80A058B4 00000000 */  nop
/* 03C88 80A058B8 468021A0 */  cvt.s.w $f6, $f4
/* 03C8C 80A058BC 4608303C */  c.lt.s  $f6, $f8
/* 03C90 80A058C0 00000000 */  nop
/* 03C94 80A058C4 4502018E */  bc1fl   .L80A05F00
/* 03C98 80A058C8 8FBF002C */  lw      $ra, 0x002C($sp)
/* 03C9C 80A058CC 8CC40000 */  lw      $a0, 0x0000($a2)           ## 00000000
.L80A058D0:
/* 03CA0 80A058D0 0C031A73 */  jal     Graph_Alloc

/* 03CA4 80A058D4 AFA60084 */  sw      $a2, 0x0084($sp)
/* 03CA8 80A058D8 8FAD0084 */  lw      $t5, 0x0084($sp)
/* 03CAC 80A058DC 3C0680A0 */  lui     $a2, %hi(D_80A06124)       ## $a2 = 80A00000
/* 03CB0 80A058E0 24C66124 */  addiu   $a2, $a2, %lo(D_80A06124)  ## $a2 = 80A06124
/* 03CB4 80A058E4 8DA50000 */  lw      $a1, 0x0000($t5)           ## 00000000
/* 03CB8 80A058E8 AFA20068 */  sw      $v0, 0x0068($sp)
/* 03CBC 80A058EC 27A40050 */  addiu   $a0, $sp, 0x0050           ## $a0 = FFFFFFD0
/* 03CC0 80A058F0 24070AAA */  addiu   $a3, $zero, 0x0AAA         ## $a3 = 00000AAA
/* 03CC4 80A058F4 0C031AB1 */  jal     Graph_OpenDisp
/* 03CC8 80A058F8 AFA50060 */  sw      $a1, 0x0060($sp)
/* 03CCC 80A058FC 8FAE0084 */  lw      $t6, 0x0084($sp)
/* 03CD0 80A05900 0C0252D6 */  jal     func_80094B58
/* 03CD4 80A05904 8DC40000 */  lw      $a0, 0x0000($t6)           ## 00000000
/* 03CD8 80A05908 960402BE */  lhu     $a0, 0x02BE($s0)           ## 000002BE
/* 03CDC 80A0590C 8FA80068 */  lw      $t0, 0x0068($sp)
/* 03CE0 80A05910 8FA90060 */  lw      $t1, 0x0060($sp)
/* 03CE4 80A05914 00800821 */  addu    $at, $a0, $zero
/* 03CE8 80A05918 00042080 */  sll     $a0, $a0,  2
/* 03CEC 80A0591C 00812023 */  subu    $a0, $a0, $at
/* 03CF0 80A05920 000420C0 */  sll     $a0, $a0,  3
/* 03CF4 80A05924 00812021 */  addu    $a0, $a0, $at
/* 03CF8 80A05928 00042040 */  sll     $a0, $a0,  1
/* 03CFC 80A0592C 308401FF */  andi    $a0, $a0, 0x01FF           ## $a0 = 00000000
/* 03D00 80A05930 28810100 */  slti    $at, $a0, 0x0100
/* 03D04 80A05934 14200003 */  bne     $at, $zero, .L80A05944
/* 03D08 80A05938 3C19DB06 */  lui     $t9, 0xDB06                ## $t9 = DB060000
/* 03D0C 80A0593C 240F01FF */  addiu   $t7, $zero, 0x01FF         ## $t7 = 000001FF
/* 03D10 80A05940 01E42023 */  subu    $a0, $t7, $a0
.L80A05944:
/* 03D14 80A05944 860202C2 */  lh      $v0, 0x02C2($s0)           ## 000002C2
/* 03D18 80A05948 37390020 */  ori     $t9, $t9, 0x0020           ## $t9 = DB060020
/* 03D1C 80A0594C 3C0AE700 */  lui     $t2, 0xE700                ## $t2 = E7000000
/* 03D20 80A05950 0441000A */  bgez    $v0, .L80A0597C
/* 03D24 80A05954 3C0BFA00 */  lui     $t3, 0xFA00                ## $t3 = FA000000
/* 03D28 80A05958 44825000 */  mtc1    $v0, $f10                  ## $f10 = 0.00
/* 03D2C 80A0595C 3C0180A0 */  lui     $at, %hi(D_80A0623C)       ## $at = 80A00000
/* 03D30 80A05960 C432623C */  lwc1    $f18, %lo(D_80A0623C)($at)
/* 03D34 80A05964 46805420 */  cvt.s.w $f16, $f10
/* 03D38 80A05968 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 03D3C 80A0596C 44813000 */  mtc1    $at, $f6                   ## $f6 = 1.00
/* 03D40 80A05970 46128102 */  mul.s   $f4, $f16, $f18
/* 03D44 80A05974 10000004 */  beq     $zero, $zero, .L80A05988
/* 03D48 80A05978 46062000 */  add.s   $f0, $f4, $f6
.L80A0597C:
/* 03D4C 80A0597C 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 03D50 80A05980 44810000 */  mtc1    $at, $f0                   ## $f0 = 1.00
/* 03D54 80A05984 00000000 */  nop
.L80A05988:
/* 03D58 80A05988 8D2302D0 */  lw      $v1, 0x02D0($t1)           ## 000002D0
/* 03D5C 80A0598C 01001025 */  or      $v0, $t0, $zero            ## $v0 = 00000000
/* 03D60 80A05990 356B0001 */  ori     $t3, $t3, 0x0001           ## $t3 = FA000001
/* 03D64 80A05994 24780008 */  addiu   $t8, $v1, 0x0008           ## $t8 = 00000008
/* 03D68 80A05998 AD3802D0 */  sw      $t8, 0x02D0($t1)           ## 000002D0
/* 03D6C 80A0599C AC680004 */  sw      $t0, 0x0004($v1)           ## 00000004
/* 03D70 80A059A0 AC790000 */  sw      $t9, 0x0000($v1)           ## 00000000
/* 03D74 80A059A4 AC4A0000 */  sw      $t2, 0x0000($v0)           ## 00000000
/* 03D78 80A059A8 AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 03D7C 80A059AC 25080008 */  addiu   $t0, $t0, 0x0008           ## $t0 = 00000008
/* 03D80 80A059B0 01001025 */  or      $v0, $t0, $zero            ## $v0 = 00000008
/* 03D84 80A059B4 AC4B0000 */  sw      $t3, 0x0000($v0)           ## 00000008
/* 03D88 80A059B8 444CF800 */  cfc1    $t4, $31
/* 03D8C 80A059BC 240D0001 */  addiu   $t5, $zero, 0x0001         ## $t5 = 00000001
/* 03D90 80A059C0 44CDF800 */  ctc1    $t5, $31
/* 03D94 80A059C4 C6080244 */  lwc1    $f8, 0x0244($s0)           ## 00000244
/* 03D98 80A059C8 3C0780A0 */  lui     $a3, %hi(func_80A05734)    ## $a3 = 80A00000
/* 03D9C 80A059CC 25080008 */  addiu   $t0, $t0, 0x0008           ## $t0 = 00000010
/* 03DA0 80A059D0 460042A4 */  cvt.w.s $f10, $f8
/* 03DA4 80A059D4 24E75734 */  addiu   $a3, $a3, %lo(func_80A05734) ## $a3 = 80A05734
/* 03DA8 80A059D8 444DF800 */  cfc1    $t5, $31
/* 03DAC 80A059DC 00000000 */  nop
/* 03DB0 80A059E0 31AD0078 */  andi    $t5, $t5, 0x0078           ## $t5 = 00000000
/* 03DB4 80A059E4 11A00012 */  beq     $t5, $zero, .L80A05A30
/* 03DB8 80A059E8 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 03DBC 80A059EC 44815000 */  mtc1    $at, $f10                  ## $f10 = 2147483648.00
/* 03DC0 80A059F0 240D0001 */  addiu   $t5, $zero, 0x0001         ## $t5 = 00000001
/* 03DC4 80A059F4 460A4281 */  sub.s   $f10, $f8, $f10
/* 03DC8 80A059F8 44CDF800 */  ctc1    $t5, $31
/* 03DCC 80A059FC 00000000 */  nop
/* 03DD0 80A05A00 460052A4 */  cvt.w.s $f10, $f10
/* 03DD4 80A05A04 444DF800 */  cfc1    $t5, $31
/* 03DD8 80A05A08 00000000 */  nop
/* 03DDC 80A05A0C 31AD0078 */  andi    $t5, $t5, 0x0078           ## $t5 = 00000000
/* 03DE0 80A05A10 15A00005 */  bne     $t5, $zero, .L80A05A28
/* 03DE4 80A05A14 00000000 */  nop
/* 03DE8 80A05A18 440D5000 */  mfc1    $t5, $f10
/* 03DEC 80A05A1C 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 03DF0 80A05A20 10000007 */  beq     $zero, $zero, .L80A05A40
/* 03DF4 80A05A24 01A16825 */  or      $t5, $t5, $at              ## $t5 = 80000000
.L80A05A28:
/* 03DF8 80A05A28 10000005 */  beq     $zero, $zero, .L80A05A40
/* 03DFC 80A05A2C 240DFFFF */  addiu   $t5, $zero, 0xFFFF         ## $t5 = FFFFFFFF
.L80A05A30:
/* 03E00 80A05A30 440D5000 */  mfc1    $t5, $f10
/* 03E04 80A05A34 00000000 */  nop
/* 03E08 80A05A38 05A0FFFB */  bltz    $t5, .L80A05A28
/* 03E0C 80A05A3C 00000000 */  nop
.L80A05A40:
/* 03E10 80A05A40 44CCF800 */  ctc1    $t4, $31
/* 03E14 80A05A44 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 03E18 80A05A48 C6100248 */  lwc1    $f16, 0x0248($s0)          ## 00000248
/* 03E1C 80A05A4C 000DC600 */  sll     $t8, $t5, 24
/* 03E20 80A05A50 4459F800 */  cfc1    $t9, $31
/* 03E24 80A05A54 44CAF800 */  ctc1    $t2, $31
/* 03E28 80A05A58 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 03E2C 80A05A5C 460084A4 */  cvt.w.s $f18, $f16
/* 03E30 80A05A60 444AF800 */  cfc1    $t2, $31
/* 03E34 80A05A64 00000000 */  nop
/* 03E38 80A05A68 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 03E3C 80A05A6C 51400013 */  beql    $t2, $zero, .L80A05ABC
/* 03E40 80A05A70 440A9000 */  mfc1    $t2, $f18
/* 03E44 80A05A74 44819000 */  mtc1    $at, $f18                  ## $f18 = 2147483648.00
/* 03E48 80A05A78 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 03E4C 80A05A7C 46128481 */  sub.s   $f18, $f16, $f18
/* 03E50 80A05A80 44CAF800 */  ctc1    $t2, $31
/* 03E54 80A05A84 00000000 */  nop
/* 03E58 80A05A88 460094A4 */  cvt.w.s $f18, $f18
/* 03E5C 80A05A8C 444AF800 */  cfc1    $t2, $31
/* 03E60 80A05A90 00000000 */  nop
/* 03E64 80A05A94 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 03E68 80A05A98 15400005 */  bne     $t2, $zero, .L80A05AB0
/* 03E6C 80A05A9C 00000000 */  nop
/* 03E70 80A05AA0 440A9000 */  mfc1    $t2, $f18
/* 03E74 80A05AA4 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 03E78 80A05AA8 10000007 */  beq     $zero, $zero, .L80A05AC8
/* 03E7C 80A05AAC 01415025 */  or      $t2, $t2, $at              ## $t2 = 80000000
.L80A05AB0:
/* 03E80 80A05AB0 10000005 */  beq     $zero, $zero, .L80A05AC8
/* 03E84 80A05AB4 240AFFFF */  addiu   $t2, $zero, 0xFFFF         ## $t2 = FFFFFFFF
/* 03E88 80A05AB8 440A9000 */  mfc1    $t2, $f18
.L80A05ABC:
/* 03E8C 80A05ABC 00000000 */  nop
/* 03E90 80A05AC0 0540FFFB */  bltz    $t2, .L80A05AB0
/* 03E94 80A05AC4 00000000 */  nop
.L80A05AC8:
/* 03E98 80A05AC8 44D9F800 */  ctc1    $t9, $31
/* 03E9C 80A05ACC 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 03EA0 80A05AD0 C604024C */  lwc1    $f4, 0x024C($s0)           ## 0000024C
/* 03EA4 80A05AD4 314C00FF */  andi    $t4, $t2, 0x00FF           ## $t4 = 000000FF
/* 03EA8 80A05AD8 444FF800 */  cfc1    $t7, $31
/* 03EAC 80A05ADC 44D9F800 */  ctc1    $t9, $31
/* 03EB0 80A05AE0 000C6C00 */  sll     $t5, $t4, 16
/* 03EB4 80A05AE4 030D7025 */  or      $t6, $t8, $t5              ## $t6 = FFFFFFFF
/* 03EB8 80A05AE8 460021A4 */  cvt.w.s $f6, $f4
/* 03EBC 80A05AEC 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 03EC0 80A05AF0 4459F800 */  cfc1    $t9, $31
/* 03EC4 80A05AF4 00000000 */  nop
/* 03EC8 80A05AF8 33390078 */  andi    $t9, $t9, 0x0078           ## $t9 = 00000000
/* 03ECC 80A05AFC 53200013 */  beql    $t9, $zero, .L80A05B4C
/* 03ED0 80A05B00 44193000 */  mfc1    $t9, $f6
/* 03ED4 80A05B04 44813000 */  mtc1    $at, $f6                   ## $f6 = 2147483648.00
/* 03ED8 80A05B08 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 03EDC 80A05B0C 46062181 */  sub.s   $f6, $f4, $f6
/* 03EE0 80A05B10 44D9F800 */  ctc1    $t9, $31
/* 03EE4 80A05B14 00000000 */  nop
/* 03EE8 80A05B18 460031A4 */  cvt.w.s $f6, $f6
/* 03EEC 80A05B1C 4459F800 */  cfc1    $t9, $31
/* 03EF0 80A05B20 00000000 */  nop
/* 03EF4 80A05B24 33390078 */  andi    $t9, $t9, 0x0078           ## $t9 = 00000000
/* 03EF8 80A05B28 17200005 */  bne     $t9, $zero, .L80A05B40
/* 03EFC 80A05B2C 00000000 */  nop
/* 03F00 80A05B30 44193000 */  mfc1    $t9, $f6
/* 03F04 80A05B34 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 03F08 80A05B38 10000007 */  beq     $zero, $zero, .L80A05B58
/* 03F0C 80A05B3C 0321C825 */  or      $t9, $t9, $at              ## $t9 = 80000000
.L80A05B40:
/* 03F10 80A05B40 10000005 */  beq     $zero, $zero, .L80A05B58
/* 03F14 80A05B44 2419FFFF */  addiu   $t9, $zero, 0xFFFF         ## $t9 = FFFFFFFF
/* 03F18 80A05B48 44193000 */  mfc1    $t9, $f6
.L80A05B4C:
/* 03F1C 80A05B4C 00000000 */  nop
/* 03F20 80A05B50 0720FFFB */  bltz    $t9, .L80A05B40
/* 03F24 80A05B54 00000000 */  nop
.L80A05B58:
/* 03F28 80A05B58 44CFF800 */  ctc1    $t7, $31
/* 03F2C 80A05B5C C6080250 */  lwc1    $f8, 0x0250($s0)           ## 00000250
/* 03F30 80A05B60 240F0001 */  addiu   $t7, $zero, 0x0001         ## $t7 = 00000001
/* 03F34 80A05B64 332B00FF */  andi    $t3, $t9, 0x00FF           ## $t3 = 000000FF
/* 03F38 80A05B68 46004282 */  mul.s   $f10, $f8, $f0
/* 03F3C 80A05B6C 000B6200 */  sll     $t4, $t3,  8
/* 03F40 80A05B70 3C19E200 */  lui     $t9, 0xE200                ## $t9 = E2000000
/* 03F44 80A05B74 01CCC025 */  or      $t8, $t6, $t4              ## $t8 = FFFFFFFF
/* 03F48 80A05B78 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 03F4C 80A05B7C 3739001C */  ori     $t9, $t9, 0x001C           ## $t9 = E200001C
/* 03F50 80A05B80 444DF800 */  cfc1    $t5, $31
/* 03F54 80A05B84 44CFF800 */  ctc1    $t7, $31
/* 03F58 80A05B88 00000000 */  nop
/* 03F5C 80A05B8C 46005424 */  cvt.w.s $f16, $f10
/* 03F60 80A05B90 444FF800 */  cfc1    $t7, $31
/* 03F64 80A05B94 00000000 */  nop
/* 03F68 80A05B98 31EF0078 */  andi    $t7, $t7, 0x0078           ## $t7 = 00000000
/* 03F6C 80A05B9C 51E00013 */  beql    $t7, $zero, .L80A05BEC
/* 03F70 80A05BA0 440F8000 */  mfc1    $t7, $f16
/* 03F74 80A05BA4 44818000 */  mtc1    $at, $f16                  ## $f16 = 2147483648.00
/* 03F78 80A05BA8 240F0001 */  addiu   $t7, $zero, 0x0001         ## $t7 = 00000001
/* 03F7C 80A05BAC 46105401 */  sub.s   $f16, $f10, $f16
/* 03F80 80A05BB0 44CFF800 */  ctc1    $t7, $31
/* 03F84 80A05BB4 00000000 */  nop
/* 03F88 80A05BB8 46008424 */  cvt.w.s $f16, $f16
/* 03F8C 80A05BBC 444FF800 */  cfc1    $t7, $31
/* 03F90 80A05BC0 00000000 */  nop
/* 03F94 80A05BC4 31EF0078 */  andi    $t7, $t7, 0x0078           ## $t7 = 00000000
/* 03F98 80A05BC8 15E00005 */  bne     $t7, $zero, .L80A05BE0
/* 03F9C 80A05BCC 00000000 */  nop
/* 03FA0 80A05BD0 440F8000 */  mfc1    $t7, $f16
/* 03FA4 80A05BD4 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 03FA8 80A05BD8 10000007 */  beq     $zero, $zero, .L80A05BF8
/* 03FAC 80A05BDC 01E17825 */  or      $t7, $t7, $at              ## $t7 = 80000000
.L80A05BE0:
/* 03FB0 80A05BE0 10000005 */  beq     $zero, $zero, .L80A05BF8
/* 03FB4 80A05BE4 240FFFFF */  addiu   $t7, $zero, 0xFFFF         ## $t7 = FFFFFFFF
/* 03FB8 80A05BE8 440F8000 */  mfc1    $t7, $f16
.L80A05BEC:
/* 03FBC 80A05BEC 00000000 */  nop
/* 03FC0 80A05BF0 05E0FFFB */  bltz    $t7, .L80A05BE0
/* 03FC4 80A05BF4 00000000 */  nop
.L80A05BF8:
/* 03FC8 80A05BF8 31EA00FF */  andi    $t2, $t7, 0x00FF           ## $t2 = 000000FF
/* 03FCC 80A05BFC 030A5825 */  or      $t3, $t8, $t2              ## $t3 = FFFFFFFF
/* 03FD0 80A05C00 AC4B0004 */  sw      $t3, 0x0004($v0)           ## 0000000C
/* 03FD4 80A05C04 960E02C4 */  lhu     $t6, 0x02C4($s0)           ## 000002C4
/* 03FD8 80A05C08 44CDF800 */  ctc1    $t5, $31
/* 03FDC 80A05C0C 3C0ADF00 */  lui     $t2, 0xDF00                ## $t2 = DF000000
/* 03FE0 80A05C10 31CC0004 */  andi    $t4, $t6, 0x0004           ## $t4 = 00000004
/* 03FE4 80A05C14 1180000A */  beq     $t4, $zero, .L80A05C40
/* 03FE8 80A05C18 3C0EFB00 */  lui     $t6, 0xFB00                ## $t6 = FB000000
/* 03FEC 80A05C1C 01001025 */  or      $v0, $t0, $zero            ## $v0 = 00000010
/* 03FF0 80A05C20 3C0DE200 */  lui     $t5, 0xE200                ## $t5 = E2000000
/* 03FF4 80A05C24 3C0F0C18 */  lui     $t7, 0x0C18                ## $t7 = 0C180000
/* 03FF8 80A05C28 35EF4340 */  ori     $t7, $t7, 0x4340           ## $t7 = 0C184340
/* 03FFC 80A05C2C 35AD001C */  ori     $t5, $t5, 0x001C           ## $t5 = E200001C
/* 04000 80A05C30 AC4D0000 */  sw      $t5, 0x0000($v0)           ## 00000010
/* 04004 80A05C34 AC4F0004 */  sw      $t7, 0x0004($v0)           ## 00000014
/* 04008 80A05C38 10000007 */  beq     $zero, $zero, .L80A05C58
/* 0400C 80A05C3C 25080008 */  addiu   $t0, $t0, 0x0008           ## $t0 = 00000018
.L80A05C40:
/* 04010 80A05C40 01001025 */  or      $v0, $t0, $zero            ## $v0 = 00000018
/* 04014 80A05C44 3C180C18 */  lui     $t8, 0x0C18                ## $t8 = 0C180000
/* 04018 80A05C48 37184B50 */  ori     $t8, $t8, 0x4B50           ## $t8 = 0C184B50
/* 0401C 80A05C4C AC580004 */  sw      $t8, 0x0004($v0)           ## 0000001C
/* 04020 80A05C50 AC590000 */  sw      $t9, 0x0000($v0)           ## 00000018
/* 04024 80A05C54 25080008 */  addiu   $t0, $t0, 0x0008           ## $t0 = 00000020
.L80A05C58:
/* 04028 80A05C58 AD0A0000 */  sw      $t2, 0x0000($t0)           ## 00000020
/* 0402C 80A05C5C AD000004 */  sw      $zero, 0x0004($t0)         ## 00000024
/* 04030 80A05C60 8D2202D0 */  lw      $v0, 0x02D0($t1)           ## 000002D0
/* 04034 80A05C64 444CF800 */  cfc1    $t4, $31
/* 04038 80A05C68 240D0001 */  addiu   $t5, $zero, 0x0001         ## $t5 = 00000001
/* 0403C 80A05C6C 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000020
/* 04040 80A05C70 AD2B02D0 */  sw      $t3, 0x02D0($t1)           ## 000002D0
/* 04044 80A05C74 AC4E0000 */  sw      $t6, 0x0000($v0)           ## 00000018
/* 04048 80A05C78 44CDF800 */  ctc1    $t5, $31
/* 0404C 80A05C7C C6120254 */  lwc1    $f18, 0x0254($s0)          ## 00000254
/* 04050 80A05C80 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 04054 80A05C84 46009124 */  cvt.w.s $f4, $f18
/* 04058 80A05C88 444DF800 */  cfc1    $t5, $31
/* 0405C 80A05C8C 00000000 */  nop
/* 04060 80A05C90 31AD0078 */  andi    $t5, $t5, 0x0078           ## $t5 = 00000000
/* 04064 80A05C94 51A00013 */  beql    $t5, $zero, .L80A05CE4
/* 04068 80A05C98 440D2000 */  mfc1    $t5, $f4
/* 0406C 80A05C9C 44812000 */  mtc1    $at, $f4                   ## $f4 = 2147483648.00
/* 04070 80A05CA0 240D0001 */  addiu   $t5, $zero, 0x0001         ## $t5 = 00000001
/* 04074 80A05CA4 46049101 */  sub.s   $f4, $f18, $f4
/* 04078 80A05CA8 44CDF800 */  ctc1    $t5, $31
/* 0407C 80A05CAC 00000000 */  nop
/* 04080 80A05CB0 46002124 */  cvt.w.s $f4, $f4
/* 04084 80A05CB4 444DF800 */  cfc1    $t5, $31
/* 04088 80A05CB8 00000000 */  nop
/* 0408C 80A05CBC 31AD0078 */  andi    $t5, $t5, 0x0078           ## $t5 = 00000000
/* 04090 80A05CC0 15A00005 */  bne     $t5, $zero, .L80A05CD8
/* 04094 80A05CC4 00000000 */  nop
/* 04098 80A05CC8 440D2000 */  mfc1    $t5, $f4
/* 0409C 80A05CCC 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 040A0 80A05CD0 10000007 */  beq     $zero, $zero, .L80A05CF0
/* 040A4 80A05CD4 01A16825 */  or      $t5, $t5, $at              ## $t5 = 80000000
.L80A05CD8:
/* 040A8 80A05CD8 10000005 */  beq     $zero, $zero, .L80A05CF0
/* 040AC 80A05CDC 240DFFFF */  addiu   $t5, $zero, 0xFFFF         ## $t5 = FFFFFFFF
/* 040B0 80A05CE0 440D2000 */  mfc1    $t5, $f4
.L80A05CE4:
/* 040B4 80A05CE4 00000000 */  nop
/* 040B8 80A05CE8 05A0FFFB */  bltz    $t5, .L80A05CD8
/* 040BC 80A05CEC 00000000 */  nop
.L80A05CF0:
/* 040C0 80A05CF0 44CCF800 */  ctc1    $t4, $31
/* 040C4 80A05CF4 240B0001 */  addiu   $t3, $zero, 0x0001         ## $t3 = 00000001
/* 040C8 80A05CF8 C6060258 */  lwc1    $f6, 0x0258($s0)           ## 00000258
/* 040CC 80A05CFC 000DC600 */  sll     $t8, $t5, 24
/* 040D0 80A05D00 444AF800 */  cfc1    $t2, $31
/* 040D4 80A05D04 44CBF800 */  ctc1    $t3, $31
/* 040D8 80A05D08 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 040DC 80A05D0C 46003224 */  cvt.w.s $f8, $f6
/* 040E0 80A05D10 444BF800 */  cfc1    $t3, $31
/* 040E4 80A05D14 00000000 */  nop
/* 040E8 80A05D18 316B0078 */  andi    $t3, $t3, 0x0078           ## $t3 = 00000000
/* 040EC 80A05D1C 51600013 */  beql    $t3, $zero, .L80A05D6C
/* 040F0 80A05D20 440B4000 */  mfc1    $t3, $f8
/* 040F4 80A05D24 44814000 */  mtc1    $at, $f8                   ## $f8 = 2147483648.00
/* 040F8 80A05D28 240B0001 */  addiu   $t3, $zero, 0x0001         ## $t3 = 00000001
/* 040FC 80A05D2C 46083201 */  sub.s   $f8, $f6, $f8
/* 04100 80A05D30 44CBF800 */  ctc1    $t3, $31
/* 04104 80A05D34 00000000 */  nop
/* 04108 80A05D38 46004224 */  cvt.w.s $f8, $f8
/* 0410C 80A05D3C 444BF800 */  cfc1    $t3, $31
/* 04110 80A05D40 00000000 */  nop
/* 04114 80A05D44 316B0078 */  andi    $t3, $t3, 0x0078           ## $t3 = 00000000
/* 04118 80A05D48 15600005 */  bne     $t3, $zero, .L80A05D60
/* 0411C 80A05D4C 00000000 */  nop
/* 04120 80A05D50 440B4000 */  mfc1    $t3, $f8
/* 04124 80A05D54 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 04128 80A05D58 10000007 */  beq     $zero, $zero, .L80A05D78
/* 0412C 80A05D5C 01615825 */  or      $t3, $t3, $at              ## $t3 = 80000000
.L80A05D60:
/* 04130 80A05D60 10000005 */  beq     $zero, $zero, .L80A05D78
/* 04134 80A05D64 240BFFFF */  addiu   $t3, $zero, 0xFFFF         ## $t3 = FFFFFFFF
/* 04138 80A05D68 440B4000 */  mfc1    $t3, $f8
.L80A05D6C:
/* 0413C 80A05D6C 00000000 */  nop
/* 04140 80A05D70 0560FFFB */  bltz    $t3, .L80A05D60
/* 04144 80A05D74 00000000 */  nop
.L80A05D78:
/* 04148 80A05D78 44CAF800 */  ctc1    $t2, $31
/* 0414C 80A05D7C 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 04150 80A05D80 C60A025C */  lwc1    $f10, 0x025C($s0)          ## 0000025C
/* 04154 80A05D84 316C00FF */  andi    $t4, $t3, 0x00FF           ## $t4 = 000000FF
/* 04158 80A05D88 4459F800 */  cfc1    $t9, $31
/* 0415C 80A05D8C 44CAF800 */  ctc1    $t2, $31
/* 04160 80A05D90 000C6C00 */  sll     $t5, $t4, 16
/* 04164 80A05D94 030D7825 */  or      $t7, $t8, $t5              ## $t7 = FFFFFFFF
/* 04168 80A05D98 46005424 */  cvt.w.s $f16, $f10
/* 0416C 80A05D9C 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 04170 80A05DA0 444AF800 */  cfc1    $t2, $31
/* 04174 80A05DA4 00000000 */  nop
/* 04178 80A05DA8 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 0417C 80A05DAC 51400013 */  beql    $t2, $zero, .L80A05DFC
/* 04180 80A05DB0 440A8000 */  mfc1    $t2, $f16
/* 04184 80A05DB4 44818000 */  mtc1    $at, $f16                  ## $f16 = 2147483648.00
/* 04188 80A05DB8 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 0418C 80A05DBC 46105401 */  sub.s   $f16, $f10, $f16
/* 04190 80A05DC0 44CAF800 */  ctc1    $t2, $31
/* 04194 80A05DC4 00000000 */  nop
/* 04198 80A05DC8 46008424 */  cvt.w.s $f16, $f16
/* 0419C 80A05DCC 444AF800 */  cfc1    $t2, $31
/* 041A0 80A05DD0 00000000 */  nop
/* 041A4 80A05DD4 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 041A8 80A05DD8 15400005 */  bne     $t2, $zero, .L80A05DF0
/* 041AC 80A05DDC 00000000 */  nop
/* 041B0 80A05DE0 440A8000 */  mfc1    $t2, $f16
/* 041B4 80A05DE4 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 041B8 80A05DE8 10000007 */  beq     $zero, $zero, .L80A05E08
/* 041BC 80A05DEC 01415025 */  or      $t2, $t2, $at              ## $t2 = 80000000
.L80A05DF0:
/* 041C0 80A05DF0 10000005 */  beq     $zero, $zero, .L80A05E08
/* 041C4 80A05DF4 240AFFFF */  addiu   $t2, $zero, 0xFFFF         ## $t2 = FFFFFFFF
/* 041C8 80A05DF8 440A8000 */  mfc1    $t2, $f16
.L80A05DFC:
/* 041CC 80A05DFC 00000000 */  nop
/* 041D0 80A05E00 0540FFFB */  bltz    $t2, .L80A05DF0
/* 041D4 80A05E04 00000000 */  nop
.L80A05E08:
/* 041D8 80A05E08 44D9F800 */  ctc1    $t9, $31
/* 041DC 80A05E0C 44849000 */  mtc1    $a0, $f18                  ## $f18 = 0.00
/* 041E0 80A05E10 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 041E4 80A05E14 314E00FF */  andi    $t6, $t2, 0x00FF           ## $t6 = 000000FF
/* 041E8 80A05E18 46809120 */  cvt.s.w $f4, $f18
/* 041EC 80A05E1C 000E6200 */  sll     $t4, $t6,  8
/* 041F0 80A05E20 01ECC025 */  or      $t8, $t7, $t4              ## $t8 = FFFFFFFF
/* 041F4 80A05E24 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 041F8 80A05E28 46002182 */  mul.s   $f6, $f4, $f0
/* 041FC 80A05E2C 444DF800 */  cfc1    $t5, $31
/* 04200 80A05E30 44D9F800 */  ctc1    $t9, $31
/* 04204 80A05E34 00000000 */  nop
/* 04208 80A05E38 46003224 */  cvt.w.s $f8, $f6
/* 0420C 80A05E3C 4459F800 */  cfc1    $t9, $31
/* 04210 80A05E40 00000000 */  nop
/* 04214 80A05E44 33390078 */  andi    $t9, $t9, 0x0078           ## $t9 = 00000000
/* 04218 80A05E48 53200013 */  beql    $t9, $zero, .L80A05E98
/* 0421C 80A05E4C 44194000 */  mfc1    $t9, $f8
/* 04220 80A05E50 44814000 */  mtc1    $at, $f8                   ## $f8 = 2147483648.00
/* 04224 80A05E54 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 04228 80A05E58 46083201 */  sub.s   $f8, $f6, $f8
/* 0422C 80A05E5C 44D9F800 */  ctc1    $t9, $31
/* 04230 80A05E60 00000000 */  nop
/* 04234 80A05E64 46004224 */  cvt.w.s $f8, $f8
/* 04238 80A05E68 4459F800 */  cfc1    $t9, $31
/* 0423C 80A05E6C 00000000 */  nop
/* 04240 80A05E70 33390078 */  andi    $t9, $t9, 0x0078           ## $t9 = 00000000
/* 04244 80A05E74 17200005 */  bne     $t9, $zero, .L80A05E8C
/* 04248 80A05E78 00000000 */  nop
/* 0424C 80A05E7C 44194000 */  mfc1    $t9, $f8
/* 04250 80A05E80 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 04254 80A05E84 10000007 */  beq     $zero, $zero, .L80A05EA4
/* 04258 80A05E88 0321C825 */  or      $t9, $t9, $at              ## $t9 = 80000000
.L80A05E8C:
/* 0425C 80A05E8C 10000005 */  beq     $zero, $zero, .L80A05EA4
/* 04260 80A05E90 2419FFFF */  addiu   $t9, $zero, 0xFFFF         ## $t9 = FFFFFFFF
/* 04264 80A05E94 44194000 */  mfc1    $t9, $f8
.L80A05E98:
/* 04268 80A05E98 00000000 */  nop
/* 0426C 80A05E9C 0720FFFB */  bltz    $t9, .L80A05E8C
/* 04270 80A05EA0 00000000 */  nop
.L80A05EA4:
/* 04274 80A05EA4 332B00FF */  andi    $t3, $t9, 0x00FF           ## $t3 = 000000FF
/* 04278 80A05EA8 030B7025 */  or      $t6, $t8, $t3              ## $t6 = FFFFFFFF
/* 0427C 80A05EAC AC4E0004 */  sw      $t6, 0x0004($v0)           ## 0000001C
/* 04280 80A05EB0 8E06016C */  lw      $a2, 0x016C($s0)           ## 0000016C
/* 04284 80A05EB4 8E050150 */  lw      $a1, 0x0150($s0)           ## 00000150
/* 04288 80A05EB8 AFB00014 */  sw      $s0, 0x0014($sp)
/* 0428C 80A05EBC AFA00010 */  sw      $zero, 0x0010($sp)
/* 04290 80A05EC0 8D2F02D0 */  lw      $t7, 0x02D0($t1)           ## 000002D0
/* 04294 80A05EC4 44CDF800 */  ctc1    $t5, $31
/* 04298 80A05EC8 AFA90060 */  sw      $t1, 0x0060($sp)
/* 0429C 80A05ECC 8FA40084 */  lw      $a0, 0x0084($sp)
/* 042A0 80A05ED0 0C0288A2 */  jal     SkelAnime_Draw2
/* 042A4 80A05ED4 AFAF0018 */  sw      $t7, 0x0018($sp)
/* 042A8 80A05ED8 8FA90060 */  lw      $t1, 0x0060($sp)
/* 042AC 80A05EDC 3C0680A0 */  lui     $a2, %hi(D_80A06134)       ## $a2 = 80A00000
/* 042B0 80A05EE0 24C66134 */  addiu   $a2, $a2, %lo(D_80A06134)  ## $a2 = 80A06134
/* 042B4 80A05EE4 AD2202D0 */  sw      $v0, 0x02D0($t1)           ## 000002D0
/* 042B8 80A05EE8 8FAC0084 */  lw      $t4, 0x0084($sp)
/* 042BC 80A05EEC 27A40050 */  addiu   $a0, $sp, 0x0050           ## $a0 = FFFFFFD0
/* 042C0 80A05EF0 24070AE9 */  addiu   $a3, $zero, 0x0AE9         ## $a3 = 00000AE9
/* 042C4 80A05EF4 0C031AD5 */  jal     Graph_CloseDisp
/* 042C8 80A05EF8 8D850000 */  lw      $a1, 0x0000($t4)           ## 00000000
.L80A05EFC:
/* 042CC 80A05EFC 8FBF002C */  lw      $ra, 0x002C($sp)
.L80A05F00:
/* 042D0 80A05F00 8FB00028 */  lw      $s0, 0x0028($sp)
/* 042D4 80A05F04 27BD0080 */  addiu   $sp, $sp, 0x0080           ## $sp = 00000000
/* 042D8 80A05F08 03E00008 */  jr      $ra
/* 042DC 80A05F0C 00000000 */  nop


