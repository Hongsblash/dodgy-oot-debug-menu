glabel ArrowLight_Draw
/* 00580 8086A230 27BDFF68 */  addiu   $sp, $sp, 0xFF68           ## $sp = FFFFFF68
/* 00584 8086A234 AFBF003C */  sw      $ra, 0x003C($sp)
/* 00588 8086A238 AFB00038 */  sw      $s0, 0x0038($sp)
/* 0058C 8086A23C AFA40098 */  sw      $a0, 0x0098($sp)
/* 00590 8086A240 AFA5009C */  sw      $a1, 0x009C($sp)
/* 00594 8086A244 8CAF009C */  lw      $t7, 0x009C($a1)           ## 0000009C
/* 00598 8086A248 AFAF008C */  sw      $t7, 0x008C($sp)
/* 0059C 8086A24C 8C830118 */  lw      $v1, 0x0118($a0)           ## 00000118
/* 005A0 8086A250 50600112 */  beql    $v1, $zero, .L8086A69C
/* 005A4 8086A254 8FBF003C */  lw      $ra, 0x003C($sp)
/* 005A8 8086A258 8C790130 */  lw      $t9, 0x0130($v1)           ## 00000130
/* 005AC 8086A25C 5320010F */  beql    $t9, $zero, .L8086A69C
/* 005B0 8086A260 8FBF003C */  lw      $ra, 0x003C($sp)
/* 005B4 8086A264 948C014E */  lhu     $t4, 0x014E($a0)           ## 0000014E
/* 005B8 8086A268 3C068087 */  lui     $a2, %hi(D_8086BB2C)       ## $a2 = 80870000
/* 005BC 8086A26C 24C6BB2C */  addiu   $a2, $a2, %lo(D_8086BB2C)  ## $a2 = 8086BB2C
/* 005C0 8086A270 298100FF */  slti    $at, $t4, 0x00FF
/* 005C4 8086A274 10200108 */  beq     $at, $zero, .L8086A698
/* 005C8 8086A278 8FAF009C */  lw      $t7, 0x009C($sp)
/* 005CC 8086A27C 906D0249 */  lbu     $t5, 0x0249($v1)           ## 00000249
/* 005D0 8086A280 24070256 */  addiu   $a3, $zero, 0x0256         ## $a3 = 00000256
/* 005D4 8086A284 31AE0002 */  andi    $t6, $t5, 0x0002           ## $t6 = 00000000
/* 005D8 8086A288 51C00004 */  beql    $t6, $zero, .L8086A29C
/* 005DC 8086A28C 00601025 */  or      $v0, $v1, $zero            ## $v0 = 00000000
/* 005E0 8086A290 10000002 */  beq     $zero, $zero, .L8086A29C
/* 005E4 8086A294 00801025 */  or      $v0, $a0, $zero            ## $v0 = 00000000
/* 005E8 8086A298 00601025 */  or      $v0, $v1, $zero            ## $v0 = 00000000
.L8086A29C:
/* 005EC 8086A29C 8DE50000 */  lw      $a1, 0x0000($t7)           ## 00000000
/* 005F0 8086A2A0 AFA20084 */  sw      $v0, 0x0084($sp)
/* 005F4 8086A2A4 27A40070 */  addiu   $a0, $sp, 0x0070           ## $a0 = FFFFFFD8
/* 005F8 8086A2A8 0C031AB1 */  jal     Graph_OpenDisps
/* 005FC 8086A2AC 00A08025 */  or      $s0, $a1, $zero            ## $s0 = 00000000
/* 00600 8086A2B0 8FA20084 */  lw      $v0, 0x0084($sp)
/* 00604 8086A2B4 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 00608 8086A2B8 C44C0024 */  lwc1    $f12, 0x0024($v0)          ## 00000024
/* 0060C 8086A2BC C44E0028 */  lwc1    $f14, 0x0028($v0)          ## 00000028
/* 00610 8086A2C0 0C034261 */  jal     Matrix_Translate
/* 00614 8086A2C4 8C46002C */  lw      $a2, 0x002C($v0)           ## 0000002C
/* 00618 8086A2C8 8FA20084 */  lw      $v0, 0x0084($sp)
/* 0061C 8086A2CC 3C018087 */  lui     $at, %hi(D_8086BB7C)       ## $at = 80870000
/* 00620 8086A2D0 C428BB7C */  lwc1    $f8, %lo(D_8086BB7C)($at)
/* 00624 8086A2D4 845800B6 */  lh      $t8, 0x00B6($v0)           ## 000000B6
/* 00628 8086A2D8 24050001 */  addiu   $a1, $zero, 0x0001         ## $a1 = 00000001
/* 0062C 8086A2DC 44982000 */  mtc1    $t8, $f4                   ## $f4 = 0.00
/* 00630 8086A2E0 00000000 */  nop
/* 00634 8086A2E4 468021A0 */  cvt.s.w $f6, $f4
/* 00638 8086A2E8 46083302 */  mul.s   $f12, $f6, $f8
/* 0063C 8086A2EC 0C034348 */  jal     Matrix_RotateY
/* 00640 8086A2F0 00000000 */  nop
/* 00644 8086A2F4 8FA20084 */  lw      $v0, 0x0084($sp)
/* 00648 8086A2F8 3C018087 */  lui     $at, %hi(D_8086BB80)       ## $at = 80870000
/* 0064C 8086A2FC C432BB80 */  lwc1    $f18, %lo(D_8086BB80)($at)
/* 00650 8086A300 845900B4 */  lh      $t9, 0x00B4($v0)           ## 000000B4
/* 00654 8086A304 24050001 */  addiu   $a1, $zero, 0x0001         ## $a1 = 00000001
/* 00658 8086A308 44995000 */  mtc1    $t9, $f10                  ## $f10 = 0.00
/* 0065C 8086A30C 00000000 */  nop
/* 00660 8086A310 46805420 */  cvt.s.w $f16, $f10
/* 00664 8086A314 46128302 */  mul.s   $f12, $f16, $f18
/* 00668 8086A318 0C0342DC */  jal     Matrix_RotateX
/* 0066C 8086A31C 00000000 */  nop
/* 00670 8086A320 8FA20084 */  lw      $v0, 0x0084($sp)
/* 00674 8086A324 3C018087 */  lui     $at, %hi(D_8086BB84)       ## $at = 80870000
/* 00678 8086A328 C428BB84 */  lwc1    $f8, %lo(D_8086BB84)($at)
/* 0067C 8086A32C 844B00B8 */  lh      $t3, 0x00B8($v0)           ## 000000B8
/* 00680 8086A330 24050001 */  addiu   $a1, $zero, 0x0001         ## $a1 = 00000001
/* 00684 8086A334 448B2000 */  mtc1    $t3, $f4                   ## $f4 = 0.00
/* 00688 8086A338 00000000 */  nop
/* 0068C 8086A33C 468021A0 */  cvt.s.w $f6, $f4
/* 00690 8086A340 46083302 */  mul.s   $f12, $f6, $f8
/* 00694 8086A344 0C0343B5 */  jal     Matrix_RotateZ
/* 00698 8086A348 00000000 */  nop
/* 0069C 8086A34C 3C018087 */  lui     $at, %hi(D_8086BB88)       ## $at = 80870000
/* 006A0 8086A350 C42CBB88 */  lwc1    $f12, %lo(D_8086BB88)($at)
/* 006A4 8086A354 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 006A8 8086A358 44066000 */  mfc1    $a2, $f12
/* 006AC 8086A35C 0C0342A3 */  jal     Matrix_Scale
/* 006B0 8086A360 46006386 */  mov.s   $f14, $f12
/* 006B4 8086A364 8FAC0098 */  lw      $t4, 0x0098($sp)
/* 006B8 8086A368 44805000 */  mtc1    $zero, $f10                ## $f10 = 0.00
/* 006BC 8086A36C C5900164 */  lwc1    $f16, 0x0164($t4)          ## 00000164
/* 006C0 8086A370 4610503C */  c.lt.s  $f10, $f16
/* 006C4 8086A374 00000000 */  nop
/* 006C8 8086A378 45020039 */  bc1fl   .L8086A460
/* 006CC 8086A37C 8FAE009C */  lw      $t6, 0x009C($sp)
/* 006D0 8086A380 0C024DF0 */  jal     func_800937C0
/* 006D4 8086A384 8E0402D0 */  lw      $a0, 0x02D0($s0)           ## 000002D0
/* 006D8 8086A388 AE0202D0 */  sw      $v0, 0x02D0($s0)           ## 000002D0
/* 006DC 8086A38C 244D0008 */  addiu   $t5, $v0, 0x0008           ## $t5 = 00000008
/* 006E0 8086A390 AE0D02D0 */  sw      $t5, 0x02D0($s0)           ## 000002D0
/* 006E4 8086A394 3C0EFA00 */  lui     $t6, 0xFA00                ## $t6 = FA000000
/* 006E8 8086A398 AC4E0000 */  sw      $t6, 0x0000($v0)           ## 00000000
/* 006EC 8086A39C 8FAF0098 */  lw      $t7, 0x0098($sp)
/* 006F0 8086A3A0 3C0141F0 */  lui     $at, 0x41F0                ## $at = 41F00000
/* 006F4 8086A3A4 44819000 */  mtc1    $at, $f18                  ## $f18 = 30.00
/* 006F8 8086A3A8 C5E00164 */  lwc1    $f0, 0x0164($t7)           ## 00000164
/* 006FC 8086A3AC 3C014220 */  lui     $at, 0x4220                ## $at = 42200000
/* 00700 8086A3B0 44814000 */  mtc1    $at, $f8                   ## $f8 = 40.00
/* 00704 8086A3B4 46120102 */  mul.s   $f4, $f0, $f18
/* 00708 8086A3B8 3C014316 */  lui     $at, 0x4316                ## $at = 43160000
/* 0070C 8086A3BC 44819000 */  mtc1    $at, $f18                  ## $f18 = 150.00
/* 00710 8086A3C0 46004282 */  mul.s   $f10, $f8, $f0
/* 00714 8086A3C4 3C0EE300 */  lui     $t6, 0xE300                ## $t6 = E3000000
/* 00718 8086A3C8 35CE1A01 */  ori     $t6, $t6, 0x1A01           ## $t6 = E3001A01
/* 0071C 8086A3CC 240F0030 */  addiu   $t7, $zero, 0x0030         ## $t7 = 00000030
/* 00720 8086A3D0 4600218D */  trunc.w.s $f6, $f4
/* 00724 8086A3D4 46009102 */  mul.s   $f4, $f18, $f0
/* 00728 8086A3D8 440C3000 */  mfc1    $t4, $f6
/* 0072C 8086A3DC 4600540D */  trunc.w.s $f16, $f10
/* 00730 8086A3E0 000C6E00 */  sll     $t5, $t4, 24
/* 00734 8086A3E4 4600218D */  trunc.w.s $f6, $f4
/* 00738 8086A3E8 44188000 */  mfc1    $t8, $f16
/* 0073C 8086A3EC 00000000 */  nop
/* 00740 8086A3F0 331900FF */  andi    $t9, $t8, 0x00FF           ## $t9 = 00000000
/* 00744 8086A3F4 44183000 */  mfc1    $t8, $f6
/* 00748 8086A3F8 00195C00 */  sll     $t3, $t9, 16
/* 0074C 8086A3FC 01AB6025 */  or      $t4, $t5, $t3              ## $t4 = 00000008
/* 00750 8086A400 331900FF */  andi    $t9, $t8, 0x00FF           ## $t9 = 00000000
/* 00754 8086A404 01996825 */  or      $t5, $t4, $t9              ## $t5 = 00000008
/* 00758 8086A408 AC4D0004 */  sw      $t5, 0x0004($v0)           ## 00000004
/* 0075C 8086A40C 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 00760 8086A410 3C0CE300 */  lui     $t4, 0xE300                ## $t4 = E3000000
/* 00764 8086A414 358C1801 */  ori     $t4, $t4, 0x1801           ## $t4 = E3001801
/* 00768 8086A418 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 0076C 8086A41C AE0B02D0 */  sw      $t3, 0x02D0($s0)           ## 000002D0
/* 00770 8086A420 AC4F0004 */  sw      $t7, 0x0004($v0)           ## 00000004
/* 00774 8086A424 AC4E0000 */  sw      $t6, 0x0000($v0)           ## 00000000
/* 00778 8086A428 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 0077C 8086A42C 241900C0 */  addiu   $t9, $zero, 0x00C0         ## $t9 = 000000C0
/* 00780 8086A430 3C0BF64F */  lui     $t3, 0xF64F                ## $t3 = F64F0000
/* 00784 8086A434 24580008 */  addiu   $t8, $v0, 0x0008           ## $t8 = 00000008
/* 00788 8086A438 AE1802D0 */  sw      $t8, 0x02D0($s0)           ## 000002D0
/* 0078C 8086A43C AC590004 */  sw      $t9, 0x0004($v0)           ## 00000004
/* 00790 8086A440 AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
/* 00794 8086A444 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 00798 8086A448 356BC3BC */  ori     $t3, $t3, 0xC3BC           ## $t3 = F64FC3BC
/* 0079C 8086A44C 244D0008 */  addiu   $t5, $v0, 0x0008           ## $t5 = 00000008
/* 007A0 8086A450 AE0D02D0 */  sw      $t5, 0x02D0($s0)           ## 000002D0
/* 007A4 8086A454 AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 007A8 8086A458 AC4B0000 */  sw      $t3, 0x0000($v0)           ## 00000000
/* 007AC 8086A45C 8FAE009C */  lw      $t6, 0x009C($sp)
.L8086A460:
/* 007B0 8086A460 0C024F61 */  jal     func_80093D84
/* 007B4 8086A464 8DC40000 */  lw      $a0, 0x0000($t6)           ## E3001A01
/* 007B8 8086A468 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 007BC 8086A46C 3C18FA00 */  lui     $t8, 0xFA00                ## $t8 = FA000000
/* 007C0 8086A470 37188080 */  ori     $t8, $t8, 0x8080           ## $t8 = FA008080
/* 007C4 8086A474 244F0008 */  addiu   $t7, $v0, 0x0008           ## $t7 = 00000008
/* 007C8 8086A478 AE0F02D0 */  sw      $t7, 0x02D0($s0)           ## 000002D0
/* 007CC 8086A47C AC580000 */  sw      $t8, 0x0000($v0)           ## 00000000
/* 007D0 8086A480 8FAC0098 */  lw      $t4, 0x0098($sp)
/* 007D4 8086A484 2401AA00 */  addiu   $at, $zero, 0xAA00         ## $at = FFFFAA00
/* 007D8 8086A488 3C18FFFF */  lui     $t8, 0xFFFF                ## $t8 = FFFF0000
/* 007DC 8086A48C 918D0150 */  lbu     $t5, 0x0150($t4)           ## 00000150
/* 007E0 8086A490 37180080 */  ori     $t8, $t8, 0x0080           ## $t8 = FFFF0080
/* 007E4 8086A494 3C0FFB00 */  lui     $t7, 0xFB00                ## $t7 = FB000000
/* 007E8 8086A498 01A15825 */  or      $t3, $t5, $at              ## $t3 = FFFFAA00
/* 007EC 8086A49C AC4B0004 */  sw      $t3, 0x0004($v0)           ## 00000004
/* 007F0 8086A4A0 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 007F4 8086A4A4 24044000 */  addiu   $a0, $zero, 0x4000         ## $a0 = 00004000
/* 007F8 8086A4A8 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 007FC 8086A4AC 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 00800 8086A4B0 AE0E02D0 */  sw      $t6, 0x02D0($s0)           ## 000002D0
/* 00804 8086A4B4 00003025 */  or      $a2, $zero, $zero          ## $a2 = 00000000
/* 00808 8086A4B8 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 0080C 8086A4BC AC580004 */  sw      $t8, 0x0004($v0)           ## 00000004
/* 00810 8086A4C0 0C034421 */  jal     Matrix_RotateZYX
/* 00814 8086A4C4 AC4F0000 */  sw      $t7, 0x0000($v0)           ## 00000000
/* 00818 8086A4C8 8FAC0098 */  lw      $t4, 0x0098($sp)
/* 0081C 8086A4CC 3C018087 */  lui     $at, %hi(D_8086BB8C)       ## $at = 80870000
/* 00820 8086A4D0 24060000 */  addiu   $a2, $zero, 0x0000         ## $a2 = 00000000
/* 00824 8086A4D4 9599014E */  lhu     $t9, 0x014E($t4)           ## 0000014E
/* 00828 8086A4D8 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 0082C 8086A4DC 53200009 */  beql    $t9, $zero, .L8086A504
/* 00830 8086A4E0 44806000 */  mtc1    $zero, $f12                ## $f12 = 0.00
/* 00834 8086A4E4 44806000 */  mtc1    $zero, $f12                ## $f12 = 0.00
/* 00838 8086A4E8 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 0083C 8086A4EC 44066000 */  mfc1    $a2, $f12
/* 00840 8086A4F0 0C034261 */  jal     Matrix_Translate
/* 00844 8086A4F4 46006386 */  mov.s   $f14, $f12
/* 00848 8086A4F8 10000005 */  beq     $zero, $zero, .L8086A510
/* 0084C 8086A4FC 8FAD0098 */  lw      $t5, 0x0098($sp)
/* 00850 8086A500 44806000 */  mtc1    $zero, $f12                ## $f12 = 0.00
.L8086A504:
/* 00854 8086A504 0C034261 */  jal     Matrix_Translate
/* 00858 8086A508 C42EBB8C */  lwc1    $f14, %lo(D_8086BB8C)($at)
/* 0085C 8086A50C 8FAD0098 */  lw      $t5, 0x0098($sp)
.L8086A510:
/* 00860 8086A510 3C018087 */  lui     $at, %hi(D_8086BB90)       ## $at = 80870000
/* 00864 8086A514 C430BB90 */  lwc1    $f16, %lo(D_8086BB90)($at)
/* 00868 8086A518 85AB014C */  lh      $t3, 0x014C($t5)           ## 0000014C
/* 0086C 8086A51C 3C014080 */  lui     $at, 0x4080                ## $at = 40800000
/* 00870 8086A520 44812000 */  mtc1    $at, $f4                   ## $f4 = 4.00
/* 00874 8086A524 448B4000 */  mtc1    $t3, $f8                   ## $f8 = 0.00
/* 00878 8086A528 C5B20160 */  lwc1    $f18, 0x0160($t5)          ## 00000160
/* 0087C 8086A52C 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 00880 8086A530 468042A0 */  cvt.s.w $f10, $f8
/* 00884 8086A534 46105302 */  mul.s   $f12, $f10, $f16
/* 00888 8086A538 44066000 */  mfc1    $a2, $f12
/* 0088C 8086A53C 46049382 */  mul.s   $f14, $f18, $f4
/* 00890 8086A540 0C0342A3 */  jal     Matrix_Scale
/* 00894 8086A544 00000000 */  nop
/* 00898 8086A548 44806000 */  mtc1    $zero, $f12                ## $f12 = 0.00
/* 0089C 8086A54C 3C01C42F */  lui     $at, 0xC42F                ## $at = C42F0000
/* 008A0 8086A550 44817000 */  mtc1    $at, $f14                  ## $f14 = -700.00
/* 008A4 8086A554 44066000 */  mfc1    $a2, $f12
/* 008A8 8086A558 0C034261 */  jal     Matrix_Translate
/* 008AC 8086A55C 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 008B0 8086A560 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 008B4 8086A564 3C0FDA38 */  lui     $t7, 0xDA38                ## $t7 = DA380000
/* 008B8 8086A568 35EF0003 */  ori     $t7, $t7, 0x0003           ## $t7 = DA380003
/* 008BC 8086A56C 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 008C0 8086A570 AE0E02D0 */  sw      $t6, 0x02D0($s0)           ## 000002D0
/* 008C4 8086A574 AC4F0000 */  sw      $t7, 0x0000($v0)           ## 00000000
/* 008C8 8086A578 8FB8009C */  lw      $t8, 0x009C($sp)
/* 008CC 8086A57C 3C058087 */  lui     $a1, %hi(D_8086BB40)       ## $a1 = 80870000
/* 008D0 8086A580 24A5BB40 */  addiu   $a1, $a1, %lo(D_8086BB40)  ## $a1 = 8086BB40
/* 008D4 8086A584 8F040000 */  lw      $a0, 0x0000($t8)           ## 00000000
/* 008D8 8086A588 24060288 */  addiu   $a2, $zero, 0x0288         ## $a2 = 00000288
/* 008DC 8086A58C 0C0346A2 */  jal     Matrix_NewMtx
/* 008E0 8086A590 AFA20054 */  sw      $v0, 0x0054($sp)
/* 008E4 8086A594 8FA30054 */  lw      $v1, 0x0054($sp)
/* 008E8 8086A598 3C198087 */  lui     $t9, %hi(D_8086B960)       ## $t9 = 80870000
/* 008EC 8086A59C 2739B960 */  addiu   $t9, $t9, %lo(D_8086B960)  ## $t9 = 8086B960
/* 008F0 8086A5A0 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 008F4 8086A5A4 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 008F8 8086A5A8 3C09DE00 */  lui     $t1, 0xDE00                ## $t1 = DE000000
/* 008FC 8086A5AC 240301FF */  addiu   $v1, $zero, 0x01FF         ## $v1 = 000001FF
/* 00900 8086A5B0 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 00000008
/* 00904 8086A5B4 AE0C02D0 */  sw      $t4, 0x02D0($s0)           ## 000002D0
/* 00908 8086A5B8 AC490000 */  sw      $t1, 0x0000($v0)           ## 00000000
/* 0090C 8086A5BC AC590004 */  sw      $t9, 0x0004($v0)           ## 00000004
/* 00910 8086A5C0 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 00914 8086A5C4 8FAA008C */  lw      $t2, 0x008C($sp)
/* 00918 8086A5C8 24180004 */  addiu   $t8, $zero, 0x0004         ## $t8 = 00000004
/* 0091C 8086A5CC 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 00920 8086A5D0 AE0B02D0 */  sw      $t3, 0x02D0($s0)           ## 000002D0
/* 00924 8086A5D4 AC490000 */  sw      $t1, 0x0000($v0)           ## 00000000
/* 00928 8086A5D8 8FAD009C */  lw      $t5, 0x009C($sp)
/* 0092C 8086A5DC 000A7080 */  sll     $t6, $t2,  2
/* 00930 8086A5E0 01CA7021 */  addu    $t6, $t6, $t2
/* 00934 8086A5E4 8DA40000 */  lw      $a0, 0x0000($t5)           ## 00000000
/* 00938 8086A5E8 31CF01FF */  andi    $t7, $t6, 0x01FF           ## $t7 = 00000000
/* 0093C 8086A5EC 006F3023 */  subu    $a2, $v1, $t7
/* 00940 8086A5F0 000A5880 */  sll     $t3, $t2,  2
/* 00944 8086A5F4 016A5821 */  addu    $t3, $t3, $t2
/* 00948 8086A5F8 000A7900 */  sll     $t7, $t2,  4
/* 0094C 8086A5FC 000B5840 */  sll     $t3, $t3,  1
/* 00950 8086A600 01EA7823 */  subu    $t7, $t7, $t2
/* 00954 8086A604 000F7840 */  sll     $t7, $t7,  1
/* 00958 8086A608 316D01FF */  andi    $t5, $t3, 0x01FF           ## $t5 = 00000008
/* 0095C 8086A60C AFB80010 */  sw      $t8, 0x0010($sp)
/* 00960 8086A610 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 00964 8086A614 240C0020 */  addiu   $t4, $zero, 0x0020         ## $t4 = 00000020
/* 00968 8086A618 AFAC0014 */  sw      $t4, 0x0014($sp)
/* 0096C 8086A61C AFB90018 */  sw      $t9, 0x0018($sp)
/* 00970 8086A620 31F801FF */  andi    $t8, $t7, 0x01FF           ## $t8 = 00000000
/* 00974 8086A624 00786023 */  subu    $t4, $v1, $t8
/* 00978 8086A628 24190008 */  addiu   $t9, $zero, 0x0008         ## $t9 = 00000008
/* 0097C 8086A62C 006D7023 */  subu    $t6, $v1, $t5
/* 00980 8086A630 240B0010 */  addiu   $t3, $zero, 0x0010         ## $t3 = 00000010
/* 00984 8086A634 AFAB0028 */  sw      $t3, 0x0028($sp)
/* 00988 8086A638 AFAE001C */  sw      $t6, 0x001C($sp)
/* 0098C 8086A63C AFB90024 */  sw      $t9, 0x0024($sp)
/* 00990 8086A640 AFAC0020 */  sw      $t4, 0x0020($sp)
/* 00994 8086A644 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 00998 8086A648 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 0099C 8086A64C 0C0253D0 */  jal     Gfx_TwoTexScroll
/* 009A0 8086A650 AFA2004C */  sw      $v0, 0x004C($sp)
/* 009A4 8086A654 8FA8004C */  lw      $t0, 0x004C($sp)
/* 009A8 8086A658 3C0F8087 */  lui     $t7, %hi(D_8086BA10)       ## $t7 = 80870000
/* 009AC 8086A65C 25EFBA10 */  addiu   $t7, $t7, %lo(D_8086BA10)  ## $t7 = 8086BA10
/* 009B0 8086A660 AD020004 */  sw      $v0, 0x0004($t0)           ## 00000004
/* 009B4 8086A664 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 009B8 8086A668 3C0EDE00 */  lui     $t6, 0xDE00                ## $t6 = DE000000
/* 009BC 8086A66C 3C068087 */  lui     $a2, %hi(D_8086BB54)       ## $a2 = 80870000
/* 009C0 8086A670 244D0008 */  addiu   $t5, $v0, 0x0008           ## $t5 = 00000008
/* 009C4 8086A674 AE0D02D0 */  sw      $t5, 0x02D0($s0)           ## 000002D0
/* 009C8 8086A678 AC4F0004 */  sw      $t7, 0x0004($v0)           ## 00000004
/* 009CC 8086A67C AC4E0000 */  sw      $t6, 0x0000($v0)           ## 00000000
/* 009D0 8086A680 8FB8009C */  lw      $t8, 0x009C($sp)
/* 009D4 8086A684 24C6BB54 */  addiu   $a2, $a2, %lo(D_8086BB54)  ## $a2 = 8086BB54
/* 009D8 8086A688 27A40070 */  addiu   $a0, $sp, 0x0070           ## $a0 = FFFFFFD8
/* 009DC 8086A68C 24070298 */  addiu   $a3, $zero, 0x0298         ## $a3 = 00000298
/* 009E0 8086A690 0C031AD5 */  jal     Graph_CloseDisps
/* 009E4 8086A694 8F050000 */  lw      $a1, 0x0000($t8)           ## 00000000
.L8086A698:
/* 009E8 8086A698 8FBF003C */  lw      $ra, 0x003C($sp)
.L8086A69C:
/* 009EC 8086A69C 8FB00038 */  lw      $s0, 0x0038($sp)
/* 009F0 8086A6A0 27BD0098 */  addiu   $sp, $sp, 0x0098           ## $sp = 00000000
/* 009F4 8086A6A4 03E00008 */  jr      $ra
/* 009F8 8086A6A8 00000000 */  nop
/* 009FC 8086A6AC 00000000 */  nop

