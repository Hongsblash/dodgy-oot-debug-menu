.rdata
glabel D_80B37B54
    .asciz "../z_en_wf.c"
    .balign 4

glabel D_80B37B64
    .asciz "../z_en_wf.c"
    .balign 4

.text
glabel EnWf_Draw
/* 0398C 80B3763C 27BDFFA8 */  addiu   $sp, $sp, 0xFFA8           ## $sp = FFFFFFA8
/* 03990 80B37640 AFBF002C */  sw      $ra, 0x002C($sp)
/* 03994 80B37644 AFB00028 */  sw      $s0, 0x0028($sp)
/* 03998 80B37648 AFA5005C */  sw      $a1, 0x005C($sp)
/* 0399C 80B3764C 8CA50000 */  lw      $a1, 0x0000($a1)           ## 00000000
/* 039A0 80B37650 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 039A4 80B37654 3C0680B3 */  lui     $a2, %hi(D_80B37B54)       ## $a2 = 80B30000
/* 039A8 80B37658 24C67B54 */  addiu   $a2, $a2, %lo(D_80B37B54)  ## $a2 = 80B37B54
/* 039AC 80B3765C 27A40040 */  addiu   $a0, $sp, 0x0040           ## $a0 = FFFFFFE8
/* 039B0 80B37660 2407086D */  addiu   $a3, $zero, 0x086D         ## $a3 = 0000086D
/* 039B4 80B37664 0C031AB1 */  jal     Graph_OpenDisps
/* 039B8 80B37668 AFA50050 */  sw      $a1, 0x0050($sp)
/* 039BC 80B3766C 8E0F02D4 */  lw      $t7, 0x02D4($s0)           ## 000002D4
/* 039C0 80B37670 8FA80050 */  lw      $t0, 0x0050($sp)
/* 039C4 80B37674 8FB9005C */  lw      $t9, 0x005C($sp)
/* 039C8 80B37678 55E00005 */  bnel    $t7, $zero, .L80B37690
/* 039CC 80B3767C 8F240000 */  lw      $a0, 0x0000($t9)           ## 00000000
/* 039D0 80B37680 86180300 */  lh      $t8, 0x0300($s0)           ## 00000300
/* 039D4 80B37684 5700005F */  bnel    $t8, $zero, .L80B37804
/* 039D8 80B37688 8FAF005C */  lw      $t7, 0x005C($sp)
/* 039DC 80B3768C 8F240000 */  lw      $a0, 0x0000($t9)           ## 00000000
.L80B37690:
/* 039E0 80B37690 0C024F46 */  jal     func_80093D18
/* 039E4 80B37694 AFA80050 */  sw      $t0, 0x0050($sp)
/* 039E8 80B37698 8609001C */  lh      $t1, 0x001C($s0)           ## 0000001C
/* 039EC 80B3769C 8FA80050 */  lw      $t0, 0x0050($sp)
/* 039F0 80B376A0 5520001B */  bnel    $t1, $zero, .L80B37710
/* 039F4 80B376A4 8D0302C0 */  lw      $v1, 0x02C0($t0)           ## 000002C0
/* 039F8 80B376A8 8D0302C0 */  lw      $v1, 0x02C0($t0)           ## 000002C0
/* 039FC 80B376AC 3C0BDB06 */  lui     $t3, 0xDB06                ## $t3 = DB060000
/* 03A00 80B376B0 356B0020 */  ori     $t3, $t3, 0x0020           ## $t3 = DB060020
/* 03A04 80B376B4 246A0008 */  addiu   $t2, $v1, 0x0008           ## $t2 = 00000008
/* 03A08 80B376B8 AD0A02C0 */  sw      $t2, 0x02C0($t0)           ## 000002C0
/* 03A0C 80B376BC AC6B0000 */  sw      $t3, 0x0000($v1)           ## 00000000
/* 03A10 80B376C0 920C0302 */  lbu     $t4, 0x0302($s0)           ## 00000302
/* 03A14 80B376C4 3C0480B3 */  lui     $a0, %hi(D_80B37AF4)       ## $a0 = 80B30000
/* 03A18 80B376C8 3C098016 */  lui     $t1, %hi(gSegments)
/* 03A1C 80B376CC 000C6880 */  sll     $t5, $t4,  2
/* 03A20 80B376D0 008D2021 */  addu    $a0, $a0, $t5
/* 03A24 80B376D4 8C847AF4 */  lw      $a0, %lo(D_80B37AF4)($a0)
/* 03A28 80B376D8 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 03A2C 80B376DC 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 03A30 80B376E0 00047900 */  sll     $t7, $a0,  4
/* 03A34 80B376E4 000FC702 */  srl     $t8, $t7, 28
/* 03A38 80B376E8 0018C880 */  sll     $t9, $t8,  2
/* 03A3C 80B376EC 01394821 */  addu    $t1, $t1, $t9
/* 03A40 80B376F0 8D296FA8 */  lw      $t1, %lo(gSegments)($t1)
/* 03A44 80B376F4 00817024 */  and     $t6, $a0, $at
/* 03A48 80B376F8 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 03A4C 80B376FC 01C95021 */  addu    $t2, $t6, $t1
/* 03A50 80B37700 01415821 */  addu    $t3, $t2, $at
/* 03A54 80B37704 10000019 */  beq     $zero, $zero, .L80B3776C
/* 03A58 80B37708 AC6B0004 */  sw      $t3, 0x0004($v1)           ## 00000004
/* 03A5C 80B3770C 8D0302C0 */  lw      $v1, 0x02C0($t0)           ## 000002C0
.L80B37710:
/* 03A60 80B37710 3C0DDB06 */  lui     $t5, 0xDB06                ## $t5 = DB060000
/* 03A64 80B37714 35AD0020 */  ori     $t5, $t5, 0x0020           ## $t5 = DB060020
/* 03A68 80B37718 246C0008 */  addiu   $t4, $v1, 0x0008           ## $t4 = 00000008
/* 03A6C 80B3771C AD0C02C0 */  sw      $t4, 0x02C0($t0)           ## 000002C0
/* 03A70 80B37720 AC6D0000 */  sw      $t5, 0x0000($v1)           ## 00000000
/* 03A74 80B37724 920F0302 */  lbu     $t7, 0x0302($s0)           ## 00000302
/* 03A78 80B37728 3C0480B3 */  lui     $a0, %hi(D_80B37B04)       ## $a0 = 80B30000
/* 03A7C 80B3772C 3C0B8016 */  lui     $t3, %hi(gSegments)
/* 03A80 80B37730 000FC080 */  sll     $t8, $t7,  2
/* 03A84 80B37734 00982021 */  addu    $a0, $a0, $t8
/* 03A88 80B37738 8C847B04 */  lw      $a0, %lo(D_80B37B04)($a0)
/* 03A8C 80B3773C 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 03A90 80B37740 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 03A94 80B37744 00047100 */  sll     $t6, $a0,  4
/* 03A98 80B37748 000E4F02 */  srl     $t1, $t6, 28
/* 03A9C 80B3774C 00095080 */  sll     $t2, $t1,  2
/* 03AA0 80B37750 016A5821 */  addu    $t3, $t3, $t2
/* 03AA4 80B37754 8D6B6FA8 */  lw      $t3, %lo(gSegments)($t3)
/* 03AA8 80B37758 0081C824 */  and     $t9, $a0, $at
/* 03AAC 80B3775C 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 03AB0 80B37760 032B6021 */  addu    $t4, $t9, $t3
/* 03AB4 80B37764 01816821 */  addu    $t5, $t4, $at
/* 03AB8 80B37768 AC6D0004 */  sw      $t5, 0x0004($v1)           ## 00000004
.L80B3776C:
/* 03ABC 80B3776C 8E05018C */  lw      $a1, 0x018C($s0)           ## 0000018C
/* 03AC0 80B37770 8E0601A8 */  lw      $a2, 0x01A8($s0)           ## 000001A8
/* 03AC4 80B37774 9207018A */  lbu     $a3, 0x018A($s0)           ## 0000018A
/* 03AC8 80B37778 3C0F80B3 */  lui     $t7, %hi(func_80B37454)    ## $t7 = 80B30000
/* 03ACC 80B3777C 3C1880B3 */  lui     $t8, %hi(func_80B37494)    ## $t8 = 80B30000
/* 03AD0 80B37780 27187494 */  addiu   $t8, $t8, %lo(func_80B37494) ## $t8 = 80B37494
/* 03AD4 80B37784 25EF7454 */  addiu   $t7, $t7, %lo(func_80B37454) ## $t7 = 80B37454
/* 03AD8 80B37788 AFAF0010 */  sw      $t7, 0x0010($sp)
/* 03ADC 80B3778C AFB80014 */  sw      $t8, 0x0014($sp)
/* 03AE0 80B37790 AFB00018 */  sw      $s0, 0x0018($sp)
/* 03AE4 80B37794 0C0286B2 */  jal     SkelAnime_DrawSV
/* 03AE8 80B37798 8FA4005C */  lw      $a0, 0x005C($sp)
/* 03AEC 80B3779C 860E02E4 */  lh      $t6, 0x02E4($s0)           ## 000002E4
/* 03AF0 80B377A0 51C00018 */  beql    $t6, $zero, .L80B37804
/* 03AF4 80B377A4 8FAF005C */  lw      $t7, 0x005C($sp)
/* 03AF8 80B377A8 861902E4 */  lh      $t9, 0x02E4($s0)           ## 000002E4
/* 03AFC 80B377AC 92090114 */  lbu     $t1, 0x0114($s0)           ## 00000114
/* 03B00 80B377B0 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 03B04 80B377B4 272BFFFF */  addiu   $t3, $t9, 0xFFFF           ## $t3 = FFFFFFFF
/* 03B08 80B377B8 A60B02E4 */  sh      $t3, 0x02E4($s0)           ## 000002E4
/* 03B0C 80B377BC 860302E4 */  lh      $v1, 0x02E4($s0)           ## 000002E4
/* 03B10 80B377C0 252A0001 */  addiu   $t2, $t1, 0x0001           ## $t2 = 00000001
/* 03B14 80B377C4 A20A0114 */  sb      $t2, 0x0114($s0)           ## 00000114
/* 03B18 80B377C8 306C0003 */  andi    $t4, $v1, 0x0003           ## $t4 = 00000000
/* 03B1C 80B377CC 1580000C */  bne     $t4, $zero, .L80B37800
/* 03B20 80B377D0 00031083 */  sra     $v0, $v1,  2
/* 03B24 80B377D4 00026880 */  sll     $t5, $v0,  2
/* 03B28 80B377D8 01A26823 */  subu    $t5, $t5, $v0
/* 03B2C 80B377DC 000D6840 */  sll     $t5, $t5,  1
/* 03B30 80B377E0 020D3021 */  addu    $a2, $s0, $t5
/* 03B34 80B377E4 24C6014C */  addiu   $a2, $a2, 0x014C           ## $a2 = 0000014C
/* 03B38 80B377E8 8FA4005C */  lw      $a0, 0x005C($sp)
/* 03B3C 80B377EC 2407004B */  addiu   $a3, $zero, 0x004B         ## $a3 = 0000004B
/* 03B40 80B377F0 AFA00010 */  sw      $zero, 0x0010($sp)
/* 03B44 80B377F4 AFA00014 */  sw      $zero, 0x0014($sp)
/* 03B48 80B377F8 0C00A953 */  jal     func_8002A54C
/* 03B4C 80B377FC AFA20018 */  sw      $v0, 0x0018($sp)
.L80B37800:
/* 03B50 80B37800 8FAF005C */  lw      $t7, 0x005C($sp)
.L80B37804:
/* 03B54 80B37804 3C0680B3 */  lui     $a2, %hi(D_80B37B64)       ## $a2 = 80B30000
/* 03B58 80B37808 24C67B64 */  addiu   $a2, $a2, %lo(D_80B37B64)  ## $a2 = 80B37B64
/* 03B5C 80B3780C 27A40040 */  addiu   $a0, $sp, 0x0040           ## $a0 = FFFFFFE8
/* 03B60 80B37810 2407088E */  addiu   $a3, $zero, 0x088E         ## $a3 = 0000088E
/* 03B64 80B37814 0C031AD5 */  jal     Graph_CloseDisps
/* 03B68 80B37818 8DE50000 */  lw      $a1, 0x0000($t7)           ## 00000000
/* 03B6C 80B3781C 8FBF002C */  lw      $ra, 0x002C($sp)
/* 03B70 80B37820 8FB00028 */  lw      $s0, 0x0028($sp)
/* 03B74 80B37824 27BD0058 */  addiu   $sp, $sp, 0x0058           ## $sp = 00000000
/* 03B78 80B37828 03E00008 */  jr      $ra
/* 03B7C 80B3782C 00000000 */  nop
