glabel EnBili_Draw
/* 01BF8 809C1498 27BDFFA8 */  addiu   $sp, $sp, 0xFFA8           ## $sp = FFFFFFA8
/* 01BFC 809C149C AFBF002C */  sw      $ra, 0x002C($sp)
/* 01C00 809C14A0 AFB10028 */  sw      $s1, 0x0028($sp)
/* 01C04 809C14A4 AFB00024 */  sw      $s0, 0x0024($sp)
/* 01C08 809C14A8 AFA5005C */  sw      $a1, 0x005C($sp)
/* 01C0C 809C14AC 8CA50000 */  lw      $a1, 0x0000($a1)           ## 00000000
/* 01C10 809C14B0 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 01C14 809C14B4 3C06809C */  lui     $a2, %hi(D_809C1710)       ## $a2 = 809C0000
/* 01C18 809C14B8 24C61710 */  addiu   $a2, $a2, %lo(D_809C1710)  ## $a2 = 809C1710
/* 01C1C 809C14BC 27A40040 */  addiu   $a0, $sp, 0x0040           ## $a0 = FFFFFFE8
/* 01C20 809C14C0 240705F1 */  addiu   $a3, $zero, 0x05F1         ## $a3 = 000005F1
/* 01C24 809C14C4 0C031AB1 */  jal     func_800C6AC4
/* 01C28 809C14C8 00A08825 */  or      $s1, $a1, $zero            ## $s1 = 00000000
/* 01C2C 809C14CC 8FAF005C */  lw      $t7, 0x005C($sp)
/* 01C30 809C14D0 0C024F61 */  jal     func_80093D84
/* 01C34 809C14D4 8DE40000 */  lw      $a0, 0x0000($t7)           ## 00000000
/* 01C38 809C14D8 92020194 */  lbu     $v0, 0x0194($s0)           ## 00000194
/* 01C3C 809C14DC 3C07809C */  lui     $a3, %hi(func_809C13A8)    ## $a3 = 809C0000
/* 01C40 809C14E0 24180007 */  addiu   $t8, $zero, 0x0007         ## $t8 = 00000007
/* 01C44 809C14E4 28410008 */  slti    $at, $v0, 0x0008
/* 01C48 809C14E8 14200003 */  bne     $at, $zero, .L809C14F8
/* 01C4C 809C14EC 24E713A8 */  addiu   $a3, $a3, %lo(func_809C13A8) ## $a3 = 809C13A8
/* 01C50 809C14F0 10000002 */  beq     $zero, $zero, .L809C14FC
/* 01C54 809C14F4 A2180194 */  sb      $t8, 0x0194($s0)           ## 00000194
.L809C14F8:
/* 01C58 809C14F8 A2020194 */  sb      $v0, 0x0194($s0)           ## 00000194
.L809C14FC:
/* 01C5C 809C14FC 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
/* 01C60 809C1500 3C08DB06 */  lui     $t0, 0xDB06                ## $t0 = DB060000
/* 01C64 809C1504 35080020 */  ori     $t0, $t0, 0x0020           ## $t0 = DB060020
/* 01C68 809C1508 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 01C6C 809C150C AE3902D0 */  sw      $t9, 0x02D0($s1)           ## 000002D0
/* 01C70 809C1510 AC480000 */  sw      $t0, 0x0000($v0)           ## 00000000
/* 01C74 809C1514 92090194 */  lbu     $t1, 0x0194($s0)           ## 00000194
/* 01C78 809C1518 3C04809C */  lui     $a0, %hi(D_809C16CC)       ## $a0 = 809C0000
/* 01C7C 809C151C 3C0F8016 */  lui     $t7, 0x8016                ## $t7 = 80160000
/* 01C80 809C1520 00095080 */  sll     $t2, $t1,  2
/* 01C84 809C1524 008A2021 */  addu    $a0, $a0, $t2
/* 01C88 809C1528 8C8416CC */  lw      $a0, %lo(D_809C16CC)($a0)
/* 01C8C 809C152C 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 01C90 809C1530 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 01C94 809C1534 00046100 */  sll     $t4, $a0,  4
/* 01C98 809C1538 000C6F02 */  srl     $t5, $t4, 28
/* 01C9C 809C153C 000D7080 */  sll     $t6, $t5,  2
/* 01CA0 809C1540 01EE7821 */  addu    $t7, $t7, $t6
/* 01CA4 809C1544 8DEF6FA8 */  lw      $t7, 0x6FA8($t7)           ## 80166FA8
/* 01CA8 809C1548 00815824 */  and     $t3, $a0, $at
/* 01CAC 809C154C 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 01CB0 809C1550 016FC021 */  addu    $t8, $t3, $t7
/* 01CB4 809C1554 0301C821 */  addu    $t9, $t8, $at
/* 01CB8 809C1558 AC590004 */  sw      $t9, 0x0004($v0)           ## 00000004
/* 01CBC 809C155C 8E090190 */  lw      $t1, 0x0190($s0)           ## 00000190
/* 01CC0 809C1560 3C08809C */  lui     $t0, %hi(func_809C02B8)    ## $t0 = 809C0000
/* 01CC4 809C1564 250802B8 */  addiu   $t0, $t0, %lo(func_809C02B8) ## $t0 = 809C02B8
/* 01CC8 809C1568 1509000F */  bne     $t0, $t1, .L809C15A8
/* 01CCC 809C156C 3C18DB06 */  lui     $t8, 0xDB06                ## $t8 = DB060000
/* 01CD0 809C1570 860A0196 */  lh      $t2, 0x0196($s0)           ## 00000196
/* 01CD4 809C1574 3C0EDB06 */  lui     $t6, 0xDB06                ## $t6 = DB060000
/* 01CD8 809C1578 35CE0024 */  ori     $t6, $t6, 0x0024           ## $t6 = DB060024
/* 01CDC 809C157C 314C0001 */  andi    $t4, $t2, 0x0001           ## $t4 = 00000000
/* 01CE0 809C1580 5180000A */  beql    $t4, $zero, .L809C15AC
/* 01CE4 809C1584 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
/* 01CE8 809C1588 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
/* 01CEC 809C158C 3C0B809C */  lui     $t3, %hi(D_809C16F0)       ## $t3 = 809C0000
/* 01CF0 809C1590 256B16F0 */  addiu   $t3, $t3, %lo(D_809C16F0)  ## $t3 = 809C16F0
/* 01CF4 809C1594 244D0008 */  addiu   $t5, $v0, 0x0008           ## $t5 = 00000008
/* 01CF8 809C1598 AE2D02D0 */  sw      $t5, 0x02D0($s1)           ## 000002D0
/* 01CFC 809C159C AC4B0004 */  sw      $t3, 0x0004($v0)           ## 00000004
/* 01D00 809C15A0 10000009 */  beq     $zero, $zero, .L809C15C8
/* 01D04 809C15A4 AC4E0000 */  sw      $t6, 0x0000($v0)           ## 00000000
.L809C15A8:
/* 01D08 809C15A8 8E2202D0 */  lw      $v0, 0x02D0($s1)           ## 000002D0
.L809C15AC:
/* 01D0C 809C15AC 3C19809C */  lui     $t9, %hi(D_809C1700)       ## $t9 = 809C0000
/* 01D10 809C15B0 27391700 */  addiu   $t9, $t9, %lo(D_809C1700)  ## $t9 = 809C1700
/* 01D14 809C15B4 244F0008 */  addiu   $t7, $v0, 0x0008           ## $t7 = 00000008
/* 01D18 809C15B8 AE2F02D0 */  sw      $t7, 0x02D0($s1)           ## 000002D0
/* 01D1C 809C15BC 37180024 */  ori     $t8, $t8, 0x0024           ## $t8 = DB060024
/* 01D20 809C15C0 AC580000 */  sw      $t8, 0x0000($v0)           ## 00000000
/* 01D24 809C15C4 AC590004 */  sw      $t9, 0x0004($v0)           ## 00000004
.L809C15C8:
/* 01D28 809C15C8 8E050150 */  lw      $a1, 0x0150($s0)           ## 00000150
/* 01D2C 809C15CC 8E06016C */  lw      $a2, 0x016C($s0)           ## 0000016C
/* 01D30 809C15D0 AFB00014 */  sw      $s0, 0x0014($sp)
/* 01D34 809C15D4 AFA00010 */  sw      $zero, 0x0010($sp)
/* 01D38 809C15D8 8E2802D0 */  lw      $t0, 0x02D0($s1)           ## 000002D0
/* 01D3C 809C15DC 8FA4005C */  lw      $a0, 0x005C($sp)
/* 01D40 809C15E0 0C0288A2 */  jal     SkelAnime_Draw2
/* 01D44 809C15E4 AFA80018 */  sw      $t0, 0x0018($sp)
/* 01D48 809C15E8 AE2202D0 */  sw      $v0, 0x02D0($s1)           ## 000002D0
/* 01D4C 809C15EC 8FA9005C */  lw      $t1, 0x005C($sp)
/* 01D50 809C15F0 3C06809C */  lui     $a2, %hi(D_809C1720)       ## $a2 = 809C0000
/* 01D54 809C15F4 24C61720 */  addiu   $a2, $a2, %lo(D_809C1720)  ## $a2 = 809C1720
/* 01D58 809C15F8 27A40040 */  addiu   $a0, $sp, 0x0040           ## $a0 = FFFFFFE8
/* 01D5C 809C15FC 24070610 */  addiu   $a3, $zero, 0x0610         ## $a3 = 00000610
/* 01D60 809C1600 0C031AD5 */  jal     func_800C6B54
/* 01D64 809C1604 8D250000 */  lw      $a1, 0x0000($t1)           ## 00000000
/* 01D68 809C1608 8FBF002C */  lw      $ra, 0x002C($sp)
/* 01D6C 809C160C 8FB00024 */  lw      $s0, 0x0024($sp)
/* 01D70 809C1610 8FB10028 */  lw      $s1, 0x0028($sp)
/* 01D74 809C1614 03E00008 */  jr      $ra
/* 01D78 809C1618 27BD0058 */  addiu   $sp, $sp, 0x0058           ## $sp = 00000000
/* 01D7C 809C161C 00000000 */  nop

