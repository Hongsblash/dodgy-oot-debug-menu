glabel EnClearTag_Init
/* 0040C 809D39BC 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 00410 809D39C0 AFB10018 */  sw      $s1, 0x0018($sp)           
/* 00414 809D39C4 00808825 */  or      $s1, $a0, $zero            ## $s1 = 00000000
/* 00418 809D39C8 AFA50034 */  sw      $a1, 0x0034($sp)           
/* 0041C 809D39CC 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 00420 809D39D0 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 00424 809D39D4 26250194 */  addiu   $a1, $s1, 0x0194           ## $a1 = 00000194
/* 00428 809D39D8 AFB00014 */  sw      $s0, 0x0014($sp)           
/* 0042C 809D39DC 0C0170D9 */  jal     ActorCollider_AllocCylinder
              
/* 00430 809D39E0 AFA50020 */  sw      $a1, 0x0020($sp)           
/* 00434 809D39E4 862E001C */  lh      $t6, 0x001C($s1)           ## 0000001C
/* 00438 809D39E8 24100064 */  addiu   $s0, $zero, 0x0064         ## $s0 = 00000064
/* 0043C 809D39EC 240B0005 */  addiu   $t3, $zero, 0x0005         ## $t3 = 00000005
/* 00440 809D39F0 160E002C */  bne     $s0, $t6, .L809D3AA4       
/* 00444 809D39F4 02203025 */  or      $a2, $s1, $zero            ## $a2 = 00000000
/* 00448 809D39F8 3C01420C */  lui     $at, 0x420C                ## $at = 420C0000
/* 0044C 809D39FC 44812000 */  mtc1    $at, $f4                   ## $f4 = 35.00
/* 00450 809D3A00 240F0064 */  addiu   $t7, $zero, 0x0064         ## $t7 = 00000064
/* 00454 809D3A04 24180046 */  addiu   $t8, $zero, 0x0046         ## $t8 = 00000046
/* 00458 809D3A08 A22F014E */  sb      $t7, 0x014E($s1)           ## 0000014E
/* 0045C 809D3A0C A6380150 */  sh      $t8, 0x0150($s1)           ## 00000150
/* 00460 809D3A10 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 00464 809D3A14 0C00B642 */  jal     func_8002D908              
/* 00468 809D3A18 E6240068 */  swc1    $f4, 0x0068($s1)           ## 00000068
/* 0046C 809D3A1C 00008025 */  or      $s0, $zero, $zero          ## $s0 = 00000000
.L809D3A20:
/* 00470 809D3A20 0C00B5FB */  jal     func_8002D7EC              
/* 00474 809D3A24 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 00478 809D3A28 26100001 */  addiu   $s0, $s0, 0x0001           ## $s0 = 00000001
/* 0047C 809D3A2C 00108400 */  sll     $s0, $s0, 16               
/* 00480 809D3A30 00108403 */  sra     $s0, $s0, 16               
/* 00484 809D3A34 1A00FFFA */  blez    $s0, .L809D3A20            
/* 00488 809D3A38 00000000 */  nop
/* 0048C 809D3A3C 3C01809E */  lui     $at, %hi(D_809DC0E4)       ## $at = 809E0000
/* 00490 809D3A40 C420C0E4 */  lwc1    $f0, %lo(D_809DC0E4)($at)  
/* 00494 809D3A44 3C014000 */  lui     $at, 0x4000                ## $at = 40000000
/* 00498 809D3A48 44813000 */  mtc1    $at, $f6                   ## $f6 = 2.00
/* 0049C 809D3A4C 863900B4 */  lh      $t9, 0x00B4($s1)           ## 000000B4
/* 004A0 809D3A50 3C01428C */  lui     $at, 0x428C                ## $at = 428C0000
/* 004A4 809D3A54 44814000 */  mtc1    $at, $f8                   ## $f8 = 70.00
/* 004A8 809D3A58 00194023 */  subu    $t0, $zero, $t9            
/* 004AC 809D3A5C A62800B4 */  sh      $t0, 0x00B4($s1)           ## 000000B4
/* 004B0 809D3A60 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 004B4 809D3A64 E6200050 */  swc1    $f0, 0x0050($s1)           ## 00000050
/* 004B8 809D3A68 E6200054 */  swc1    $f0, 0x0054($s1)           ## 00000054
/* 004BC 809D3A6C E6260058 */  swc1    $f6, 0x0058($s1)           ## 00000058
/* 004C0 809D3A70 0C00B642 */  jal     func_8002D908              
/* 004C4 809D3A74 E6280068 */  swc1    $f8, 0x0068($s1)           ## 00000068
/* 004C8 809D3A78 3C07809D */  lui     $a3, %hi(D_809D5C6C)       ## $a3 = 809D0000
/* 004CC 809D3A7C 24E75C6C */  addiu   $a3, $a3, %lo(D_809D5C6C)  ## $a3 = 809D5C6C
/* 004D0 809D3A80 8FA40034 */  lw      $a0, 0x0034($sp)           
/* 004D4 809D3A84 8FA50020 */  lw      $a1, 0x0020($sp)           
/* 004D8 809D3A88 0C01712B */  jal     ActorCollider_InitCylinder
              
/* 004DC 809D3A8C 02203025 */  or      $a2, $s1, $zero            ## $a2 = 00000000
/* 004E0 809D3A90 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 004E4 809D3A94 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 004E8 809D3A98 2405182A */  addiu   $a1, $zero, 0x182A         ## $a1 = 0000182A
/* 004EC 809D3A9C 10000034 */  beq     $zero, $zero, .L809D3B70   
/* 004F0 809D3AA0 8FBF001C */  lw      $ra, 0x001C($sp)           
.L809D3AA4:
/* 004F4 809D3AA4 8E290004 */  lw      $t1, 0x0004($s1)           ## 00000004
/* 004F8 809D3AA8 A22B001F */  sb      $t3, 0x001F($s1)           ## 0000001F
/* 004FC 809D3AAC 3C07809D */  lui     $a3, %hi(D_809D5C40)       ## $a3 = 809D0000
/* 00500 809D3AB0 352A0001 */  ori     $t2, $t1, 0x0001           ## $t2 = 00000001
/* 00504 809D3AB4 AE2A0004 */  sw      $t2, 0x0004($s1)           ## 00000004
/* 00508 809D3AB8 8FA50020 */  lw      $a1, 0x0020($sp)           
/* 0050C 809D3ABC 8FA40034 */  lw      $a0, 0x0034($sp)           
/* 00510 809D3AC0 0C01712B */  jal     ActorCollider_InitCylinder
              
/* 00514 809D3AC4 24E75C40 */  addiu   $a3, $a3, %lo(D_809D5C40)  ## $a3 = 809D5C40
/* 00518 809D3AC8 862D001C */  lh      $t5, 0x001C($s1)           ## 0000001C
/* 0051C 809D3ACC 240C0003 */  addiu   $t4, $zero, 0x0003         ## $t4 = 00000003
/* 00520 809D3AD0 A22C00AF */  sb      $t4, 0x00AF($s1)           ## 000000AF
/* 00524 809D3AD4 15A0000E */  bne     $t5, $zero, .L809D3B10     
/* 00528 809D3AD8 3C03809D */  lui     $v1, %hi(D_809D5C30)       ## $v1 = 809D0000
/* 0052C 809D3ADC 24050001 */  addiu   $a1, $zero, 0x0001         ## $a1 = 00000001
/* 00530 809D3AE0 240E0046 */  addiu   $t6, $zero, 0x0046         ## $t6 = 00000046
/* 00534 809D3AE4 240F00FA */  addiu   $t7, $zero, 0x00FA         ## $t7 = 000000FA
/* 00538 809D3AE8 24180002 */  addiu   $t8, $zero, 0x0002         ## $t8 = 00000002
/* 0053C 809D3AEC 24194000 */  addiu   $t9, $zero, 0x4000         ## $t9 = 00004000
/* 00540 809D3AF0 24080014 */  addiu   $t0, $zero, 0x0014         ## $t0 = 00000014
/* 00544 809D3AF4 A62E0150 */  sh      $t6, 0x0150($s1)           ## 00000150
/* 00548 809D3AF8 A62F0152 */  sh      $t7, 0x0152($s1)           ## 00000152
/* 0054C 809D3AFC A238014E */  sb      $t8, 0x014E($s1)           ## 0000014E
/* 00550 809D3B00 A6390030 */  sh      $t9, 0x0030($s1)           ## 00000030
/* 00554 809D3B04 A22501E0 */  sb      $a1, 0x01E0($s1)           ## 000001E0
/* 00558 809D3B08 A63001FC */  sh      $s0, 0x01FC($s1)           ## 000001FC
/* 0055C 809D3B0C A6280154 */  sh      $t0, 0x0154($s1)           ## 00000154
.L809D3B10:
/* 00560 809D3B10 24635C30 */  addiu   $v1, $v1, %lo(D_809D5C30)  ## $v1 = 809D5C30
/* 00564 809D3B14 90690000 */  lbu     $t1, 0x0000($v1)           ## 809D5C30
/* 00568 809D3B18 3C04809E */  lui     $a0, %hi(D_809DC3D0)       ## $a0 = 809E0000
/* 0056C 809D3B1C 24050001 */  addiu   $a1, $zero, 0x0001         ## $a1 = 00000001
/* 00570 809D3B20 15200012 */  bne     $t1, $zero, .L809D3B6C     
/* 00574 809D3B24 2484C3D0 */  addiu   $a0, $a0, %lo(D_809DC3D0)  ## $a0 = 809DC3D0
/* 00578 809D3B28 8FAA0034 */  lw      $t2, 0x0034($sp)           
/* 0057C 809D3B2C A0650000 */  sb      $a1, 0x0000($v1)           ## 809D5C30
/* 00580 809D3B30 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00584 809D3B34 002A0821 */  addu    $at, $at, $t2              
/* 00588 809D3B38 AC241E10 */  sw      $a0, 0x1E10($at)           ## 00011E10
/* 0058C 809D3B3C 2403006C */  addiu   $v1, $zero, 0x006C         ## $v1 = 0000006C
/* 00590 809D3B40 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L809D3B44:
/* 00594 809D3B44 00430019 */  multu   $v0, $v1                   
/* 00598 809D3B48 24420001 */  addiu   $v0, $v0, 0x0001           ## $v0 = 00000001
/* 0059C 809D3B4C 00021400 */  sll     $v0, $v0, 16               
/* 005A0 809D3B50 00021403 */  sra     $v0, $v0, 16               
/* 005A4 809D3B54 28410064 */  slti    $at, $v0, 0x0064           
/* 005A8 809D3B58 00005812 */  mflo    $t3                        
/* 005AC 809D3B5C 008B6021 */  addu    $t4, $a0, $t3              
/* 005B0 809D3B60 1420FFF8 */  bne     $at, $zero, .L809D3B44     
/* 005B4 809D3B64 A1800000 */  sb      $zero, 0x0000($t4)         ## 00000003
/* 005B8 809D3B68 A225014D */  sb      $a1, 0x014D($s1)           ## 0000014D
.L809D3B6C:
/* 005BC 809D3B6C 8FBF001C */  lw      $ra, 0x001C($sp)           
.L809D3B70:
/* 005C0 809D3B70 8FB00014 */  lw      $s0, 0x0014($sp)           
/* 005C4 809D3B74 8FB10018 */  lw      $s1, 0x0018($sp)           
/* 005C8 809D3B78 03E00008 */  jr      $ra                        
/* 005CC 809D3B7C 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000


