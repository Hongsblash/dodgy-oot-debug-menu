glabel func_80A5E39C
/* 030AC 80A5E39C 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 030B0 80A5E3A0 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 030B4 80A5E3A4 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 030B8 80A5E3A8 AFA5002C */  sw      $a1, 0x002C($sp)           
/* 030BC 80A5E3AC 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 030C0 80A5E3B0 C4800068 */  lwc1    $f0, 0x0068($a0)           ## 00000068
/* 030C4 80A5E3B4 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 030C8 80A5E3B8 3C0180A6 */  lui     $at, %hi(D_80A668CC)       ## $at = 80A60000
/* 030CC 80A5E3BC 4600103C */  c.lt.s  $f2, $f0                   
/* 030D0 80A5E3C0 00000000 */  nop
/* 030D4 80A5E3C4 4502000B */  bc1fl   .L80A5E3F4                 
/* 030D8 80A5E3C8 8E0E01F0 */  lw      $t6, 0x01F0($s0)           ## 000001F0
/* 030DC 80A5E3CC C42468CC */  lwc1    $f4, %lo(D_80A668CC)($at)  
/* 030E0 80A5E3D0 46040181 */  sub.s   $f6, $f0, $f4              
/* 030E4 80A5E3D4 E4860068 */  swc1    $f6, 0x0068($a0)           ## 00000068
/* 030E8 80A5E3D8 C4880068 */  lwc1    $f8, 0x0068($a0)           ## 00000068
/* 030EC 80A5E3DC 4602403C */  c.lt.s  $f8, $f2                   
/* 030F0 80A5E3E0 00000000 */  nop
/* 030F4 80A5E3E4 45020003 */  bc1fl   .L80A5E3F4                 
/* 030F8 80A5E3E8 8E0E01F0 */  lw      $t6, 0x01F0($s0)           ## 000001F0
/* 030FC 80A5E3EC E4820068 */  swc1    $f2, 0x0068($a0)           ## 00000068
/* 03100 80A5E3F0 8E0E01F0 */  lw      $t6, 0x01F0($s0)           ## 000001F0
.L80A5E3F4:
/* 03104 80A5E3F4 3C0141E8 */  lui     $at, 0x41E8                ## $at = 41E80000
/* 03108 80A5E3F8 31CF0400 */  andi    $t7, $t6, 0x0400           ## $t7 = 00000000
/* 0310C 80A5E3FC 51E00035 */  beql    $t7, $zero, .L80A5E4D4     
/* 03110 80A5E400 3C0141E8 */  lui     $at, 0x41E8                ## $at = 41E80000
/* 03114 80A5E404 44815000 */  mtc1    $at, $f10                  ## $f10 = 29.00
/* 03118 80A5E408 C61001C4 */  lwc1    $f16, 0x01C4($s0)          ## 000001C4
/* 0311C 80A5E40C 4610503C */  c.lt.s  $f10, $f16                 
/* 03120 80A5E410 00000000 */  nop
/* 03124 80A5E414 4502002F */  bc1fl   .L80A5E4D4                 
/* 03128 80A5E418 3C0141E8 */  lui     $at, 0x41E8                ## $at = 41E80000
/* 0312C 80A5E41C 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 03130 80A5E420 E6020068 */  swc1    $f2, 0x0068($s0)           ## 00000068
/* 03134 80A5E424 3C013FE0 */  lui     $at, 0x3FE0                ## $at = 3FE00000
/* 03138 80A5E428 44819800 */  mtc1    $at, $f19                  ## $f19 = 1.75
/* 0313C 80A5E42C 44809000 */  mtc1    $zero, $f18                ## $f18 = 0.00
/* 03140 80A5E430 46000121 */  cvt.d.s $f4, $f0                   
/* 03144 80A5E434 4624903C */  c.lt.d  $f18, $f4                  
/* 03148 80A5E438 00000000 */  nop
/* 0314C 80A5E43C 45000020 */  bc1f    .L80A5E4C0                 
/* 03150 80A5E440 00000000 */  nop
/* 03154 80A5E444 8E190228 */  lw      $t9, 0x0228($s0)           ## 00000228
/* 03158 80A5E448 2605021C */  addiu   $a1, $s0, 0x021C           ## $a1 = 0000021C
/* 0315C 80A5E44C 3C078013 */  lui     $a3, 0x8013                ## $a3 = 80130000
/* 03160 80A5E450 ACB90000 */  sw      $t9, 0x0000($a1)           ## 0000021C
/* 03164 80A5E454 8E18022C */  lw      $t8, 0x022C($s0)           ## 0000022C
/* 03168 80A5E458 3C0A8013 */  lui     $t2, 0x8013                ## $t2 = 80130000
/* 0316C 80A5E45C 24E733E0 */  addiu   $a3, $a3, 0x33E0           ## $a3 = 801333E0
/* 03170 80A5E460 ACB80004 */  sw      $t8, 0x0004($a1)           ## 00000220
/* 03174 80A5E464 8E190230 */  lw      $t9, 0x0230($s0)           ## 00000230
/* 03178 80A5E468 254A33E8 */  addiu   $t2, $t2, 0x33E8           ## $t2 = 801333E8
/* 0317C 80A5E46C 24042805 */  addiu   $a0, $zero, 0x2805         ## $a0 = 00002805
/* 03180 80A5E470 ACB90008 */  sw      $t9, 0x0008($a1)           ## 00000224
/* 03184 80A5E474 8E0801F0 */  lw      $t0, 0x01F0($s0)           ## 000001F0
/* 03188 80A5E478 24060004 */  addiu   $a2, $zero, 0x0004         ## $a2 = 00000004
/* 0318C 80A5E47C 00084900 */  sll     $t1, $t0,  4               
/* 03190 80A5E480 05230005 */  bgezl   $t1, .L80A5E498            
/* 03194 80A5E484 44806000 */  mtc1    $zero, $f12                ## $f12 = 0.00
/* 03198 80A5E488 AFA70010 */  sw      $a3, 0x0010($sp)           
/* 0319C 80A5E48C 0C03DCE3 */  jal     Audio_PlaySoundGeneral
              
/* 031A0 80A5E490 AFAA0014 */  sw      $t2, 0x0014($sp)           
/* 031A4 80A5E494 44806000 */  mtc1    $zero, $f12                ## $f12 = 0.00
.L80A5E498:
/* 031A8 80A5E498 240500B4 */  addiu   $a1, $zero, 0x00B4         ## $a1 = 000000B4
/* 031AC 80A5E49C 24060014 */  addiu   $a2, $zero, 0x0014         ## $a2 = 00000014
/* 031B0 80A5E4A0 0C02A800 */  jal     func_800AA000              
/* 031B4 80A5E4A4 24070064 */  addiu   $a3, $zero, 0x0064         ## $a3 = 00000064
/* 031B8 80A5E4A8 8E0B01F0 */  lw      $t3, 0x01F0($s0)           ## 000001F0
/* 031BC 80A5E4AC 2401FBFF */  addiu   $at, $zero, 0xFBFF         ## $at = FFFFFBFF
/* 031C0 80A5E4B0 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 031C4 80A5E4B4 01616024 */  and     $t4, $t3, $at              
/* 031C8 80A5E4B8 10000005 */  beq     $zero, $zero, .L80A5E4D0   
/* 031CC 80A5E4BC AE0C01F0 */  sw      $t4, 0x01F0($s0)           ## 000001F0
.L80A5E4C0:
/* 031D0 80A5E4C0 0C2973CA */  jal     func_80A5CF28              
/* 031D4 80A5E4C4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 031D8 80A5E4C8 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 031DC 80A5E4CC 00000000 */  nop
.L80A5E4D0:
/* 031E0 80A5E4D0 3C0141E8 */  lui     $at, 0x41E8                ## $at = 41E80000
.L80A5E4D4:
/* 031E4 80A5E4D4 44813000 */  mtc1    $at, $f6                   ## $f6 = 29.00
/* 031E8 80A5E4D8 C60801C4 */  lwc1    $f8, 0x01C4($s0)           ## 000001C4
/* 031EC 80A5E4DC 3C014040 */  lui     $at, 0x4040                ## $at = 40400000
/* 031F0 80A5E4E0 4608303C */  c.lt.s  $f6, $f8                   
/* 031F4 80A5E4E4 00000000 */  nop
/* 031F8 80A5E4E8 45020004 */  bc1fl   .L80A5E4FC                 
/* 031FC 80A5E4EC 44810000 */  mtc1    $at, $f0                   ## $f0 = 3.00
/* 03200 80A5E4F0 1000000C */  beq     $zero, $zero, .L80A5E524   
/* 03204 80A5E4F4 E6020068 */  swc1    $f2, 0x0068($s0)           ## 00000068
/* 03208 80A5E4F8 44810000 */  mtc1    $at, $f0                   ## $f0 = 3.00
.L80A5E4FC:
/* 0320C 80A5E4FC C60A0068 */  lwc1    $f10, 0x0068($s0)          ## 00000068
/* 03210 80A5E500 460A003C */  c.lt.s  $f0, $f10                  
/* 03214 80A5E504 00000000 */  nop
/* 03218 80A5E508 45000006 */  bc1f    .L80A5E524                 
/* 0321C 80A5E50C 00000000 */  nop
/* 03220 80A5E510 8E0D01F0 */  lw      $t5, 0x01F0($s0)           ## 000001F0
/* 03224 80A5E514 31AE0010 */  andi    $t6, $t5, 0x0010           ## $t6 = 00000000
/* 03228 80A5E518 11C00002 */  beq     $t6, $zero, .L80A5E524     
/* 0322C 80A5E51C 00000000 */  nop
/* 03230 80A5E520 E6000068 */  swc1    $f0, 0x0068($s0)           ## 00000068
.L80A5E524:
/* 03234 80A5E524 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 03238 80A5E528 260401AC */  addiu   $a0, $s0, 0x01AC           ## $a0 = 000001AC
/* 0323C 80A5E52C 50400012 */  beql    $v0, $zero, .L80A5E578     
/* 03240 80A5E530 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 03244 80A5E534 8E0F01F0 */  lw      $t7, 0x01F0($s0)           ## 000001F0
/* 03248 80A5E538 24020064 */  addiu   $v0, $zero, 0x0064         ## $v0 = 00000064
/* 0324C 80A5E53C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03250 80A5E540 31F80010 */  andi    $t8, $t7, 0x0010           ## $t8 = 00000000
/* 03254 80A5E544 13000009 */  beq     $t8, $zero, .L80A5E56C     
/* 03258 80A5E548 00000000 */  nop
/* 0325C 80A5E54C AE020150 */  sw      $v0, 0x0150($s0)           ## 00000150
/* 03260 80A5E550 0C29796B */  jal     func_80A5E5AC              
/* 03264 80A5E554 AE020154 */  sw      $v0, 0x0154($s0)           ## 00000154
/* 03268 80A5E558 8E1901F0 */  lw      $t9, 0x01F0($s0)           ## 000001F0
/* 0326C 80A5E55C 2401FFEF */  addiu   $at, $zero, 0xFFEF         ## $at = FFFFFFEF
/* 03270 80A5E560 03214024 */  and     $t0, $t9, $at              
/* 03274 80A5E564 10000003 */  beq     $zero, $zero, .L80A5E574   
/* 03278 80A5E568 AE0801F0 */  sw      $t0, 0x01F0($s0)           ## 000001F0
.L80A5E56C:
/* 0327C 80A5E56C 0C2973CA */  jal     func_80A5CF28              
/* 03280 80A5E570 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
.L80A5E574:
/* 03284 80A5E574 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80A5E578:
/* 03288 80A5E578 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 0328C 80A5E57C 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 03290 80A5E580 03E00008 */  jr      $ra                        
/* 03294 80A5E584 00000000 */  nop


