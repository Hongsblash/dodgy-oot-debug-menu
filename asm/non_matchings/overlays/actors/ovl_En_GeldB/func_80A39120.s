glabel func_80A39120
/* 03E10 80A39120 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 03E14 80A39124 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 03E18 80A39128 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 03E1C 80A3912C AFA5002C */  sw      $a1, 0x002C($sp)           
/* 03E20 80A39130 94820088 */  lhu     $v0, 0x0088($a0)           ## 00000088
/* 03E24 80A39134 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 03E28 80A39138 3C063F80 */  lui     $a2, 0x3F80                ## $a2 = 3F800000
/* 03E2C 80A3913C 304E0002 */  andi    $t6, $v0, 0x0002           ## $t6 = 00000000
/* 03E30 80A39140 11C00004 */  beq     $t6, $zero, .L80A39154     
/* 03E34 80A39144 3C073F00 */  lui     $a3, 0x3F00                ## $a3 = 3F000000
/* 03E38 80A39148 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 03E3C 80A3914C 94820088 */  lhu     $v0, 0x0088($a0)           ## 00000088
/* 03E40 80A39150 E4800068 */  swc1    $f0, 0x0068($a0)           ## 00000068
.L80A39154:
/* 03E44 80A39154 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 03E48 80A39158 304F0001 */  andi    $t7, $v0, 0x0001           ## $t7 = 00000000
/* 03E4C 80A3915C 11E00005 */  beq     $t7, $zero, .L80A39174     
/* 03E50 80A39160 26040068 */  addiu   $a0, $s0, 0x0068           ## $a0 = 00000068
/* 03E54 80A39164 44050000 */  mfc1    $a1, $f0                   
/* 03E58 80A39168 0C01E0C4 */  jal     Math_SmoothScaleMaxMinF
              
/* 03E5C 80A3916C E7A00010 */  swc1    $f0, 0x0010($sp)           
/* 03E60 80A39170 A6000318 */  sh      $zero, 0x0318($s0)         ## 00000318
.L80A39174:
/* 03E64 80A39174 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 03E68 80A39178 26040188 */  addiu   $a0, $s0, 0x0188           ## $a0 = 00000188
/* 03E6C 80A3917C 50400006 */  beql    $v0, $zero, .L80A39198     
/* 03E70 80A39180 C60401A0 */  lwc1    $f4, 0x01A0($s0)           ## 000001A0
/* 03E74 80A39184 0C28D6E3 */  jal     func_80A35B8C              
/* 03E78 80A39188 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03E7C 80A3918C 1000000E */  beq     $zero, $zero, .L80A391C8   
/* 03E80 80A39190 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 03E84 80A39194 C60401A0 */  lwc1    $f4, 0x01A0($s0)           ## 000001A0
.L80A39198:
/* 03E88 80A39198 2401000A */  addiu   $at, $zero, 0x000A         ## $at = 0000000A
/* 03E8C 80A3919C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 03E90 80A391A0 4600218D */  trunc.w.s $f6, $f4                   
/* 03E94 80A391A4 44193000 */  mfc1    $t9, $f6                   
/* 03E98 80A391A8 00000000 */  nop
/* 03E9C 80A391AC 57210006 */  bnel    $t9, $at, .L80A391C8       
/* 03EA0 80A391B0 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 03EA4 80A391B4 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 03EA8 80A391B8 2405387A */  addiu   $a1, $zero, 0x387A         ## $a1 = 0000387A
/* 03EAC 80A391BC 0C03D6D6 */  jal     func_800F5B58              
/* 03EB0 80A391C0 00000000 */  nop
/* 03EB4 80A391C4 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80A391C8:
/* 03EB8 80A391C8 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 03EBC 80A391CC 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 03EC0 80A391D0 03E00008 */  jr      $ra                        
/* 03EC4 80A391D4 00000000 */  nop


