glabel func_80A38430
/* 03120 80A38430 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 03124 80A38434 AFB00028 */  sw      $s0, 0x0028($sp)           
/* 03128 80A38438 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0312C 80A3843C AFBF002C */  sw      $ra, 0x002C($sp)           
/* 03130 80A38440 3C040600 */  lui     $a0, 0x0600                ## $a0 = 06000000
/* 03134 80A38444 0C028800 */  jal     SkelAnime_GetFrameCount
              
/* 03138 80A38448 24841578 */  addiu   $a0, $a0, 0x1578           ## $a0 = 06001578
/* 0313C 80A3844C 860E0310 */  lh      $t6, 0x0310($s0)           ## 00000310
/* 03140 80A38450 44822000 */  mtc1    $v0, $f4                   ## $f4 = 0.00
/* 03144 80A38454 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 03148 80A38458 11C00003 */  beq     $t6, $zero, .L80A38468     
/* 0314C 80A3845C 468023A0 */  cvt.s.w $f14, $f4                  
/* 03150 80A38460 240FFFFF */  addiu   $t7, $zero, 0xFFFF         ## $t7 = FFFFFFFF
/* 03154 80A38464 A60F0310 */  sh      $t7, 0x0310($s0)           ## 00000310
.L80A38468:
/* 03158 80A38468 24180006 */  addiu   $t8, $zero, 0x0006         ## $t8 = 00000006
/* 0315C 80A3846C 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 03160 80A38470 E6020068 */  swc1    $f2, 0x0068($s0)           ## 00000068
/* 03164 80A38474 AE1802EC */  sw      $t8, 0x02EC($s0)           ## 000002EC
/* 03168 80A38478 44816000 */  mtc1    $at, $f12                  ## $f12 = 10.00
/* 0316C 80A3847C 0C00CFC8 */  jal     Math_Rand_CenteredFloat
              
/* 03170 80A38480 E7AE0034 */  swc1    $f14, 0x0034($sp)          
/* 03174 80A38484 4600018D */  trunc.w.s $f6, $f0                   
/* 03178 80A38488 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 0317C 80A3848C C7AE0034 */  lwc1    $f14, 0x0034($sp)          
/* 03180 80A38490 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 03184 80A38494 44083000 */  mfc1    $t0, $f6                   
/* 03188 80A38498 240A0002 */  addiu   $t2, $zero, 0x0002         ## $t2 = 00000002
/* 0318C 80A3849C 44061000 */  mfc1    $a2, $f2                   
/* 03190 80A384A0 2509000A */  addiu   $t1, $t0, 0x000A           ## $t1 = 0000000A
/* 03194 80A384A4 AE090300 */  sw      $t1, 0x0300($s0)           ## 00000300
/* 03198 80A384A8 44071000 */  mfc1    $a3, $f2                   
/* 0319C 80A384AC AFAA0014 */  sw      $t2, 0x0014($sp)           
/* 031A0 80A384B0 24A51578 */  addiu   $a1, $a1, 0x1578           ## $a1 = 06001578
/* 031A4 80A384B4 26040188 */  addiu   $a0, $s0, 0x0188           ## $a0 = 00000188
/* 031A8 80A384B8 E7A20018 */  swc1    $f2, 0x0018($sp)           
/* 031AC 80A384BC 0C029468 */  jal     SkelAnime_ChangeAnimation
              
/* 031B0 80A384C0 E7AE0010 */  swc1    $f14, 0x0010($sp)          
/* 031B4 80A384C4 3C0580A4 */  lui     $a1, %hi(func_80A384E8)    ## $a1 = 80A40000
/* 031B8 80A384C8 24A584E8 */  addiu   $a1, $a1, %lo(func_80A384E8) ## $a1 = 80A384E8
/* 031BC 80A384CC 0C28D4C4 */  jal     func_80A35310              
/* 031C0 80A384D0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 031C4 80A384D4 8FBF002C */  lw      $ra, 0x002C($sp)           
/* 031C8 80A384D8 8FB00028 */  lw      $s0, 0x0028($sp)           
/* 031CC 80A384DC 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 031D0 80A384E0 03E00008 */  jr      $ra                        
/* 031D4 80A384E4 00000000 */  nop


