glabel func_80B559C4
/* 02614 80B559C4 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 02618 80B559C8 AFA40020 */  sw      $a0, 0x0020($sp)           
/* 0261C 80B559CC AFBF001C */  sw      $ra, 0x001C($sp)           
/* 02620 80B559D0 3C040600 */  lui     $a0, 0x0600                ## $a0 = 06000000
/* 02624 80B559D4 0C028800 */  jal     SkelAnime_GetFrameCount
              
/* 02628 80B559D8 24845248 */  addiu   $a0, $a0, 0x5248           ## $a0 = 06005248
/* 0262C 80B559DC 8FA80020 */  lw      $t0, 0x0020($sp)           
/* 02630 80B559E0 240F0003 */  addiu   $t7, $zero, 0x0003         ## $t7 = 00000003
/* 02634 80B559E4 3044FFFF */  andi    $a0, $v0, 0xFFFF           ## $a0 = 00000000
/* 02638 80B559E8 C5040164 */  lwc1    $f4, 0x0164($t0)           ## 00000164
/* 0263C 80B559EC AFAF0010 */  sw      $t7, 0x0010($sp)           
/* 02640 80B559F0 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 02644 80B559F4 4600218D */  trunc.w.s $f6, $f4                   
/* 02648 80B559F8 24070003 */  addiu   $a3, $zero, 0x0003         ## $a3 = 00000003
/* 0264C 80B559FC 44063000 */  mfc1    $a2, $f6                   
/* 02650 80B55A00 0C01BE6F */  jal     func_8006F9BC              
/* 02654 80B55A04 30C6FFFF */  andi    $a2, $a2, 0xFFFF           ## $a2 = 00000000
/* 02658 80B55A08 8FA80020 */  lw      $t0, 0x0020($sp)           
/* 0265C 80B55A0C 2503032C */  addiu   $v1, $t0, 0x032C           ## $v1 = 0000032C
/* 02660 80B55A10 25040338 */  addiu   $a0, $t0, 0x0338           ## $a0 = 00000338
/* 02664 80B55A14 C4880000 */  lwc1    $f8, 0x0000($a0)           ## 00000338
/* 02668 80B55A18 C4620000 */  lwc1    $f2, 0x0000($v1)           ## 0000032C
/* 0266C 80B55A1C 25020024 */  addiu   $v0, $t0, 0x0024           ## $v0 = 00000024
/* 02670 80B55A20 46024281 */  sub.s   $f10, $f8, $f2             
/* 02674 80B55A24 460A0402 */  mul.s   $f16, $f0, $f10            
/* 02678 80B55A28 46101480 */  add.s   $f18, $f2, $f16            
/* 0267C 80B55A2C E4520000 */  swc1    $f18, 0x0000($v0)          ## 00000024
/* 02680 80B55A30 C4840008 */  lwc1    $f4, 0x0008($a0)           ## 00000340
/* 02684 80B55A34 C46C0008 */  lwc1    $f12, 0x0008($v1)          ## 00000334
/* 02688 80B55A38 460C2181 */  sub.s   $f6, $f4, $f12             
/* 0268C 80B55A3C 46060202 */  mul.s   $f8, $f0, $f6              
/* 02690 80B55A40 46086280 */  add.s   $f10, $f12, $f8            
/* 02694 80B55A44 E44A0008 */  swc1    $f10, 0x0008($v0)          ## 0000002C
/* 02698 80B55A48 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 0269C 80B55A4C 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 026A0 80B55A50 03E00008 */  jr      $ra                        
/* 026A4 80B55A54 00000000 */  nop


