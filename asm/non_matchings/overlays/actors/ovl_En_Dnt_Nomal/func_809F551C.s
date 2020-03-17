glabel func_809F551C
/* 0225C 809F551C 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 02260 809F5520 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 02264 809F5524 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 02268 809F5528 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 0226C 809F552C AFA50034 */  sw      $a1, 0x0034($sp)           
/* 02270 809F5530 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 02274 809F5534 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 02278 809F5538 C604021C */  lwc1    $f4, 0x021C($s0)           ## 0000021C
/* 0227C 809F553C C6060024 */  lwc1    $f6, 0x0024($s0)           ## 00000024
/* 02280 809F5540 C6080224 */  lwc1    $f8, 0x0224($s0)           ## 00000224
/* 02284 809F5544 C60A002C */  lwc1    $f10, 0x002C($s0)          ## 0000002C
/* 02288 809F5548 46062301 */  sub.s   $f12, $f4, $f6             
/* 0228C 809F554C 460A4381 */  sub.s   $f14, $f8, $f10            
/* 02290 809F5550 E7AC002C */  swc1    $f12, 0x002C($sp)          
/* 02294 809F5554 0C03F494 */  jal     func_800FD250              
/* 02298 809F5558 E7AE0028 */  swc1    $f14, 0x0028($sp)          
/* 0229C 809F555C 3C01809F */  lui     $at, %hi(D_809F6068)       ## $at = 809F0000
/* 022A0 809F5560 C4306068 */  lwc1    $f16, %lo(D_809F6068)($at) 
/* 022A4 809F5564 260400B6 */  addiu   $a0, $s0, 0x00B6           ## $a0 = 000000B6
/* 022A8 809F5568 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 022AC 809F556C 46100482 */  mul.s   $f18, $f0, $f16            
/* 022B0 809F5570 24070BB8 */  addiu   $a3, $zero, 0x0BB8         ## $a3 = 00000BB8
/* 022B4 809F5574 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 022B8 809F5578 4600910D */  trunc.w.s $f4, $f18                  
/* 022BC 809F557C 44052000 */  mfc1    $a1, $f4                   
/* 022C0 809F5580 00000000 */  nop
/* 022C4 809F5584 00052C00 */  sll     $a1, $a1, 16               
/* 022C8 809F5588 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 022CC 809F558C 00052C03 */  sra     $a1, $a1, 16               
/* 022D0 809F5590 8602025E */  lh      $v0, 0x025E($s0)           ## 0000025E
/* 022D4 809F5594 240F000A */  addiu   $t7, $zero, 0x000A         ## $t7 = 0000000A
/* 022D8 809F5598 14400003 */  bne     $v0, $zero, .L809F55A8     
/* 022DC 809F559C 30580001 */  andi    $t8, $v0, 0x0001           ## $t8 = 00000000
/* 022E0 809F55A0 10000005 */  beq     $zero, $zero, .L809F55B8   
/* 022E4 809F55A4 A60F025E */  sh      $t7, 0x025E($s0)           ## 0000025E
.L809F55A8:
/* 022E8 809F55A8 17000003 */  bne     $t8, $zero, .L809F55B8     
/* 022EC 809F55AC 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 022F0 809F55B0 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 022F4 809F55B4 2405387F */  addiu   $a1, $zero, 0x387F         ## $a1 = 0000387F
.L809F55B8:
/* 022F8 809F55B8 C7A0002C */  lwc1    $f0, 0x002C($sp)           
/* 022FC 809F55BC 3C0140E0 */  lui     $at, 0x40E0                ## $at = 40E00000
/* 02300 809F55C0 44811000 */  mtc1    $at, $f2                   ## $f2 = 7.00
/* 02304 809F55C4 46000005 */  abs.s   $f0, $f0                   
/* 02308 809F55C8 4602003C */  c.lt.s  $f0, $f2                   
/* 0230C 809F55CC C7A00028 */  lwc1    $f0, 0x0028($sp)           
/* 02310 809F55D0 4502000F */  bc1fl   .L809F5610                 
/* 02314 809F55D4 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 02318 809F55D8 46000005 */  abs.s   $f0, $f0                   
/* 0231C 809F55DC 3C19809F */  lui     $t9, %hi(func_809F4E18)    ## $t9 = 809F0000
/* 02320 809F55E0 4602003C */  c.lt.s  $f0, $f2                   
/* 02324 809F55E4 27394E18 */  addiu   $t9, $t9, %lo(func_809F4E18) ## $t9 = 809F4E18
/* 02328 809F55E8 45020009 */  bc1fl   .L809F5610                 
/* 0232C 809F55EC 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 02330 809F55F0 C606021C */  lwc1    $f6, 0x021C($s0)           ## 0000021C
/* 02334 809F55F4 C6080224 */  lwc1    $f8, 0x0224($s0)           ## 00000224
/* 02338 809F55F8 44805000 */  mtc1    $zero, $f10                ## $f10 = 0.00
/* 0233C 809F55FC AE190214 */  sw      $t9, 0x0214($s0)           ## 00000214
/* 02340 809F5600 E6060024 */  swc1    $f6, 0x0024($s0)           ## 00000024
/* 02344 809F5604 E608002C */  swc1    $f8, 0x002C($s0)           ## 0000002C
/* 02348 809F5608 E60A0068 */  swc1    $f10, 0x0068($s0)          ## 00000068
/* 0234C 809F560C 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L809F5610:
/* 02350 809F5610 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 02354 809F5614 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000
/* 02358 809F5618 03E00008 */  jr      $ra                        
/* 0235C 809F561C 00000000 */  nop


