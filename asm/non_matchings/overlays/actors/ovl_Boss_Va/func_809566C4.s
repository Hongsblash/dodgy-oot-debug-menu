glabel func_809566C4
/* 07404 809566C4 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 07408 809566C8 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 0740C 809566CC 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 07410 809566D0 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 07414 809566D4 AFA5002C */  sw      $a1, 0x002C($sp)           
/* 07418 809566D8 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 0741C 809566DC 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 07420 809566E0 0C253CB2 */  jal     func_8094F2C8              
/* 07424 809566E4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 07428 809566E8 260401E6 */  addiu   $a0, $s0, 0x01E6           ## $a0 = 000001E6
/* 0742C 809566EC 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 07430 809566F0 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 07434 809566F4 24071770 */  addiu   $a3, $zero, 0x1770         ## $a3 = 00001770
/* 07438 809566F8 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 0743C 809566FC AFA00010 */  sw      $zero, 0x0010($sp)         
/* 07440 80956700 260401E4 */  addiu   $a0, $s0, 0x01E4           ## $a0 = 000001E4
/* 07444 80956704 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 07448 80956708 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 0744C 8095670C 24071770 */  addiu   $a3, $zero, 0x1770         ## $a3 = 00001770
/* 07450 80956710 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 07454 80956714 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 07458 80956718 260401EC */  addiu   $a0, $s0, 0x01EC           ## $a0 = 000001EC
/* 0745C 8095671C 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 07460 80956720 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 07464 80956724 24071770 */  addiu   $a3, $zero, 0x1770         ## $a3 = 00001770
/* 07468 80956728 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 0746C 8095672C AFA00010 */  sw      $zero, 0x0010($sp)         
/* 07470 80956730 260401EA */  addiu   $a0, $s0, 0x01EA           ## $a0 = 000001EA
/* 07474 80956734 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 07478 80956738 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 0747C 8095673C 24071770 */  addiu   $a3, $zero, 0x1770         ## $a3 = 00001770
/* 07480 80956740 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 07484 80956744 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 07488 80956748 860500B6 */  lh      $a1, 0x00B6($s0)           ## 000000B6
/* 0748C 8095674C AFA00010 */  sw      $zero, 0x0010($sp)         
/* 07490 80956750 260401F2 */  addiu   $a0, $s0, 0x01F2           ## $a0 = 000001F2
/* 07494 80956754 24A5C000 */  addiu   $a1, $a1, 0xC000           ## $a1 = FFFFC000
/* 07498 80956758 00052C00 */  sll     $a1, $a1, 16               
/* 0749C 8095675C 00052C03 */  sra     $a1, $a1, 16               
/* 074A0 80956760 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 074A4 80956764 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 074A8 80956768 24072710 */  addiu   $a3, $zero, 0x2710         ## $a3 = 00002710
/* 074AC 8095676C 8E0E016C */  lw      $t6, 0x016C($s0)           ## 0000016C
/* 074B0 80956770 260401F0 */  addiu   $a0, $s0, 0x01F0           ## $a0 = 000001F0
/* 074B4 80956774 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
/* 074B8 80956778 85C5002E */  lh      $a1, 0x002E($t6)           ## 0000002E
/* 074BC 8095677C AFA00010 */  sw      $zero, 0x0010($sp)         
/* 074C0 80956780 24071770 */  addiu   $a3, $zero, 0x1770         ## $a3 = 00001770
/* 074C4 80956784 24A5EC78 */  addiu   $a1, $a1, 0xEC78           ## $a1 = FFFFEC78
/* 074C8 80956788 00052C00 */  sll     $a1, $a1, 16               
/* 074CC 8095678C 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 074D0 80956790 00052C03 */  sra     $a1, $a1, 16               
/* 074D4 80956794 8E0F0118 */  lw      $t7, 0x0118($s0)           ## 00000118
/* 074D8 80956798 44802000 */  mtc1    $zero, $f4                 ## $f4 = 0.00
/* 074DC 8095679C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 074E0 809567A0 C5E60068 */  lwc1    $f6, 0x0068($t7)           ## 00000068
/* 074E4 809567A4 46062032 */  c.eq.s  $f4, $f6                   
/* 074E8 809567A8 00000000 */  nop
/* 074EC 809567AC 45020004 */  bc1fl   .L809567C0                 
/* 074F0 809567B0 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 074F4 809567B4 0C2552E3 */  jal     func_80954B8C              
/* 074F8 809567B8 8FA5002C */  lw      $a1, 0x002C($sp)           
/* 074FC 809567BC 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L809567C0:
/* 07500 809567C0 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 07504 809567C4 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 07508 809567C8 03E00008 */  jr      $ra                        
/* 0750C 809567CC 00000000 */  nop


